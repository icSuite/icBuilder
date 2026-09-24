"""Collect CS-covered DMSP samples with detector precipitation values."""

#%% Imports

import argparse
from multiprocessing import get_context
from pathlib import Path

from icreader import open_product
from icbuilder.footprints import (
    infer_footprints,
    overlap_mapping,
    point_in_quadrilaterals,
)
import numpy as np
from tqdm import tqdm
import xarray as xr


#%% Defaults

DEFAULT_DMSP_PATH = Path(
    "/home/bing/Dropbox/work/data/dmsp/dmsp_ssj_yearly"
)
DEFAULT_IMAGE_PATH = Path(
    "/home/bing/Dropbox/work/data/IMAGE_FUV/"
    "precipitation_detector/IR_hardy"
)
DEFAULT_CS_PATH = Path(
    "/home/bing/Dropbox/work/data/IMAGE_FUV/precipitation_cs/"
    "image_apex_130km_46x46_v1/IR_hardy"
)
DEFAULT_OUTPUT_PATH = Path(
    "/home/bing/Dropbox/work/data/IMAGE_FUV/"
    "ratio_relation_validation/detector_crossings"
)

SATELLITES = ("f12", "f13", "f14", "f15")
SUPPORT_SECONDS = 60

IMAGE_FIELDS = {
    "wic_corrected": "img_wic",
    "dwic_corrected": "img_wic_std",
    "si13_corrected": "img_si13",
    "dsi13_corrected": "img_si13_std",
    "R": "img_ratio",
    "dR": "img_ratio_std",
    "E0": "img_energy",
    "dE0": "img_energy_std",
    "dza": "wic_dza",
    "sza": "wic_sza",
    "method_quality_weight": "quality_weight",
    "method_valid": "method_valid",
}

CS_GATE_FIELDS = {
    "cs_row", "cs_column", "cs_wic_coverage", "cs_si13_coverage",
    "detector_footprint_match_count", "detector_pair_valid",
}
OUTPUT_FIELDS = (
    set(IMAGE_FIELDS.values())
    | CS_GATE_FIELDS
    | {"orbit", "wic_source_index"}
)

SPATIAL_FILTER = (
    "DMSP location inside the source precipitation_cs grid with positive "
    "WIC and SI13 coverage and finite binned WIC/SI13 data at the IMAGE "
    "frame time, followed by DMSP-point containment in an inferred WIC "
    "detector-pixel footprint"
)

DMSP_FIELDS = (
    "electron_mean_energy",
    "electron_mean_energy_fractional_std",
    "electron_total_energy_flux",
    "electron_total_energy_flux_fractional_std",
)


#%% DMSP data

def index_dmsp_files(path):
    """Index yearly northern-hemisphere files by satellite and year."""

    files = {}
    for filename in path.glob("dmsp_*_ssj_*_north.nc"):
        parts = filename.stem.split("_")
        if len(parts) == 5 and parts[1].lower() in SATELLITES:
            files[parts[1].lower(), int(parts[3])] = filename
    return files


def load_dmsp(files, cache, satellite, start, stop):
    """Load the short interval needed for one IMAGE orbit."""

    start_year = int(np.datetime_as_string(start, unit="Y"))
    stop_year = int(np.datetime_as_string(stop, unit="Y"))
    pieces = []

    for year in range(start_year, stop_year + 1):
        filename = files.get((satellite, year))
        if filename is None:
            continue
        if filename not in cache:
            cache[filename] = xr.open_dataset(filename)
        pieces.append(cache[filename].sel(time=slice(start, stop)).load())

    if not pieces:
        return None
    if len(pieces) == 1:
        return pieces[0]
    return xr.concat(pieces, dim="time").sortby("time")


#%% CS spatial gate and detector matching

def unit_vectors(mlat, mlt):
    """Convert Modified Apex latitude and MLT to unit-sphere vectors."""

    latitude = np.deg2rad(np.asarray(mlat, dtype=float))
    longitude = np.deg2rad(np.asarray(mlt, dtype=float) * 15)
    cos_latitude = np.cos(latitude)
    return np.column_stack((
        cos_latitude * np.cos(longitude),
        cos_latitude * np.sin(longitude),
        np.sin(latitude),
    ))


def cs_supported_samples(
    image, frame, dmsp, wic_coverage, si13_coverage, binned_wic, binned_si13
):
    """Keep DMSP samples inside CS cells covered by both IMAGE sensors."""

    dmsp_mlat = np.asarray(dmsp.mlat.values, dtype=float)
    dmsp_mlt = np.asarray(dmsp.mlt.values, dtype=float)
    finite = np.isfinite(dmsp_mlat) & np.isfinite(dmsp_mlt)
    finite_indices = np.flatnonzero(finite)
    if not finite_indices.size:
        return dmsp.isel(time=slice(0, 0))

    longitude = dmsp_mlt[finite] * 15
    latitude = dmsp_mlat[finite]
    inside = np.asarray(image.grid.ingrid(longitude, latitude), dtype=bool)
    inside_indices = finite_indices[inside]
    if not inside_indices.size:
        return dmsp.isel(time=slice(0, 0))

    row, column = image.grid.bin_index(longitude[inside], latitude[inside])
    row = np.asarray(row, dtype=np.int16)
    column = np.asarray(column, dtype=np.int16)
    wic = np.asarray(wic_coverage[frame, row, column], dtype=float)
    si13 = np.asarray(si13_coverage[frame, row, column], dtype=float)
    wic_values = np.asarray(binned_wic[frame, row, column], dtype=float)
    si13_values = np.asarray(binned_si13[frame, row, column], dtype=float)
    supported = (
        np.isfinite(wic) & (wic > 0)
        & np.isfinite(si13) & (si13 > 0)
        & np.isfinite(wic_values) & np.isfinite(si13_values)
    )

    selected = dmsp.isel(time=inside_indices[supported]).copy()
    selected["cs_row"] = ("time", row[supported])
    selected["cs_column"] = ("time", column[supported])
    selected["cs_wic_coverage"] = ("time", wic[supported].astype(np.float32))
    selected["cs_si13_coverage"] = (
        "time", si13[supported].astype(np.float32)
    )
    return selected


def make_detector_footprint_matcher(image_mlat, image_mlt, grid):
    """Prepare projected detector footprints and their CS-cell candidates."""

    corners, _ = infer_footprints(image_mlat, image_mlt, grid)
    mapping, _ = overlap_mapping(image_mlat, image_mlt, grid)
    corners = corners.reshape(-1, 4, 2)
    centres = np.mean(corners, axis=1)
    return corners, centres, mapping


def candidate_footprints(mapping, grid_shape, row, column):
    """Return footprints intersecting the point's CS cell or a neighbour."""

    ny, nx = grid_shape
    candidates = []
    for nearby_row in range(max(row - 1, 0), min(row + 2, ny)):
        for nearby_column in range(max(column - 1, 0), min(column + 2, nx)):
            target = nearby_row * nx + nearby_column
            start, stop = mapping.indptr[target:target + 2]
            candidates.append(mapping.indices[start:stop])
    if not candidates:
        return np.empty(0, dtype=int)
    return np.unique(np.concatenate(candidates))


def match_detector_footprints(
    matcher, grid, cs_row, cs_column, dmsp_mlat, dmsp_mlt,
    image_mlat, image_mlt,
):
    """Find the detector quadrilateral containing each DMSP location."""

    corners, centres, mapping = matcher
    longitude = np.asarray(dmsp_mlt, dtype=float) * 15
    latitude = np.asarray(dmsp_mlat, dtype=float)
    point_x, point_y = grid.projection.geo2cube(
        longitude, latitude, set_points_off_cube_to_nan=True
    )
    point_x = np.asarray(point_x, dtype=float)
    point_y = np.asarray(point_y, dtype=float)

    flat_pixel = np.full(point_x.size, -1, dtype=int)
    match_count = np.zeros(point_x.size, dtype=np.int16)
    for sample in range(point_x.size):
        candidates = candidate_footprints(
            mapping, grid.shape, int(cs_row[sample]), int(cs_column[sample])
        )
        inside = point_in_quadrilaterals(
            point_x[sample], point_y[sample], corners[candidates]
        )
        containing = candidates[inside]
        match_count[sample] = containing.size
        if not containing.size:
            continue
        if containing.size == 1:
            flat_pixel[sample] = containing[0]
            continue

        distance_squared = (
            (centres[containing, 0] - point_x[sample]) ** 2
            + (centres[containing, 1] - point_y[sample]) ** 2
        )
        flat_pixel[sample] = containing[np.argmin(distance_squared)]

    contained = flat_pixel >= 0
    row = np.full(point_x.size, -1, dtype=np.int16)
    column = np.full(point_x.size, -1, dtype=np.int16)
    separation = np.full(point_x.size, np.nan, dtype=np.float32)
    if np.any(contained):
        matched_row, matched_column = np.unravel_index(
            flat_pixel[contained], image_mlat.shape
        )
        row[contained] = matched_row
        column[contained] = matched_column

        dmsp_vectors = unit_vectors(
            latitude[contained], np.asarray(dmsp_mlt)[contained]
        )
        pixel_vectors = unit_vectors(
            image_mlat[matched_row, matched_column],
            image_mlt[matched_row, matched_column],
        )
        cosine = np.sum(dmsp_vectors * pixel_vectors, axis=1)
        separation[contained] = np.rad2deg(
            np.arccos(np.clip(cosine, -1, 1))
        ).astype(np.float32)

    return row, column, separation, match_count


def append_frame_samples(
    samples, image, frame, satellite, dmsp, mlat, mlt, matcher, grid
):
    """Append one frame/satellite match to ordinary NumPy lists."""

    dmsp_mlat = np.asarray(dmsp.mlat.values, dtype=float)
    dmsp_mlt = np.asarray(dmsp.mlt.values, dtype=float)
    finite = np.isfinite(dmsp_mlat) & np.isfinite(dmsp_mlt)
    if not np.any(finite):
        return 0

    dmsp_time = np.asarray(dmsp.time.values, dtype="datetime64[ns]")[finite]
    dmsp_mlat = dmsp_mlat[finite]
    dmsp_mlt = dmsp_mlt[finite]
    row, column, separation, footprint_count = match_detector_footprints(
        matcher,
        grid,
        np.asarray(dmsp.cs_row.values)[finite],
        np.asarray(dmsp.cs_column.values)[finite],
        dmsp_mlat,
        dmsp_mlt,
        mlat[frame],
        mlt[frame],
    )
    contained = footprint_count > 0
    missed = int(np.count_nonzero(~contained))
    if not np.any(contained):
        return missed

    dmsp_time = dmsp_time[contained]
    dmsp_mlat = dmsp_mlat[contained]
    dmsp_mlt = dmsp_mlt[contained]
    row = row[contained]
    column = column[contained]
    separation = separation[contained]
    footprint_count = footprint_count[contained]
    count = dmsp_time.size
    image_time = np.datetime64(image.time[frame], "ns")

    samples["image_frame"].append(np.full(count, frame, dtype=np.int32))
    samples["image_time"].append(np.full(count, image_time, dtype="datetime64[ns]"))
    samples["dmsp_sat"].append(np.full(count, satellite, dtype="U3"))
    samples["dmsp_time"].append(dmsp_time)
    samples["time_offset_seconds"].append(
        ((dmsp_time - image_time) / np.timedelta64(1, "s")).astype(np.float32)
    )
    samples["dmsp_mlat"].append(dmsp_mlat.astype(np.float32))
    samples["dmsp_mlt"].append(dmsp_mlt.astype(np.float32))
    samples["detector_row"].append(row.astype(np.int16))
    samples["detector_column"].append(column.astype(np.int16))
    samples["detector_separation_deg"].append(separation.astype(np.float32))
    samples["detector_footprint_match_count"].append(footprint_count)
    for name in (
        "cs_row", "cs_column", "cs_wic_coverage", "cs_si13_coverage"
    ):
        samples[name].append(np.asarray(dmsp[name].values)[finite][contained])

    for name in DMSP_FIELDS:
        samples[f"dmsp_{name}"].append(
            np.asarray(dmsp[name].values, dtype=np.float32)[finite][contained]
        )
    return missed


def make_samples():
    names = [
        "image_frame", "image_time", "dmsp_sat", "dmsp_time",
        "time_offset_seconds", "dmsp_mlat", "dmsp_mlt",
        "detector_row", "detector_column", "detector_separation_deg",
        "detector_footprint_match_count",
        "cs_row", "cs_column", "cs_wic_coverage", "cs_si13_coverage",
    ]
    names += [f"dmsp_{name}" for name in DMSP_FIELDS]
    return {name: [] for name in names}


#%% One output file per IMAGE orbit

def process_orbit(image_file, cs_file, dmsp_files, dmsp_cache):
    """Collect CS-gated DMSP matches and sample the detector product."""

    orbit = int(image_file.stem.split("_")[-1])
    samples = make_samples()
    candidates = []

    # The compact CS product establishes temporal and spatial overlap before
    # the much larger detector product is opened.
    with open_product(cs_file) as cs_image:
        if (
            cs_image.product_type != "precipitation_cs"
            or cs_image.representation != "cs"
        ):
            raise ValueError(f"{cs_file} is not CS Product 2")

        cs_times = np.asarray(cs_image.time, dtype="datetime64[ns]")
        start = cs_times[0] - np.timedelta64(SUPPORT_SECONDS, "s")
        stop = cs_times[-1] + np.timedelta64(SUPPORT_SECONDS, "s")
        dmsp = {
            satellite: load_dmsp(
                dmsp_files, dmsp_cache, satellite, start, stop
            )
            for satellite in SATELLITES
        }
        wic_coverage = np.asarray(cs_image.read("wic_coverage"), dtype=float)
        si13_coverage = np.asarray(cs_image.read("si13_coverage"), dtype=float)
        binned_wic = np.asarray(cs_image.read("wic_corrected"), dtype=float)
        binned_si13 = np.asarray(cs_image.read("si13_corrected"), dtype=float)

        for frame, image_time in enumerate(cs_times):
            frame_start = image_time - np.timedelta64(SUPPORT_SECONDS, "s")
            frame_stop = image_time + np.timedelta64(SUPPORT_SECONDS, "s")
            for satellite, data in dmsp.items():
                if data is None:
                    continue
                subset = data.sel(time=slice(frame_start, frame_stop))
                if subset.sizes["time"]:
                    supported = cs_supported_samples(
                        cs_image, frame, subset, wic_coverage, si13_coverage,
                        binned_wic, binned_si13
                    )
                    if supported.sizes["time"]:
                        candidates.append((frame, satellite, supported))

        cs_grid_id = cs_image.grid_id
        cs_grid = cs_image.grid

    cs_gated_count = sum(
        candidate.sizes["time"] for _, _, candidate in candidates
    )

    empty_attrs = {
        "product_type": "dmsp_image_crossings",
        "schema_version": 1,
        "representation": "detector",
        "frame_half_support_seconds": SUPPORT_SECONDS,
        "spatial_filter": SPATIAL_FILTER,
        "source_precipitation_detector": str(image_file.resolve()),
        "source_precipitation_cs": str(cs_file.resolve()),
        "source_precipitation_cs_grid_id": cs_grid_id,
        "cs_gated_sample_count": cs_gated_count,
        "detector_footprint_miss_count": 0,
    }
    if not candidates:
        return xr.Dataset(
            coords={"sample": np.arange(0)}, attrs=empty_attrs
        )

    with open_product(image_file) as image:
        if image.product_type != "precipitation_detector":
            raise ValueError(f"{image_file} is not detector Product 2")

        image_times = np.asarray(image.time, dtype="datetime64[ns]")
        if not np.array_equal(image_times, cs_times):
            raise ValueError(
                f"Detector and CS frame times do not match for orbit {orbit:04d}"
            )

        mlat = np.asarray(image.read("mlat"), dtype=float)
        mlt = np.asarray(image.read("mlt"), dtype=float)
        footprint_frame = None
        footprint_matcher = None
        footprint_misses = 0
        for frame, satellite, supported in candidates:
            if frame != footprint_frame:
                footprint_matcher = make_detector_footprint_matcher(
                    mlat[frame], mlt[frame], cs_grid
                )
                footprint_frame = frame
            footprint_misses += append_frame_samples(
                samples, image, frame, satellite, supported, mlat, mlt,
                footprint_matcher, cs_grid
            )

        if not samples["image_frame"]:
            empty_attrs["detector_footprint_miss_count"] = footprint_misses
            return xr.Dataset(
                coords={"sample": np.arange(0)}, attrs=empty_attrs
            )

        values = {name: np.concatenate(parts) for name, parts in samples.items()}
        frame = values["image_frame"]
        row = values["detector_row"]
        column = values["detector_column"]

        for source_name, output_name in IMAGE_FIELDS.items():
            cube = np.asarray(image.read(source_name))
            values[output_name] = cube[frame, row, column]

        values["detector_pair_valid"] = (
            np.isfinite(values["img_wic"])
            & np.isfinite(values["img_si13"])
        )

        count = frame.size
        values["orbit"] = np.full(count, orbit, dtype=np.int32)
        values["wic_source_index"] = np.asarray(
            image.wic_source_index, dtype=np.int32
        )[frame]

        result = xr.Dataset({
            name: ("sample", value) for name, value in values.items()
        })
        result.attrs.update({
            "product_type": "dmsp_image_crossings",
            "schema_version": 1,
            "representation": "detector",
            "frame_half_support_seconds": SUPPORT_SECONDS,
            "matching_method": (
                "CS coverage gate followed by DMSP-point containment in an "
                "inferred WIC detector-pixel footprint"
            ),
            "spatial_filter": SPATIAL_FILTER,
            "source_precipitation_detector": str(image_file.resolve()),
            "source_precipitation_cs": str(cs_file.resolve()),
            "source_precipitation_cs_grid_id": cs_grid_id,
            "source_fuv_detector_sha256": image.source_fuv_detector_sha256,
            "source_preprocessing_label": image.source_preprocessing_label,
            "proton_energy_model": image.proton_energy_model,
            "source_count_source": image.attrs.get(
                "count_source", "background_subtracted"
            ),
            "cs_gated_sample_count": cs_gated_count,
            "detector_footprint_miss_count": footprint_misses,
            "detector_footprint_ambiguous_count": int(np.count_nonzero(
                values["detector_footprint_match_count"] > 1
            )),
            "detector_pair_missing_count": int(
                np.count_nonzero(~values["detector_pair_valid"])
            ),
        })
        return result


def save_orbit(data, filename):
    """Write one orbit atomically so interrupted runs cannot look complete."""

    filename.parent.mkdir(parents=True, exist_ok=True)
    temporary = filename.with_suffix(filename.suffix + ".partial")
    encoding = {
        name: {"zlib": True, "complevel": 4}
        for name in data.data_vars
        if data[name].dtype.kind not in ("U", "O", "M")
    }
    data.to_netcdf(temporary, engine="netcdf4", encoding=encoding)
    temporary.replace(filename)


def output_is_complete(
    filename, required_fields=OUTPUT_FIELDS, required_spatial_filter=None
):
    """Return whether an existing orbit contains the completed output fields."""

    if not filename.exists():
        return False
    try:
        with xr.open_dataset(filename) as data:
            if (
                required_spatial_filter is not None
                and data.attrs.get("spatial_filter") != required_spatial_filter
            ):
                return False
            if data.sizes.get("sample", 0) == 0:
                return True
            return required_fields.issubset(data.variables)
    except (OSError, ValueError):
        return False


def process_and_save(image_file, cs_file, output_path, dmsp_files, dmsp_cache):
    """Process one orbit and return its orbit number and match count."""

    orbit = int(image_file.stem.split("_")[-1])
    result = process_orbit(image_file, cs_file, dmsp_files, dmsp_cache)
    save_orbit(result, output_path / f"or_{orbit:04d}.nc")
    count = result.sizes["sample"]
    footprint_misses = int(result.attrs.get("detector_footprint_miss_count", 0))
    missing = int(result.attrs.get("detector_pair_missing_count", 0))
    result.close()
    return orbit, count, footprint_misses, missing


_worker_dmsp_files = None
_worker_dmsp_cache = None


def initialize_worker(dmsp_path):
    """Give each process independent DMSP file handles and a local cache."""

    global _worker_dmsp_files, _worker_dmsp_cache
    _worker_dmsp_files = index_dmsp_files(dmsp_path)
    _worker_dmsp_cache = {}


def process_in_worker(task):
    """Multiprocessing entry point; all netCDF files are opened here."""

    image_file, cs_file, output_path = task
    return process_and_save(
        image_file,
        cs_file,
        output_path,
        _worker_dmsp_files,
        _worker_dmsp_cache,
    )


#%% Run

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dmsp-path", type=Path, default=DEFAULT_DMSP_PATH)
    parser.add_argument("--image-path", type=Path, default=DEFAULT_IMAGE_PATH)
    parser.add_argument(
        "--cs-path", type=Path, default=DEFAULT_CS_PATH,
        help="directory containing compact precipitation_cs orbit files",
    )
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--orbit", type=int, nargs="+")
    parser.add_argument(
        "--workers", type=int, default=1,
        help="number of independent orbit workers; 1 runs serially",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")

    image_files = sorted(args.image_path.glob("or_*.nc"))
    if args.orbit:
        selected = set(args.orbit)
        image_files = [
            file for file in image_files
            if int(file.stem.split("_")[-1]) in selected
        ]
    if not image_files:
        raise FileNotFoundError(f"No orbit files found in {args.image_path}")

    dmsp_files = index_dmsp_files(args.dmsp_path)
    if not dmsp_files:
        raise FileNotFoundError(f"No yearly DMSP files found in {args.dmsp_path}")

    tasks = []
    for image_file in image_files:
        orbit = int(image_file.stem.split("_")[-1])
        cs_file = args.cs_path / image_file.name
        if not cs_file.exists():
            raise FileNotFoundError(f"Missing CS gate product: {cs_file}")
        output_file = args.output_path / f"or_{orbit:04d}.nc"
        if (
            output_is_complete(
                output_file, required_spatial_filter=SPATIAL_FILTER
            )
            and not args.overwrite
        ):
            continue
        tasks.append((image_file, cs_file, args.output_path))

    if args.workers == 1:
        dmsp_cache = {}
        try:
            for image_file, cs_file, output_path in tqdm(
                tasks, desc="DMSP crossings"
            ):
                orbit, count, footprint_misses, missing = process_and_save(
                    image_file, cs_file, output_path, dmsp_files, dmsp_cache
                )
                print(f"Orbit {orbit:04d}: {count} matches")
                if footprint_misses:
                    print(
                        f" Orbit {orbit:04d}: {footprint_misses} CS-gated "
                        "samples fall outside every detector footprint"
                    )
                if missing:
                    print(
                        f" Orbit {orbit:04d}: {missing} CS-supported samples "
                        "lack a detector WIC/SI13 pair"
                    )
        finally:
            for data in dmsp_cache.values():
                data.close()
    else:
        context = get_context("spawn")
        with context.Pool(
            args.workers,
            initializer=initialize_worker,
            initargs=(args.dmsp_path,),
        ) as pool:
            results = pool.imap_unordered(process_in_worker, tasks)
            for orbit, count, footprint_misses, missing in tqdm(
                results, total=len(tasks), desc="DMSP crossings"
            ):
                print(f"Orbit {orbit:04d}: {count} matches")
                if footprint_misses:
                    print(
                        f" Orbit {orbit:04d}: {footprint_misses} CS-gated "
                        "samples fall outside every detector footprint"
                    )
                if missing:
                    print(
                        f" Orbit {orbit:04d}: {missing} CS-supported samples "
                        "lack a detector WIC/SI13 pair"
                    )


if __name__ == "__main__":
    main()
