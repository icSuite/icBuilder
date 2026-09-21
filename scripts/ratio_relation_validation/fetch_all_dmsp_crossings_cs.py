"""Collect DMSP samples matched to fixed-grid CS precipitation products."""

#%% Imports

import argparse
from multiprocessing import get_context
from pathlib import Path

from icreader import open_product
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm
import xarray as xr

from fetch_all_dmsp_crossings import (
    DMSP_FIELDS,
    SATELLITES,
    SUPPORT_SECONDS,
    index_dmsp_files,
    load_dmsp,
    save_orbit,
    unit_vectors,
)


#%% Defaults

DEFAULT_DMSP_PATH = Path(
    "/home/bing/Dropbox/work/data/dmsp/dmsp_ssj_yearly"
)
DEFAULT_IMAGE_PATH = Path(
    "/home/bing/Dropbox/work/data/IMAGE_FUV/precipitation_cs/"
    "image_apex_130km_46x46_v1/IR_hardy"
)
DEFAULT_OUTPUT_PATH = Path(
    "/home/bing/Dropbox/work/data/IMAGE_FUV/"
    "ratio_relation_validation/cs_crossings"
)

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
    "method_quality_weight": "quality_weight",
    "method_valid": "method_valid",
}


#%% Fixed-grid matching

def make_grid_matcher(image):
    """Prepare the fixed CS grid lookup used by every frame in an orbit."""

    mlat = np.asarray(image.mlat, dtype=float)
    mlt = np.asarray(image.mlt, dtype=float)
    finite = np.isfinite(mlat) & np.isfinite(mlt)
    flat_cells = np.flatnonzero(finite)
    tree = cKDTree(unit_vectors(
        mlat.ravel()[flat_cells], mlt.ravel()[flat_cells]
    ))
    return mlat, mlt, flat_cells, tree


def nearest_grid_cells(image, matcher, dmsp_mlat, dmsp_mlt):
    """Find the nearest CS cell and identify footprints inside the grid."""

    mlat, mlt, flat_cells, tree = matcher
    chord, nearest = tree.query(unit_vectors(dmsp_mlat, dmsp_mlt))
    flat_nearest = flat_cells[nearest]
    row, column = np.unravel_index(flat_nearest, mlat.shape)
    separation = np.rad2deg(2 * np.arcsin(np.clip(chord / 2, 0, 1)))

    track_xi, track_eta = image.grid.projection.geo2cube(
        dmsp_mlt * 15,
        dmsp_mlat,
    )
    inside = (
        np.isfinite(track_xi)
        & np.isfinite(track_eta)
        & (track_xi >= np.nanmin(image.grid.xi))
        & (track_xi <= np.nanmax(image.grid.xi))
        & (track_eta >= np.nanmin(image.grid.eta))
        & (track_eta <= np.nanmax(image.grid.eta))
    )
    return row, column, separation, inside


def make_samples():
    names = [
        "image_frame", "image_time", "dmsp_sat", "dmsp_time",
        "time_offset_seconds", "dmsp_mlat", "dmsp_mlt",
        "cs_row", "cs_column", "cs_separation_deg", "cs_inside",
    ]
    names += [f"dmsp_{name}" for name in DMSP_FIELDS]
    return {name: [] for name in names}


def append_frame_samples(samples, image, matcher, frame, satellite, dmsp):
    """Append one frame/satellite match to ordinary NumPy lists."""

    dmsp_mlat = np.asarray(dmsp.mlat.values, dtype=float)
    dmsp_mlt = np.asarray(dmsp.mlt.values, dtype=float)
    finite = np.isfinite(dmsp_mlat) & np.isfinite(dmsp_mlt)
    if not np.any(finite):
        return

    dmsp_time = np.asarray(dmsp.time.values, dtype="datetime64[ns]")[finite]
    dmsp_mlat = dmsp_mlat[finite]
    dmsp_mlt = dmsp_mlt[finite]
    row, column, separation, inside = nearest_grid_cells(
        image, matcher, dmsp_mlat, dmsp_mlt
    )
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
    samples["cs_row"].append(row.astype(np.int16))
    samples["cs_column"].append(column.astype(np.int16))
    samples["cs_separation_deg"].append(separation.astype(np.float32))
    samples["cs_inside"].append(inside)

    for name in DMSP_FIELDS:
        samples[f"dmsp_{name}"].append(
            np.asarray(dmsp[name].values, dtype=np.float32)[finite]
        )


#%% One output file per IMAGE orbit

def process_orbit(image_file, dmsp_files, dmsp_cache):
    """Collect all temporal DMSP matches for one CS precipitation orbit."""

    orbit = int(image_file.stem.split("_")[-1])
    samples = make_samples()

    with open_product(image_file) as image:
        if image.product_type != "precipitation_cs" or image.representation != "cs":
            raise ValueError(f"{image_file} is not CS Product 2")

        image_times = np.asarray(image.time, dtype="datetime64[ns]")
        start = image_times[0] - np.timedelta64(SUPPORT_SECONDS, "s")
        stop = image_times[-1] + np.timedelta64(SUPPORT_SECONDS, "s")
        dmsp = {
            satellite: load_dmsp(dmsp_files, dmsp_cache, satellite, start, stop)
            for satellite in SATELLITES
        }
        matcher = make_grid_matcher(image)

        for frame, image_time in enumerate(image_times):
            frame_start = image_time - np.timedelta64(SUPPORT_SECONDS, "s")
            frame_stop = image_time + np.timedelta64(SUPPORT_SECONDS, "s")
            for satellite, data in dmsp.items():
                if data is None:
                    continue
                subset = data.sel(time=slice(frame_start, frame_stop))
                if subset.sizes["time"]:
                    append_frame_samples(
                        samples, image, matcher, frame, satellite, subset
                    )

        if not samples["image_frame"]:
            return xr.Dataset(coords={"sample": np.arange(0)})

        values = {name: np.concatenate(parts) for name, parts in samples.items()}
        frame = values["image_frame"]
        row = values["cs_row"]
        column = values["cs_column"]

        for source_name, output_name in IMAGE_FIELDS.items():
            cube = np.asarray(image.read(source_name))
            values[output_name] = cube[frame, row, column]

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
            "representation": "cs",
            "frame_half_support_seconds": SUPPORT_SECONDS,
            "matching_method": "nearest fixed Cubed-Sphere cell in Modified Apex coordinates",
            "spatial_filter": "none; cs_inside is stored for downstream selection",
            "source_precipitation_cs": str(image_file.resolve()),
            "source_precipitation_detector_sha256": image.attrs[
                "source_precipitation_detector_sha256"
            ],
            "source_fuv_detector_sha256": image.source_fuv_detector_sha256,
            "source_preprocessing_label": image.attrs["source_preprocessing_label"],
            "proton_energy_model": image.proton_energy_model,
            "grid_id": image.grid_id,
            "grid_coordinate_sha256": image.grid_coordinate_sha256,
        })
        return result


def process_and_save(image_file, output_path, dmsp_files, dmsp_cache):
    """Process one orbit and return its orbit number and match count."""

    orbit = int(image_file.stem.split("_")[-1])
    result = process_orbit(image_file, dmsp_files, dmsp_cache)
    save_orbit(result, output_path / f"or_{orbit:04d}.nc")
    count = result.sizes["sample"]
    result.close()
    return orbit, count


_worker_dmsp_files = None
_worker_dmsp_cache = None


def initialize_worker(dmsp_path):
    """Give each process independent DMSP file handles and a local cache."""

    global _worker_dmsp_files, _worker_dmsp_cache
    _worker_dmsp_files = index_dmsp_files(dmsp_path)
    _worker_dmsp_cache = {}


def process_in_worker(task):
    """Multiprocessing entry point; all NetCDF files are opened here."""

    image_file, output_path = task
    return process_and_save(
        image_file, output_path, _worker_dmsp_files, _worker_dmsp_cache
    )


#%% Run

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dmsp-path", type=Path, default=DEFAULT_DMSP_PATH)
    parser.add_argument("--image-path", type=Path, default=DEFAULT_IMAGE_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--orbit", type=int, nargs="+")
    parser.add_argument(
        "--workers", type=int, default=1,
        help="number of independent orbit workers; 1 runs serially"
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
        output_file = args.output_path / f"or_{orbit:04d}.nc"
        if output_file.exists() and not args.overwrite:
            continue
        tasks.append((image_file, args.output_path))

    if args.workers == 1:
        dmsp_cache = {}
        try:
            for image_file, output_path in tqdm(tasks, desc="DMSP/CS crossings"):
                orbit, count = process_and_save(
                    image_file, output_path, dmsp_files, dmsp_cache
                )
                print(f"Orbit {orbit:04d}: {count} matches")
        finally:
            for data in dmsp_cache.values():
                data.close()
    else:
        context = get_context("spawn")
        with context.Pool(
            args.workers,
            initializer=initialize_worker,
            initargs=(args.dmsp_path,)
        ) as pool:
            results = pool.imap_unordered(process_in_worker, tasks)
            for orbit, count in tqdm(
                results, total=len(tasks), desc="DMSP/CS crossings"
            ):
                print(f"Orbit {orbit:04d}: {count} matches")


if __name__ == "__main__":
    main()
