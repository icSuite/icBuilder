"""Collect DMSP samples matched to detector precipitation products."""

#%% Imports

import argparse
from multiprocessing import get_context
from pathlib import Path

from icreader import open_product
import numpy as np
from scipy.spatial import cKDTree
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


#%% Detector matching

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


def nearest_pixels(image_mlat, image_mlt, dmsp_mlat, dmsp_mlt):
    """Match DMSP positions to the nearest finite detector pixels."""

    finite = np.isfinite(image_mlat) & np.isfinite(image_mlt)
    flat_pixels = np.flatnonzero(finite)
    tree = cKDTree(unit_vectors(
        image_mlat.ravel()[flat_pixels],
        image_mlt.ravel()[flat_pixels],
    ))
    chord, nearest = tree.query(unit_vectors(dmsp_mlat, dmsp_mlt))
    flat_nearest = flat_pixels[nearest]
    row, column = np.unravel_index(flat_nearest, image_mlat.shape)
    separation = np.rad2deg(2 * np.arcsin(np.clip(chord / 2, 0, 1)))
    return row, column, separation


def append_frame_samples(samples, image, frame, satellite, dmsp, mlat, mlt):
    """Append one frame/satellite match to ordinary NumPy lists."""

    dmsp_mlat = np.asarray(dmsp.mlat.values, dtype=float)
    dmsp_mlt = np.asarray(dmsp.mlt.values, dtype=float)
    finite = np.isfinite(dmsp_mlat) & np.isfinite(dmsp_mlt)
    if not np.any(finite):
        return

    dmsp_time = np.asarray(dmsp.time.values, dtype="datetime64[ns]")[finite]
    dmsp_mlat = dmsp_mlat[finite]
    dmsp_mlt = dmsp_mlt[finite]
    row, column, separation = nearest_pixels(
        mlat[frame], mlt[frame], dmsp_mlat, dmsp_mlt
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
    samples["detector_row"].append(row.astype(np.int16))
    samples["detector_column"].append(column.astype(np.int16))
    samples["detector_separation_deg"].append(separation.astype(np.float32))

    for name in DMSP_FIELDS:
        samples[f"dmsp_{name}"].append(
            np.asarray(dmsp[name].values, dtype=np.float32)[finite]
        )


def make_samples():
    names = [
        "image_frame", "image_time", "dmsp_sat", "dmsp_time",
        "time_offset_seconds", "dmsp_mlat", "dmsp_mlt",
        "detector_row", "detector_column", "detector_separation_deg",
    ]
    names += [f"dmsp_{name}" for name in DMSP_FIELDS]
    return {name: [] for name in names}


#%% One output file per IMAGE orbit

def process_orbit(image_file, dmsp_files, dmsp_cache):
    """Collect all temporal DMSP matches for one detector orbit."""

    orbit = int(image_file.stem.split("_")[-1])
    samples = make_samples()

    with open_product(image_file) as image:
        if image.product_type != "precipitation_detector":
            raise ValueError(f"{image_file} is not detector Product 2")

        image_times = np.asarray(image.time, dtype="datetime64[ns]")
        start = image_times[0] - np.timedelta64(SUPPORT_SECONDS, "s")
        stop = image_times[-1] + np.timedelta64(SUPPORT_SECONDS, "s")
        dmsp = {
            satellite: load_dmsp(
                dmsp_files, dmsp_cache, satellite, start, stop
            )
            for satellite in SATELLITES
        }

        mlat = np.asarray(image.read("mlat"), dtype=float)
        mlt = np.asarray(image.read("mlt"), dtype=float)

        for frame, image_time in enumerate(image_times):
            frame_start = image_time - np.timedelta64(SUPPORT_SECONDS, "s")
            frame_stop = image_time + np.timedelta64(SUPPORT_SECONDS, "s")
            for satellite, data in dmsp.items():
                if data is None:
                    continue
                subset = data.sel(time=slice(frame_start, frame_stop))
                if subset.sizes["time"]:
                    append_frame_samples(
                        samples, image, frame, satellite, subset, mlat, mlt
                    )

        if not samples["image_frame"]:
            return xr.Dataset(coords={"sample": np.arange(0)})

        values = {name: np.concatenate(parts) for name, parts in samples.items()}
        frame = values["image_frame"]
        row = values["detector_row"]
        column = values["detector_column"]

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
            "representation": "detector",
            "frame_half_support_seconds": SUPPORT_SECONDS,
            "matching_method": "nearest detector pixel in Modified Apex latitude/MLT",
            "spatial_filter": "none",
            "source_precipitation_detector": str(image_file.resolve()),
            "source_fuv_detector_sha256": image.source_fuv_detector_sha256,
            "source_preprocessing_label": image.source_preprocessing_label,
            "proton_energy_model": image.proton_energy_model,
            "source_count_source": image.attrs.get(
                "count_source", "background_subtracted"
            ),
        })
        return result


def save_orbit(data, filename):
    filename.parent.mkdir(parents=True, exist_ok=True)
    encoding = {
        name: {"zlib": True, "complevel": 4}
        for name in data.data_vars
        if data[name].dtype.kind not in ("U", "O", "M")
    }
    data.to_netcdf(filename, engine="netcdf4", encoding=encoding)


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
    """Multiprocessing entry point; all netCDF files are opened here."""

    image_file, output_path = task
    return process_and_save(
        image_file,
        output_path,
        _worker_dmsp_files,
        _worker_dmsp_cache,
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
        output_file = args.output_path / f"or_{orbit:04d}.nc"
        if output_file.exists() and not args.overwrite:
            continue
        tasks.append((image_file, args.output_path))

    if args.workers == 1:
        dmsp_cache = {}
        try:
            for image_file, output_path in tqdm(tasks, desc="DMSP crossings"):
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
            initargs=(args.dmsp_path,),
        ) as pool:
            results = pool.imap_unordered(process_in_worker, tasks)
            for orbit, count in tqdm(
                results, total=len(tasks), desc="DMSP crossings"
            ):
                print(f"Orbit {orbit:04d}: {count} matches")


if __name__ == "__main__":
    main()
