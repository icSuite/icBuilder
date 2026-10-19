"""Build detector Product 1 from raw IDL with shared pointing and 2D backgrounds.

Sigma freezes after temporal fit 2; the temporal model then converges before
each frame is fitted independently. Both results are saved in fuv_native/
fuvpy_bs2d_sigma2_pointing_v1/{wic,s12,s13}. Product 1 uses independent counts
and weights and is written to fuv_detector/fuvpy_bs2d_sigma2_pointing_v1.
"""

#%% Imports and fixed scientific settings

import argparse
import hashlib
import os
from functools import partial
from pathlib import Path

import fuvpy as fuv
import numpy as np
import pandas as pd
from icreader import open_product
from netCDF4 import Dataset
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

from icbuilder.fuvdetector import (
    BS2D_PREPROCESSING_LABEL,
    SCHEMA_VERSION,
    TIME_TOLERANCE_SECONDS,
    FUVDetector,
    load_detector_source
)


PREPROCESSING_LABEL = BS2D_PREPROCESSING_LABEL
SENSOR_PREFIXES = {"WIC": "wic", "SI12": "s12", "SI13": "s13"}
BACKGROUND_MODEL_SETTINGS = {"dgsigma_freeze_iteration": 2}


#%% Product checks and atomic writing

def native_paths(output_base, orbit):
    directory = output_base / "fuv_native" / PREPROCESSING_LABEL
    return {sensor: directory / prefix / f"{prefix}_or{orbit:04d}.nc" for sensor, prefix in SENSOR_PREFIXES.items()}


def product_status(filename, sources, input_digest, software_version=None):
    """Check the reader contract, selected fit and native source identity."""

    if not filename.is_file():
        return "missing"
    try:
        with open_product(filename) as product:
            if product.product_type != "fuv_detector" or product.schema_version != SCHEMA_VERSION:
                return "invalid"
            if product.preprocessing_label != PREPROCESSING_LABEL:
                return "mismatch"
            if product.attrs.get("background_fit_stage") != "independent" or product.attrs.get("dgsigma_freeze_iteration") != 2:
                return "mismatch"
            if product.attrs.get("raw_input_digest") != input_digest:
                return "mismatch"
            if software_version is not None and product.software_version != software_version:
                return "mismatch"
            if not np.isclose(product.time_match_tolerance_seconds, TIME_TOLERANCE_SECONDS):
                return "mismatch"
            if product.attrs.get("pointing_correction_method") != "clock_phase" or product.attrs.get("pointing_reference_sensor") != "WIC":
                return "mismatch"
            for sensor, source in sources.items():
                name = sensor.lower()
                if not source.is_file() or product.source_products[name] != str(source):
                    return "mismatch"
                if product.attrs.get(f"{name}_image_field") != "dgimg_independent":
                    return "mismatch"
                if product.attrs.get(f"source_{name}_size_bytes") != source.stat().st_size or product.attrs.get(f"source_{name}_mtime_ns") != source.stat().st_mtime_ns:
                    return "mismatch"
    except (OSError, ValueError, KeyError, AttributeError, RuntimeError):
        return "invalid"
    return "complete"


def save_native(data, filename, sensor):
    """Retain temporal, independent and pointing fields in a native orbit."""

    filename.parent.mkdir(parents=True, exist_ok=True)
    partial_file = filename.with_name(filename.name + ".partial")
    data = data.assign(t_start=np.datetime_as_string(data.date.values[0], unit="s"))
    data.attrs["preprocessing_label"] = PREPROCESSING_LABEL
    encoding = {name: {"zlib": True, "complevel": 4} for name in data.data_vars}
    try:
        data.to_netcdf(partial_file, format="NETCDF4", encoding=encoding)
        load_detector_source(partial_file, sensor, preprocessing_label=PREPROCESSING_LABEL)
        os.replace(partial_file, filename)
    finally:
        partial_file.unlink(missing_ok=True)


def raw_input_digest(file_lists):
    """Identify the selected raw files without rereading every IDL payload."""

    digest = hashlib.sha256()
    for sensor, files in file_lists.items():
        digest.update(sensor.encode())
        for filename in files:
            info = filename.stat()
            digest.update(f"{filename}\0{info.st_size}\0{info.st_mtime_ns}\n".encode())
    return digest.hexdigest()


def implementation_version():
    """Fingerprint the runner, coregistration, FUVPY calculations and lookup."""

    repository = Path(__file__).resolve().parents[2]
    package = Path(fuv.__file__).resolve().parent
    files = [Path(__file__), *(repository / "icbuilder" / name for name in ("fuvdetector.py", "detector_coregistration.py", "footprints.py"))]
    files += sorted(package.rglob("*.py"))
    files.append(package / "data" / "pointing" / "mlatlon.1997a.xdr")
    digest = hashlib.sha256()
    for filename in files:
        digest.update(filename.name.encode())
        digest.update(filename.read_bytes())
    return "icBuilder+fuvpy.bs2d-" + digest.hexdigest()[:16]


#%% One orbit: common correction, two fits, native archives, coregistration

def process_orbit(task, output_base, output_folder, overwrite=False):
    orbit, file_lists = task
    sources = native_paths(output_base, orbit)
    output_file = output_base / output_folder / PREPROCESSING_LABEL / f"or_{orbit:04d}.nc"
    input_digest = raw_input_digest(file_lists)
    software_version = implementation_version()
    if not overwrite:
        status = product_status(output_file, sources, input_digest, software_version)
        if status == "complete":
            return orbit, "complete"
        if status == "mismatch":
            raise ValueError(f"{output_file} has different inputs/configuration; use --overwrite")

    # All three full streams must be present when estimating the common motion.
    sensors = {}
    for sensor, filenames in file_lists.items():
        sensors[sensor] = fuv.read_idl(
            [str(filename) for filename in filenames], dzalim=75,
            reflat=(sensor == "WIC"), measurement_variance=True,
            poisson_offset=0.5, wic_variance_model="reference", pointing_metadata=True
        ).sortby("pointing_time")
        precise_times = sensors[sensor].pointing_time.values
        if np.any(np.isnat(precise_times)) or np.any(np.diff(precise_times) <= np.timedelta64(0, "ns")):
            raise ValueError(f"orbit {orbit} {sensor} has invalid or duplicate precise times")

    print(f"Orbit {orbit}: shared WIC/SI12/SI13 metadata pointing correction", flush=True)
    sensors = fuv.correct_pointing(sensors, method="clock_phase")

    # Select the northern science interval only after correcting shared metadata.
    for sensor, data in sensors.items():
        data = data.isel(date=data.hemisphere.values == "north")
        if not data.sizes["date"] or not np.any(data.frame_quality.values >= 2):
            raise ValueError(f"orbit {orbit} {sensor} has no northern science-ready frames")
        print(f"Orbit {orbit} {sensor}: temporal + independent 2D, sigma after fit 2", flush=True)
        fitted = fuv.backgroundmodel_BS2D(data, **BACKGROUND_MODEL_SETTINGS)
        fitted.attrs["raw_input_digest"] = input_digest
        fitted.attrs["software_version"] = software_version
        fitted.attrs["raw_source_files"] = "\n".join(str(path) for path in file_lists[sensor])
        save_native(fitted, sources[sensor], sensor)
    del sensors, data, fitted

    # Select the independent fields explicitly; keep both fits in native archives.
    print(f"Orbit {orbit}: coregistering independent backgrounds on WIC pixels", flush=True)
    product = FUVDetector.from_files(
        sources["WIC"], sources["SI12"], sources["SI13"],
        preprocessing_label=PREPROCESSING_LABEL,
        software_version=software_version
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    partial_file = output_file.with_name(output_file.name + ".partial")
    try:
        product.to_nc(partial_file)
        # Record raw provenance alongside the native-file hashes in Product 1.
        with Dataset(partial_file, "r+") as nc:
            nc.raw_input_digest = input_digest
        if product_status(partial_file, sources, input_digest, software_version) != "complete":
            raise ValueError(f"{partial_file} failed Product-1 validation")
        os.replace(partial_file, output_file)
    finally:
        partial_file.unlink(missing_ok=True)
    return orbit, product.shape[0]


#%% Command line and orbit selection

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-input", "--base", type=Path, default=Path(__file__).resolve().parents[2] / "example_data", help="Directory with *files.h5 indices and raw *_data folders")
    parser.add_argument("--base-output", type=Path, help="Output base; defaults to the input base, under a separate 2D label")
    parser.add_argument("--wic-folder", default="wic_data", help="Raw WIC IDL folder relative to the input base")
    parser.add_argument("--s12-folder", default="s12_data", help="Raw SI12 IDL folder relative to the input base")
    parser.add_argument("--s13-folder", default="s13_data", help="Raw SI13 IDL folder relative to the input base")
    parser.add_argument("--output-folder", default="fuv_detector")
    parser.add_argument("--orbit", type=int, action="append", help="Orbit to process; repeat for multiple orbits")
    parser.add_argument("--workers", type=int, default=1, help="Concurrent orbits; each worker holds all three sensor streams")
    parser.add_argument("--overwrite", action="store_true", help="Rebuild completed 2D products and native archives")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.workers < 1:
        raise ValueError("workers must be at least 1")
    base_input = args.base_input.expanduser().resolve()
    base_output = args.base_output.expanduser().resolve() if args.base_output else base_input
    folders = {"WIC": args.wic_folder, "SI12": args.s12_folder, "SI13": args.s13_folder}
    indices = {sensor: pd.read_hdf(base_input / f"{prefix}files.h5", key="data") for sensor, prefix in SENSOR_PREFIXES.items()}
    available = set.intersection(*(set(index.orbit.astype(int)) for index in indices.values()))
    available.discard(-1)
    selected = sorted(set(args.orbit)) if args.orbit else sorted(available)
    missing = set(selected) - available
    if missing:
        raise ValueError(f"orbits missing from one or more sensor indices: {sorted(missing)}")
    if not selected:
        raise ValueError("no orbits are available in all three sensor indices")

    tasks = []
    for orbit in selected:
        files = {}
        for sensor, index in indices.items():
            names = index.loc[index.orbit == orbit].sort_index().filename
            files[sensor] = [(base_input / folders[sensor] / str(name)).resolve() for name in names]
            if len(set(files[sensor])) != len(files[sensor]):
                raise ValueError(f"orbit {orbit} {sensor} index contains duplicate files")
        tasks.append((orbit, files))

    print(f"2D model, sigma frozen after fit 2, shared metadata pointing: {len(tasks)} orbits", flush=True)
    print(f"Native fits: {base_output / 'fuv_native' / PREPROCESSING_LABEL}", flush=True)
    print(f"Product 1: {base_output / args.output_folder / PREPROCESSING_LABEL}", flush=True)
    function = partial(process_orbit, output_base=base_output, output_folder=args.output_folder, overwrite=args.overwrite)
    if args.workers > 1:
        return process_map(function, tasks, max_workers=args.workers, chunksize=1, desc="2D detector orbits")
    return [function(task) for task in tqdm(tasks, desc="2D detector orbits")]


if __name__ == "__main__":
    main()
