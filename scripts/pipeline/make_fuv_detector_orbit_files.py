"""Build experimental detector-space Product-1 orbit files."""

#%% Imports

import argparse
import hashlib
import os
import re
import subprocess
from functools import partial
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

from icbuilder.fuvdetector import (
    COREGISTRATION_FLOAT_FIELDS,
    COREGISTRATION_INTEGER_FIELDS,
    IMAGE_FIELDS,
    PREPROCESSING_LABEL,
    SCHEMA_VERSION,
    SOURCE_TIME_DECODING,
    TIME_TOLERANCE_SECONDS,
    FUVDetector,
)


#%% Product validation and atomic publication

REQUIRED_FRAME_FIELDS = (
    "wic_counts", "si12_counts", "si13_counts",
    "wic_variance", "si12_variance", "si13_variance",
    "wic_quality_weight", "si12_quality_weight", "si13_quality_weight",
    "wic_coverage", "si12_coverage", "si13_coverage",
    "wic_valid", "si12_valid", "si13_valid",
    "si12_source_count", "si13_source_count",
    "glat", "glon", "mlat", "mlon", "mlt", "sza", "dza",
)

BASE_REQUIRED_TIME_FIELDS = (
    "time", "wic_source_time", "si12_source_time", "si13_source_time",
    "wic_source_index", "si12_source_index", "si13_source_index",
    "wic_frame_quality", "si12_frame_quality", "si13_frame_quality",
    "ssalon",
)

REQUIRED_TIME_FIELDS = BASE_REQUIRED_TIME_FIELDS + tuple(
    f"{sensor}_coreg_{name}"
    for sensor in ("si12", "si13")
    for name in (*COREGISTRATION_INTEGER_FIELDS, *COREGISTRATION_FLOAT_FIELDS)
)

def validate_label(label):
    """Require the one preprocessing configuration implemented in this slice."""

    if not re.fullmatch(r"[A-Za-z0-9_.-]+", label):
        raise ValueError(
            "preprocessing label may contain only letters, numbers, _, -, and ."
        )
    if label != PREPROCESSING_LABEL:
        raise ValueError(
            f"this implementation supports only {PREPROCESSING_LABEL}; "
            "another label requires an explicit preprocessing branch"
        )
    return label

def source_paths(base, orbit, wic_folder, si12_folder, si13_folder):
    """Resolve one WIC-led orbit and optional SI source files."""

    paths = {
        "wic": base / wic_folder / f"wic_or{orbit:04d}.nc",
        "si12": base / si12_folder / f"s12_or{orbit:04d}.nc",
        "si13": base / si13_folder / f"s13_or{orbit:04d}.nc",
    }
    if not paths["wic"].is_file():
        raise FileNotFoundError(paths["wic"])
    for sensor in ("si12", "si13"):
        if not paths[sensor].is_file():
            paths[sensor] = None
    return paths

def current_software_version(repository):
    """Describe Git plus the exact Product-1 implementation source."""

    try:
        revision = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=repository,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=repository,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        revision = "unknown"
        status = "unknown"

    digest = hashlib.sha256()
    source_files = (
        "icbuilder/fuvdetector.py",
        "icbuilder/detector_coregistration.py",
        "icbuilder/footprints.py",
        "scripts/pipeline/make_fuv_detector_orbit_files.py",
    )
    for relative_path in source_files:
        path = repository / relative_path
        digest.update(relative_path.encode())
        digest.update(path.read_bytes())
    worktree = "+worktree" if status else ""
    return f"{revision}{worktree}.product1-{digest.hexdigest()[:12]}"


def fuv_detector_file_status(
    filename,
    preprocessing_label,
    source_files,
    time_tolerance_seconds=TIME_TOLERANCE_SECONDS,
):
    """Return missing, invalid, mismatch, or complete for one Product-1 file."""

    filename = Path(filename)
    if not filename.is_file():
        return "missing"

    try:
        with Dataset(filename) as nc:
            if (
                nc.product_type != "fuv_detector"
                or nc.representation != "detector"
                or int(nc.schema_version) != SCHEMA_VERSION
            ):
                return "invalid"
            if (
                nc.preprocessing_label != preprocessing_label
                or not np.isclose(
                    nc.time_match_tolerance_seconds, time_tolerance_seconds
                )
            ):
                return "mismatch"
            if getattr(nc, "source_time_decoding", None) != SOURCE_TIME_DECODING:
                return "invalid"

            for sensor in ("wic", "si12", "si13"):
                source = source_files[sensor]
                expected_path = "" if source is None else str(source)
                if (
                    nc.getncattr(f"source_{sensor}") != expected_path
                    or nc.getncattr(f"{sensor}_image_field")
                    != IMAGE_FIELDS[sensor.upper()]
                ):
                    return "mismatch"
            shape = (
                len(nc.dimensions["time"]),
                len(nc.dimensions["row"]),
                len(nc.dimensions["column"]),
            )
            if any(length == 0 for length in shape):
                return "invalid"
            for name in REQUIRED_FRAME_FIELDS:
                if nc.variables[name].shape != shape:
                    return "invalid"
            expected_units = {
                "counts": "counts",
                "variance": "counts^2",
                "quality_weight": "1",
                "coverage": "1",
            }
            for sensor in ("wic", "si12", "si13"):
                for field, units in expected_units.items():
                    variable = nc.variables[f"{sensor}_{field}"]
                    if getattr(variable, "units", None) != units:
                        return "invalid"
                frame_quality = nc.variables[f"{sensor}_frame_quality"]
                if (
                    getattr(frame_quality, "flag_meanings", None)
                    != "rejected usable science_ready"
                ):
                    return "invalid"
            for name in REQUIRED_TIME_FIELDS:
                if nc.variables[name].shape != (shape[0],):
                    return "invalid"
            if nc.variables["detector_row"].shape != (shape[1],):
                return "invalid"
            if nc.variables["detector_column"].shape != (shape[2],):
                return "invalid"
    except (OSError, RuntimeError, KeyError, AttributeError, ValueError):
        return "invalid"

    return "complete"


def save_fuv_detector_file(product, filename, source_files):
    """Write beside the final path and atomically publish after validation."""

    filename = Path(filename)
    partial = Path(str(filename) + ".partial")
    try:
        product.to_nc(partial)
        status = fuv_detector_file_status(
            partial,
            product.preprocessing_label,
            source_files,
            product.time_tolerance_seconds,
        )
        if status != "complete":
            raise RuntimeError(f"incomplete fuv_detector file: {partial}")
        os.replace(partial, filename)
    finally:
        if partial.exists():
            partial.unlink()


#%% Orbit processing

def get_orbits(wic_directory):
    """Discover orbit numbers from WIC source NetCDF files."""

    orbits = []
    for filename in sorted(wic_directory.glob("wic_or*.nc")):
        orbit_text = filename.stem[-4:]
        if orbit_text.isdigit():
            orbits.append(int(orbit_text))
    return np.unique(orbits)


def process_orbit(
    orbit,
    base,
    output_directory,
    preprocessing_label,
    software_version,
    wic_folder="wic",
    si12_folder="s12",
    si13_folder="s13",
):
    """Coregister one WIC-led orbit and publish its detector Product 1."""

    paths = source_paths(
        base, orbit, wic_folder, si12_folder, si13_folder
    )
    try:
        product = FUVDetector.from_files(
            paths["wic"],
            paths["si12"],
            paths["si13"],
            preprocessing_label=preprocessing_label,
            software_version=software_version,
        )
        output_file = output_directory / f"or_{orbit:04d}.nc"
        save_fuv_detector_file(product, output_file, paths)
    except Exception as error:
        raise RuntimeError(f"fuv_detector orbit {orbit:04d} failed") from error

    return orbit, product.shape[0]


#%% Command line

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Create WIC-detector coregistered IMAGE-FUV orbit files."
    )
    parser.add_argument(
        "--base-input", "--base",
        dest="base_input",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "example_data",
        help="Base directory containing corrected WIC, SI12, and SI13 orbits.",
    )
    parser.add_argument(
        "--base-output",
        type=Path,
        help="Base directory for Product 1 (default: --base-input).",
    )
    parser.add_argument("--wic-folder", default="wic")
    parser.add_argument("--s12-folder", default="s12")
    parser.add_argument("--s13-folder", default="s13")
    parser.add_argument("--output-folder", default="fuv_detector")
    parser.add_argument(
        "--preprocessing-label", default=PREPROCESSING_LABEL
    )
    parser.add_argument(
        "--orbit", action="append", type=int,
        help="Process only this orbit; repeat the option for more than one.",
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Number of orbit workers; 1 runs serially (default: 1).",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.workers < 1:
        raise ValueError("workers must be at least 1")

    label = validate_label(args.preprocessing_label)
    base_input = args.base_input.expanduser()
    base_output = (
        args.base_output.expanduser()
        if args.base_output is not None else base_input
    )
    output_directory = base_output / args.output_folder / label
    output_directory.mkdir(parents=True, exist_ok=True)

    available = get_orbits(base_input / args.wic_folder)
    if args.orbit is None:
        selected = available
    else:
        selected = np.unique(args.orbit)
        missing = selected[~np.isin(selected, available)]
        if missing.size:
            raise ValueError(f"WIC source orbit is missing: {missing.tolist()}")

    repository = Path(__file__).resolve().parents[2]
    software_version = current_software_version(repository)
    pending = []
    for orbit in selected:
        paths = source_paths(
            base_input,
            int(orbit),
            args.wic_folder,
            args.s12_folder,
            args.s13_folder,
        )
        output_file = output_directory / f"or_{int(orbit):04d}.nc"
        if args.overwrite:
            pending.append(int(orbit))
            continue
        status = fuv_detector_file_status(
            output_file,
            label,
            paths,
        )
        if status == "mismatch":
            raise ValueError(
                f"{output_file} does not match the requested Product-1 "
                "configuration; choose another label or use --overwrite"
            )
        if status != "complete":
            pending.append(int(orbit))

    print(f"{args.output_folder}/{label}: {len(selected) - len(pending)} complete, "
          f"{len(pending)} pending")
    if not pending:
        return []

    function = partial(
        process_orbit,
        base=base_input,
        output_directory=output_directory,
        preprocessing_label=label,
        software_version=software_version,
        wic_folder=args.wic_folder,
        si12_folder=args.s12_folder,
        si13_folder=args.s13_folder,
    )
    if args.workers > 1:
        return process_map(
            function,
            pending,
            max_workers=args.workers,
            chunksize=1,
            desc="Create detector Product 1 orbits",
        )

    return [
        function(orbit)
        for orbit in tqdm(pending, desc="Create detector Product 1 orbits")
    ]


if __name__ == "__main__":
    main()
