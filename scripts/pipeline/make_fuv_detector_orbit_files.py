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
from icreader import open_product
from netCDF4 import Dataset
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

from icbuilder.fuvdetector import (
    BS2D_IMAGE_FIELDS,
    BS2D_POINTING_METHODS,
    BS2D_PREPROCESSING_LABEL,
    IMAGE_FIELDS,
    PREPROCESSING_LABEL,
    SCHEMA_VERSION,
    SOURCE_TIME_DECODING,
    TIME_TOLERANCE_SECONDS,
    UNSUBTRACTED_IMAGE_FIELDS,
    FUVDetector,
)


#%% Product validation and atomic publication

def validate_label(label):
    """Require an implemented temporal1D or independent2D configuration."""

    if not re.fullmatch(r"[A-Za-z0-9_.-]+", label):
        raise ValueError(
            "preprocessing label may contain only letters, numbers, _, -, and ."
        )
    if label not in (PREPROCESSING_LABEL, BS2D_PREPROCESSING_LABEL):
        raise ValueError(
            f"this implementation supports {PREPROCESSING_LABEL} and {BS2D_PREPROCESSING_LABEL}; "
            "another label requires an explicit preprocessing branch"
        )
    return label

def load_successful_orbits(base):
    """Read the completed native run's per-sensor success arrays."""

    successful = {}
    filenames = {
        "wic": "wic_avail_orbit.npy",
        "si12": "s12_avail_orbit.npy",
        "si13": "s13_avail_orbit.npy",
    }
    for sensor, name in filenames.items():
        filename = base / name
        if not filename.is_file():
            raise FileNotFoundError(
                f"{filename} is missing; finish the native orbit run first"
            )
        status = np.load(filename, allow_pickle=False)
        if status.ndim != 2 or status.shape[1] != 2:
            raise ValueError(f"{filename} must contain orbit/status pairs")
        successful[sensor] = set(status[status[:, 1] == 1, 0].astype(int))
    return successful


def source_paths(
    base, orbit, successful_orbits, wic_folder, si12_folder, si13_folder
):
    """Resolve sources marked successful by the completed native run."""

    paths = {
        "wic": base / wic_folder / f"wic_or{orbit:04d}.nc",
        "si12": base / si12_folder / f"s12_or{orbit:04d}.nc",
        "si13": base / si13_folder / f"s13_or{orbit:04d}.nc",
    }
    if orbit not in successful_orbits["wic"]:
        raise ValueError(f"WIC orbit {orbit} is not marked successful")
    for sensor in ("si12", "si13"):
        if orbit not in successful_orbits[sensor]:
            paths[sensor] = None
    for filename in paths.values():
        if filename is not None and not filename.is_file():
            raise FileNotFoundError(filename)
    return paths


def source_pointing_method(source_files):
    """Read and cross-check the shared pointing method in native 2D files."""

    methods = set()
    for filename in source_files.values():
        if filename is None:
            continue
        with Dataset(filename) as source:
            methods.add(getattr(source, "pointing_correction_method", ""))
    if len(methods) != 1 or not methods.issubset(BS2D_POINTING_METHODS):
        raise ValueError("native sensor files must use one supported pointing method")
    return next(iter(methods))


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
    pointing_method=None,
):
    """Return missing, invalid, mismatch, or complete for one Product-1 file."""

    filename = Path(filename)
    if not filename.is_file():
        return "missing"

    try:
        with open_product(filename) as product:
            if (
                product.product_type != "fuv_detector"
                or product.representation != "detector"
                or int(product.schema_version) != SCHEMA_VERSION
            ):
                return "invalid"
            if (
                product.preprocessing_label != preprocessing_label
                or not np.isclose(
                    product.time_match_tolerance_seconds,
                    time_tolerance_seconds,
                )
            ):
                return "mismatch"
            if product.source_time_decoding != SOURCE_TIME_DECODING:
                return "invalid"
            if int(product.attrs.get("unsubtracted_counts_stored", 0)) != 1:
                return "invalid"

            image_fields = IMAGE_FIELDS
            if preprocessing_label == BS2D_PREPROCESSING_LABEL:
                image_fields = BS2D_IMAGE_FIELDS
                if pointing_method is None:
                    pointing_method = product.attrs.get("pointing_correction_method")
                expected = {
                    "background_fit_stage": "independent",
                    "dgsigma_freeze_iteration": 2,
                    "pointing_correction_method": pointing_method,
                    "pointing_reference_sensor": "WIC"
                }
                if pointing_method not in BS2D_POINTING_METHODS:
                    return "mismatch"
                for name, value in expected.items():
                    if product.attrs.get(name) != value:
                        return "mismatch"

            for sensor in ("wic", "si12", "si13"):
                source = source_files[sensor]
                expected_path = "" if source is None else str(source)
                if (
                    product.source_products[sensor] != expected_path
                    or product.attrs[f"{sensor}_image_field"]
                    != image_fields[sensor.upper()]
                ):
                    return "mismatch"
                if (
                    product.attrs.get(f"{sensor}_unsubtracted_image_field")
                    != UNSUBTRACTED_IMAGE_FIELDS[sensor.upper()]
                ):
                    return "invalid"
                for field in (
                    f"{sensor}_unsubtracted_counts",
                    f"{sensor}_unsubtracted_valid",
                ):
                    if product.read(field, 0).shape != product.shape[1:]:
                        return "invalid"
            expected_units = {
                "counts": "counts",
                "variance": "counts^2",
                "quality_weight": "1",
                "coverage": "1",
            }
            for sensor in ("wic", "si12", "si13"):
                for field, units in expected_units.items():
                    metadata = product.variable_attrs[f"{sensor}_{field}"]
                    if metadata.get("units") != units:
                        return "invalid"
                if (
                    product.variable_attrs[f"{sensor}_frame_quality"].get(
                        "flag_meanings"
                    )
                    != "rejected usable science_ready"
                ):
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
            product.pointing_correction_method,
        )
        if status != "complete":
            raise RuntimeError(f"incomplete fuv_detector file: {partial}")
        os.replace(partial, filename)
    finally:
        if partial.exists():
            partial.unlink()


#%% Orbit processing

def process_orbit(
    orbit,
    base,
    output_directory,
    preprocessing_label,
    software_version,
    successful_orbits,
    wic_folder="wic",
    si12_folder="s12",
    si13_folder="s13",
):
    """Coregister one WIC-led orbit and publish its detector Product 1."""

    paths = source_paths(
        base, orbit, successful_orbits, wic_folder, si12_folder, si13_folder
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
        "--preprocessing-label", default=BS2D_PREPROCESSING_LABEL,
        help="Input background configuration (default: independent2D with shared pointing)."
    )
    parser.add_argument(
        "--orbit", action="append", type=int,
        help="Process only this orbit; repeat the option for more than one.",
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Number of orbit workers; 1 runs serially (default: 1).",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Rebuild existing Product 1 files after native inputs were regenerated."
    )
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

    successful_orbits = load_successful_orbits(base_input)
    available = np.asarray(sorted(successful_orbits["wic"]), dtype=int)
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
            successful_orbits,
            args.wic_folder,
            args.s12_folder,
            args.s13_folder,
        )
        output_file = output_directory / f"or_{int(orbit):04d}.nc"
        if args.overwrite:
            pending.append(int(orbit))
            continue
        pointing_method = None
        if label == BS2D_PREPROCESSING_LABEL and output_file.is_file():
            pointing_method = source_pointing_method(paths)
        status = fuv_detector_file_status(
            output_file,
            label,
            paths,
            pointing_method=pointing_method,
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
        successful_orbits=successful_orbits,
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
