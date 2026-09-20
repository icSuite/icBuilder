"""Build detector-space conductance from precipitation Product 2."""

#%% Imports

import argparse
import os
import subprocess
from functools import partial
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

from icbuilder.conductancedetector import (
    CONDUCTANCE_MODEL,
    GEOMETRY_FIELDS,
    PRECIPITATION_FIELDS,
    SCHEMA_VERSION,
    TIME_FIELDS,
    TIME_VALUE_FIELDS,
    ConductanceDetector,
)


#%% Product validation and atomic publication

RESULT_FIELDS = (
    "P", "H", "dP", "dH", "method_quality_weight",
    "method_valid", "conductance_valid", "conductance_uncertainty_valid",
    "Ep_clipping_flag",
)


def conductance_detector_file_status(
    filename,
    source_precipitation,
    conductance_model=CONDUCTANCE_MODEL,
):
    """Return missing, invalid, mismatch, or complete for Product 3."""

    filename = Path(filename)
    source_precipitation = Path(source_precipitation)
    if not filename.is_file():
        return "missing"
    if not source_precipitation.is_file():
        return "mismatch"

    try:
        source_stat = source_precipitation.stat()
        with Dataset(filename) as nc, Dataset(source_precipitation) as source:
            if (
                nc.product_type != "conductance_detector"
                or nc.representation != "detector"
                or int(nc.schema_version) != SCHEMA_VERSION
            ):
                return "invalid"
            if (
                nc.conductance_model != conductance_model
                or nc.source_precipitation_detector
                != str(source_precipitation)
                or int(nc.source_precipitation_detector_size_bytes)
                != source_stat.st_size
                or int(nc.source_precipitation_detector_mtime_ns)
                != source_stat.st_mtime_ns
            ):
                return "mismatch"
            if (
                source.product_type != "precipitation_detector"
                or source.representation != "detector"
                or nc.precipitation_method != source.method
                or nc.proton_flux_source != source.proton_flux_source
                or nc.proton_energy_model != source.proton_energy_model
                or nc.count_uncertainty_mode
                != source.count_uncertainty_mode
                or int(nc.source_precipitation_detector_schema_version)
                != int(source.schema_version)
            ):
                return "mismatch"
            if source.proton_energy_model == "constant" and (
                not np.isclose(
                    nc.proton_energy_constant,
                    source.proton_energy_constant,
                )
                or not np.isclose(
                    nc.proton_energy_uncertainty_constant,
                    source.proton_energy_uncertainty_constant,
                )
            ):
                return "mismatch"

            shape = (
                len(nc.dimensions["time"]),
                len(nc.dimensions["row"]),
                len(nc.dimensions["column"]),
            )
            if any(length == 0 for length in shape):
                return "invalid"
            for name in GEOMETRY_FIELDS + PRECIPITATION_FIELDS + RESULT_FIELDS:
                if nc.variables[name].shape != shape:
                    return "invalid"
            for name in TIME_FIELDS + TIME_VALUE_FIELDS:
                if nc.variables[name].shape != (shape[0],):
                    return "invalid"
            if nc.variables["detector_row"].shape != (shape[1],):
                return "invalid"
            if nc.variables["detector_column"].shape != (shape[2],):
                return "invalid"
    except (OSError, RuntimeError, KeyError, AttributeError, ValueError):
        return "invalid"

    return "complete"


def save_conductance_detector(product, filename):
    """Write beside the final path and atomically publish after validation."""

    filename = Path(filename)
    partial = Path(str(filename) + ".partial")
    try:
        product.to_nc(partial)
        status = conductance_detector_file_status(
            partial,
            product.source_precipitation_detector,
            product.conductance_model,
        )
        if status != "complete":
            raise RuntimeError(f"incomplete conductance_detector file: {partial}")
        os.replace(partial, filename)
    finally:
        if partial.exists():
            partial.unlink()


#%% Orbit processing

def get_orbits(input_directory):
    """Discover orbit numbers from detector Product-2 files."""

    orbits = []
    for filename in sorted(input_directory.glob("or_*.nc")):
        orbit_text = filename.stem[-4:]
        if orbit_text.isdigit():
            orbits.append(int(orbit_text))
    return np.unique(orbits)


def current_revision(repository):
    """Return the Git revision recorded in experimental Product 3."""

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
        return "unknown"
    return revision + ("+worktree" if status else "")


def process_orbit(
    orbit,
    input_directory,
    output_directory,
    software_version,
):
    """Calculate and save one detector-space conductance orbit."""

    source = input_directory / f"or_{orbit:04d}.nc"
    try:
        product = ConductanceDetector(
            source,
            software_version=software_version,
        )
        output = output_directory / f"or_{orbit:04d}.nc"
        save_conductance_detector(product, output)
    except Exception as error:
        raise RuntimeError(
            f"conductance_detector orbit {orbit:04d} failed"
        ) from error

    return orbit, product.shape[0]


#%% Command line

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Create Robinson conductance on WIC detector pixels."
    )
    parser.add_argument(
        "--base-input", "--base",
        dest="base_input",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "example_data",
        help="Base directory containing detector Product-2 orbit files.",
    )
    parser.add_argument(
        "--base-output",
        type=Path,
        help="Base directory for Product 3 (default: --base-input).",
    )
    parser.add_argument(
        "--input-folder", default="precipitation_detector"
    )
    parser.add_argument(
        "--output-folder", default="conductance_detector"
    )
    parser.add_argument(
        "--retrieval-label", default="IR_hardy",
        help="Product-2 retrieval subfolder (default: IR_hardy).",
    )
    parser.add_argument(
        "--conductance-model",
        choices=(CONDUCTANCE_MODEL,),
        default=CONDUCTANCE_MODEL,
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

    base_input = args.base_input.expanduser()
    base_output = (
        args.base_output.expanduser()
        if args.base_output is not None else base_input
    )
    input_directory = (
        base_input / args.input_folder / args.retrieval_label
    )
    output_directory = (
        base_output / args.output_folder / args.retrieval_label
        / args.conductance_model
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    available = get_orbits(input_directory)
    selected = available if args.orbit is None else np.unique(args.orbit)
    missing = selected[~np.isin(selected, available)]
    if missing.size:
        raise ValueError(
            f"precipitation_detector orbit is missing: {missing.tolist()}"
        )

    pending = []
    for orbit in selected:
        source = input_directory / f"or_{int(orbit):04d}.nc"
        output = output_directory / f"or_{int(orbit):04d}.nc"
        if args.overwrite:
            pending.append(int(orbit))
            continue
        status = conductance_detector_file_status(
            output, source, args.conductance_model
        )
        if status == "mismatch":
            raise ValueError(
                f"{output} does not match the requested detector "
                "conductance configuration or Product-2 source; use another "
                "retrieval/model folder or --overwrite"
            )
        if status != "complete":
            pending.append(int(orbit))

    print(
        f"conductance_detector/{args.retrieval_label}/"
        f"{args.conductance_model}: "
        f"{len(selected) - len(pending)} complete, {len(pending)} pending"
    )
    if not pending:
        return []

    repository = Path(__file__).resolve().parents[2]
    software_version = current_revision(repository)
    function = partial(
        process_orbit,
        input_directory=input_directory,
        output_directory=output_directory,
        software_version=software_version,
    )
    if args.workers > 1:
        return process_map(
            function,
            pending,
            max_workers=args.workers,
            chunksize=1,
            desc="Create detector conductance orbits",
        )

    return [
        function(orbit)
        for orbit in tqdm(pending, desc="Create detector conductance orbits")
    ]


if __name__ == "__main__":
    main()
