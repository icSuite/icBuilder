"""Build paired fixed-grid precipitation and conductance orbit products."""

#%% Imports

import argparse
import os
import subprocess
from functools import partial
from pathlib import Path

import numpy as np
from icreader import open_product
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

from icbuilder.detectorcs import (
    BINNING_METHOD,
    build_detector_cs_products,
)
from icbuilder.grids import DETECTOR_CS_GRID_ID
from icbuilder.precipitationcs import SCHEMA_VERSION as PRECIPITATION_CS_SCHEMA_VERSION
from icbuilder.conductancecs import SCHEMA_VERSION as CONDUCTANCE_CS_SCHEMA_VERSION


#%% Restart validation and atomic publication

def _common_cs_status(product, product_type, schema_version):
    if (
        product.product_type != product_type
        or product.representation != "cs"
        or int(product.schema_version) != schema_version
        or product.grid_id != DETECTOR_CS_GRID_ID
        or product.shape[1:] != (46, 46)
        or product.binning_method != BINNING_METHOD
    ):
        return "invalid"
    return "complete"


def precipitation_cs_file_status(filename, source_precipitation):
    """Return missing, invalid, mismatch, or complete for Product-2 CS."""

    filename = Path(filename)
    source_precipitation = Path(source_precipitation)
    if not filename.is_file():
        return "missing"
    if not source_precipitation.is_file():
        return "mismatch"

    try:
        source_stat = source_precipitation.stat()
        with open_product(filename) as product:
            status = _common_cs_status(
                product,
                "precipitation_cs",
                PRECIPITATION_CS_SCHEMA_VERSION,
            )
            if status != "complete":
                return status
            if (
                product.source_precipitation_detector
                != str(source_precipitation)
                or int(
                    product.attrs["source_precipitation_detector_size_bytes"]
                )
                != source_stat.st_size
                or int(
                    product.attrs["source_precipitation_detector_mtime_ns"]
                )
                != source_stat.st_mtime_ns
            ):
                return "mismatch"
    except (OSError, RuntimeError, KeyError, AttributeError, ValueError):
        return "invalid"
    return "complete"


def conductance_cs_file_status(
    filename,
    source_conductance,
    companion_precipitation_cs,
):
    """Return missing, invalid, mismatch, or complete for Product-3 CS."""

    filename = Path(filename)
    source_conductance = Path(source_conductance)
    companion_precipitation_cs = Path(companion_precipitation_cs)
    if not filename.is_file():
        return "missing"
    if not source_conductance.is_file():
        return "mismatch"

    try:
        source_stat = source_conductance.stat()
        with open_product(filename) as product:
            status = _common_cs_status(
                product,
                "conductance_cs",
                CONDUCTANCE_CS_SCHEMA_VERSION,
            )
            if status != "complete":
                return status
            if (
                product.source_conductance_detector != str(source_conductance)
                or int(product.attrs["source_conductance_detector_size_bytes"])
                != source_stat.st_size
                or int(product.attrs["source_conductance_detector_mtime_ns"])
                != source_stat.st_mtime_ns
                or product.companion_precipitation_cs
                != str(companion_precipitation_cs)
            ):
                return "mismatch"
    except (OSError, RuntimeError, KeyError, AttributeError, ValueError):
        return "invalid"
    return "complete"


def _atomic_save(product, filename, status_function, *status_arguments):
    filename = Path(filename)
    partial_file = Path(str(filename) + ".partial")
    try:
        product.to_nc(partial_file)
        status = status_function(partial_file, *status_arguments)
        if status != "complete":
            raise RuntimeError(f"incomplete CS product: {partial_file}")
        os.replace(partial_file, filename)
    finally:
        if partial_file.exists():
            partial_file.unlink()


#%% Orbit processing

def get_orbits(directory):
    """Discover four-digit orbit files in one detector-product directory."""

    orbits = []
    for filename in sorted(directory.glob("or_*.nc")):
        orbit_text = filename.stem[-4:]
        if orbit_text.isdigit():
            orbits.append(int(orbit_text))
    return np.unique(orbits)


def current_revision(repository):
    """Return the Git revision recorded in the derived products."""

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
    task,
    precipitation_input_directory,
    conductance_input_directory,
    precipitation_output_directory,
    conductance_output_directory,
    software_version,
):
    """Reduce and atomically publish one paired detector orbit."""

    orbit, write_precipitation, write_conductance = task
    source_precipitation = (
        precipitation_input_directory / f"or_{orbit:04d}.nc"
    )
    source_conductance = conductance_input_directory / f"or_{orbit:04d}.nc"
    precipitation_output = (
        precipitation_output_directory / f"or_{orbit:04d}.nc"
    )
    conductance_output = conductance_output_directory / f"or_{orbit:04d}.nc"

    try:
        precipitation, conductance = build_detector_cs_products(
            source_precipitation,
            source_conductance,
            software_version=software_version,
        )
        conductance.companion_precipitation_cs = str(precipitation_output)

        if write_precipitation:
            _atomic_save(
                precipitation,
                precipitation_output,
                precipitation_cs_file_status,
                source_precipitation,
            )
        if write_conductance:
            _atomic_save(
                conductance,
                conductance_output,
                conductance_cs_file_status,
                source_conductance,
                precipitation_output,
            )
    except Exception as error:
        raise RuntimeError(f"detector CS orbit {orbit:04d} failed") from error

    return orbit, precipitation.shape[0]


#%% Command line

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Reduce paired detector precipitation and conductance onto the "
            "frozen 46-by-46 Cubed-Sphere grid."
        )
    )
    parser.add_argument(
        "--base-input", "--base",
        dest="base_input",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "example_data",
    )
    parser.add_argument(
        "--base-output",
        type=Path,
        help="Output base (default: --base-input).",
    )
    parser.add_argument("--retrieval-label", default="IR_hardy")
    parser.add_argument("--conductance-model", default="robinson")
    parser.add_argument(
        "--precipitation-input-folder", default="precipitation_detector"
    )
    parser.add_argument(
        "--conductance-input-folder", default="conductance_detector"
    )
    parser.add_argument(
        "--precipitation-output-folder", default="precipitation_cs"
    )
    parser.add_argument(
        "--conductance-output-folder", default="conductance_cs"
    )
    parser.add_argument("--orbit", action="append", type=int)
    parser.add_argument("--workers", type=int, default=1)
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
    precipitation_input_directory = (
        base_input / args.precipitation_input_folder / args.retrieval_label
    )
    conductance_input_directory = (
        base_input / args.conductance_input_folder / args.retrieval_label
        / args.conductance_model
    )
    precipitation_output_directory = (
        base_output / args.precipitation_output_folder / DETECTOR_CS_GRID_ID
        / args.retrieval_label
    )
    conductance_output_directory = (
        base_output / args.conductance_output_folder / DETECTOR_CS_GRID_ID
        / args.retrieval_label / args.conductance_model
    )
    precipitation_output_directory.mkdir(parents=True, exist_ok=True)
    conductance_output_directory.mkdir(parents=True, exist_ok=True)

    precipitation_orbits = get_orbits(precipitation_input_directory)
    conductance_orbits = get_orbits(conductance_input_directory)
    available = np.intersect1d(precipitation_orbits, conductance_orbits)
    selected = available if args.orbit is None else np.unique(args.orbit)
    missing = selected[~np.isin(selected, available)]
    if missing.size:
        raise ValueError(
            "paired detector Product 2/Product 3 orbit is missing: "
            f"{missing.tolist()}"
        )

    tasks = []
    for orbit_value in selected:
        orbit = int(orbit_value)
        source_precipitation = (
            precipitation_input_directory / f"or_{orbit:04d}.nc"
        )
        source_conductance = (
            conductance_input_directory / f"or_{orbit:04d}.nc"
        )
        precipitation_output = (
            precipitation_output_directory / f"or_{orbit:04d}.nc"
        )
        conductance_output = conductance_output_directory / f"or_{orbit:04d}.nc"

        if args.overwrite:
            precipitation_status = conductance_status = "invalid"
        else:
            precipitation_status = precipitation_cs_file_status(
                precipitation_output, source_precipitation
            )
            conductance_status = conductance_cs_file_status(
                conductance_output, source_conductance, precipitation_output
            )
        if precipitation_status == "mismatch" or conductance_status == "mismatch":
            raise ValueError(
                f"orbit {orbit:04d} CS output does not match its source or "
                "configuration; use another output folder or --overwrite"
            )
        write_precipitation = precipitation_status != "complete"
        write_conductance = conductance_status != "complete"
        if write_precipitation or write_conductance:
            tasks.append((orbit, write_precipitation, write_conductance))

    print(
        f"detector CS/{DETECTOR_CS_GRID_ID}/{args.retrieval_label}/"
        f"{args.conductance_model}: "
        f"{len(selected) - len(tasks)} complete, {len(tasks)} pending"
    )
    if not tasks:
        return []

    repository = Path(__file__).resolve().parents[2]
    function = partial(
        process_orbit,
        precipitation_input_directory=precipitation_input_directory,
        conductance_input_directory=conductance_input_directory,
        precipitation_output_directory=precipitation_output_directory,
        conductance_output_directory=conductance_output_directory,
        software_version=current_revision(repository),
    )
    if args.workers > 1:
        return process_map(
            function,
            tasks,
            max_workers=args.workers,
            chunksize=1,
            desc="Create detector CS orbit products",
        )
    return [
        function(task)
        for task in tqdm(tasks, desc="Create detector CS orbit products")
    ]


if __name__ == "__main__":
    main()
