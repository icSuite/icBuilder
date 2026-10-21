"""Build fixed-grid precipitation and optional conductance orbit products."""

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
    build_conductance_cs_product,
    build_detector_cs_products,
    build_precipitation_cs_product,
)
from icbuilder.grids import DETECTOR_CS_GRID_ID
from icbuilder.precipitationdetector import ELECTRON_METADATA_ATTRIBUTES
from icbuilder.precipitationcs import SCHEMA_VERSION as PRECIPITATION_CS_SCHEMA_VERSION
from icbuilder.conductancecs import SCHEMA_VERSION as CONDUCTANCE_CS_SCHEMA_VERSION


#%% Pipeline defaults

DEFAULT_RETRIEVAL_LABEL = "IR_hardy_smoothed_wic0p8_si120_si131p2"


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
    method = product.attrs.get("method", product.attrs.get("precipitation_method"))
    if method == "hardy":
        for name in ELECTRON_METADATA_ATTRIBUTES:
            if name not in product.attrs:
                return "invalid"
        for name in ("E0_response_range_fraction", "E0_response_range_any"):
            if not hasattr(product, name) or getattr(product, name).shape != product.shape:
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
        with open_product(filename) as product, open_product(source_precipitation) as source:
            status = _common_cs_status(
                product,
                "precipitation_cs",
                PRECIPITATION_CS_SCHEMA_VERSION,
            )
            if status != "complete":
                return status
            if source.method == "hardy":
                if product.method != source.method:
                    return "mismatch"
                for name in ELECTRON_METADATA_ATTRIBUTES:
                    if product.attrs[name] != source.attrs[name]:
                        return "mismatch"
            if not hasattr(product, "si12") or not hasattr(product, "dsi12"):
                return "invalid"
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
):
    """Return missing, invalid, mismatch, or complete for Product-3 CS."""

    filename = Path(filename)
    source_conductance = Path(source_conductance)
    if not filename.is_file():
        return "missing"
    if not source_conductance.is_file():
        return "mismatch"

    try:
        source_stat = source_conductance.stat()
        with open_product(filename) as product, open_product(source_conductance) as source:
            status = _common_cs_status(
                product,
                "conductance_cs",
                CONDUCTANCE_CS_SCHEMA_VERSION,
            )
            if status != "complete":
                return status
            if source.precipitation_method == "hardy":
                if product.precipitation_method != source.precipitation_method:
                    return "mismatch"
                for name in ELECTRON_METADATA_ATTRIBUTES:
                    if product.attrs[name] != source.attrs[name]:
                        return "mismatch"
            if (
                product.source_conductance_detector != str(source_conductance)
                or int(product.attrs["source_conductance_detector_size_bytes"])
                != source_stat.st_size
                or int(product.attrs["source_conductance_detector_mtime_ns"])
                != source_stat.st_mtime_ns
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
    """Reduce and atomically publish one detector orbit."""

    orbit, write_precipitation, write_conductance = task
    source_precipitation = (
        precipitation_input_directory / f"or_{orbit:04d}.nc"
    )
    precipitation_output = (
        precipitation_output_directory / f"or_{orbit:04d}.nc"
    )

    try:
        if write_precipitation and write_conductance:
            source_conductance = (
                conductance_input_directory / f"or_{orbit:04d}.nc"
            )
            conductance_output = (
                conductance_output_directory / f"or_{orbit:04d}.nc"
            )
            precipitation, conductance = build_detector_cs_products(
                source_precipitation,
                source_conductance,
                software_version=software_version,
            )
            _atomic_save(
                precipitation,
                precipitation_output,
                precipitation_cs_file_status,
                source_precipitation,
            )
            _atomic_save(
                conductance,
                conductance_output,
                conductance_cs_file_status,
                source_conductance,
            )
            return orbit, precipitation.shape[0]

        if write_precipitation:
            precipitation = build_precipitation_cs_product(
                source_precipitation,
                software_version=software_version,
            )
            _atomic_save(
                precipitation,
                precipitation_output,
                precipitation_cs_file_status,
                source_precipitation,
            )
            return orbit, precipitation.shape[0]

        source_conductance = (
            conductance_input_directory / f"or_{orbit:04d}.nc"
        )
        conductance_output = (
            conductance_output_directory / f"or_{orbit:04d}.nc"
        )
        conductance = build_conductance_cs_product(
            source_conductance,
            software_version=software_version,
        )
        _atomic_save(
            conductance,
            conductance_output,
            conductance_cs_file_status,
            source_conductance,
        )
    except Exception as error:
        raise RuntimeError(f"detector CS orbit {orbit:04d} failed") from error

    return orbit, conductance.shape[0]


#%% Command line

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Reduce detector precipitation and conductance independently "
            "onto the frozen 46-by-46 Cubed-Sphere grid."
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
    parser.add_argument(
        "--retrieval-label", default=DEFAULT_RETRIEVAL_LABEL
    )
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
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--precipitation-only",
        action="store_true",
        help=(
            "Write precipitation_cs directly from Product 2 without "
            "requiring or writing conductance products."
        ),
    )
    mode.add_argument(
        "--conductance-only",
        action="store_true",
        help=(
            "Write conductance_cs from Product 3 while requiring an existing "
            "matching precipitation_cs product."
        ),
    )
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
    if not args.conductance_only:
        precipitation_output_directory.mkdir(parents=True, exist_ok=True)
    if not args.precipitation_only:
        conductance_output_directory.mkdir(parents=True, exist_ok=True)

    precipitation_orbits = get_orbits(precipitation_input_directory)
    conductance_orbits = get_orbits(conductance_input_directory)
    if args.precipitation_only:
        available = precipitation_orbits
    elif args.conductance_only:
        available = conductance_orbits
    else:
        available = np.union1d(precipitation_orbits, conductance_orbits)
    selected = available if args.orbit is None else np.unique(args.orbit)
    missing = selected[~np.isin(selected, available)]
    if missing.size:
        source_description = (
            "detector Product 2"
            if args.precipitation_only
            else (
                "detector Product 3"
                if args.conductance_only
                else "detector Product 2 or Product 3"
            )
        )
        raise ValueError(
            f"{source_description} orbit is missing: {missing.tolist()}"
        )

    tasks = []
    for orbit_value in selected:
        orbit = int(orbit_value)
        precipitation_available = np.isin(orbit, precipitation_orbits)
        conductance_available = np.isin(orbit, conductance_orbits)
        write_precipitation = False
        write_conductance = False

        if not args.conductance_only and precipitation_available:
            source_precipitation = (
                precipitation_input_directory / f"or_{orbit:04d}.nc"
            )
            precipitation_output = (
                precipitation_output_directory / f"or_{orbit:04d}.nc"
            )
            precipitation_status = precipitation_cs_file_status(
                precipitation_output, source_precipitation
            ) if not args.overwrite else "invalid"
            if precipitation_status == "mismatch":
                raise ValueError(
                    f"orbit {orbit:04d} precipitation CS output does not "
                    "match its Product-2 source; use another output folder "
                    "or --overwrite"
                )
            write_precipitation = precipitation_status != "complete"

        if not args.precipitation_only and conductance_available:
            source_conductance = (
                conductance_input_directory / f"or_{orbit:04d}.nc"
            )
            conductance_output = (
                conductance_output_directory / f"or_{orbit:04d}.nc"
            )
            conductance_status = (
                "invalid"
                if args.overwrite
                else conductance_cs_file_status(
                    conductance_output,
                    source_conductance,
                )
            )
            if conductance_status == "mismatch":
                raise ValueError(
                    f"orbit {orbit:04d} conductance CS output does not match "
                    "its Product-3 source; use another output folder or "
                    "--overwrite"
                )
            write_conductance = conductance_status != "complete"

        if write_precipitation or write_conductance:
            tasks.append((orbit, write_precipitation, write_conductance))

    product_label = (
        "precipitation CS"
        if args.precipitation_only
        else (
            "conductance CS"
            if args.conductance_only
            else f"detector CS/{args.conductance_model}"
        )
    )
    print(
        f"{product_label}/{DETECTOR_CS_GRID_ID}/{args.retrieval_label}: "
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
