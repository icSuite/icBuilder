"""Build image-ratio precipitation on the WIC detector geometry."""

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

from icbuilder.kp import load_gfz_kp
from icbuilder.precipitationdetector import (
    COUNT_SOURCES,
    COUNT_UNCERTAINTY_MODE,
    PRECIPITATION_METHOD,
    PROTON_ENERGY_MODELS,
    SCHEMA_VERSION,
    SOURCE_TIME_DECODING,
    SPATIAL_SMOOTHING_KERNELS,
    PrecipitationDetector,
    resolve_smoothing_configuration,
)


#%% Product validation and atomic publication

def precipitation_detector_file_status(
    filename,
    source_fuv_detector,
    proton_energy_model,
    proton_energy,
    proton_energy_uncertainty,
    count_source="background_subtracted",
    spatial_smoothing_kernel="none",
    wic_smoothing_width=0.0,
    si13_smoothing_width=0.0,
):
    """Return missing, invalid, mismatch, or complete for one Product-2 file."""

    filename = Path(filename)
    if not filename.is_file():
        return "missing"

    try:
        with open_product(filename) as product:
            if (
                product.product_type != "precipitation_detector"
                or product.representation != "detector"
                or int(product.schema_version) != SCHEMA_VERSION
            ):
                return "invalid"
            if (
                product.method != PRECIPITATION_METHOD
                or product.proton_flux_source != "SI12"
                or product.proton_energy_model != proton_energy_model
                or product.source_fuv_detector != str(source_fuv_detector)
                or product.attrs.get(
                    "count_source", "background_subtracted"
                ) != count_source
                or product.attrs.get(
                    "spatial_smoothing_kernel", "none"
                ) != spatial_smoothing_kernel
            ):
                return "mismatch"
            product_wic_width = float(
                product.attrs.get("wic_smoothing_width_pixels", 0.0)
            )
            product_si13_width = float(
                product.attrs.get("si13_smoothing_width_pixels", 0.0)
            )
            if (
                not np.isclose(product_wic_width, wic_smoothing_width)
                or not np.isclose(product_si13_width, si13_smoothing_width)
            ):
                return "mismatch"
            if (
                product.source_fuv_detector_time_decoding != SOURCE_TIME_DECODING
                or product.count_uncertainty_mode != COUNT_UNCERTAINTY_MODE
                or not hasattr(product, "si12")
                or not hasattr(product, "dsi12")
            ):
                return "invalid"
            if proton_energy_model == "constant" and (
                not np.isclose(product.proton_energy_constant, proton_energy)
                or not np.isclose(
                    product.proton_energy_uncertainty_constant,
                    proton_energy_uncertainty,
                )
            ):
                return "mismatch"
    except (OSError, RuntimeError, KeyError, AttributeError, ValueError):
        return "invalid"

    return "complete"


def save_precipitation_detector(product, filename):
    """Write beside the final path and atomically publish after validation."""

    filename = Path(filename)
    partial = Path(str(filename) + ".partial")
    try:
        product.to_nc(partial)
        status = precipitation_detector_file_status(
            partial,
            product.source_fuv_detector,
            product.proton_energy_model,
            product.proton_energy_constant,
            product.proton_energy_uncertainty_constant,
            product.count_source,
            product.spatial_smoothing_kernel,
            product.wic_smoothing_width_pixels,
            product.si13_smoothing_width_pixels,
        )
        if status != "complete":
            raise RuntimeError(f"incomplete precipitation_detector file: {partial}")
        os.replace(partial, filename)
    finally:
        if partial.exists():
            partial.unlink()


#%% Orbit processing

def get_orbits(input_directory):
    """Discover orbit numbers from detector Product-1 files."""

    orbits = []
    for filename in sorted(input_directory.glob("or_*.nc")):
        orbit_text = filename.stem[-4:]
        if orbit_text.isdigit():
            orbits.append(int(orbit_text))
    return np.unique(orbits)


def current_revision(repository):
    """Return the Git revision recorded in experimental Product 2."""

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
    kp_series,
    proton_energy_model,
    proton_energy,
    proton_energy_uncertainty,
    count_source,
    software_version,
    spatial_smoothing_kernel="none",
    wic_smoothing_width=0.0,
    si13_smoothing_width=0.0,
):
    """Calculate and save one detector-space precipitation orbit."""

    source = input_directory / f"or_{orbit:04d}.nc"
    try:
        product = PrecipitationDetector(
            source,
            kp_series=kp_series,
            proton_energy_model=proton_energy_model,
            proton_energy=proton_energy,
            proton_energy_uncertainty=proton_energy_uncertainty,
            count_source=count_source,
            spatial_smoothing_kernel=spatial_smoothing_kernel,
            wic_smoothing_width=wic_smoothing_width,
            si13_smoothing_width=si13_smoothing_width,
            software_version=software_version,
        )
        output = output_directory / f"or_{orbit:04d}.nc"
        save_precipitation_detector(product, output)
    except Exception as error:
        raise RuntimeError(
            f"precipitation_detector orbit {orbit:04d} failed"
        ) from error

    return orbit, product.shape[0]


#%% Command line

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Create image-ratio precipitation on WIC detector pixels."
    )
    parser.add_argument(
        "--base-input", "--base",
        dest="base_input",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "example_data",
        help="Base directory containing detector Product-1 orbit files.",
    )
    parser.add_argument(
        "--base-output",
        type=Path,
        help="Base directory for Product 2 (default: --base-input).",
    )
    parser.add_argument(
        "--input-folder", default="fuv_detector/fuvpy_bs_directional_v1"
    )
    parser.add_argument("--output-folder", default="precipitation_detector")
    parser.add_argument(
        "--retrieval-label",
        help=(
            "output subfolder; otherwise derived from proton model, count "
            "source, smoothing kernel, and widths"
        ),
    )
    parser.add_argument("--orbit", action="append", type=int)
    parser.add_argument(
        "--proton-energy-model",
        choices=PROTON_ENERGY_MODELS,
        default="hardy",
        help="Proton mean-energy model (default: hardy).",
    )
    parser.add_argument("--proton-energy", type=float, default=2.0)
    parser.add_argument(
        "--proton-energy-uncertainty", type=float, default=0.0
    )
    parser.add_argument(
        "--count-source",
        choices=COUNT_SOURCES,
        default="background_subtracted",
        help="Product-1 count stage (default: background_subtracted).",
    )
    parser.add_argument(
        "--spatial-smoothing-kernel",
        choices=SPATIAL_SMOOTHING_KERNELS,
        default="none",
        help=(
            "Diagnostic WIC/SI13 smoothing before proton correction "
            "(default: none)."
        ),
    )
    parser.add_argument(
        "--wic-smoothing-width",
        type=float,
        help=(
            "WIC Gaussian sigma or odd boxcar side width in WIC pixels; "
            "defaults to 1 for Gaussian and 3 for boxcar."
        ),
    )
    parser.add_argument(
        "--si13-smoothing-width",
        type=float,
        help=(
            "SI13 Gaussian sigma or odd boxcar side width in WIC pixels; "
            "defaults to 1 for Gaussian and 3 for boxcar."
        ),
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
    smoothing_widths = resolve_smoothing_configuration(
        args.spatial_smoothing_kernel,
        args.wic_smoothing_width,
        args.si13_smoothing_width,
    )
    wic_smoothing_width, si13_smoothing_width = smoothing_widths

    base_input = args.base_input.expanduser()
    base_output = (
        args.base_output.expanduser()
        if args.base_output is not None else base_input
    )
    input_directory = base_input / args.input_folder
    if args.retrieval_label is None:
        count_suffix = (
            ""
            if args.count_source == "background_subtracted"
            else "_unsubtracted"
        )
        if args.spatial_smoothing_kernel == "none":
            smoothing_suffix = ""
        else:
            wic_label = str(wic_smoothing_width).replace(".", "p")
            si13_label = str(si13_smoothing_width).replace(".", "p")
            smoothing_suffix = (
                f"_smooth_{args.spatial_smoothing_kernel}"
                f"_wic{wic_label}_si13{si13_label}"
            )
        retrieval_label = (
            f"IR_{args.proton_energy_model}{count_suffix}{smoothing_suffix}"
        )
    else:
        retrieval_label = args.retrieval_label
    output_directory = base_output / args.output_folder / retrieval_label
    output_directory.mkdir(parents=True, exist_ok=True)

    available = get_orbits(input_directory)
    selected = available if args.orbit is None else np.unique(args.orbit)
    missing = selected[~np.isin(selected, available)]
    if missing.size:
        raise ValueError(f"fuv_detector orbit is missing: {missing.tolist()}")

    pending = []
    for orbit in selected:
        source = input_directory / f"or_{int(orbit):04d}.nc"
        output = output_directory / f"or_{int(orbit):04d}.nc"
        if args.overwrite:
            pending.append(int(orbit))
            continue
        status = precipitation_detector_file_status(
            output,
            source,
            args.proton_energy_model,
            args.proton_energy,
            args.proton_energy_uncertainty,
            args.count_source,
            args.spatial_smoothing_kernel,
            wic_smoothing_width,
            si13_smoothing_width,
        )
        if status == "mismatch":
            raise ValueError(
                f"{output} does not match the requested detector "
                "precipitation configuration; use another retrieval label "
                "or --overwrite"
            )
        if status != "complete":
            pending.append(int(orbit))

    print(
        f"precipitation_detector/{retrieval_label}: "
        f"{len(selected) - len(pending)} complete, {len(pending)} pending"
    )
    if not pending:
        return []

    kp_series = load_gfz_kp()
    repository = Path(__file__).resolve().parents[2]
    software_version = current_revision(repository)
    function = partial(
        process_orbit,
        input_directory=input_directory,
        output_directory=output_directory,
        kp_series=kp_series,
        proton_energy_model=args.proton_energy_model,
        proton_energy=args.proton_energy,
        proton_energy_uncertainty=args.proton_energy_uncertainty,
        count_source=args.count_source,
        spatial_smoothing_kernel=args.spatial_smoothing_kernel,
        wic_smoothing_width=wic_smoothing_width,
        si13_smoothing_width=si13_smoothing_width,
        software_version=software_version,
    )
    if args.workers > 1:
        return process_map(
            function,
            pending,
            max_workers=args.workers,
            chunksize=1,
            desc="Create detector precipitation orbits",
        )

    return [
        function(orbit)
        for orbit in tqdm(pending, desc="Create detector precipitation orbits")
    ]


if __name__ == "__main__":
    main()
