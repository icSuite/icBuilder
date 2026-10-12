"""Build method-specific precipitation orbit files from binned FUV images."""

import argparse
import os
from functools import partial
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

from icbuilder import PrecipitationImage
from icbuilder.kp import load_gfz_kp
from icbuilder.precipitationimage import REGRID_METHOD, REGRID_UNCERTAINTY
from icbuilder.precipitationdetector import hardy_electron_metadata
from icbuilder.precipitationimage import HARDY_COORDINATE_NOTE


#%%

PRECIPITATION_FIELDS = (
    "wic", "dwic", "si12", "dsi12", "si13", "dsi13",
    "wic_weight", "si12_weight", "si13_weight", "w",
    "wic_corrected", "dwic_corrected", "si13_corrected",
    "dsi13_corrected", "Ep_model", "Ep", "dEp", "Ep_clipping_flag",
    "Fp", "dFp", "E0", "dE0", "Fe", "dFe", "varE0Fe",
)

def save_precipitation_file(precipitation, filename):
    """Write beside the final file, then rename after a successful close."""

    filename = Path(filename)
    partial_file = Path(str(filename) + ".partial")

    try:
        precipitation.to_nc(partial_file)
        status = precipitation_file_status(
            partial_file,
            precipitation.method,
            precipitation.proton_energy_model,
            precipitation.proton_energy_constant,
            precipitation.proton_energy_uncertainty_constant,
        )
        if status != "complete":
            raise RuntimeError(f"incomplete precipitation file: {partial_file}")
        os.replace(partial_file, filename)
    finally:
        if partial_file.exists():
            partial_file.unlink()


def precipitation_file_status(
    filename,
    method,
    proton_energy_model,
    proton_energy,
    proton_energy_uncertainty,
):
    """Return missing, invalid, mismatch, or complete for one Product-2 file."""

    filename = Path(filename)
    if not filename.is_file():
        return "missing"

    try:
        with Dataset(filename) as nc:
            if nc.product_type != "precipitation" or int(nc.schema_version) != 2:
                return "invalid"
            if (
                nc.method != method
                or nc.proton_flux_source != "SI12"
                or nc.proton_energy_model != proton_energy_model
                or nc.regrid_method != REGRID_METHOD
                or nc.regrid_uncertainty != REGRID_UNCERTAINTY
            ):
                return "mismatch"
            if proton_energy_model == "constant" and (
                not np.isclose(nc.proton_energy_constant, proton_energy)
                or not np.isclose(
                    nc.proton_energy_uncertainty_constant,
                    proton_energy_uncertainty,
                )
            ):
                return "mismatch"

            shape = (
                len(nc.dimensions["time"]),
                len(nc.dimensions["dim1"]),
                len(nc.dimensions["dim2"]),
            )
            if shape[0] == 0:
                return "invalid"
            for name in PRECIPITATION_FIELDS:
                if nc.variables[name].shape != shape:
                    return "invalid"
            if method == "hardy":
                for name, value in hardy_electron_metadata(HARDY_COORDINATE_NOTE).items():
                    if getattr(nc, name, None) != value:
                        return "mismatch"
                if nc.variables["E0_response_range_flag"].shape != shape:
                    return "invalid"
            if method == "image_ratio":
                for name in ("R", "dR"):
                    if nc.variables[name].shape != shape:
                        return "invalid"
            for name in (
                "time", "wic_source_index", "si12_source_index",
                "si13_source_index", "Kp", "Kp_interval_start", "ssalon",
            ):
                if nc.variables[name].shape != (shape[0],):
                    return "invalid"
            grid = nc.groups["grid"]
            for name in ("xi", "eta", "mlat", "mlt"):
                if grid.variables[name].shape != shape[1:]:
                    return "invalid"
    except (OSError, RuntimeError, KeyError, AttributeError, ValueError):
        return "invalid"

    return "complete"


def process_orbit(
    orbit,
    input_paths,
    output_dir,
    kp_series,
    method,
    proton_energy_model,
    proton_energy,
    proton_energy_uncertainty,
):
    """Create one precipitation product from the required binned sensors."""

    orbit_file = f"or_{orbit:04d}.nc"

    wic_file = input_paths["wic"] / orbit_file
    si12_file = input_paths["si12"] / orbit_file

    si13_file = input_paths["si13"] / orbit_file
    if not si13_file.exists():
        si13_file = None

    precipitation = PrecipitationImage(
        wic=wic_file,
        si12=si12_file,
        si13=si13_file,
        method=method,
        kp_series=kp_series,
        proton_energy_model=proton_energy_model,
        proton_energy=proton_energy,
        proton_energy_uncertainty=proton_energy_uncertainty,
    )

    output_file = output_dir / orbit_file
    save_precipitation_file(precipitation, output_file)

    return orbit, precipitation.time.size

def get_orbits(input_dir):
    """Return orbit numbers discovered from the NetCDF files present."""

    orbits = []
    for filename in sorted(input_dir.glob("*.nc")):
        orbit_text = filename.stem[-4:]
        if orbit_text.isdigit():
            orbits.append(int(orbit_text))
    return np.unique(orbits)


def str2bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in ("yes", "true", "t", "1"):
        return True
    if value.lower() in ("no", "false", "f", "0"):
        return False
    raise argparse.ArgumentTypeError("Boolean value expected")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Create precipitation orbit files from binned IMAGE-FUV data.")
    parser.add_argument("--parallel", type=str2bool, default=False)
    parser.add_argument("--pool_size", type=int, default=4)
    parser.add_argument("--base", type=Path, default=Path(__file__).resolve().parents[1] / "example_data")
    parser.add_argument("--wic-folder", default="binned/wic")
    parser.add_argument("--s12-folder", default="binned/si12")
    parser.add_argument("--s13-folder", default="binned/si13")
    parser.add_argument("--output-folder", default="precipitation")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--precipitation-method", "--precipitation_method",
        dest="precipitation_method",
        choices=("zhang_paxton", "image_ratio", "hardy"),
        default="zhang_paxton",
    )
    parser.add_argument(
        "--proton-energy-model",
        choices=("hardy", "constant"),
        default="hardy",
        help="proton mean-energy model (default: hardy)",
    )
    parser.add_argument(
        "--proton-energy", type=float, default=2.0,
        help="constant proton mean energy in keV; used only with --proton-energy-model constant",
    )
    parser.add_argument(
        "--proton-energy-uncertainty", type=float, default=0.0,
        help="constant dEp in keV; used only with --proton-energy-model constant",
    )
    return parser.parse_args(argv)


#%% Build the available orbits

def main(argv=None):
    args = parse_args(argv)

    base = args.base.expanduser()
    input_paths = {"wic": base / args.wic_folder, 
                   "si12": base / args.s12_folder, 
                   "si13": base / args.s13_folder}
    output_dir = base / args.output_folder
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Precipitation method: {args.precipitation_method}")
    print("Proton flux source: SI12")
    print(f"Proton energy model: {args.proton_energy_model}")
    if args.proton_energy_model == "constant":
        print(
            f"Constant Ep = {args.proton_energy:g} keV, "
            f"dEp = {args.proton_energy_uncertainty:g} keV"
        )
    print(f"Output: {output_dir}\n")

    # Zhang-Paxton needs WIC and SI12. The image-ratio method also needs SI13.
    common_orbits = np.intersect1d(get_orbits(input_paths["wic"]), get_orbits(input_paths["si12"]))
    if args.precipitation_method == "image_ratio":
        common_orbits = np.intersect1d(common_orbits, get_orbits(input_paths["si13"]))

    pending = []
    for orbit in common_orbits:
        if args.overwrite:
            pending.append(int(orbit))
            continue

        filename = output_dir / f"or_{orbit:04d}.nc"
        status = precipitation_file_status(
            filename,
            args.precipitation_method,
            args.proton_energy_model,
            args.proton_energy,
            args.proton_energy_uncertainty,
        )
        if status == "mismatch":
            raise ValueError(
                f"{filename} does not match the requested precipitation "
                "configuration; choose another output folder or use --overwrite"
            )
        if status != "complete":
            pending.append(int(orbit))

    print(
        f"{args.precipitation_method}: "
        f"{len(common_orbits) - len(pending)} complete, {len(pending)} pending"
    )
    if not pending:
        print("All precipitation files already exist.")
        return []

    # Load once before multiprocessing rather than once for every orbit.
    kp_series = load_gfz_kp()

    function = partial(
        process_orbit,
        input_paths=input_paths,
        output_dir=output_dir,
        kp_series=kp_series,
        method=args.precipitation_method,
        proton_energy_model=args.proton_energy_model,
        proton_energy=args.proton_energy,
        proton_energy_uncertainty=args.proton_energy_uncertainty,
    )

    if args.parallel:
        return process_map(function, pending, max_workers=args.pool_size, chunksize=1, desc="Create precipitation orbits")

    return [function(orbit) for orbit in tqdm(pending, desc="Create precipitation orbits")]


#%%

if __name__ == "__main__":
    main()
