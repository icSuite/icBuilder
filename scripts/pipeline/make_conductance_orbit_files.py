"""Build conductance orbit files from completed precipitation products."""

#%% Imports

import argparse
import os
from functools import partial
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

from icbuilder import ConductanceImage
from icbuilder.precipitationdetector import ELECTRON_METADATA_ATTRIBUTES


#%% File handling

REQUIRED_FIELDS = (
    "Ep_model", "Ep", "dEp", "Ep_clipping_flag", "Fp", "dFp",
    "E0", "dE0", "Fe", "dFe", "varE0Fe", "P", "H", "dP", "dH", "w",
)


def get_orbits(input_dir):
    """Return orbit numbers discovered from precipitation NetCDF files."""

    orbits = []
    for filename in sorted(input_dir.glob("or_*.nc")):
        orbit_text = filename.stem[-4:]
        if orbit_text.isdigit():
            orbits.append(int(orbit_text))
    return orbits


def conductance_file_is_complete(filename):
    """Check the small Product-3 schema used for restart decisions."""

    filename = Path(filename)
    if not filename.is_file():
        return False

    try:
        with Dataset(filename) as nc:
            if nc.product_type != "conductance" or int(nc.schema_version) != 2:
                return False
            shape = (
                len(nc.dimensions["time"]),
                len(nc.dimensions["dim1"]),
                len(nc.dimensions["dim2"]),
            )
            if shape[0] == 0:
                return False
            for name in REQUIRED_FIELDS:
                if nc.variables[name].shape != shape:
                    return False
            for name in ("time", "Kp", "Kp_interval_start", "ssalon"):
                if nc.variables[name].shape != (shape[0],):
                    return False
            if nc.precipitation_method == "hardy":
                for name in ELECTRON_METADATA_ATTRIBUTES:
                    nc.getncattr(name)
                if nc.variables["E0_response_range_flag"].shape != shape:
                    return False
            nc.getncattr("precipitation_method")
            nc.getncattr("proton_flux_source")
            nc.getncattr("proton_energy_model")
            nc.getncattr("conductance_model")
            nc.groups["grid"]
    except (OSError, RuntimeError, KeyError, AttributeError, ValueError):
        return False

    return True


def conductance_matches_precipitation(filename, precipitation_file):
    """Check that an existing Product-3 file matches its requested Product 2."""

    if not conductance_file_is_complete(filename):
        return False

    try:
        with Dataset(filename) as conductance, Dataset(precipitation_file) as precipitation:
            matches = (
                conductance.precipitation_method == precipitation.method
                and conductance.proton_flux_source == precipitation.proton_flux_source
                and conductance.proton_energy_model == precipitation.proton_energy_model
            )
            if not matches:
                return False
            if precipitation.method == "hardy":
                if any(conductance.getncattr(name) != precipitation.getncattr(name)
                       for name in ELECTRON_METADATA_ATTRIBUTES):
                    return False
            if precipitation.proton_energy_model == "constant":
                return (
                    np.isclose(
                        conductance.proton_energy_constant,
                        precipitation.proton_energy_constant,
                    )
                    and np.isclose(
                        conductance.proton_energy_uncertainty_constant,
                        precipitation.proton_energy_uncertainty_constant,
                    )
                )
            return True
    except (OSError, RuntimeError, KeyError, AttributeError, ValueError):
        return False


def save_conductance_file(conductance, filename):
    """Write beside the final file and publish it after validation."""

    filename = Path(filename)
    partial_file = Path(str(filename) + ".partial")

    try:
        conductance.to_nc(partial_file)
        if not conductance_file_is_complete(partial_file):
            raise RuntimeError(f"incomplete conductance file: {partial_file}")
        os.replace(partial_file, filename)
    finally:
        if partial_file.exists():
            partial_file.unlink()


def process_orbit(orbit, input_dir, output_dir):
    """Convert one precipitation orbit into conductance."""

    orbit_file = f"or_{orbit:04d}.nc"
    precipitation_file = input_dir / orbit_file
    output_file = output_dir / orbit_file

    conductance = ConductanceImage(precipitation_file)
    save_conductance_file(conductance, output_file)

    return orbit, conductance.shape[0]


#%% Command line

def str2bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in ("yes", "true", "t", "1"):
        return True
    if value.lower() in ("no", "false", "f", "0"):
        return False
    raise argparse.ArgumentTypeError("Boolean value expected")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Create conductance orbit files from precipitation files.")
    parser.add_argument("--parallel", type=str2bool, default=False)
    parser.add_argument("--pool_size", type=int, default=4)
    parser.add_argument("--base", type=Path,
                        default=Path(__file__).resolve().parents[1] / "example_data")
    parser.add_argument("--input-folder", default="precipitation")
    parser.add_argument("--output-folder", default="conductance")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    base = args.base.expanduser()
    input_dir = base / args.input_folder
    output_dir = base / args.output_folder
    output_dir.mkdir(parents=True, exist_ok=True)

    orbits = get_orbits(input_dir)
    pending = []
    for orbit in orbits:
        output_file = output_dir / f"or_{orbit:04d}.nc"
        precipitation_file = input_dir / f"or_{orbit:04d}.nc"
        if args.overwrite or not output_file.exists():
            pending.append(orbit)
        elif not conductance_matches_precipitation(
            output_file, precipitation_file
        ):
            raise ValueError(
                f"{output_file} is invalid or does not match "
                f"{precipitation_file}; choose another output folder or "
                "use --overwrite"
            )

    print(f"Conductance: {len(orbits) - len(pending)} complete, {len(pending)} pending")
    if not pending:
        return []

    function = partial(process_orbit, input_dir=input_dir, output_dir=output_dir)
    if args.parallel:
        return process_map(function, pending, max_workers=args.pool_size, chunksize=1, desc="Create conductance orbits")

    return [function(orbit) for orbit in tqdm(pending, desc="Create conductance orbits")]


if __name__ == "__main__":
    main()
