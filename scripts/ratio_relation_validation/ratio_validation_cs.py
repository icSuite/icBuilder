"""Explore CS-grid DMSP matches and the IMAGE ratio relation."""

#%% Imports and paths

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from ratio_validation import (
    MAX_ENERGY_FRACTIONAL_UNCERTAINTY,
    MAX_FLUX_FRACTIONAL_UNCERTAINTY,
    apply_annotations,
    bin_ratio_energy,
    finite_measurements,
    plot_quality_histograms,
    plot_ratio_energy,
    plot_ratio_energy_binned_by_energy,
    plot_selected_distributions,
    quality_selection,
)


REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_MATCHES = Path(
    "/home/bing/Dropbox/work/data/IMAGE_FUV/"
    "ratio_relation_validation/cs_crossings"
)
DEFAULT_ANNOTATIONS = REPOSITORY / "data/dmsp_frame_annotations.csv"
DEFAULT_FIGURES = REPOSITORY / "figures/ratio_relation_validation_cs"
DEFAULT_BINNED = REPOSITORY / "data/ratio_energy_bins_cs.nc"

# This is a provisional centre-to-footprint threshold, matching the detector
# analysis numerically. Review it together with the full CS separation
# histogram because fixed-grid cells have different spatial support.
MAX_CS_SEPARATION_DEG = 1.0

MATCH_COLUMNS = [
    "orbit", "image_frame", "image_time", "dmsp_sat", "dmsp_time",
    "dmsp_mlat", "dmsp_mlt", "cs_row", "cs_column",
    "cs_separation_deg", "cs_inside",
    "dmsp_electron_mean_energy",
    "dmsp_electron_mean_energy_fractional_std",
    "dmsp_electron_total_energy_flux",
    "dmsp_electron_total_energy_flux_fractional_std",
    "img_wic", "img_wic_std", "img_si13", "img_si13_std",
    "img_ratio", "img_ratio_std", "wic_dza", "quality_weight",
    "method_valid",
]


#%% Load per-orbit CS matches

def load_matches(path):
    """Load selected columns from every CS-crossing orbit file."""

    frames = []
    files = sorted(path.glob("or_*.nc"))
    if not files:
        raise FileNotFoundError(f"No CS-crossing files found in {path}")

    for filename in files:
        with xr.open_dataset(filename) as data:
            if data.sizes.get("sample", 0) == 0:
                continue
            if data.attrs.get("representation") != "cs":
                raise ValueError(f"{filename} is not a CS crossing file")

            missing = set(MATCH_COLUMNS).difference(data.variables)
            if missing:
                raise ValueError(f"{filename} is missing {sorted(missing)}")
            frames.append(
                data[MATCH_COLUMNS].to_dataframe().reset_index(drop=True)
            )

    if not frames:
        raise ValueError(f"No CS matches found in {path}")

    matches = pd.concat(frames, ignore_index=True)
    matches["dmsp_sat"] = matches["dmsp_sat"].str.lower().astype("category")
    return matches


def cs_quality_selection(data):
    """Apply DMSP cuts and require an adequately close cell inside the grid."""

    return (
        quality_selection(
            data,
            separation_column="cs_separation_deg",
            maximum_separation=MAX_CS_SEPARATION_DEG,
        )
        & data["cs_inside"]
    )


def add_cs_provenance(binned):
    """Record the selections specific to the fixed-grid analysis."""

    binned.attrs["representation"] = "cs"
    binned.attrs["spatial_selection"] = "cs_inside and nearest cell separation"
    binned.attrs["cs_separation_maximum_degree"] = MAX_CS_SEPARATION_DEG
    binned.attrs["energy_fractional_uncertainty_maximum"] = (
        MAX_ENERGY_FRACTIONAL_UNCERTAINTY
    )
    binned.attrs["flux_fractional_uncertainty_maximum"] = (
        MAX_FLUX_FRACTIONAL_UNCERTAINTY
    )
    # Remove the detector-specific default inherited from the shared binning.
    binned.attrs.pop("detector_separation_maximum_degree", None)
    return binned


#%% Run

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matches", type=Path, default=DEFAULT_MATCHES)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--figures", type=Path, default=DEFAULT_FIGURES)
    parser.add_argument("--binned-output", type=Path, default=DEFAULT_BINNED)
    parser.add_argument("--ignore-annotations", action="store_true")
    args = parser.parse_args()

    matches = load_matches(args.matches)
    if not args.ignore_annotations:
        matches = apply_annotations(matches, args.annotations)
        if matches.empty:
            raise ValueError(
                "No CS matches correspond to accepted annotations"
            )

    finite = matches.loc[
        finite_measurements(matches, "cs_separation_deg")
    ].copy()
    plot_quality_histograms(
        finite, args.figures, separation_column="cs_separation_deg"
    )

    selected = matches.loc[cs_quality_selection(matches)].copy()
    if selected.empty:
        raise ValueError("No CS matches remain after the provisional quality cuts")
    plot_selected_distributions(selected, args.figures)

    binned = add_cs_provenance(bin_ratio_energy(selected))
    args.binned_output.parent.mkdir(parents=True, exist_ok=True)
    binned.to_netcdf(args.binned_output)
    plot_ratio_energy(binned, args.figures, representation="CS-grid")
    plot_ratio_energy_binned_by_energy(
        binned, args.figures, representation="CS-grid"
    )

    print(f"Loaded CS matches: {len(matches):,}")
    print(f"Finite ratio-energy pairs: {len(finite):,}")
    print(f"Selected pairs: {len(selected):,}")
    print(f"Saved {args.binned_output}")
    print(f"Saved figures under {args.figures}")


if __name__ == "__main__":
    main()
