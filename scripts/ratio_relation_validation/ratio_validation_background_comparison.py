"""Compare background-subtracted and unsubtracted IMAGE ratio retrievals."""

#%% Imports and paths

import argparse
from pathlib import Path

from icphysics.image import EE_ratio, wic_to_s13
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import xarray as xr

from ratio_validation import (
    DEFAULT_ANNOTATIONS,
    ENERGY_EDGES,
    MATCH_COLUMNS,
    QUANTILES,
    RATIO_EDGES,
    apply_annotations,
    bin_ratio_energy,
    load_matches,
    quality_selection,
)


REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_BACKGROUND = Path(
    "/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test/"
    "matched_dmsp_image_data"
)
DEFAULT_UNSUBTRACTED = Path(
    "/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test/"
    "matched_dmsp_image_data_unsubtracted"
)
DEFAULT_FIGURES = REPOSITORY / "figures/ratio_relation_validation_background_comparison"
DEFAULT_BINNED = REPOSITORY / "data/ratio_energy_bins_background_comparison.nc"

PAIR_KEYS = ["orbit", "image_time", "dmsp_sat", "dmsp_time"]
SZA_MINIMUM = 105.0


#%% Paired selection

def add_wic_sza(matches, path):
    """Add SZA in the same sorted orbit/sample order used by load_matches."""

    pieces = []
    for filename in sorted(path.glob("or_*.nc")):
        with xr.open_dataset(filename) as data:
            if data.sizes.get("sample", 0):
                pieces.append(np.asarray(data["wic_sza"]))
    sza = np.concatenate(pieces)
    if sza.size != len(matches):
        raise ValueError(f"SZA sample count does not match crossings in {path}")
    matches["wic_sza"] = sza
    return matches


def load_paired_matches(background_path, unsubtracted_path, annotations=None):
    """Load both retrievals and retain samples with identical match keys."""

    background = add_wic_sza(load_matches(background_path), background_path)
    unsubtracted = add_wic_sza(
        load_matches(unsubtracted_path), unsubtracted_path
    )

    if annotations is not None:
        background = apply_annotations(background, annotations)
        unsubtracted = apply_annotations(unsubtracted, annotations)

    if background.duplicated(PAIR_KEYS).any():
        raise ValueError("Background crossing data contain duplicate paired keys")
    if unsubtracted.duplicated(PAIR_KEYS).any():
        raise ValueError("Unsubtracted crossing data contain duplicate paired keys")

    return background.merge(
        unsubtracted,
        on=PAIR_KEYS,
        how="inner",
        suffixes=("_background", "_unsubtracted"),
        validate="one_to_one",
    )


def branch_frame(paired, suffix):
    """Recover one ordinary validation table from the paired wide table."""

    data = paired[PAIR_KEYS].copy()
    for name in MATCH_COLUMNS + ["wic_sza"]:
        if name not in PAIR_KEYS:
            data[name] = paired[f"{name}_{suffix}"].to_numpy()
    return data


def select_common_support(paired, minimum_sza=SZA_MINIMUM):
    """Apply the ordinary quality rules to both sides and require darkness."""

    background = branch_frame(paired, "background")
    unsubtracted = branch_frame(paired, "unsubtracted")
    selected = (
        quality_selection(background)
        & quality_selection(unsubtracted)
        & (background["wic_sza"] > minimum_sza)
        & (unsubtracted["wic_sza"] > minimum_sza)
    )
    return paired.loc[selected].copy()


#%% Common-support histogram

def plot_comparison(
    background, unsubtracted, output, crossings, orbits, minimum_sza
):
    """Plot both ratio-energy histograms with identical bins and scaling."""

    products = [bin_ratio_energy(background), bin_ratio_energy(unsubtracted)]
    maximum = max(float(product["count"].max()) for product in products)
    if maximum < 1:
        raise ValueError("No selected samples fall inside the plotted bins")
    norm = LogNorm(vmin=1, vmax=max(2, maximum))
    figure, axes = plt.subplots(1, 2, figsize=(15, 6), sharex=True, sharey=True,
                               constrained_layout=True)

    plotted = None
    titles = ("Background-subtracted counts", "Unsubtracted counts")
    colours = plt.cm.coolwarm(np.linspace(0, 1, QUANTILES.size))
    for axis, product, title in zip(axes, products, titles):
        count = np.asarray(product["count"])
        plotted = axis.pcolormesh(
            product.ratio_edge, product.energy_edge,
            np.ma.masked_equal(count.T, 0), norm=norm, cmap="viridis",
            shading="flat"
        )
        for index, (quantile, colour) in enumerate(zip(QUANTILES, colours)):
            axis.plot(product.ratio_bin, product.energy_quantile[:, index],
                      color=colour, linewidth=1.5, label=f"DMSP {quantile:.0%}")
        axis.plot(wic_to_s13, EE_ratio, color="black", linewidth=4)
        axis.plot(wic_to_s13, EE_ratio, color="white", linewidth=2.5,
                  label="Frey et al. (2003)")
        axis.set(xlim=(RATIO_EDGES[0], RATIO_EDGES[-1]),
                 ylim=(ENERGY_EDGES[0], ENERGY_EDGES[-1]),
                 xlabel="IMAGE corrected WIC/SI13 ratio", title=title)
        axis.legend(framealpha=.5, fontsize=8, ncol=2)

    axes[0].set_ylabel("DMSP electron mean energy [keV]")
    figure.suptitle(
        f"Paired common support, WIC SZA > {minimum_sza:g}°\n"
        f"N={len(background):,} samples, {crossings:,} crossings, {orbits:,} orbits"
    )
    figure.colorbar(plotted, ax=axes, label="Matched sample count")
    output.mkdir(parents=True, exist_ok=True)
    figure.savefig(output / "ratio_energy_background_comparison.png", dpi=220)
    figure.savefig(output / "ratio_energy_background_comparison.pdf")
    plt.close(figure)
    return products


#%% Run

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--background-matches", type=Path, default=DEFAULT_BACKGROUND)
    parser.add_argument("--unsubtracted-matches", type=Path, default=DEFAULT_UNSUBTRACTED)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--figures", type=Path, default=DEFAULT_FIGURES)
    parser.add_argument("--binned-output", type=Path, default=DEFAULT_BINNED)
    parser.add_argument("--minimum-sza", type=float, default=SZA_MINIMUM)
    parser.add_argument("--ignore-annotations", action="store_true")
    args = parser.parse_args()

    annotations = None if args.ignore_annotations else args.annotations
    paired = load_paired_matches(
        args.background_matches, args.unsubtracted_matches, annotations
    )
    selected = select_common_support(paired, args.minimum_sza)
    if selected.empty:
        raise ValueError("No paired samples remain after quality and SZA selection")

    background = branch_frame(selected, "background")
    unsubtracted = branch_frame(selected, "unsubtracted")
    crossing_keys = ["orbit", "image_time", "dmsp_sat"]
    crossings = selected[crossing_keys].drop_duplicates().shape[0]
    orbits = selected["orbit"].nunique()
    products = plot_comparison(
        background, unsubtracted, args.figures, crossings, orbits,
        args.minimum_sza
    )

    combined = xr.concat(products, pd.Index(
        ["background_subtracted", "unsubtracted"], name="count_source"
    ))
    combined.attrs.update({
        "paired_sample_count": len(selected),
        "paired_crossing_count": crossings,
        "paired_orbit_count": orbits,
        "minimum_wic_sza_degree": args.minimum_sza,
        "selection": "ordinary quality rules on both count sources; exact paired keys",
    })
    args.binned_output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_netcdf(args.binned_output)

    print(f"Exact paired matches: {len(paired):,}")
    print(f"Selected common-support samples: {len(selected):,}")
    print(f"Selected crossings: {crossings:,}")
    print(f"Selected orbits: {orbits:,}")
    print(f"Saved {args.binned_output}")
    print(f"Saved figures under {args.figures}")


if __name__ == "__main__":
    main()
