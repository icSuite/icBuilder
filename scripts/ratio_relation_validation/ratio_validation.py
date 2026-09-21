"""Explore detector-matched DMSP quality and the IMAGE ratio relation."""

#%% Imports and paths

import argparse
from pathlib import Path

from icphysics.image import EE_ratio, wic_to_s13
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import xarray as xr


REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_MATCHES = Path("/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test/matched_dmsp_image_data/")
DEFAULT_ANNOTATIONS = REPOSITORY / "data/dmsp_frame_annotations.csv"
DEFAULT_FIGURES = REPOSITORY / "figures/ratio_relation_validation"
DEFAULT_BINNED = REPOSITORY / "data/ratio_energy_bins.nc"

# These are provisional analysis choices, not extraction rules. Revisit them
# after inspecting the unfiltered quality histograms.
MAX_ENERGY_FRACTIONAL_UNCERTAINTY = 0.25
MAX_FLUX_FRACTIONAL_UNCERTAINTY = 0.20
MAX_DETECTOR_SEPARATION_DEG = 1.0

RATIO_EDGES = np.linspace(0, 150, 51)
ENERGY_EDGES = np.linspace(0, 8, 41)
QUANTILES = np.array([0.10, 0.25, 0.50, 0.75, 0.90])
MINIMUM_QUANTILE_COUNT = 20

MATCH_COLUMNS = [
    "orbit", "image_frame", "image_time", "dmsp_sat", "dmsp_time",
    "dmsp_mlat", "dmsp_mlt", "detector_separation_deg",
    "dmsp_electron_mean_energy",
    "dmsp_electron_mean_energy_fractional_std",
    "dmsp_electron_total_energy_flux",
    "dmsp_electron_total_energy_flux_fractional_std",
    "img_wic", "img_wic_std", "img_si13", "img_si13_std",
    "img_ratio", "img_ratio_std", "wic_dza", "quality_weight",
    "method_valid",
]


#%% Load detector matches and migrate the old annotation keys

def load_matches(path):
    """Load selected columns from every detector-crossing orbit file."""

    frames = []
    files = sorted(path.glob("or_*.nc"))
    if not files:
        raise FileNotFoundError(f"No detector-crossing files found in {path}")

    for filename in files:
        with xr.open_dataset(filename) as data:
            # The extractor writes an empty ``sample`` dataset for an orbit
            # with no temporally overlapping DMSP observations.  That is a
            # completed orbit, not a malformed crossing file.
            if data.sizes.get("sample", 0) == 0:
                continue

            missing = set(MATCH_COLUMNS).difference(data.variables)
            if missing:
                raise ValueError(f"{filename} is missing {sorted(missing)}")
            frames.append(
                data[MATCH_COLUMNS].to_dataframe().reset_index(drop=True)
            )

    if not frames:
        raise ValueError(f"No detector matches found in {path}")

    matches = pd.concat(frames, ignore_index=True)
    matches["dmsp_sat"] = matches["dmsp_sat"].str.lower().astype("category")
    return matches


def apply_annotations(matches, filename):
    """Select accepted crossings by orbit, absolute IMAGE time, and satellite."""

    annotations = pd.read_csv(filename)
    annotations = annotations.loc[
        annotations["accepted"].eq(1),
        ["orbit", "img_time", "satellites"],
    ].copy()
    annotations["image_time"] = pd.to_datetime(annotations.pop("img_time"))
    annotations["dmsp_sat"] = annotations.pop("satellites").str.lower()
    annotations = annotations.drop_duplicates(
        ["orbit", "image_time", "dmsp_sat"]
    )

    return matches.merge(
        annotations,
        on=["orbit", "image_time", "dmsp_sat"],
        how="inner",
        validate="many_to_one",
    )


#%% Quality selections

def finite_measurements(data, separation_column="detector_separation_deg"):
    """Require the observables used in the ratio-energy comparison."""

    return (
        np.isfinite(data["img_ratio"])
        & (data["img_ratio"] > 0)
        & np.isfinite(data["dmsp_electron_mean_energy"])
        & (data["dmsp_electron_mean_energy"] > 0)
        & np.isfinite(data["wic_dza"])
        & np.isfinite(data[separation_column])
    )


def quality_selection(
    data,
    separation_column="detector_separation_deg",
    maximum_separation=MAX_DETECTOR_SEPARATION_DEG,
):
    """Apply the provisional DMSP and detector-quality thresholds."""

    return (
        finite_measurements(data, separation_column)
        & data["method_valid"]
        & (
            data["dmsp_electron_mean_energy_fractional_std"]
            <= MAX_ENERGY_FRACTIONAL_UNCERTAINTY
        )
        & (
            data["dmsp_electron_total_energy_flux_fractional_std"]
            <= MAX_FLUX_FRACTIONAL_UNCERTAINTY
        )
        & (
            data[separation_column]
            <= maximum_separation
        )
    )


#%% Initial one-dimensional diagnostics

def histogram_panel(axis, values, bins, title, label):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        axis.set(title=f"{title}\nN=0", xlabel=label, ylabel="Count")
        return
    axis.hist(values, bins=bins, color="0.35", edgecolor="white")
    for quantile in (0.10, 0.50, 0.90):
        value = np.quantile(values, quantile)
        axis.axvline(value, color="tab:orange", linewidth=1.5)
        axis.text(
            value, 0.96, f"{quantile:.0%}: {value:.2g}",
            rotation=90, va="top", ha="right",
            transform=axis.get_xaxis_transform(), fontsize=8
        )
    axis.set(title=f"{title}\nN={values.size:,}", xlabel=label, ylabel="Count")
    axis.grid(alpha=0.2)


def plot_quality_histograms(
    data, output, separation_column="detector_separation_deg"
):
    """Plot the distributions used to choose DMSP and matching cutoffs."""

    figure, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    panels = [
        ("dmsp_electron_mean_energy_fractional_std", np.linspace(0, 1, 101),
         "DMSP energy uncertainty", "Fractional uncertainty"),
        ("dmsp_electron_total_energy_flux_fractional_std", np.linspace(0, 1, 101),
         "DMSP flux uncertainty", "Fractional uncertainty"),
        (separation_column, np.linspace(0, 5, 101),
         "Spatial matching", "Nearest-pixel separation [degree]"),
        ("dmsp_electron_mean_energy", np.linspace(0, 20, 101),
         "DMSP mean energy", "Energy [keV]"),
        ("dmsp_electron_total_energy_flux", 100,
         "DMSP total energy flux", "eV cm$^{-2}$ sr$^{-1}$ s$^{-1}$"),
        ("img_ratio", np.linspace(0, 300, 101),
         "IMAGE WIC/SI13", "Ratio"),
    ]
    for axis, (name, bins, title, label) in zip(axes.flat, panels):
        histogram_panel(axis, data[name], bins, title, label)

    output.mkdir(parents=True, exist_ok=True)
    figure.savefig(output / "quality_histograms.png", dpi=200)
    plt.close(figure)


def plot_selected_distributions(data, output):
    """Show the primary distributions after the provisional cuts."""

    figure, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    panels = [
        ("dmsp_electron_mean_energy", np.linspace(0, 20, 101),
         "DMSP mean energy", "Energy [keV]"),
        ("dmsp_electron_total_energy_flux", 100,
         "DMSP total energy flux", "eV cm$^{-2}$ sr$^{-1}$ s$^{-1}$"),
        ("img_wic", np.linspace(0, 4000, 101),
         "Corrected WIC", "Counts"),
        ("img_si13", np.linspace(0, 80, 101),
         "Corrected SI13", "Counts"),
        ("img_ratio", np.linspace(0, 300, 101),
         "IMAGE WIC/SI13", "Ratio"),
        ("wic_dza", np.linspace(0, 75, 76),
         "WIC detector zenith angle", "Degree"),
    ]
    for axis, (name, bins, title, label) in zip(axes.flat, panels):
        histogram_panel(axis, data[name], bins, title, label)

    figure.savefig(output / "selected_distributions.png", dpi=200)
    plt.close(figure)


#%% Ratio and energy bins

def bin_ratio_energy(data):
    """Bin samples and calculate conditional quantiles along both axes."""

    ratio = np.asarray(data["img_ratio"], dtype=float)
    energy = np.asarray(data["dmsp_electron_mean_energy"], dtype=float)
    in_domain = (
        (ratio >= RATIO_EDGES[0]) & (ratio < RATIO_EDGES[-1])
        & (energy >= ENERGY_EDGES[0]) & (energy < ENERGY_EDGES[-1])
    )
    ratio = ratio[in_domain]
    energy = energy[in_domain]
    counts, _, _ = np.histogram2d(
        ratio, energy, bins=(RATIO_EDGES, ENERGY_EDGES)
    )

    ratio_bin = np.digitize(ratio, RATIO_EDGES) - 1
    energy_bin = np.digitize(energy, ENERGY_EDGES) - 1

    energy_quantiles = np.full((RATIO_EDGES.size - 1, QUANTILES.size), np.nan)
    ratio_bin_sample_count = np.zeros(RATIO_EDGES.size - 1, dtype=int)
    for index in range(RATIO_EDGES.size - 1):
        values = energy[ratio_bin == index]
        values = values[np.isfinite(values)]
        ratio_bin_sample_count[index] = values.size
        if values.size >= MINIMUM_QUANTILE_COUNT:
            energy_quantiles[index] = np.quantile(values, QUANTILES)

    ratio_quantiles = np.full((ENERGY_EDGES.size - 1, QUANTILES.size), np.nan)
    energy_bin_sample_count = np.zeros(ENERGY_EDGES.size - 1, dtype=int)
    for index in range(ENERGY_EDGES.size - 1):
        values = ratio[energy_bin == index]
        values = values[np.isfinite(values)]
        energy_bin_sample_count[index] = values.size
        if values.size >= MINIMUM_QUANTILE_COUNT:
            ratio_quantiles[index] = np.quantile(values, QUANTILES)

    binned = xr.Dataset(
        data_vars={
            "count": (("ratio_bin", "energy_bin"), counts.astype(np.int64)),
            "ratio_bin_sample_count": (
                "ratio_bin", ratio_bin_sample_count
            ),
            "energy_bin_sample_count": (
                "energy_bin", energy_bin_sample_count
            ),
            "energy_quantile": (
                ("ratio_bin", "quantile"), energy_quantiles
            ),
            "ratio_quantile": (
                ("energy_bin", "quantile"), ratio_quantiles
            ),
        },
        coords={
            "ratio_bin": (RATIO_EDGES[:-1] + RATIO_EDGES[1:]) / 2,
            "energy_bin": (ENERGY_EDGES[:-1] + ENERGY_EDGES[1:]) / 2,
            "ratio_edge": RATIO_EDGES,
            "energy_edge": ENERGY_EDGES,
            "quantile": QUANTILES,
        },
        attrs={
            "minimum_quantile_count": MINIMUM_QUANTILE_COUNT,
            "selected_sample_count": len(data),
            "binned_sample_count": int(counts.sum()),
            "energy_fractional_uncertainty_maximum": MAX_ENERGY_FRACTIONAL_UNCERTAINTY,
            "flux_fractional_uncertainty_maximum": MAX_FLUX_FRACTIONAL_UNCERTAINTY,
            "detector_separation_maximum_degree": MAX_DETECTOR_SEPARATION_DEG,
        },
    )
    binned["count"].attrs["units"] = "1"
    binned["ratio_bin_sample_count"].attrs["units"] = "1"
    binned["energy_bin_sample_count"].attrs["units"] = "1"
    binned["energy_quantile"].attrs["units"] = "keV"
    binned["ratio_quantile"].attrs["units"] = "1"
    binned["ratio_bin"].attrs["units"] = "1"
    binned["energy_bin"].attrs["units"] = "keV"
    binned["ratio_edge"].attrs["units"] = "1"
    binned["energy_edge"].attrs["units"] = "keV"
    return binned


def plot_ratio_energy(binned, output, representation="Detector"):
    """Plot the 2-D count distribution, conditional quantiles, and Frey curve."""

    figure, axis = plt.subplots(figsize=(10, 7), constrained_layout=True)
    count = np.asarray(binned["count"])
    plotted = axis.pcolormesh(
        binned.ratio_edge,
        binned.energy_edge,
        np.ma.masked_equal(count.T, 0),
        norm=LogNorm(vmin=1),
        cmap="viridis",
        shading="flat",
    )

    colours = plt.cm.coolwarm(np.linspace(0, 1, QUANTILES.size))
    for index, (quantile, colour) in enumerate(zip(QUANTILES, colours)):
        axis.plot(
            binned.ratio_bin,
            binned.energy_quantile[:, index],
            color=colour,
            linewidth=2,
            label=f"DMSP {quantile:.0%}",
        )

    axis.plot(wic_to_s13, EE_ratio, color="black", linewidth=4)
    axis.plot(
        wic_to_s13, EE_ratio, color="white", linewidth=2.5,
        label="Frey et al. (2003)"
    )
    axis.set(
        xlim=(RATIO_EDGES[0], RATIO_EDGES[-1]),
        ylim=(ENERGY_EDGES[0], ENERGY_EDGES[-1]),
        xlabel="IMAGE corrected WIC/SI13 ratio",
        ylabel="DMSP electron mean energy [keV]",
        title=f"{representation} IMAGE ratio versus DMSP energy",
    )
    axis.legend(frameon=False, ncol=2)
    figure.colorbar(plotted, ax=axis, label="Matched sample count")
    figure.savefig(output / "ratio_energy_relation.png", dpi=220)
    figure.savefig(output / "ratio_energy_relation.pdf")
    plt.close(figure)


def plot_ratio_energy_binned_by_energy(
    binned, output, representation="Detector"
):
    """Plot IMAGE-ratio quantiles calculated within DMSP-energy bins."""

    figure, axis = plt.subplots(figsize=(10, 7), constrained_layout=True)
    count = np.asarray(binned["count"])
    plotted = axis.pcolormesh(
        binned.ratio_edge,
        binned.energy_edge,
        np.ma.masked_equal(count.T, 0),
        norm=LogNorm(vmin=1),
        cmap="viridis",
        shading="flat",
    )

    colours = plt.cm.coolwarm(np.linspace(0, 1, QUANTILES.size))
    for index, (quantile, colour) in enumerate(zip(QUANTILES, colours)):
        axis.plot(
            binned.ratio_quantile[:, index],
            binned.energy_bin,
            color=colour,
            linewidth=2,
            label=f"IMAGE ratio {quantile:.0%}",
        )

    axis.plot(wic_to_s13, EE_ratio, color="black", linewidth=4)
    axis.plot(
        wic_to_s13, EE_ratio, color="white", linewidth=2.5,
        label="Frey et al. (2003)"
    )
    axis.set(
        xlim=(RATIO_EDGES[0], RATIO_EDGES[-1]),
        ylim=(ENERGY_EDGES[0], ENERGY_EDGES[-1]),
        xlabel="IMAGE corrected WIC/SI13 ratio",
        ylabel="DMSP electron mean energy [keV]",
        title=(
            f"{representation} IMAGE-ratio quantiles within "
            "DMSP-energy bins"
        ),
    )
    axis.legend(frameon=False, ncol=2)
    figure.colorbar(plotted, ax=axis, label="Matched sample count")
    figure.savefig(output / "ratio_energy_relation_energy_binned.png", dpi=220)
    figure.savefig(output / "ratio_energy_relation_energy_binned.pdf")
    plt.close(figure)


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
                "No detector matches correspond to accepted annotations"
            )
    finite = matches.loc[finite_measurements(matches)].copy()

    plot_quality_histograms(finite, args.figures)
    selected = matches.loc[quality_selection(matches)].copy()
    if selected.empty:
        raise ValueError("No matches remain after the provisional quality cuts")
    plot_selected_distributions(selected, args.figures)

    binned = bin_ratio_energy(selected)
    args.binned_output.parent.mkdir(parents=True, exist_ok=True)
    binned.to_netcdf(args.binned_output)
    plot_ratio_energy(binned, args.figures)
    plot_ratio_energy_binned_by_energy(binned, args.figures)

    print(f"Loaded matches: {len(matches):,}")
    print(f"Finite ratio-energy pairs: {len(finite):,}")
    print(f"Selected pairs: {len(selected):,}")
    print(f"Saved {args.binned_output}")
    print(f"Saved figures under {args.figures}")


if __name__ == "__main__":
    main()
