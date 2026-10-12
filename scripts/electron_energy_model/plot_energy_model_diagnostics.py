"""Plot DMSP and Zhang--Paxton electron-energy model diagnostics.

The input is one result directory written by ``evaluate_energy_models.py``.
Each invocation therefore describes one explicit flux selection.  Run the
script once for the no-floor results and once for the strong-signal results
rather than mixing differently selected samples in one opaque figure.
"""

#%% Imports and plotting choices

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd


REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = REPOSITORY / "data/electron_energy_model/results"
DEFAULT_FIGURES = REPOSITORY / "figures/electron_energy_model"
DEFAULT_SUPPORT_SECONDS = 120

COHORTS = ("all", "accepted")
MODEL_ORDER = (
    "constant_2kev",
    "constant_physical",
    "constant_log_median",
    "hardy_pointwise",
    "zhang_paxton_pointwise",
    "zhang_paxton_area_mean",
    "zhang_paxton_area_median",
    "zhang_paxton_flux_weighted",
    "zhang_paxton_scaled",
    "zhang_paxton_calibrated"
)
MODEL_LABELS = {
    "constant_2kev": "hard-coded 2 keV",
    "constant_physical": "physical constant",
    "constant_log_median": "log-median constant",
    "hardy_pointwise": "Hardy pointwise",
    "zhang_paxton_pointwise": "ZP pointwise",
    "zhang_paxton_area_mean": "ZP area mean",
    "zhang_paxton_area_median": "ZP area median",
    "zhang_paxton_flux_weighted": "ZP flux weighted",
    "zhang_paxton_scaled": "ZP area mean, scaled",
    "zhang_paxton_calibrated": "ZP area mean, calibrated"
}
MODEL_COLOURS = {
    "constant_2kev": "tab:brown",
    "constant_physical": "0.60",
    "constant_log_median": "0.20",
    "hardy_pointwise": "tab:pink",
    "zhang_paxton_pointwise": "tab:cyan",
    "zhang_paxton_area_mean": "tab:orange",
    "zhang_paxton_area_median": "tab:purple",
    "zhang_paxton_flux_weighted": "tab:red",
    "zhang_paxton_scaled": "tab:green",
    "zhang_paxton_calibrated": "tab:blue"
}
ZP_DISTRIBUTIONS = {
    "hardy_pointwise_energy": ("Hardy pointwise", "tab:pink"),
    "zp_pointwise_energy": ("ZP pointwise", "tab:cyan"),
    "zp_area_mean_energy": ("ZP area mean", "tab:orange"),
    "zp_area_median_energy": ("ZP area median", "tab:purple"),
    "zp_flux_weighted_energy": ("ZP flux weighted", "tab:red")
}
RESIDUAL_MODELS = (
    "constant_2kev",
    "constant_log_median",
    "hardy_pointwise",
    "zhang_paxton_pointwise",
    "zhang_paxton_area_mean",
    "zhang_paxton_scaled",
    "zhang_paxton_calibrated"
)

ENERGY_BINS = np.geomspace(0.03, 30, 55)
FLUX_BINS = np.geomspace(5e8, 2e15, 55)
RESIDUAL_BINS = np.linspace(-2, 2, 65)


#%% Load one explicit evaluation

def load_results(path, support_seconds):
    """Load and validate the tables needed for plotting."""

    predictions = pd.read_csv(path / "predictions.csv")
    scores = pd.read_csv(path / "scores.csv")
    counts = pd.read_csv(path / "sample_counts.csv")

    prediction_columns = {
        "cohort", "support_seconds", "dmsp_energy", "dmsp_flux",
        *(name for name in ZP_DISTRIBUTIONS if name != "hardy_pointwise_energy" or name in predictions),
        *(f"prediction_{model}" for model in MODEL_ORDER if model != "hardy_pointwise" or "hardy_pointwise_energy" in predictions)
    }
    missing = prediction_columns.difference(predictions.columns)
    if missing:
        raise ValueError(f"{path / 'predictions.csv'} is missing {sorted(missing)}")

    predictions = predictions.loc[
        predictions["support_seconds"].eq(support_seconds)
    ].copy()
    scores = scores.loc[scores["support_seconds"].eq(support_seconds)].copy()
    counts = counts.loc[counts["support_seconds"].eq(support_seconds)].copy()
    if predictions.empty or scores.empty or counts.empty:
        raise ValueError(f"No {support_seconds}-s results found in {path}")

    available = set(predictions["cohort"])
    cohorts = [cohort for cohort in COHORTS if cohort in available]
    return predictions, scores, counts, cohorts


def selection_description(counts, support_seconds):
    """Describe the temporal support and absolute flux selection."""

    minimum_flux = counts["minimum_flux"].unique()
    maximum_flux = counts["maximum_flux"].unique()
    if len(minimum_flux) != 1 or len(maximum_flux) != 1:
        raise ValueError("result directory mixes multiple flux selections")

    lower = minimum_flux[0]
    upper = maximum_flux[0]
    if lower <= 0 and np.isinf(upper):
        flux_text = "no absolute flux floor"
    elif np.isinf(upper):
        flux_text = f"Q >= {lower:.1e} eV cm$^{{-2}}$ sr$^{{-1}}$ s$^{{-1}}$"
    else:
        flux_text = f"{lower:.1e} <= Q <= {upper:.1e}"
    return f"{support_seconds}-s support; {flux_text}"


def save_figure(figure, output, stem):
    """Save matching high-resolution raster and vector figures."""

    output.mkdir(parents=True, exist_ok=True)
    figure.savefig(output / f"{stem}.png", dpi=220)
    figure.savefig(output / f"{stem}.pdf")
    plt.close(figure)


#%% Histograms of the selected samples

def plot_sample_histograms(predictions, cohorts, description, output):
    """Compare DMSP and ZP distributions and show selected DMSP flux."""

    figure, axes = plt.subplots(
        len(cohorts), 2, figsize=(12, 4.2 * len(cohorts)),
        constrained_layout=True, squeeze=False
    )
    for row, cohort in enumerate(cohorts):
        selected = predictions.loc[predictions["cohort"].eq(cohort)]

        axes[row, 0].hist(
            selected["dmsp_energy"], bins=ENERGY_BINS,
            weights=np.full(len(selected), 1 / len(selected)),
            histtype="step", linewidth=2.5, color="black", label="DMSP"
        )
        for name, (label, colour) in ZP_DISTRIBUTIONS.items():
            if name not in selected:
                continue
            axes[row, 0].hist(
                selected[name], bins=ENERGY_BINS,
                weights=np.full(len(selected), 1 / len(selected)),
                histtype="step", linewidth=1.8, color=colour, label=label
            )
        axes[row, 0].set(
            xscale="log", xlim=(ENERGY_BINS[0], ENERGY_BINS[-1]),
            xlabel="Electron mean energy [keV]", ylabel="Fraction per log bin",
            title=f"{cohort.capitalize()} crossings: energy, N={len(selected):,}"
        )
        axes[row, 0].legend(frameon=False, fontsize=9)

        axes[row, 1].hist(
            selected["dmsp_flux"], bins=FLUX_BINS,
            color="0.35", edgecolor="white"
        )
        median_flux = np.median(selected["dmsp_flux"])
        axes[row, 1].axvline(
            median_flux, color="tab:orange", linewidth=2,
            label=f"median {median_flux:.2e}"
        )
        axes[row, 1].set(
            xscale="log", xlim=(FLUX_BINS[0], FLUX_BINS[-1]),
            xlabel="DMSP electron energy flux [eV cm$^{-2}$ sr$^{-1}$ s$^{-1}$]",
            ylabel="Crossings",
            title=f"{cohort.capitalize()} crossings: DMSP flux"
        )
        axes[row, 1].legend(frameon=False)
        for axis in axes[row]:
            axis.grid(alpha=0.2)

    figure.suptitle(f"Electron-energy validation samples\n{description}")
    save_figure(figure, output, "sample_histograms")


#%% Direct Zhang--Paxton comparison

def plot_zhang_paxton_density(predictions, cohorts, description, output):
    """Plot each raw statistical-model energy against matched DMSP energy."""

    distributions = {name: style for name, style in ZP_DISTRIBUTIONS.items() if name in predictions}

    figure, axes = plt.subplots(
        len(cohorts), len(distributions),
        figsize=(4.8 * len(distributions), 4.3 * len(cohorts)),
        constrained_layout=True, squeeze=False, sharex=True, sharey=True
    )
    limits = (0.03, 30)
    for row, cohort in enumerate(cohorts):
        selected = predictions.loc[predictions["cohort"].eq(cohort)]
        for column, (name, (label, _)) in enumerate(distributions.items()):
            axis = axes[row, column]
            plotted = axis.hexbin(
                selected[name], selected["dmsp_energy"],
                xscale="log", yscale="log", gridsize=42, mincnt=1,
                norm=LogNorm(), cmap="viridis"
            )
            axis.plot(limits, limits, color="white", linewidth=3)
            axis.plot(limits, limits, color="black", linewidth=1.4)
            axis.set(
                xlim=limits, ylim=limits,
                xlabel=f"{label} [keV]",
                ylabel="DMSP mean energy [keV]" if column == 0 else None,
                title=f"{cohort.capitalize()}, N={len(selected):,}"
            )
            figure.colorbar(plotted, ax=axis, label="Crossings")

    figure.suptitle(f"Raw electron-energy model predictions versus DMSP\n{description}")
    save_figure(figure, output, "zhang_paxton_vs_dmsp")


#%% Held-out residual histograms

def plot_residual_histograms(predictions, cohorts, description, output):
    """Show signed held-out log errors for the principal model variants."""

    figure, axes = plt.subplots(
        len(cohorts), 1, figsize=(10, 4.2 * len(cohorts)),
        constrained_layout=True, squeeze=False
    )
    for row, cohort in enumerate(cohorts):
        selected = predictions.loc[predictions["cohort"].eq(cohort)]
        axis = axes[row, 0]
        truth = np.log10(selected["dmsp_energy"])
        for model in RESIDUAL_MODELS:
            if model == "hardy_pointwise" and f"prediction_{model}" not in selected:
                continue
            residual = truth - np.log10(selected[f"prediction_{model}"])
            axis.hist(
                residual, bins=RESIDUAL_BINS, density=True,
                histtype="step", linewidth=2, color=MODEL_COLOURS[model],
                label=(f"{MODEL_LABELS[model]} "
                       f"(median {np.median(residual):+.2f} dex)")
            )
        axis.axvline(0, color="black", linewidth=1)
        axis.set(
            xlim=(RESIDUAL_BINS[0], RESIDUAL_BINS[-1]),
            xlabel="log$_{10}$(DMSP / held-out prediction)",
            ylabel="Probability density",
            title=f"{cohort.capitalize()} crossings, N={len(selected):,}"
        )
        axis.legend(frameon=False, fontsize=9)
        axis.grid(alpha=0.2)

    figure.suptitle(f"Held-out electron-energy residuals\n{description}")
    save_figure(figure, output, "held_out_residual_histograms")


#%% Equal-orbit score and improvement intervals

def plot_model_scores(scores, cohorts, description, output):
    """Plot primary error and paired improvement over the log constant."""

    figure, axes = plt.subplots(
        len(cohorts), 2, figsize=(13, 5.2 * len(cohorts)),
        constrained_layout=True, squeeze=False
    )
    model_order = tuple(model for model in MODEL_ORDER if model != "hardy_pointwise" or model in set(scores["model"]))
    positions = np.arange(len(model_order))
    labels = [MODEL_LABELS[model] for model in model_order]

    for row, cohort in enumerate(cohorts):
        selected = scores.loc[scores["cohort"].eq(cohort)].set_index("model")
        selected = selected.loc[list(model_order)]
        colours = [MODEL_COLOURS[model] for model in model_order]

        axes[row, 0].barh(
            positions, selected["grouped_mean_absolute_log_error"],
            color=colours
        )
        axes[row, 0].set_yticks(positions, labels)
        axes[row, 0].invert_yaxis()
        axes[row, 0].set(
            xlabel="Equal-orbit mean absolute log error [dex]",
            title=f"{cohort.capitalize()} crossings"
        )

        improvement = selected["improvement_over_log_constant"].to_numpy()
        low = selected["improvement_ci_low"].to_numpy()
        high = selected["improvement_ci_high"].to_numpy()
        baseline = np.array(model_order) == "constant_log_median"
        low[baseline] = 0
        high[baseline] = 0
        uncertainty = np.vstack((improvement - low, high - improvement))
        for position, value, error, colour in zip(
            positions, improvement, uncertainty.T, colours
        ):
            axes[row, 1].errorbar(
                value, position, xerr=error[:, None],
                fmt="o", color=colour, ecolor=colour, capsize=3
            )
        axes[row, 1].axvline(0, color="0.4", linewidth=1)
        axes[row, 1].set_yticks(positions, labels)
        axes[row, 1].invert_yaxis()
        axes[row, 1].set(
            xlabel="Improvement over log-median constant [dex]",
            title="Positive is better; 95% paired bootstrap interval"
        )
        for axis in axes[row]:
            axis.grid(axis="x", alpha=0.2)

    figure.suptitle(f"Held-out electron-energy model scores\n{description}")
    save_figure(figure, output, "model_scores")


#%% Command line

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--figures", type=Path)
    parser.add_argument(
        "--support-seconds", type=int, default=DEFAULT_SUPPORT_SECONDS
    )
    args = parser.parse_args()

    output = args.figures
    if output is None:
        output = (
            DEFAULT_FIGURES / args.results.name
            / f"{args.support_seconds:d}s"
        )

    predictions, scores, counts, cohorts = load_results(
        args.results, args.support_seconds
    )
    description = selection_description(counts, args.support_seconds)
    plot_sample_histograms(predictions, cohorts, description, output)
    plot_zhang_paxton_density(predictions, cohorts, description, output)
    plot_residual_histograms(predictions, cohorts, description, output)
    plot_model_scores(scores, cohorts, description, output)
    print(f"Saved four figure pairs to {output}")


if __name__ == "__main__":
    main()
