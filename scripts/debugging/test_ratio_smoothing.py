"""Test whether simple spatial smoothing repairs the WIC/SI13 ratio.

This is a sensitivity test on the common-grid, proton-corrected fields. It is
not a reconstruction of Meurant et al. (2003)'s detector-space registration.
"""

#%% Imports

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from icreader import load as icload
import numpy as np
from scipy.ndimage import gaussian_filter

from icphysics.image import wic_to_s13


#%% Paths and test widths

REPOSITORY = Path(__file__).resolve().parents[2]
INPUT = REPOSITORY / "example_data" / "precipitation" / "or_0968.nc"
OUTPUT = REPOSITORY / "figures" / "debugging" / "ratio_smoothing"

SIGMAS = [0.0, 0.5, 1.0, 1.5, 2.0]
RATIO_MIN = float(np.min(wic_to_s13))
RATIO_MAX = float(np.max(wic_to_s13))


#%% NaN-aware smoothing

def smooth_image(values, sigma):
    """Gaussian-smooth each frame without spreading NaNs into valid data."""

    values = np.asarray(values, dtype=float)
    if sigma == 0:
        return values.copy()

    output = np.full(values.shape, np.nan)

    for frame in range(values.shape[0]):
        valid = np.isfinite(values[frame])
        numerator = gaussian_filter(
            np.where(valid, values[frame], 0.0),
            sigma=sigma,
            mode="constant",
            cval=0.0,
        )
        denominator = gaussian_filter(
            valid.astype(float),
            sigma=sigma,
            mode="constant",
            cval=0.0,
        )

        enough_support = denominator >= 0.5
        output[frame, enough_support] = (
            numerator[enough_support] / denominator[enough_support]
        )

    return output


def ratio_summary(wic, si13, label, sigma):
    """Calculate the ratio statistics used by the active retrieval."""

    supported = (
        np.isfinite(wic)
        & np.isfinite(si13)
        & (wic > 0)
        & (si13 >= 3)
    )
    ratio = np.full(wic.shape, np.nan)
    ratio[supported] = wic[supported] / si13[supported]
    values = ratio[supported]

    inside = (values >= RATIO_MIN) & (values <= RATIO_MAX)
    above = values > RATIO_MAX

    summary = {
        "case": label,
        "sigma_grid_cells": sigma,
        "supported_pixels": values.size,
        "median_ratio": float(np.median(values)),
        "ratio_p90": float(np.percentile(values, 90)),
        "ratio_p95": float(np.percentile(values, 95)),
        "inside_table_percent": float(100 * np.mean(inside)),
        "above_table_percent": float(100 * np.mean(above)),
    }
    return ratio, summary


#%% Run both smoothing interpretations

def calculate_cases(wic, si13):
    """Test smoothing both cameras and adding smoothing only to WIC."""

    ratio, summary = ratio_summary(wic, si13, "no smoothing", 0.0)
    cases = [(ratio, summary)]

    for sigma in SIGMAS[1:]:
        wic_smooth = smooth_image(wic, sigma)
        si13_smooth = smooth_image(si13, sigma)
        ratio, summary = ratio_summary(
            wic_smooth,
            si13_smooth,
            "smooth both",
            sigma,
        )
        cases.append((ratio, summary))

    for sigma in SIGMAS[1:]:
        wic_smooth = smooth_image(wic, sigma)
        ratio, summary = ratio_summary(
            wic_smooth,
            si13,
            "smooth WIC only",
            sigma,
        )
        cases.append((ratio, summary))

    return cases


#%% Output

def save_summary(cases):
    """Write one compact CSV and summary figure."""

    summaries = [summary for _, summary in cases]
    csv_file = OUTPUT / "ratio_smoothing_summary.csv"

    with csv_file.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=summaries[0])
        writer.writeheader()
        writer.writerows(summaries)

    figure, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    baseline = summaries[0]
    metrics = ["median_ratio", "ratio_p90", "above_table_percent"]
    for axis, metric in zip(axes, metrics):
        axis.plot(
            baseline["sigma_grid_cells"],
            baseline[metric],
            marker="D",
            color="black",
            linestyle="none",
            label="no smoothing",
        )

    colours = {"smooth both": "tab:blue", "smooth WIC only": "tab:orange"}

    for label in colours:
        selected = [row for row in summaries if row["case"] == label]
        sigma = [row["sigma_grid_cells"] for row in selected]
        axes[0].plot(sigma, [row["median_ratio"] for row in selected], "o-", label=label, color=colours[label])
        axes[1].plot(sigma, [row["ratio_p90"] for row in selected], "o-", label=label, color=colours[label])
        axes[2].plot(sigma, [row["above_table_percent"] for row in selected], "o-", label=label, color=colours[label])

    axes[0].axhline(RATIO_MAX, color="black", linestyle="--", linewidth=1)
    axes[1].axhline(RATIO_MAX, color="black", linestyle="--", linewidth=1)
    axes[0].set_ylabel("median WIC/SI13")
    axes[1].set_ylabel("90th percentile WIC/SI13")
    axes[2].set_ylabel("pixels above Frey table [%]")

    for axis in axes:
        axis.set_xlabel("Gaussian sigma [200-km grid cells]")
        axis.grid(alpha=0.2)

    axes[0].legend(frameon=False)
    figure.suptitle("Orbit 0968: smoothing does not reproduce detector-space co-registration")
    figure.savefig(OUTPUT / "ratio_smoothing_summary.png", dpi=180)
    plt.close(figure)


def save_example_maps(cases, times):
    """Show the unsmoothed result and two representative smoothing cases."""

    choices = [
        cases[0],
        next(case for case in cases if case[1]["case"] == "smooth both" and case[1]["sigma_grid_cells"] == 1.0),
        next(case for case in cases if case[1]["case"] == "smooth WIC only" and case[1]["sigma_grid_cells"] == 1.0),
    ]
    frame = 1

    figure, axes = plt.subplots(1, 3, figsize=(11, 3.8), constrained_layout=True)
    for axis, (ratio, summary) in zip(axes, choices):
        image = axis.imshow(
            ratio[frame],
            origin="lower",
            cmap="viridis",
            vmin=RATIO_MIN,
            vmax=300,
        )
        title = "no smoothing" if summary["sigma_grid_cells"] == 0 else (
            f"{summary['case']}\n$\\sigma$ = {summary['sigma_grid_cells']:.1f} cells"
        )
        axis.set_title(title)
        axis.set_xticks([])
        axis.set_yticks([])

    figure.colorbar(image, ax=axes, label="WIC/SI13", shrink=0.8)
    figure.suptitle(f"Orbit 0968, {np.datetime_as_string(times[frame], unit='s')} UTC")
    figure.savefig(OUTPUT / "ratio_smoothing_maps.png", dpi=180)
    plt.close(figure)


#%% Main calculation

def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)

    product = icload(INPUT)
    wic = np.asarray(product.wic_corrected, dtype=float)
    si13 = np.asarray(product.si13_corrected, dtype=float)
    times = np.asarray(product.time, dtype="datetime64[ns]")

    cases = calculate_cases(wic, si13)
    save_summary(cases)
    save_example_maps(cases, times)

    for _, summary in cases:
        print(summary)


if __name__ == "__main__":
    main()
