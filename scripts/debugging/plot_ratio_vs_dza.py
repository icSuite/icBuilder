"""Plot the WIC/SI13 ratio against WIC detector zenith angle.

By default this uses the same image-ratio orbit snapshot recorded in
dza_orbit_summary.csv. The external IMAGE products are only read.
"""

#%% Imports

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from icreader import load as icload
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from tqdm import tqdm


#%% Settings

DEFAULT_BASE = Path("/home/bing/dtu_server/IMAGE_FUV")
DEFAULT_SUMMARY = (
    Path(__file__).resolve().parents[2]
    / "figures" / "debugging" / "dza_sensitivity" / "dza_orbit_summary.csv"
)
DEFAULT_OUTPUT = DEFAULT_SUMMARY.parent / "ratio_vs_dza.png"

FREY_RATIO_MAX = 101 / 0.74
TAIL_LEVELS = [FREY_RATIO_MAX, 150, 200, 300]


#%% Read the ratio and corresponding WIC DZA

def orbit_number(filename):
    return int(filename.stem[-4:])


def read_ratio_pixels(precipitation_file, wic_file):
    precipitation = icload(precipitation_file)
    ratio = np.asarray(precipitation.R, dtype=float)
    wic = np.asarray(precipitation.wic_corrected, dtype=float)
    si13 = np.asarray(precipitation.si13_corrected, dtype=float)
    source_index = np.asarray(
        precipitation.wic_source_index, dtype=int
    )

    with Dataset(wic_file) as nc:
        dza = np.asarray(nc.variables["dza"][:], dtype=float)[source_index]

    keep = (
        np.isfinite(ratio)
        & np.isfinite(dza)
        & np.isfinite(wic)
        & np.isfinite(si13)
        & (wic > 0)
        & (si13 >= 3)
    )
    return dza[keep].astype(np.float32), ratio[keep].astype(np.float32)


#%% Conditional statistics within DZA bins

def binned_statistics(dza, ratio, dza_edges):
    centres = (dza_edges[:-1] + dza_edges[1:]) / 2
    quantiles = np.full((len(centres), 5), np.nan)
    tail_fractions = np.full((len(centres), len(TAIL_LEVELS)), np.nan)

    for i, (lower, upper) in enumerate(zip(dza_edges[:-1], dza_edges[1:])):
        values = ratio[(dza >= lower) & (dza < upper)]
        if values.size == 0:
            continue

        quantiles[i] = np.quantile(values, [0.25, 0.50, 0.75, 0.90, 0.95])
        tail_fractions[i] = [100 * np.mean(values > level) for level in TAIL_LEVELS]

    return centres, quantiles, tail_fractions


#%% Plot

def make_figure(dza, ratio, orbit_count, output):
    dza_edges = np.linspace(0, 75, 31)
    ratio_edges = np.linspace(0, 500, 251)
    centres, quantiles, tail_fractions = binned_statistics(dza, ratio, dza_edges)

    counts, _, _ = np.histogram2d(dza, ratio, bins=(dza_edges, ratio_edges))
    counts[counts == 0] = np.nan

    figure, axes = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1]}, constrained_layout=True,
    )

    density = axes[0].pcolormesh(
        dza_edges, ratio_edges, counts.T,
        cmap="viridis", norm=LogNorm(), shading="auto",
    )
    figure.colorbar(density, ax=axes[0], label="Pixels per 2.5° × 2 ratio bin")

    axes[0].fill_between(
        centres, quantiles[:, 0], quantiles[:, 2],
        color="white", alpha=0.20, label="25th–75th percentile",
    )
    axes[0].plot(centres, quantiles[:, 1], color="white", lw=2, label="Median")
    axes[0].plot(centres, quantiles[:, 3], color="white", ls="--", label="90th percentile")
    axes[0].plot(centres, quantiles[:, 4], color="white", ls=":", label="95th percentile")
    axes[0].axhline(FREY_RATIO_MAX, color="tab:red", lw=1.5, label="Frey table maximum")
    axes[0].set_ylim(0, 500)
    axes[0].set_ylabel("Corrected WIC/SI13 ratio")
    axes[0].legend(loc="upper right", framealpha=0.9)

    labels = ["Above Frey maximum", "Above 150", "Above 200", "Above 300"]
    for values, label in zip(tail_fractions.T, labels):
        axes[1].plot(centres, values, "o-", ms=3, label=label)

    axes[1].set_xlim(0, 75)
    axes[1].set_xlabel("WIC detector zenith angle [degrees]")
    axes[1].set_ylabel("Pixels in DZA bin [%]")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False, ncol=2)

    figure.suptitle(
        f"WIC/SI13 ratio versus WIC DZA\n"
        f"{orbit_count} image-ratio orbits, {ratio.size:,} supported pixels"
    )
    figure.savefig(output, dpi=200)
    plt.close(figure)


#%% Command line

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    summary = pd.read_csv(args.summary)
    orbits = sorted(summary.loc[summary.method == "image_ratio", "orbit"].unique())

    precipitation_files = {
        orbit_number(filename): filename
        for filename in (args.base / "=precipitation_IR_P2").glob("or_*.nc")
    }
    wic_files = {
        orbit_number(filename): filename
        for filename in (args.base / "=binned" / "wic").glob("or_*.nc")
    }

    all_dza = []
    all_ratio = []
    for orbit in tqdm(orbits, desc="Read image-ratio pixels"):
        dza, ratio = read_ratio_pixels(precipitation_files[orbit], wic_files[orbit])
        all_dza.append(dza)
        all_ratio.append(ratio)

    all_dza = np.concatenate(all_dza)
    all_ratio = np.concatenate(all_ratio)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    make_figure(all_dza, all_ratio, len(orbits), args.output)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
