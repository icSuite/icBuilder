"""Summarize how DZA masks change modular IMAGE products.

The script only reads existing binned, precipitation, and conductance files.
It does not apply a line-of-sight correction or alter any scientific product.
"""

#%% Imports

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from icreader import load as icload
import numpy as np
import pandas as pd
from tqdm import tqdm


#%% Settings

DEFAULT_BASE = Path("/home/bing/dtu_server/IMAGE_FUV")
DEFAULT_OUTPUT = Path(__file__).resolve().parents[2] / "figures" / "debugging" / "dza_sensitivity"

THRESHOLDS = np.array([30, 40, 50, 60, 65, 70, 75], dtype=float)
BAND_EDGES = np.array([0, 30, 40, 50, 60, 65, 70, 75], dtype=float)

FREY_RATIO_MIN = 446 / 12.8
FREY_RATIO_MAX = 101 / 0.74


#%% File handling

def orbit_number(filename):
    return int(filename.stem[-4:])


def orbit_files(folder):
    return {
        orbit_number(filename): filename
        for filename in sorted(folder.glob("or_*.nc"))
    }


def read_precipitation(filename, binned_wic_file, method):
    """Read the fields needed for this diagnostic and attach WIC DZA."""

    precipitation = icload(filename)
    data = {
        "time": np.asarray(precipitation.time, dtype=object),
        "wic_source_index": np.asarray(
            precipitation.wic_source_index, dtype=int
        ),
        "E0": np.asarray(precipitation.E0, dtype=float),
        "Fe": np.asarray(precipitation.Fe, dtype=float),
    }
    if method == "image_ratio":
        data["wic"] = np.asarray(precipitation.wic_corrected, dtype=float)
        data["si13"] = np.asarray(
            precipitation.si13_corrected, dtype=float
        )
        data["R"] = np.asarray(precipitation.R, dtype=float)
    else:
        shape = data["E0"].shape
        data["wic"] = np.full(shape, np.nan)
        data["si13"] = np.full(shape, np.nan)
        data["R"] = np.full(shape, np.nan)

    binned_wic = icload(binned_wic_file)
    source_dza = np.asarray(binned_wic.dza, dtype=float)

    data["dza"] = source_dza[data["wic_source_index"]]
    return data


def calculate_conductance(E0, Fe):
    """Apply the same Robinson equations used to create Product 3."""

    with np.errstate(invalid="ignore"):
        P = 40 * E0 / (16 + E0**2) * np.sqrt(Fe)
        H = 18 * E0**1.85 / (16 + E0**2) * np.sqrt(Fe)
    return P, H


#%% Orbit summaries

def safe_median(values):
    values = values[np.isfinite(values)]
    return float(np.median(values)) if values.size else np.nan


def summarize_selection(method, orbit, upper_dza, selected, data, P, H):
    """Summarize one cumulative threshold or one fixed DZA band."""

    product = (
        selected
        & np.isfinite(data["E0"])
        & np.isfinite(data["Fe"])
    )
    positive = product & (data["Fe"] > 0)

    ratio_support = (
        selected
        & np.isfinite(data["wic"])
        & np.isfinite(data["si13"])
        & (data["wic"] > 0)
        & (data["si13"] >= 3)
    )
    ratio = data["R"][ratio_support]
    inside = (ratio >= FREY_RATIO_MIN) & (ratio <= FREY_RATIO_MAX)

    row = {
        "method": method,
        "orbit": orbit,
        "upper_dza": upper_dza,
        "frames": data["time"].size,
        "product_pixels": int(product.sum()),
        "positive_pixels": int(positive.sum()),
        "ratio_pixels": int(ratio_support.sum()),
        "inside_table_pixels": int(inside.sum()),
        "above_table_pixels": int((ratio > FREY_RATIO_MAX).sum()),
        "below_table_pixels": int((ratio < FREY_RATIO_MIN).sum()),
        "median_ratio": safe_median(ratio),
        "median_E0": safe_median(data["E0"][positive]),
        "median_Fe": safe_median(data["Fe"][positive]),
        "median_P": safe_median(P[positive]),
        "median_H": safe_median(H[positive]),
        "sum_Fe": float(np.nansum(data["Fe"][product])),
        "sum_P": float(np.nansum(P[product])),
        "sum_H": float(np.nansum(H[product])),
        "E0_0p2_pixels": int((product & np.isclose(data["E0"], 0.2)).sum()),
        "E0_1_pixels": int((product & np.isclose(data["E0"], 1.0)).sum()),
        "E0_25_pixels": int((product & np.isclose(data["E0"], 25.0)).sum()),
    }
    return row


def summarize_orbit(method, orbit, precipitation_file, wic_file):
    data = read_precipitation(precipitation_file, wic_file, method)
    P, H = calculate_conductance(data["E0"], data["Fe"])

    finite_dza = np.isfinite(data["dza"])
    rows = []
    for threshold in THRESHOLDS:
        selected = finite_dza & (data["dza"] < threshold)
        row = summarize_selection(method, orbit, threshold, selected, data, P, H)
        row["selection"] = "threshold"
        row["lower_dza"] = 0.0
        rows.append(row)

    for lower, upper in zip(BAND_EDGES[:-1], BAND_EDGES[1:]):
        selected = finite_dza & (data["dza"] >= lower) & (data["dza"] < upper)
        row = summarize_selection(method, orbit, upper, selected, data, P, H)
        row["selection"] = "band"
        row["lower_dza"] = lower
        rows.append(row)

    return rows


#%% Mission summaries

def percentage(numerator, denominator):
    return 100 * numerator / denominator if denominator else np.nan


def aggregate_rows(rows, selection):
    """Combine pixel counts and summarize variation among orbit medians."""

    rows = rows.loc[rows.selection == selection]
    group_columns = ["method", "lower_dza", "upper_dza"]
    output = []

    for keys, group in rows.groupby(group_columns, sort=True):
        method, lower, upper = keys
        item = {
            "method": method,
            "lower_dza": lower,
            "upper_dza": upper,
            "orbits": group.orbit.nunique(),
            "frames": int(group.groupby("orbit").frames.first().sum()),
        }
        for name in (
            "product_pixels", "positive_pixels", "ratio_pixels",
            "inside_table_pixels", "above_table_pixels", "below_table_pixels",
            "E0_0p2_pixels", "E0_1_pixels", "E0_25_pixels",
        ):
            item[name] = int(group[name].sum())

        for name in ("sum_Fe", "sum_P", "sum_H"):
            item[name] = float(group[name].sum())

        for name in ("median_ratio", "median_E0", "median_Fe", "median_P", "median_H"):
            values = group[name].dropna()
            item[name] = values.median() if len(values) else np.nan
            item[f"{name}_q25"] = values.quantile(0.25) if len(values) else np.nan
            item[f"{name}_q75"] = values.quantile(0.75) if len(values) else np.nan

        item["inside_table_percent"] = percentage(
            item["inside_table_pixels"], item["ratio_pixels"]
        )
        item["above_table_percent"] = percentage(
            item["above_table_pixels"], item["ratio_pixels"]
        )
        item["below_table_percent"] = percentage(
            item["below_table_pixels"], item["ratio_pixels"]
        )
        output.append(item)

    output = pd.DataFrame(output)

    if selection == "threshold":
        reference = output.loc[output.upper_dza == THRESHOLDS[-1]].set_index("method")
        output["product_retained_percent"] = [
            percentage(row.product_pixels, reference.loc[row.method, "product_pixels"])
            for row in output.itertuples()
        ]
        output["positive_retained_percent"] = [
            percentage(row.positive_pixels, reference.loc[row.method, "positive_pixels"])
            for row in output.itertuples()
        ]
        output["ratio_retained_percent"] = [
            percentage(row.ratio_pixels, reference.loc[row.method, "ratio_pixels"])
            for row in output.itertuples()
        ]
        for name in ("Fe", "P", "H"):
            output[f"{name}_sum_retained_percent"] = [
                percentage(getattr(row, f"sum_{name}"), reference.loc[row.method, f"sum_{name}"])
                for row in output.itertuples()
            ]

    return output


#%% Figures

def plot_thresholds(summary, output):
    figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)

    for method, values in summary.groupby("method"):
        label = "Zhang–Paxton" if method == "zhang_paxton" else "Image ratio"
        axes[0, 0].plot(values.upper_dza, values.positive_retained_percent, "o-", label=label)
        axes[0, 0].plot(
            values.upper_dza, values.Fe_sum_retained_percent, ".--",
            label=f"{label} summed $F_e$",
        )
        if method == "image_ratio":
            axes[0, 0].plot(values.upper_dza, values.ratio_retained_percent, "s--", label="IR ratio support")

        axes[1, 0].plot(values.upper_dza, values.median_Fe, "o-", label=label)

    ir = summary.loc[summary.method == "image_ratio"]
    if len(ir):
        axes[0, 1].plot(ir.upper_dza, ir.median_ratio, "o-")
        axes[0, 1].fill_between(
            ir.upper_dza, ir.median_ratio_q25, ir.median_ratio_q75, alpha=0.2
        )
        axes[0, 1].axhspan(FREY_RATIO_MIN, FREY_RATIO_MAX, color="tab:green", alpha=0.12)
        axes[1, 1].plot(ir.upper_dza, ir.inside_table_percent, "o-", label="Inside table")
        axes[1, 1].plot(ir.upper_dza, ir.above_table_percent, "o-", label="Above table")

    axes[0, 0].set_ylabel("Retained pixels [% of DZA < 75°]")
    axes[0, 1].set_ylabel("Median orbit WIC/SI13")
    axes[1, 0].set_ylabel("Median orbit $F_e$ [mW m$^{-2}$]")
    axes[1, 1].set_ylabel("IR ratio pixels [%]")
    for axis in axes.flat:
        axis.set_xlabel("Upper DZA threshold [degrees]")
        axis.grid(alpha=0.25)
    axes[0, 0].legend(frameon=False)
    axes[1, 0].legend(frameon=False)
    axes[1, 1].legend(frameon=False)

    figure.suptitle("Sensitivity to cumulative DZA masking")
    figure.savefig(output / "dza_threshold_sensitivity.png", dpi=180)
    plt.close(figure)


def plot_bands(summary, output):
    figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)

    for method, values in summary.groupby("method"):
        x = (values.lower_dza + values.upper_dza) / 2
        label = "Zhang–Paxton" if method == "zhang_paxton" else "Image ratio"
        axes[1, 0].plot(x, values.median_Fe, "o-", label=label)
        axes[1, 1].plot(x, values.positive_pixels, "o-", label=label)

    ir = summary.loc[summary.method == "image_ratio"]
    if len(ir):
        x = (ir.lower_dza + ir.upper_dza) / 2
        axes[0, 0].plot(x, ir.median_ratio, "o-")
        axes[0, 0].fill_between(x, ir.median_ratio_q25, ir.median_ratio_q75, alpha=0.2)
        axes[0, 0].axhspan(FREY_RATIO_MIN, FREY_RATIO_MAX, color="tab:green", alpha=0.12)
        axes[0, 1].plot(x, ir.inside_table_percent, "o-", label="Inside table")
        axes[0, 1].plot(x, ir.above_table_percent, "o-", label="Above table")

    axes[0, 0].set_ylabel("Median orbit WIC/SI13")
    axes[0, 1].set_ylabel("IR ratio pixels [%]")
    axes[1, 0].set_ylabel("Median orbit $F_e$ [mW m$^{-2}$]")
    axes[1, 1].set_ylabel("Positive-flux pixels")
    for axis in axes.flat:
        axis.set_xlabel("DZA band centre [degrees]")
        axis.grid(alpha=0.25)
    axes[0, 1].legend(frameon=False)
    axes[1, 0].legend(frameon=False)
    axes[1, 1].legend(frameon=False)

    figure.suptitle("Product behavior in fixed DZA bands")
    figure.savefig(output / "dza_band_sensitivity.png", dpi=180)
    plt.close(figure)


def plot_conductance(summary, output):
    available = summary.dropna(subset=["median_P", "median_H"], how="all")
    if not len(available):
        return

    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for method, values in available.groupby("method"):
        label = "Zhang–Paxton" if method == "zhang_paxton" else "Image ratio"
        axes[0].plot(values.upper_dza, values.median_P, "o-", label=label)
        axes[1].plot(values.upper_dza, values.median_H, "o-", label=label)

    axes[0].set_ylabel("Median orbit Pedersen [S]")
    axes[1].set_ylabel("Median orbit Hall [S]")
    for axis in axes:
        axis.set_xlabel("Upper DZA threshold [degrees]")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
    figure.suptitle("Conductance sensitivity to cumulative DZA masking")
    figure.savefig(output / "dza_conductance_sensitivity.png", dpi=180)
    plt.close(figure)


#%% Command line

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--limit", type=int, help="Process at most this many orbits per method")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    base = args.base.expanduser()
    output = args.output.expanduser()
    output.mkdir(parents=True, exist_ok=True)

    wic_files = orbit_files(base / "=binned" / "wic")
    configurations = {
        "zhang_paxton": orbit_files(base / "=precipitation_ZP_P2"),
        "image_ratio": orbit_files(base / "=precipitation_IR_P2"),
    }

    all_rows = []
    for method, precipitation_files in configurations.items():
        orbits = sorted(set(precipitation_files) & set(wic_files))
        if args.limit:
            orbits = orbits[:args.limit]
        print(
            f"{method}: {len(precipitation_files)} precipitation files, "
            f"{len(orbits)} usable orbits"
        )

        for orbit in tqdm(orbits, desc=method):
            all_rows.extend(
                summarize_orbit(
                    method,
                    orbit,
                    precipitation_files[orbit],
                    wic_files[orbit],
                )
            )

    rows = pd.DataFrame(all_rows)
    rows.to_csv(output / "dza_orbit_summary.csv", index=False)

    threshold_summary = aggregate_rows(rows, "threshold")
    band_summary = aggregate_rows(rows, "band")
    threshold_summary.to_csv(output / "dza_threshold_summary.csv", index=False)
    band_summary.to_csv(output / "dza_band_summary.csv", index=False)

    plot_thresholds(threshold_summary, output)
    plot_bands(band_summary, output)
    plot_conductance(threshold_summary, output)

    print(f"Saved DZA sensitivity results in {output}")


if __name__ == "__main__":
    main()
