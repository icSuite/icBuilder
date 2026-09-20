"""Compare fixed-2-keV and Hardy proton corrections orbit-for-orbit.

The two precipitation corpora must have been generated from the same binned
images. Ratios are compared only where both products contain valid corrected
WIC and SI13 values. Pixel-frames are correlated, so the reported fractions
are descriptive rather than independent statistical samples.
"""

#%% Imports and settings

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from icreader import load as icload


DATA = Path("/home/bing/Dropbox/work/data/IMAGE_FUV")
FIXED = DATA / "precipitation_IR_P2_weighted"
HARDY = DATA / "=precipitation_IR_hardy_weighted"
OUTPUT = Path("figures/debugging/hardy_fixed_proton_comparison")

FREY_MAXIMUM = 136.48648648648648
WIC_MINIMUM = 50.0
SI13_MINIMUM = 3.0
SAMPLE_PER_ORBIT = 300
RANDOM_SEED = 20260831


#%% Read matching orbits and accumulate paired statistics

def main():
    fixed_files = {path.name: path for path in FIXED.glob("or_*.nc")}
    hardy_files = {path.name: path for path in HARDY.glob("or_*.nc")}

    if fixed_files.keys() != hardy_files.keys():
        raise ValueError("The fixed and Hardy corpora do not contain the same orbits")

    counts = {
        "common_positive": 0,
        "fixed_above_common_positive": 0,
        "hardy_above_common_positive": 0,
        "common_signal_guard": 0,
        "fixed_above_common_guard": 0,
        "hardy_above_common_guard": 0,
        "hardy_lower_common_guard": 0,
        "hardy_higher_common_guard": 0,
        "hardy_equal_common_guard": 0,
        "hardy_only_above_common_guard": 0,
        "fixed_only_above_common_guard": 0,
        "hardy_Ep_clipped_common_guard": 0,
        "fixed_positive_only": 0,
        "hardy_positive_only": 0,
    }

    rng = np.random.default_rng(RANDOM_SEED)
    orbit_rows = []
    fixed_sample = []
    hardy_sample = []
    ep_sample = []

    names = sorted(fixed_files)
    for number, name in enumerate(names, start=1):
        fixed_data = icload(fixed_files[name])
        fixed_time = np.asarray(fixed_data.time, dtype="datetime64[ns]")
        fixed_indices = np.asarray(fixed_data.wic_source_index)
        fixed_wic = np.asarray(fixed_data.wic_corrected, dtype=float)
        fixed_si13 = np.asarray(fixed_data.si13_corrected, dtype=float)

        hardy_data = icload(hardy_files[name])
        hardy_time = np.asarray(hardy_data.time, dtype="datetime64[ns]")
        hardy_indices = np.asarray(hardy_data.wic_source_index)
        hardy_wic = np.asarray(hardy_data.wic_corrected, dtype=float)
        hardy_si13 = np.asarray(hardy_data.si13_corrected, dtype=float)
        hardy_ep = np.asarray(hardy_data.Ep, dtype=float)
        hardy_ep_clipped = np.asarray(
            hardy_data.Ep_clipping_flag, dtype=bool
        )

        if not np.array_equal(fixed_time, hardy_time):
            raise ValueError(f"Time mismatch in {name}")
        if not np.array_equal(fixed_indices, hardy_indices):
            raise ValueError(f"WIC source-index mismatch in {name}")

        fixed_valid = (
            np.isfinite(fixed_wic) & np.isfinite(fixed_si13)
            & (fixed_wic > 0) & (fixed_si13 > 0)
        )
        hardy_valid = (
            np.isfinite(hardy_wic) & np.isfinite(hardy_si13)
            & (hardy_wic > 0) & (hardy_si13 > 0)
        )
        common = fixed_valid & hardy_valid

        fixed_guard = fixed_valid & (fixed_wic >= WIC_MINIMUM) & (fixed_si13 >= SI13_MINIMUM)
        hardy_guard = hardy_valid & (hardy_wic >= WIC_MINIMUM) & (hardy_si13 >= SI13_MINIMUM)
        common_guard = fixed_guard & hardy_guard

        fixed_ratio = np.full(fixed_wic.shape, np.nan)
        hardy_ratio = np.full(hardy_wic.shape, np.nan)
        fixed_ratio[common] = fixed_wic[common] / fixed_si13[common]
        hardy_ratio[common] = hardy_wic[common] / hardy_si13[common]

        fixed_above = fixed_ratio > FREY_MAXIMUM
        hardy_above = hardy_ratio > FREY_MAXIMUM

        counts["common_positive"] += int(common.sum())
        counts["fixed_above_common_positive"] += int((common & fixed_above).sum())
        counts["hardy_above_common_positive"] += int((common & hardy_above).sum())
        counts["common_signal_guard"] += int(common_guard.sum())
        counts["fixed_above_common_guard"] += int((common_guard & fixed_above).sum())
        counts["hardy_above_common_guard"] += int((common_guard & hardy_above).sum())
        counts["hardy_lower_common_guard"] += int((common_guard & (hardy_ratio < fixed_ratio)).sum())
        counts["hardy_higher_common_guard"] += int((common_guard & (hardy_ratio > fixed_ratio)).sum())
        counts["hardy_equal_common_guard"] += int((common_guard & (hardy_ratio == fixed_ratio)).sum())
        counts["hardy_only_above_common_guard"] += int(
            (common_guard & hardy_above & ~fixed_above).sum()
        )
        counts["fixed_only_above_common_guard"] += int(
            (common_guard & fixed_above & ~hardy_above).sum()
        )
        counts["hardy_Ep_clipped_common_guard"] += int((common_guard & hardy_ep_clipped).sum())
        counts["fixed_positive_only"] += int((fixed_valid & ~hardy_valid).sum())
        counts["hardy_positive_only"] += int((hardy_valid & ~fixed_valid).sum())

        guard_total = int(common_guard.sum())
        orbit_rows.append({
            "orbit": int(name[3:7]),
            "frames": fixed_time.size,
            "common_signal_guard_pixels": guard_total,
            "fixed_fraction_above_frey_maximum": (
                float((common_guard & fixed_above).sum() / guard_total) if guard_total else np.nan
            ),
            "hardy_fraction_above_frey_maximum": (
                float((common_guard & hardy_above).sum() / guard_total) if guard_total else np.nan
            ),
            "hardy_Ep_median_keV": (
                float(np.nanmedian(hardy_ep[common_guard])) if guard_total else np.nan
            ),
        })

        sample_indices = np.flatnonzero(common_guard)
        if sample_indices.size > SAMPLE_PER_ORBIT:
            sample_indices = rng.choice(sample_indices, SAMPLE_PER_ORBIT, replace=False)
        fixed_sample.append(fixed_ratio.ravel()[sample_indices])
        hardy_sample.append(hardy_ratio.ravel()[sample_indices])
        ep_sample.append(hardy_ep.ravel()[sample_indices])

        if number % 100 == 0 or number == len(names):
            print(f"Compared {number}/{len(names)} orbits")

    fixed_sample = np.concatenate(fixed_sample)
    hardy_sample = np.concatenate(hardy_sample)
    ep_sample = np.concatenate(ep_sample)

    #%% Save exact counts, paired fractions, and orbit-level results

    summary = dict(counts)
    summary["fixed_fraction_above_common_positive"] = (
        counts["fixed_above_common_positive"] / counts["common_positive"]
    )
    summary["hardy_fraction_above_common_positive"] = (
        counts["hardy_above_common_positive"] / counts["common_positive"]
    )
    summary["fixed_fraction_above_common_guard"] = (
        counts["fixed_above_common_guard"] / counts["common_signal_guard"]
    )
    summary["hardy_fraction_above_common_guard"] = (
        counts["hardy_above_common_guard"] / counts["common_signal_guard"]
    )
    summary["hardy_lower_fraction_common_guard"] = (
        counts["hardy_lower_common_guard"] / counts["common_signal_guard"]
    )
    summary["hardy_Ep_clipped_fraction_common_guard"] = (
        counts["hardy_Ep_clipped_common_guard"] / counts["common_signal_guard"]
    )
    summary["sample_median_fixed_ratio"] = float(np.median(fixed_sample))
    summary["sample_median_hardy_ratio"] = float(np.median(hardy_sample))
    summary["sample_median_hardy_to_fixed_ratio"] = float(np.median(hardy_sample / fixed_sample))
    summary["sample_median_hardy_Ep_keV"] = float(np.median(ep_sample))
    for percentile in (5, 25, 75, 95):
        summary[f"sample_hardy_to_fixed_ratio_p{percentile}"] = float(
            np.percentile(hardy_sample / fixed_sample, percentile)
        )
        summary[f"sample_hardy_Ep_p{percentile}_keV"] = float(
            np.percentile(ep_sample, percentile)
        )

    OUTPUT.mkdir(parents=True, exist_ok=True)
    with (OUTPUT / "summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("metric", "value"))
        writer.writerows(summary.items())

    with (OUTPUT / "per_orbit.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=orbit_rows[0].keys())
        writer.writeheader()
        writer.writerows(orbit_rows)

    #%% Plot the paired effect of the proton-energy change

    fixed_orbit = np.array([row["fixed_fraction_above_frey_maximum"] for row in orbit_rows])
    hardy_orbit = np.array([row["hardy_fraction_above_frey_maximum"] for row in orbit_rows])

    figure, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)

    axes[0].scatter(fixed_orbit, hardy_orbit, s=7, alpha=0.4)
    axes[0].plot((0, 1), (0, 1), "k--", lw=1)
    axes[0].set(
        xlabel="Fixed 2 keV: fraction above Frey maximum",
        ylabel="Hardy: fraction above Frey maximum",
        title="Orbit-by-orbit",
        xlim=(0, 1), ylim=(0, 1),
    )

    limits = (1, max(np.percentile(fixed_sample, 99.9), np.percentile(hardy_sample, 99.9)))
    axes[1].hexbin(fixed_sample, hardy_sample, xscale="log", yscale="log",
                   gridsize=65, bins="log", mincnt=1, cmap="viridis")
    axes[1].plot(limits, limits, "w--", lw=1)
    axes[1].set(xlim=limits, ylim=limits, xlabel="Fixed 2-keV ratio",
                ylabel="Hardy ratio", title="Orbit-balanced pixel sample")

    change = hardy_sample / fixed_sample
    axes[2].hexbin(ep_sample, change, yscale="log", gridsize=60,
                   bins="log", mincnt=1, cmap="viridis")
    axes[2].axhline(1, color="w", ls="--", lw=1)
    axes[2].set(xlabel="Hardy proton mean energy [keV]",
                ylabel="Hardy ratio / fixed-2-keV ratio",
                title="Ratio response to proton energy")

    for axis in axes:
        axis.grid(alpha=0.2)

    figure.suptitle("Image-ratio Product 2: Hardy versus fixed 2-keV proton correction")
    figure.savefig(OUTPUT / "hardy_fixed_proton_comparison.png", dpi=200)
    figure.savefig(OUTPUT / "hardy_fixed_proton_comparison.pdf")
    plt.close(figure)

    for key, value in summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
