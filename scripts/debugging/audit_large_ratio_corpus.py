"""Audit large WIC/SI13 ratios across the IMAGE ratio precipitation corpus.

The script reads one orbit at a time. Counts are exact, while the log-log
density uses a bounded, orbit-balanced random sample. The counted pixel-frames
are correlated; no statistical independence is assumed.
The WIC >= 50 and SI13 >= 3 selection is a modern diagnostic signal guard,
not a historical Coumans validity rule.
"""

#%% Imports and settings

import csv
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from icreader import load as icload
from scipy.stats import spearmanr


INPUT = Path("/home/bing/Dropbox/work/data/IMAGE_FUV/precipitation_IR_P2_weighted")
OUTPUT = Path("figures/debugging/large_ratio_corpus")

RATIO_THRESHOLDS = (120.0, 136.48648648648648, 150.0)
SI13_CUTOFFS = (0.0, 3.0, 5.0, 10.0, 20.0)
SIGNAL_GUARD_WIC = 50.0
SAMPLE_PER_ORBIT = 300
RANDOM_SEED = 20260825


#%% Small helpers

def orbit_number(path):
    return int(path.stem.split("_")[-1])


def fraction_above(ratio, threshold):
    if ratio.size == 0:
        return np.nan
    return float(np.mean(ratio > threshold))


def percentile(ratio, value):
    if ratio.size == 0:
        return np.nan
    return float(np.percentile(ratio, value))


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def positive_weight_summary(weight):
    positive = np.isfinite(weight) & (weight > 0)
    coverage = float(np.mean(positive))
    mean = float(np.mean(weight[positive])) if positive.any() else np.nan
    return coverage, mean


#%% Stream through the corpus

def audit_orbits(files):
    rng = np.random.default_rng(RANDOM_SEED)
    stages = ("pre_proton", "post_proton")
    counts = {
        (stage, cutoff, threshold): [0, 0]
        for stage in stages
        for cutoff in SI13_CUTOFFS
        for threshold in RATIO_THRESHOLDS
    }
    orbit_rows = []
    sample_wic = []
    sample_si13 = []

    for number, path in enumerate(files, start=1):
        data = icload(path)
        if data.method != "image_ratio":
            raise ValueError(f"{path.name} is not an image-ratio product")

        time = np.asarray(data.time, dtype="datetime64[ns]")
        kp = np.asarray(data.kp, dtype=float)
        wic_weight = np.asarray(data.wic_weight, dtype=float)
        si13_weight = np.asarray(data.si13_weight, dtype=float)
        fields = {
            "pre_proton": (
                np.asarray(data.wic, dtype=float),
                np.asarray(data.si13, dtype=float),
            ),
            "post_proton": (
                np.asarray(data.wic_corrected, dtype=float),
                np.asarray(data.si13_corrected, dtype=float),
            ),
        }
        stored_ratio = np.asarray(data.R, dtype=float)

        row = {
            "orbit": orbit_number(path),
            "time": np.datetime_as_string(time[len(time) // 2], unit="s"),
            "frames": time.size,
            "kp_median": float(np.nanmedian(kp)),
        }

        row["wic_weight_coverage"], row["wic_weight_mean_positive"] = (
            positive_weight_summary(wic_weight)
        )
        row["si13_weight_coverage"], row["si13_weight_mean_positive"] = (
            positive_weight_summary(si13_weight)
        )

        for stage, (wic, si13) in fields.items():
            positive = np.isfinite(wic) & np.isfinite(si13) & (wic > 0) & (si13 > 0)
            ratio = wic[positive] / si13[positive]
            signal_guard = positive & (wic >= SIGNAL_GUARD_WIC) & (si13 >= 3)
            guarded_ratio = wic[signal_guard] / si13[signal_guard]

            row[f"{stage}_positive_pixels"] = int(ratio.size)
            row[f"{stage}_signal_guard_pixels"] = int(guarded_ratio.size)
            row[f"{stage}_signal_guard_fraction"] = float(signal_guard.sum() / positive.sum()) if positive.any() else np.nan
            row[f"{stage}_ratio_median_signal_guard"] = percentile(guarded_ratio, 50)
            row[f"{stage}_ratio_p95_signal_guard"] = percentile(guarded_ratio, 95)
            row[f"{stage}_ratio_p99_signal_guard"] = percentile(guarded_ratio, 99)

            for threshold in RATIO_THRESHOLDS:
                label = f"R_gt_{threshold:g}"
                row[f"{stage}_{label}_all_positive"] = fraction_above(ratio, threshold)
                row[f"{stage}_{label}_signal_guard"] = fraction_above(guarded_ratio, threshold)

            for cutoff in SI13_CUTOFFS:
                support = positive
                if cutoff > 0:
                    support &= (wic >= SIGNAL_GUARD_WIC) & (si13 >= cutoff)
                selected_ratio = wic[support] / si13[support]
                for threshold in RATIO_THRESHOLDS:
                    total = selected_ratio.size
                    above = int(np.sum(selected_ratio > threshold))
                    counts[(stage, cutoff, threshold)][0] += total
                    counts[(stage, cutoff, threshold)][1] += above

            if stage == "post_proton":
                indices = np.flatnonzero(signal_guard.ravel())
                if indices.size > SAMPLE_PER_ORBIT:
                    indices = rng.choice(indices, SAMPLE_PER_ORBIT, replace=False)
                sample_wic.append(wic.ravel()[indices])
                sample_si13.append(si13.ravel()[indices])

                calculated = np.full(wic.shape, np.nan)
                calculated[positive] = wic[positive] / si13[positive]
                compare = positive & np.isfinite(stored_ratio)
                row["stored_R_maximum_absolute_difference"] = (
                    float(np.max(np.abs(calculated[compare] - stored_ratio[compare])))
                    if compare.any() else np.nan
                )

        orbit_rows.append(row)
        if number % 100 == 0 or number == len(files):
            print(f"Read {number}/{len(files)} orbits")

    global_rows = []
    for (stage, cutoff, threshold), (total, above) in counts.items():
        support = "all positive" if cutoff == 0 else f"WIC >= 50, SI13 >= {cutoff:g}"
        global_rows.append({
            "stage": stage,
            "support": support,
            "si13_cutoff": cutoff,
            "ratio_threshold": threshold,
            "pixels": total,
            "pixels_above_threshold": above,
            "fraction_above_threshold": above / total if total else np.nan,
        })

    return orbit_rows, global_rows, np.concatenate(sample_wic), np.concatenate(sample_si13)


#%% Orbit-level associations and calendar-month summary

def make_association_rows(orbit_rows):
    targets = (
        "pre_proton_R_gt_136.486_all_positive",
        "pre_proton_R_gt_136.486_signal_guard",
    )
    predictors = (
        "kp_median",
        "frames",
        "pre_proton_signal_guard_fraction",
        "wic_weight_coverage",
        "wic_weight_mean_positive",
        "si13_weight_coverage",
        "si13_weight_mean_positive",
    )

    rows = []
    for target in targets:
        for predictor in predictors:
            x = np.array([row[predictor] for row in orbit_rows], dtype=float)
            y = np.array([row[target] for row in orbit_rows], dtype=float)
            keep = np.isfinite(x) & np.isfinite(y)
            rho, pvalue = spearmanr(x[keep], y[keep])
            rows.append({
                "target": target,
                "predictor": predictor,
                "orbits": int(keep.sum()),
                "spearman_rho": float(rho),
                "nominal_pvalue": float(pvalue),
            })
    return rows


def make_month_rows(orbit_rows):
    rows = []
    for month in range(1, 13):
        selected = [row for row in orbit_rows if int(row["time"][5:7]) == month]
        rows.append({
            "month": month,
            "orbits": len(selected),
            "pre_proton_fraction_all_positive": float(np.nanmedian([
                row["pre_proton_R_gt_136.486_all_positive"] for row in selected
            ])),
            "pre_proton_fraction_signal_guard": float(np.nanmedian([
                row["pre_proton_R_gt_136.486_signal_guard"] for row in selected
            ])),
            "post_proton_fraction_all_positive": float(np.nanmedian([
                row["post_proton_R_gt_136.486_all_positive"] for row in selected
            ])),
            "post_proton_fraction_signal_guard": float(np.nanmedian([
                row["post_proton_R_gt_136.486_signal_guard"] for row in selected
            ])),
        })
    return rows


#%% One compact summary figure

def make_figure(orbit_rows, global_rows, sample_wic, sample_si13):
    time = np.array([row["time"] for row in orbit_rows], dtype="datetime64[s]")
    pre = np.array([row["pre_proton_R_gt_136.486_all_positive"] for row in orbit_rows])
    post = np.array([row["post_proton_R_gt_136.486_all_positive"] for row in orbit_rows])
    signal_guard = np.array([row["post_proton_R_gt_136.486_signal_guard"] for row in orbit_rows])

    figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    axes[0, 0].plot(time, post, ".", ms=2, alpha=0.5, label="all positive")
    axes[0, 0].plot(time, signal_guard, ".", ms=2, alpha=0.5,
                    label="WIC >= 50, SI13 >= 3")

    months = time.astype("datetime64[M]")
    month_values = np.unique(months)
    month_time = month_values + np.timedelta64(14, "D")
    post_monthly = [np.nanmedian(post[months == month]) for month in month_values]
    guard_monthly = [np.nanmedian(signal_guard[months == month]) for month in month_values]
    axes[0, 0].plot(month_time, post_monthly, color="tab:blue", lw=1.5)
    axes[0, 0].plot(month_time, guard_monthly, color="tab:orange", lw=1.5)

    axes[0, 0].set(ylabel="Fraction with R > 136.49", title="Post-proton ratio by orbit")
    axes[0, 0].xaxis.set_major_locator(mdates.YearLocator())
    axes[0, 0].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    axes[0, 0].legend(frameon=False)

    axes[0, 1].scatter(pre, post, s=6, alpha=0.35)
    limit = np.nanmax([pre, post])
    axes[0, 1].plot([0, limit], [0, limit], "k--", lw=1)
    axes[0, 1].set(xlabel="Pre-proton fraction", ylabel="Post-proton fraction",
                   title="Per-orbit R > 136.49")

    for stage, linestyle in (("pre_proton", "--"), ("post_proton", "-")):
        for threshold, colour in zip(RATIO_THRESHOLDS, ("tab:blue", "tab:orange", "tab:red")):
            rows = [row for row in global_rows
                    if row["stage"] == stage
                    and row["ratio_threshold"] == threshold
                    and row["si13_cutoff"] > 0]
            axes[1, 0].plot([row["si13_cutoff"] for row in rows],
                            [row["fraction_above_threshold"] for row in rows],
                            marker="o", linestyle=linestyle, color=colour,
                            label=f"{stage.replace('_', ' ')}, R>{threshold:g}")
    axes[1, 0].set(xlabel="SI13 cutoff [counts]", ylabel="Fraction above ratio threshold",
                   title="Threshold sensitivity with WIC >= 50 counts")
    axes[1, 0].legend(frameon=False, fontsize=8, ncol=2)

    axes[1, 1].hexbin(sample_si13, sample_wic, xscale="log", yscale="log",
                      gridsize=65, bins="log", mincnt=1, cmap="viridis")
    x = np.geomspace(np.min(sample_si13), np.max(sample_si13), 200)
    for threshold, colour in zip(RATIO_THRESHOLDS, ("tab:blue", "tab:orange", "tab:red")):
        axes[1, 1].plot(x, threshold * x, color=colour, lw=1.5, label=f"R={threshold:g}")
    axes[1, 1].set(xlabel="Post-proton SI13 [counts]", ylabel="Post-proton WIC [counts]",
                   title=f"Orbit-balanced sample (<= {SAMPLE_PER_ORBIT}/orbit)")
    axes[1, 1].legend(frameon=False)

    for axis in axes.ravel():
        axis.grid(alpha=0.2)
    figure.suptitle(
        "IMAGE ratio corpus: correlated pixel-frames; signal guard is diagnostic"
    )
    figure.savefig(OUTPUT / "large_ratio_corpus.png", dpi=200)
    figure.savefig(OUTPUT / "large_ratio_corpus.pdf")
    plt.close(figure)


#%% Run

def main():
    files = sorted(INPUT.glob("or_*.nc"), key=orbit_number)
    if not files:
        raise FileNotFoundError(f"No orbit files found in {INPUT}")

    OUTPUT.mkdir(parents=True, exist_ok=True)
    orbit_rows, global_rows, sample_wic, sample_si13 = audit_orbits(files)
    write_csv(OUTPUT / "large_ratio_per_orbit.csv", orbit_rows)
    write_csv(OUTPUT / "large_ratio_global.csv", global_rows)
    write_csv(OUTPUT / "large_ratio_associations.csv", make_association_rows(orbit_rows))
    write_csv(OUTPUT / "large_ratio_by_calendar_month.csv", make_month_rows(orbit_rows))
    make_figure(orbit_rows, global_rows, sample_wic, sample_si13)
    print(f"Saved audit of {len(files)} orbits in {OUTPUT}")


if __name__ == "__main__":
    main()
