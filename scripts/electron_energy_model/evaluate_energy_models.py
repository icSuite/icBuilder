"""Evaluate constant and Zhang--Paxton energies against DMSP SSJ.

Models are fitted and tested with complete IMAGE orbits held out.  Scores are
reported separately for all quality-controlled frame crossings and for the
manually accepted crossings.  The primary loss is absolute log10 energy error
averaged within orbit and then across orbits, so heavily sampled orbits do not
dominate the result.
"""

#%% Imports and defaults

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
import xarray as xr


REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = REPOSITORY / "data/electron_energy_model/comparison_data.nc"
DEFAULT_OUTPUT = REPOSITORY / "data/electron_energy_model/results"

RANDOM_SEED = 24082026
DEFAULT_FOLDS = 5
DEFAULT_BOOTSTRAP = 3000
BASELINE_MODEL = "constant_log_median"

MODEL_SPECS = {
    "constant_2kev": (None, "hardcoded_2kev"),
    "constant_physical": (None, "physical_constant"),
    "constant_log_median": (None, "log_constant"),
    "hardy_pointwise": ("hardy_pointwise_energy", "raw"),
    "zhang_paxton_pointwise": ("zp_pointwise_energy", "raw"),
    "zhang_paxton_area_mean": ("zp_area_mean_energy", "raw"),
    "zhang_paxton_area_median": ("zp_area_median_energy", "raw"),
    "zhang_paxton_flux_weighted": ("zp_flux_weighted_energy", "raw"),
    "zhang_paxton_scaled": ("zp_area_mean_energy", "scaled"),
    "zhang_paxton_calibrated": ("zp_area_mean_energy", "calibrated")
}


def available_models(data):
    """Include Hardy in new comparison files while accepting earlier ZP files."""

    return tuple(model for model in MODEL_SPECS if model != "hardy_pointwise" or "hardy_pointwise_energy" in data)


#%% Weighted fitting and orbit folds

def weighted_quantile(values, weights, quantile):
    """Return one weighted quantile after finite-value selection."""

    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values = values[good]
    weights = weights[good]
    if len(values) == 0:
        return np.nan

    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    index = np.searchsorted(cumulative, quantile * cumulative[-1])
    return float(values[order][index])


def equal_orbit_weights(frame):
    """Give every IMAGE orbit equal total influence during fitting."""

    count = frame.groupby("orbit")["orbit"].transform("size")
    return 1 / count.to_numpy(dtype=float)


def assign_folds(frame, requested_folds, seed=RANDOM_SEED):
    """Assign whole orbits to deterministic, approximately year-balanced folds."""

    orbit_year = frame.groupby("orbit")["year"].first()
    n_folds = min(requested_folds, len(orbit_year))
    if n_folds < 2:
        raise ValueError("at least two IMAGE orbits are required")

    random = np.random.default_rng(seed)
    fold_counts = np.zeros(n_folds, dtype=int)
    assignment = {}
    for year in sorted(orbit_year.unique()):
        orbits = orbit_year.index[orbit_year.eq(year)].to_numpy().copy()
        random.shuffle(orbits)
        for orbit in orbits:
            fold = int(np.argmin(fold_counts))
            assignment[int(orbit)] = fold
            fold_counts[fold] += 1

    return frame["orbit"].map(assignment).to_numpy(dtype=int), n_folds


def fit_model(frame, model):
    """Fit one deliberately small model in log10 energy space."""

    predictor, fit_type = MODEL_SPECS[model]
    energy = frame["dmsp_energy"].to_numpy(dtype=float)
    flux = frame["dmsp_flux"].to_numpy(dtype=float)
    y = np.log10(energy)
    weights = equal_orbit_weights(frame)

    if fit_type == "hardcoded_2kev":
        return {"intercept": np.log10(2.0), "slope": 0.0}
    if fit_type == "physical_constant":
        value = np.sum(weights * flux) / np.sum(weights * flux / energy)
        return {"intercept": np.log10(value), "slope": 0.0}
    if fit_type == "log_constant":
        return {
            "intercept": weighted_quantile(y, weights, 0.5),
            "slope": 0.0
        }

    x = np.log10(frame[predictor].to_numpy(dtype=float))
    if fit_type == "raw":
        return {"intercept": 0.0, "slope": 1.0}
    if fit_type == "scaled":
        return {
            "intercept": weighted_quantile(y - x, weights, 0.5),
            "slope": 1.0
        }
    if fit_type == "calibrated":
        if np.nanstd(x) < 1e-4:
            return {
                "intercept": weighted_quantile(y - x, weights, 0.5),
                "slope": 1.0
            }

        def loss(slope):
            intercept = weighted_quantile(y - slope * x, weights, 0.5)
            return np.sum(weights * np.abs(y - intercept - slope * x))

        result = minimize_scalar(loss, bounds=(0, 3), method="bounded")
        slope = float(result.x)
        intercept = weighted_quantile(y - slope * x, weights, 0.5)
        return {"intercept": intercept, "slope": slope}

    raise ValueError(f"unknown fit type: {fit_type}")


def predict_log_energy(frame, model, parameters):
    """Apply fitted model parameters and return log10 energy."""

    predictor, _ = MODEL_SPECS[model]
    if predictor is None:
        x = np.zeros(len(frame))
    else:
        x = np.log10(frame[predictor].to_numpy(dtype=float))
    return parameters["intercept"] + parameters["slope"] * x


#%% Held-out predictions and scores

def cross_validate(frame, requested_folds):
    """Return wide orbit-held-out predictions and fitted parameters."""

    frame = frame.copy()
    frame["fold"], n_folds = assign_folds(frame, requested_folds)
    prediction_parts = []
    parameter_rows = []

    keep = [
        "orbit", "image_frame", "image_time", "year", "satellite",
        "support_seconds", "annotation_status", "dmsp_energy", "dmsp_flux",
        "zp_pointwise_energy",
        "zp_area_mean_energy", "zp_area_median_energy",
        "zp_flux_weighted_energy", "fold"
    ]
    if "hardy_pointwise_energy" in frame:
        keep.append("hardy_pointwise_energy")
        if "n_hardy_pointwise" in frame:
            keep.append("n_hardy_pointwise")
    for fold in range(n_folds):
        train = frame.loc[frame["fold"].ne(fold)]
        test = frame.loc[frame["fold"].eq(fold)]
        if set(train["orbit"]).intersection(test["orbit"]):
            raise RuntimeError("an IMAGE orbit occurs in both train and test")

        part = test[keep].copy()
        for model in available_models(frame):
            parameters = fit_model(train, model)
            part[f"prediction_{model}"] = 10**predict_log_energy(
                test, model, parameters
            )
            parameter_rows.append({
                "fold": fold,
                "model": model,
                "predictor": MODEL_SPECS[model][0] or "none",
                "intercept": parameters["intercept"],
                "slope": parameters["slope"],
                "amplitude": 10**parameters["intercept"],
                "n_train": len(train),
                "n_test": len(test),
                "train_orbits": train["orbit"].nunique(),
                "test_orbits": test["orbit"].nunique()
            })
        prediction_parts.append(part)

    return pd.concat(prediction_parts, ignore_index=True), pd.DataFrame(parameter_rows)


def summarize(predictions, prediction_column):
    """Score log error after first averaging within IMAGE orbit."""

    residual = (
        np.log10(predictions["dmsp_energy"])
        - np.log10(predictions[prediction_column])
    )
    absolute = np.abs(residual)
    orbit_loss = pd.DataFrame({
        "orbit": predictions["orbit"], "loss": absolute
    }).groupby("orbit")["loss"].mean()

    return {
        "grouped_mean_absolute_log_error": float(orbit_loss.mean()),
        "row_mean_absolute_log_error": float(absolute.mean()),
        "median_absolute_log_error": float(np.median(absolute)),
        "multiplicative_error_factor": float(10**orbit_loss.mean()),
        "median_signed_log_residual": float(np.median(residual)),
        "n_rows": len(predictions),
        "n_orbits": int(predictions["orbit"].nunique())
    }


def bootstrap_improvement(predictions, model, draws, seed=RANDOM_SEED):
    """Bootstrap paired improvement over the log-constant by orbit and year."""

    truth = np.log10(predictions["dmsp_energy"])
    losses = pd.DataFrame({
        "year": predictions["year"],
        "orbit": predictions["orbit"],
        "model": np.abs(truth - np.log10(predictions[f"prediction_{model}"])),
        "baseline": np.abs(
            truth - np.log10(predictions[f"prediction_{BASELINE_MODEL}"])
        )
    })
    orbit_loss = losses.groupby(["year", "orbit"])[["model", "baseline"]].mean()
    difference = orbit_loss["baseline"] - orbit_loss["model"]

    random = np.random.default_rng(seed)
    bootstrap = np.empty(draws)
    years = difference.index.get_level_values("year").unique()
    for draw in range(draws):
        selected = []
        for year in years:
            values = difference.xs(year, level="year").to_numpy()
            selected.extend(random.choice(values, size=len(values), replace=True))
        bootstrap[draw] = np.mean(selected)
    return float(np.quantile(bootstrap, 0.025)), float(np.quantile(bootstrap, 0.975))


#%% Cohorts and complete evaluation

def select_common_sample(data, cohort, support_seconds, minimum_coverage,
                         minimum_flux, maximum_flux):
    """Select one common finite sample for every compared energy model."""

    choose = data["support_seconds"].eq(support_seconds)
    if cohort == "accepted":
        choose &= data["annotation_status"].eq(1)
    elif cohort != "all":
        raise ValueError(f"unknown cohort: {cohort}")

    choose &= (
        data["dmsp_coverage"].ge(minimum_coverage)
        & data["dmsp_energy"].gt(0)
        & data["dmsp_flux"].ge(minimum_flux)
        & data["dmsp_flux"].le(maximum_flux)
    )
    for model in available_models(data):
        predictor, _ = MODEL_SPECS[model]
        if predictor is not None:
            choose &= np.isfinite(data[predictor]) & data[predictor].gt(0)
    if "hardy_pointwise_energy" in data and {"n_hardy_pointwise", "n_dmsp"}.issubset(data.columns):
        # Do not score a partially supported model against a longer DMSP window.
        choose &= data["n_hardy_pointwise"].eq(data["n_dmsp"])
    return data.loc[choose].copy()


def evaluate(data, cohorts, supports, minimum_coverage, minimum_flux,
             maximum_flux, folds, bootstrap_draws):
    """Evaluate every requested cohort and temporal support."""

    all_predictions = []
    all_parameters = []
    all_scores = []
    sample_counts = []

    for cohort in cohorts:
        for support_seconds in supports:
            sample = select_common_sample(
                data, cohort, support_seconds, minimum_coverage,
                minimum_flux, maximum_flux
            )
            if sample["orbit"].nunique() < 2:
                raise ValueError(
                    f"{cohort}/{support_seconds} s has fewer than two orbits"
                )

            predictions, parameters = cross_validate(sample, folds)
            predictions["cohort"] = cohort
            parameters["cohort"] = cohort
            parameters["support_seconds"] = support_seconds
            all_predictions.append(predictions)
            all_parameters.append(parameters)

            sample_counts.append({
                "cohort": cohort,
                "support_seconds": support_seconds,
                "n_rows": len(sample),
                "n_orbits": sample["orbit"].nunique(),
                "n_accepted": sample["annotation_status"].eq(1).sum(),
                "hardy_included": "hardy_pointwise_energy" in sample,
                "n_support_rows": int(data["support_seconds"].eq(support_seconds).sum())
            })

            baseline_score = summarize(
                predictions, f"prediction_{BASELINE_MODEL}"
            )["grouped_mean_absolute_log_error"]
            for model in available_models(sample):
                score = summarize(predictions, f"prediction_{model}")
                score.update({
                    "cohort": cohort,
                    "support_seconds": support_seconds,
                    "model": model,
                    "improvement_over_log_constant": (
                        baseline_score
                        - score["grouped_mean_absolute_log_error"]
                    )
                })
                if model != BASELINE_MODEL and bootstrap_draws:
                    low, high = bootstrap_improvement(
                        predictions, model, bootstrap_draws
                    )
                    score["improvement_ci_low"] = low
                    score["improvement_ci_high"] = high
                all_scores.append(score)

    return (
        pd.concat(all_predictions, ignore_index=True),
        pd.concat(all_parameters, ignore_index=True),
        pd.DataFrame(all_scores),
        pd.DataFrame(sample_counts)
    )


#%% Command line

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--cohort", choices=("all", "accepted"), nargs="+",
        default=["all", "accepted"]
    )
    parser.add_argument("--support-seconds", type=int, nargs="+")
    parser.add_argument("--minimum-coverage", type=float, default=0.50)
    parser.add_argument("--minimum-flux", type=float, default=0.0)
    parser.add_argument("--maximum-flux", type=float, default=np.inf)
    parser.add_argument("--folds", type=int, default=DEFAULT_FOLDS)
    parser.add_argument("--bootstrap", type=int, default=DEFAULT_BOOTSTRAP)
    args = parser.parse_args()

    with xr.open_dataset(args.input) as dataset:
        data = dataset.to_dataframe().reset_index(drop=True)
    supports = args.support_seconds or sorted(data["support_seconds"].unique())

    predictions, parameters, scores, counts = evaluate(
        data, args.cohort, supports, args.minimum_coverage,
        args.minimum_flux, args.maximum_flux, args.folds, args.bootstrap
    )
    for frame in (predictions, parameters, scores, counts):
        frame["minimum_coverage"] = args.minimum_coverage
        frame["minimum_flux"] = args.minimum_flux
        frame["maximum_flux"] = args.maximum_flux

    args.output.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(args.output / "predictions.csv", index=False)
    parameters.to_csv(args.output / "fold_parameters.csv", index=False)
    scores.to_csv(args.output / "scores.csv", index=False)
    counts.to_csv(args.output / "sample_counts.csv", index=False)

    columns = [
        "cohort", "support_seconds", "model",
        "grouped_mean_absolute_log_error", "multiplicative_error_factor",
        "improvement_over_log_constant", "improvement_ci_low",
        "improvement_ci_high", "n_rows", "n_orbits"
    ]
    print(scores.reindex(columns=columns).round(4).to_string(index=False))
    print(f"Saved results to {args.output}")


if __name__ == "__main__":
    main()
