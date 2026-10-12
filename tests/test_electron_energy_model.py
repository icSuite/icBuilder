from pathlib import Path
import sys

import numpy as np
import pandas as pd


SCRIPT_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "scripts/electron_energy_model"
)
sys.path.insert(0, str(SCRIPT_DIRECTORY))

from build_energy_comparison import reduce_crossing
from evaluate_energy_models import (
    assign_folds,
    fit_model,
    predict_log_energy,
    select_common_sample
)
from evaluate_energy_models import cross_validate


def crossing_rows():
    time = np.datetime64("2001-01-01T00:00:00", "ns")
    offsets = np.arange(-60, 61)
    return pd.DataFrame({
        "orbit": 12,
        "image_frame": 3,
        "image_time": time,
        "dmsp_sat": "f13",
        "dmsp_time": time + offsets.astype("timedelta64[s]"),
        "time_offset_seconds": offsets,
        "dmsp_mlat": np.full(offsets.size, 70.0),
        "dmsp_mlt": np.full(offsets.size, 23.0),
        "detector_separation_deg": np.full(offsets.size, 0.2),
        "dmsp_electron_mean_energy": np.where(offsets < 0, 1.0, 2.0),
        "dmsp_electron_mean_energy_fractional_std": 0.1,
        "dmsp_electron_total_energy_flux": 4.0,
        "dmsp_electron_total_energy_flux_fractional_std": 0.1,
        "img_ratio": 20.0,
        "img_energy": 3.0,
        "img_energy_std": 0.5,
        "method_valid": True,
        "zp_pointwise_energy": np.where(offsets < 0, 2.0, 4.0),
        "zp_pointwise_flux": 6.0
    })


def test_reduce_crossing_uses_half_open_window_and_flux_number_reduction():
    rows = crossing_rows()
    key = (12, np.datetime64("2001-01-01T00:00:00", "ns"), "f13")
    result = reduce_crossing(rows, 120, {key: 1})

    selected = rows.loc[rows.time_offset_seconds.lt(60)]
    energy = selected.dmsp_electron_mean_energy.to_numpy()
    flux = selected.dmsp_electron_total_energy_flux.to_numpy()
    expected = flux.sum() / np.sum(flux / energy)

    assert result["n_available"] == 120
    assert result["n_dmsp"] == 120
    assert result["dmsp_coverage"] == 1
    assert np.isclose(result["dmsp_energy"], expected)
    pointwise_energy = selected.zp_pointwise_energy.to_numpy()
    pointwise_flux = selected.zp_pointwise_flux.to_numpy()
    expected_pointwise = (
        pointwise_flux.sum()
        / np.sum(pointwise_flux / pointwise_energy)
    )
    assert np.isclose(result["zp_pointwise_energy"], expected_pointwise)
    assert result["annotation_status"] == 1


def test_orbit_folds_never_split_an_orbit():
    frame = pd.DataFrame({
        "orbit": np.repeat(np.arange(10), 3),
        "year": np.repeat([2000, 2001], 15)
    })
    folds, n_folds = assign_folds(frame, 5)
    frame["fold"] = folds

    assert n_folds == 5
    assert frame.groupby("orbit")["fold"].nunique().max() == 1
    assert frame["fold"].nunique() == 5


def test_hardcoded_constant_predicts_exactly_two_kev():
    frame = pd.DataFrame({
        "orbit": [1, 2],
        "dmsp_energy": [0.5, 5.0],
        "dmsp_flux": [1e10, 1e12]
    })
    parameters = fit_model(frame, "constant_2kev")
    prediction = 10**predict_log_energy(
        frame, "constant_2kev", parameters
    )

    assert parameters["slope"] == 0
    assert np.allclose(prediction, 2.0)


def test_accepted_cohort_is_a_subset_of_the_same_common_sample():
    data = pd.DataFrame({
        "support_seconds": [66, 66, 66],
        "annotation_status": [1, 0, -1],
        "dmsp_coverage": [1.0, 1.0, 1.0],
        "dmsp_energy": [1.0, 2.0, 3.0],
        "dmsp_flux": [1e10, 1e10, 1e10],
        "zp_pointwise_energy": [1.0, 1.0, 1.0],
        "zp_area_mean_energy": [1.0, 1.0, 1.0],
        "zp_area_median_energy": [1.0, 1.0, 1.0],
        "zp_flux_weighted_energy": [1.0, 1.0, 1.0]
    })

    all_rows = select_common_sample(data, "all", 66, 0.5, 0, np.inf)
    accepted = select_common_sample(data, "accepted", 66, 0.5, 0, np.inf)

    assert len(all_rows) == 3
    assert len(accepted) == 1
    assert accepted.iloc[0]["annotation_status"] == 1


def test_hardy_temporal_reduction_uses_model_fluxes_on_valid_dmsp_samples():
    rows = crossing_rows()
    rows["hardy_pointwise_energy"] = np.where(rows.time_offset_seconds < 0, 1.0, 3.0)
    rows["hardy_pointwise_flux"] = np.where(rows.time_offset_seconds < 0, 2.0, 6.0)
    result = reduce_crossing(rows, 120, {})
    # Both halves supply number flux 2; total energy/number gives 2 keV.
    assert result["n_hardy_pointwise"] == 120
    assert result["hardy_pointwise_energy"] == 2.0
    assert result["hardy_pointwise_flux"] == 4.0
    rows["hardy_pointwise_energy"] = np.nan
    missing = reduce_crossing(rows, 120, {})
    assert missing["n_hardy_pointwise"] == 0
    assert np.isnan(missing["hardy_pointwise_energy"])
    assert np.isnan(missing["hardy_pointwise_flux"])


def test_hardy_common_support_and_uncalibrated_orbit_predictions():
    data = pd.DataFrame({
        "orbit": [1, 1, 2, 2, 3], "year": 2001,
        "image_frame": 0, "image_time": np.datetime64("2001-01-01"),
        "satellite": "f13", "support_seconds": 120,
        "annotation_status": 1, "dmsp_coverage": 1.0,
        "dmsp_energy": [9.0, 8.0, 7.0, 6.0, 5.0], "dmsp_flux": 1e10,
        "hardy_pointwise_energy": [1.0, 2.0, 3.0, 4.0, np.nan],
        "n_hardy_pointwise": [120, 120, 120, 120, 0], "n_dmsp": 120,
        "zp_pointwise_energy": 1.0, "zp_area_mean_energy": 1.0,
        "zp_area_median_energy": 1.0, "zp_flux_weighted_energy": 1.0
    })
    sample = select_common_sample(data, "all", 120, 0.5, 0, np.inf)
    assert len(sample) == 4
    partial = data.copy()
    partial.loc[0, "n_hardy_pointwise"] = 119
    assert len(select_common_sample(partial, "all", 120, 0.5, 0, np.inf)) == 3
    predictions, parameters = cross_validate(sample, 2)
    np.testing.assert_allclose(predictions["prediction_hardy_pointwise"], predictions["hardy_pointwise_energy"])
    hardy_parameters = parameters.loc[parameters.model.eq("hardy_pointwise")]
    assert hardy_parameters.intercept.eq(0).all()
    assert hardy_parameters.slope.eq(1).all()
    # Perturbing particle measurements cannot calibrate the raw Hardy model.
    sample["dmsp_energy"] *= 100
    changed, _ = cross_validate(sample, 2)
    np.testing.assert_array_equal(changed.prediction_hardy_pointwise, predictions.prediction_hardy_pointwise)
