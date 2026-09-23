from pathlib import Path
import sys

import numpy as np
import pandas as pd
import xarray as xr


SCRIPT_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "scripts/ratio_relation_validation"
)
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import fetch_all_dmsp_crossings as fetch
from ratio_validation import MATCH_COLUMNS
from ratio_validation_background_comparison import (
    PAIR_KEYS,
    branch_frame,
    load_paired_matches,
    select_common_support,
)


class FakeProduct:
    product_type = "precipitation_detector"
    time = np.array(["2001-01-01T00:00:00"], dtype="datetime64[ns]")
    wic_source_index = np.array([7])
    source_fuv_detector_sha256 = "abc123"
    source_preprocessing_label = "fuvpy_bs_directional_v1"
    proton_energy_model = "hardy"
    attrs = {"count_source": "unsubtracted"}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return None

    def read(self, name):
        if name == "mlat":
            return np.array([[[70.0]]])
        if name == "mlt":
            return np.array([[[12.0]]])
        if name == "method_valid":
            return np.ones((1, 1, 1), dtype=bool)
        if name == "sza":
            return np.array([[[111.0]]])
        return np.ones((1, 1, 1))


def dmsp_sample():
    time = np.array(["2001-01-01T00:00:00"], dtype="datetime64[ns]")
    values = {
        "mlat": ("time", [70.0]),
        "mlt": ("time", [12.0]),
    }
    for name in fetch.DMSP_FIELDS:
        values[name] = ("time", [1.0])
    return xr.Dataset(values, coords={"time": time})


def test_extractor_saves_wic_sza_and_count_source(monkeypatch):
    monkeypatch.setattr(fetch, "open_product", lambda filename: FakeProduct())
    monkeypatch.setattr(
        fetch, "load_dmsp",
        lambda files, cache, satellite, start, stop:
            dmsp_sample() if satellite == "f12" else None,
    )

    result = fetch.process_orbit(Path("or_0001.nc"), {}, {})

    assert result.sizes["sample"] == 1
    assert result["wic_sza"].item() == 111
    assert result.attrs["source_count_source"] == "unsubtracted"


def crossing_data(ratios, sza, valid):
    count = len(ratios)
    times = pd.date_range("2001-01-01", periods=count, freq="s")
    values = {
        "orbit": np.ones(count, dtype=int),
        "image_frame": np.zeros(count, dtype=int),
        "image_time": np.full(count, np.datetime64("2001-01-01")),
        "dmsp_sat": np.full(count, "f12"),
        "dmsp_time": times.to_numpy(),
        "dmsp_mlat": np.full(count, 70.0),
        "dmsp_mlt": np.full(count, 12.0),
        "detector_separation_deg": np.full(count, 0.1),
        "dmsp_electron_mean_energy": np.full(count, 2.0),
        "dmsp_electron_mean_energy_fractional_std": np.full(count, 0.1),
        "dmsp_electron_total_energy_flux": np.full(count, 1e12),
        "dmsp_electron_total_energy_flux_fractional_std": np.full(count, 0.1),
        "img_wic": np.full(count, 100.0),
        "img_wic_std": np.full(count, 2.0),
        "img_si13": np.full(count, 2.0),
        "img_si13_std": np.full(count, 0.1),
        "img_ratio": np.asarray(ratios),
        "img_ratio_std": np.full(count, 1.0),
        "img_energy": np.full(count, 2.0),
        "img_energy_std": np.full(count, 0.2),
        "wic_dza": np.full(count, 20.0),
        "quality_weight": np.ones(count),
        "method_valid": np.asarray(valid),
        "wic_sza": np.asarray(sza),
    }
    return pd.DataFrame(values)


def save_crossings(path, data):
    path.mkdir()
    xr.Dataset({
        name: ("sample", value.to_numpy())
        for name, value in data.items()
    }).to_netcdf(path / "or_0001.nc")


def test_paired_selection_uses_exact_keys_common_validity_and_sza(tmp_path):
    background_data = crossing_data([40.0, 50.0, 60.0], [110, 110, 100],
                                    [True, True, True])
    unsubtracted_data = crossing_data([45.0, 55.0, 65.0], [110, 110, 110],
                                      [True, False, True])
    background_path = tmp_path / "background"
    unsubtracted_path = tmp_path / "unsubtracted"
    save_crossings(background_path, background_data)
    save_crossings(unsubtracted_path, unsubtracted_data)

    paired = load_paired_matches(background_path, unsubtracted_path)
    selected = select_common_support(paired, minimum_sza=105)

    assert len(paired) == 3
    assert len(selected) == 1
    assert selected[PAIR_KEYS].iloc[0]["dmsp_time"] == pd.Timestamp(
        "2001-01-01T00:00:00"
    )
    assert branch_frame(selected, "background")["img_ratio"].item() == 40
    assert branch_frame(selected, "unsubtracted")["img_ratio"].item() == 45
    assert set(MATCH_COLUMNS).issubset(branch_frame(selected, "background"))
