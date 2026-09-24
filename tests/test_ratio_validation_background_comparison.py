from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
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
    plot_comparison,
    select_common_support,
)


class FakeGrid:
    def ingrid(self, longitude, latitude):
        return np.ones(np.asarray(latitude).shape, dtype=bool)

    def bin_index(self, longitude, latitude):
        shape = np.asarray(latitude).shape
        return np.zeros(shape, dtype=int), np.zeros(shape, dtype=int)


class IdentityProjection:
    def geo2cube(self, longitude, latitude, set_points_off_cube_to_nan=True):
        return np.asarray(longitude), np.asarray(latitude)


class TwoByTwoGrid:
    shape = (2, 2)
    size = 4
    projection = IdentityProjection()
    xi_mesh = np.array([[0.0, 1.0, 2.0]] * 3)
    eta_mesh = np.array([[0.0] * 3, [1.0] * 3, [2.0] * 3])


class FakeCSProduct:
    product_type = "precipitation_cs"
    representation = "cs"
    grid_id = "test_grid"
    grid = FakeGrid()
    time = np.array(["2001-01-01T00:00:00"], dtype="datetime64[ns]")

    def __init__(self, coverage=1.0, data=1.0):
        self.coverage = coverage
        self.data = data

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return None

    def read(self, name):
        if name in ("wic_coverage", "si13_coverage"):
            return np.full((1, 1, 1), self.coverage)
        if name in ("wic_corrected", "si13_corrected"):
            return np.full((1, 1, 1), self.data)
        raise KeyError(name)


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

    def __init__(self, detector_pair=True):
        self.detector_pair = detector_pair

    def read(self, name):
        if name == "mlat":
            return np.array([[[70.0]]])
        if name == "mlt":
            return np.array([[[12.0]]])
        if name == "method_valid":
            return np.ones((1, 1, 1), dtype=bool)
        if name == "sza":
            return np.array([[[111.0]]])
        if name in ("wic_corrected", "si13_corrected"):
            value = 1.0 if self.detector_pair else np.nan
            return np.full((1, 1, 1), value)
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


def use_fake_footprint_match(monkeypatch):
    monkeypatch.setattr(
        fetch, "make_detector_footprint_matcher",
        lambda mlat, mlt, grid: object(),
    )
    monkeypatch.setattr(
        fetch, "match_detector_footprints",
        lambda matcher, grid, cs_row, cs_column, dmsp_mlat, dmsp_mlt,
        image_mlat, image_mlt: (
            np.zeros(len(dmsp_mlat), dtype=np.int16),
            np.zeros(len(dmsp_mlat), dtype=np.int16),
            np.full(len(dmsp_mlat), 0.1, dtype=np.float32),
            np.ones(len(dmsp_mlat), dtype=np.int16),
        ),
    )


def test_detector_match_requires_point_inside_inferred_footprint():
    image_mlat = np.array([[0.5, 0.5], [1.5, 1.5]])
    image_mlt = np.array([[0.5, 1.5], [0.5, 1.5]]) / 15
    grid = TwoByTwoGrid()
    matcher = fetch.make_detector_footprint_matcher(
        image_mlat, image_mlt, grid
    )

    row, column, _, match_count = fetch.match_detector_footprints(
        matcher,
        grid,
        cs_row=np.array([0, 1, 0]),
        cs_column=np.array([0, 1, 1]),
        dmsp_mlat=np.array([0.6, 2.2, 0.5]),
        dmsp_mlt=np.array([0.6, 2.2, 1.0]) / 15,
        image_mlat=image_mlat,
        image_mlt=image_mlt,
    )

    assert (row[0], column[0], match_count[0]) == (0, 0, 1)
    assert (row[1], column[1], match_count[1]) == (-1, -1, 0)
    assert (row[2], column[2], match_count[2]) == (0, 0, 2)


def test_extractor_saves_wic_sza_and_count_source(monkeypatch):
    use_fake_footprint_match(monkeypatch)
    monkeypatch.setattr(
        fetch, "open_product",
        lambda filename: (
            FakeCSProduct() if Path(filename).parent.name == "cs"
            else FakeProduct()
        ),
    )
    monkeypatch.setattr(
        fetch, "load_dmsp",
        lambda files, cache, satellite, start, stop:
            dmsp_sample() if satellite == "f12" else None,
    )

    result = fetch.process_orbit(
        Path("or_0001.nc"), Path("cs/or_0001.nc"), {}, {}
    )

    assert result.sizes["sample"] == 1
    assert result["wic_sza"].item() == 111
    assert result["cs_wic_coverage"].item() == 1
    assert result["cs_si13_coverage"].item() == 1
    assert result["detector_footprint_match_count"].item() == 1
    assert result["detector_pair_valid"].item()
    assert result.attrs["source_count_source"] == "unsubtracted"
    assert result.attrs["spatial_filter"] == fetch.SPATIAL_FILTER


def test_cs_gate_avoids_opening_detector_without_covered_samples(monkeypatch):
    opened = []

    def open_product(filename):
        opened.append(Path(filename))
        if Path(filename).parent.name == "cs":
            return FakeCSProduct(coverage=0)
        raise AssertionError("detector product should not be opened")

    monkeypatch.setattr(fetch, "open_product", open_product)
    monkeypatch.setattr(
        fetch, "load_dmsp",
        lambda files, cache, satellite, start, stop:
            dmsp_sample() if satellite == "f12" else None,
    )

    result = fetch.process_orbit(
        Path("or_0001.nc"), Path("cs/or_0001.nc"), {}, {}
    )

    assert result.sizes["sample"] == 0
    assert opened == [Path("cs/or_0001.nc")]
    assert result.attrs["spatial_filter"] == fetch.SPATIAL_FILTER


def test_cs_gate_avoids_detector_when_binned_data_are_missing(monkeypatch):
    opened = []

    def open_product(filename):
        opened.append(Path(filename))
        if Path(filename).parent.name == "cs":
            return FakeCSProduct(coverage=1, data=np.nan)
        raise AssertionError("detector product should not be opened")

    monkeypatch.setattr(fetch, "open_product", open_product)
    monkeypatch.setattr(
        fetch, "load_dmsp",
        lambda files, cache, satellite, start, stop:
            dmsp_sample() if satellite == "f12" else None,
    )

    result = fetch.process_orbit(
        Path("or_0001.nc"), Path("cs/or_0001.nc"), {}, {}
    )

    assert result.sizes["sample"] == 0
    assert opened == [Path("cs/or_0001.nc")]


def test_cs_supported_missing_detector_pair_is_recorded(monkeypatch):
    use_fake_footprint_match(monkeypatch)
    monkeypatch.setattr(
        fetch, "open_product",
        lambda filename: (
            FakeCSProduct() if Path(filename).parent.name == "cs"
            else FakeProduct(detector_pair=False)
        ),
    )
    monkeypatch.setattr(
        fetch, "load_dmsp",
        lambda files, cache, satellite, start, stop:
            dmsp_sample() if satellite == "f12" else None,
    )

    result = fetch.process_orbit(
        Path("or_0001.nc"), Path("cs/or_0001.nc"), {}, {}
    )

    assert not result["detector_pair_valid"].item()
    assert result.attrs["detector_pair_missing_count"] == 1


def test_crossing_output_completion_detects_interrupted_file(tmp_path):
    filename = tmp_path / "or_0001.nc"
    complete = {
        name: ("sample", np.ones(2)) for name in fetch.OUTPUT_FIELDS
    }
    xr.Dataset(complete).to_netcdf(filename)
    assert fetch.output_is_complete(filename)
    assert not fetch.output_is_complete(
        filename, required_spatial_filter=fetch.SPATIAL_FILTER
    )

    xr.Dataset(complete, attrs={
        "spatial_filter": fetch.SPATIAL_FILTER
    }).to_netcdf(filename)
    assert fetch.output_is_complete(
        filename, required_spatial_filter=fetch.SPATIAL_FILTER
    )

    xr.Dataset({"img_ratio": ("sample", np.ones(2))}).to_netcdf(filename)
    assert not fetch.output_is_complete(filename)


def test_save_orbit_keeps_final_file_when_write_is_interrupted(
    tmp_path, monkeypatch
):
    filename = tmp_path / "or_0001.nc"
    filename.write_bytes(b"previous complete output")
    data = xr.Dataset({"value": ("sample", np.ones(2))})

    def interrupted_write(self, temporary, **kwargs):
        Path(temporary).write_bytes(b"partial output")
        raise RuntimeError("interrupted")

    monkeypatch.setattr(xr.Dataset, "to_netcdf", interrupted_write)
    with pytest.raises(RuntimeError, match="interrupted"):
        fetch.save_orbit(data, filename)

    assert filename.read_bytes() == b"previous complete output"
    assert filename.with_suffix(".nc.partial").exists()


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


def test_paired_comparison_plots_both_quantile_directions(tmp_path):
    ratios = np.linspace(10, 120, 40)
    data = crossing_data(ratios, np.full(40, 110.0), np.ones(40, dtype=bool))

    plot_comparison(data, data, tmp_path, crossings=1, orbits=1, minimum_sza=0)

    expected = [
        "ratio_energy_background_comparison",
        "ratio_energy_background_comparison_energy_binned",
    ]
    for filename in expected:
        assert (tmp_path / f"{filename}.png").exists()
        assert (tmp_path / f"{filename}.pdf").exists()
