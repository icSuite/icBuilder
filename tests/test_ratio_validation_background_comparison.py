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
import fetch_all_dmsp_crossings_cs as fetch_cs
import ratio_validation as validation
from ratio_validation import MATCH_COLUMNS
from ratio_validation_cs import MATCH_COLUMNS as CS_MATCH_COLUMNS
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


def test_cs_crossing_contract_carries_wic_sza():
    assert fetch_cs.IMAGE_FIELDS["sza"] == "wic_sza"
    assert "wic_sza" in fetch_cs.OUTPUT_FIELDS
    assert "wic_sza" in CS_MATCH_COLUMNS


def test_ratio_validation_quality_limits_are_configurable():
    data = pd.DataFrame({
        "img_ratio": np.ones(7),
        "dmsp_electron_mean_energy": np.ones(7),
        "wic_dza": np.ones(7),
        "detector_separation_deg": np.full(7, 0.1),
        "method_valid": np.ones(7, dtype=bool),
        "dmsp_electron_raw_counts_valid": np.ones(7, dtype=np.int8),
        "dmsp_electron_mean_energy_fractional_std": [0.1, 0.1, 0.1, 0.1, 0.1, 0.3, 0.1],
        "dmsp_electron_total_energy_flux_fractional_std": [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3],
        "dmsp_electron_total_energy_flux": [1e11, 7e10, 3e11, 1e11, 1e11, 1e11, 1e11],
        "wic_sza": [120, 120, 120, 100, 160, 120, 120],
    })

    selected = validation.quality_selection(
        data,
        min_flux=8e10,
        max_flux=2e11,
        min_sza=105,
        max_sza=150,
        max_energy_fractional_uncertainty=0.25,
        max_flux_fractional_uncertainty=0.2,
    )

    assert selected.tolist() == [True, False, False, False, False, False, False]


def test_ratio_validation_requested_defaults_and_cli_forwarding(monkeypatch):
    assert validation.MIN_FLUX == 0
    assert validation.MAX_FLUX == 1e14
    assert validation.MIN_SZA == 0
    assert validation.MAX_SZA == 180
    assert validation.MAX_ENERGY_FRACTIONAL_UNCERTAINTY == 0.25
    assert validation.MAX_FLUX_FRACTIONAL_UNCERTAINTY == 0.2

    received = {}
    monkeypatch.setattr(
        validation,
        "run",
        lambda **settings: received.update(settings),
    )
    validation.main([
        "--min-flux", "8e10",
        "--max-flux", "2e11",
        "--min-sza", "105",
        "--max-sza", "150",
        "--max-energy-fractional-uncertainty", "0.2",
        "--max-flux-fractional-uncertainty", "0.1",
    ])

    assert received["min_flux"] == 8e10
    assert received["max_flux"] == 2e11
    assert received["min_sza"] == 105
    assert received["max_sza"] == 150
    assert received["max_energy_fractional_uncertainty"] == 0.2
    assert received["max_flux_fractional_uncertainty"] == 0.1


def test_ratio_quality_rejects_bad_raw_counts_without_requiring_valid_conductance():
    data = crossing_data([40., 40., 40.], [110., 110., 110.], [True, True, True])
    data["dmsp_electron_raw_counts_valid"] = [1., 0., np.nan]
    data["dmsp_electron_conductance_valid"] = 0
    data["dmsp_electron_hall_conductance"] = np.nan
    assert validation.quality_selection(data).tolist() == [True, False, False]


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
    values["electron_pedersen_conductance_std"] = ("time", [0.2])
    values["electron_hall_conductance_std"] = ("time", [0.3])
    values["electron_pedersen_hall_conductance_covariance"] = ("time", [0.04])
    values["electron_conductance_uncertainty_channel_count"] = ("time", np.array([19], dtype=np.int16))
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


def test_crossing_restart_detects_changed_yearly_inputs_including_empty_orbits(tmp_path):
    yearly = tmp_path / "year.nc"
    yearly.write_bytes(b"initial yearly inputs")
    files = {("f13", 2000): yearly}
    initial = fetch.dmsp_source_metadata(files)
    output = tmp_path / "or_0001.nc"
    xr.Dataset(coords={"sample": np.arange(0)}, attrs=initial).to_netcdf(output)
    assert fetch.output_is_complete(output, required_dmsp_signature=initial["source_dmsp_yearly_signature"])
    yearly.write_bytes(b"updated corrected yearly inputs")
    changed = fetch.dmsp_source_metadata(files)
    assert not fetch.output_is_complete(output, required_dmsp_signature=changed["source_dmsp_yearly_signature"])


def test_cs_samples_retain_rejected_conductances_and_integer_flags(tmp_path, monkeypatch):
    first = dmsp_sample()
    second = first.assign_coords(time=first.time + np.timedelta64(1, "s")).copy(deep=True)
    second["electron_raw_counts_valid"][:] = 0
    second["electron_conductance_valid"][:] = 0
    second["electron_pedersen_conductance"][:] = np.nan
    second["electron_hall_conductance"][:] = np.nan
    second["electron_conductance_uncertainty_valid"][:] = 0
    second["electron_conductance_uncertainty_channel_count"][:] = 0
    for name in ("electron_pedersen_conductance_std", "electron_hall_conductance_std", "electron_pedersen_hall_conductance_covariance"):
        second[name][:] = np.nan
    dmsp = xr.concat([first, second], dim="time")
    monkeypatch.setattr(fetch_cs, "nearest_grid_cells", lambda *args: (np.zeros(2, dtype=int), np.zeros(2, dtype=int), np.zeros(2), np.ones(2, dtype=bool)))
    samples = fetch_cs.make_samples()
    fetch_cs.append_frame_samples(samples, FakeCSProduct(), None, 0, "f13", dmsp)
    data = xr.Dataset({name: ("sample", np.concatenate(parts)) for name, parts in samples.items()})
    output = tmp_path / "or_0001.nc"
    fetch.save_orbit(data, output)
    with xr.open_dataset(output) as saved:
        np.testing.assert_array_equal(saved.dmsp_electron_raw_counts_valid, [1, 0])
        assert saved.dmsp_electron_raw_counts_valid.dtype == np.int8
        assert saved.dmsp_electron_hall_conductance.attrs["units"] == "S"
        assert np.isnan(saved.dmsp_electron_hall_conductance.values[1])
        assert saved.sizes["sample"] == 2


@pytest.mark.parametrize("representation", ["detector", "cs"])
def test_crossings_retain_nominal_samples_with_missing_spectral_uncertainty(tmp_path, monkeypatch, representation):
    first = dmsp_sample()
    second = first.assign_coords(time=first.time + np.timedelta64(1, "s")).copy(deep=True)
    second["electron_conductance_uncertainty_valid"][:] = 0
    second["electron_conductance_uncertainty_channel_count"][:] = 0
    for name in ("electron_pedersen_conductance_std", "electron_hall_conductance_std", "electron_pedersen_hall_conductance_covariance"):
        second[name][:] = np.nan
    partial = first.assign_coords(time=first.time + np.timedelta64(2, "s")).copy(deep=True)
    partial["electron_conductance_uncertainty_channel_count"][:] = 17
    partial["electron_pedersen_conductance_std"][:] = 0.15
    partial["electron_hall_conductance_std"][:] = 0.25
    partial["electron_pedersen_hall_conductance_covariance"][:] = 0.03
    dmsp = xr.concat([first, second, partial], dim="time")
    if representation == "detector":
        use_fake_footprint_match(monkeypatch)
        monkeypatch.setattr(fetch, "open_product", lambda filename: FakeCSProduct() if Path(filename).parent.name == "cs" else FakeProduct())
        monkeypatch.setattr(fetch, "load_dmsp", lambda files, cache, satellite, start, stop: dmsp if satellite == "f13" else None)
        data = fetch.process_orbit(Path("or_0001.nc"), Path("cs/or_0001.nc"), {}, {})
        required = fetch.OUTPUT_FIELDS
    else:
        monkeypatch.setattr(fetch_cs, "nearest_grid_cells", lambda *args: (np.zeros(3, dtype=int), np.zeros(3, dtype=int), np.zeros(3), np.ones(3, dtype=bool)))
        samples = fetch_cs.make_samples()
        fetch_cs.append_frame_samples(samples, FakeCSProduct(), None, 0, "f13", dmsp)
        data = xr.Dataset({name: ("sample", np.concatenate(parts)) for name, parts in samples.items()})
        required = fetch_cs.OUTPUT_FIELDS
    output = tmp_path / "or_0001.nc"
    fetch.save_orbit(data, output)
    with xr.open_dataset(output) as saved:
        assert saved.sizes["sample"] == 3
        np.testing.assert_array_equal(saved.dmsp_electron_conductance_valid, [1, 1, 1])
        np.testing.assert_array_equal(saved.dmsp_electron_conductance_uncertainty_valid, [1, 0, 1])
        np.testing.assert_array_equal(saved.dmsp_electron_conductance_uncertainty_channel_count, [19, 0, 17])
        assert saved.dmsp_electron_conductance_uncertainty_channel_count.dtype == np.int16
        assert saved.dmsp_electron_conductance_uncertainty_channel_count.attrs["units"] == "1"
        np.testing.assert_array_equal(saved.dmsp_electron_conductance_uncertainty_channel_count.attrs["valid_range"], [0, 19])
        assert saved.dmsp_electron_conductance_uncertainty_valid.dtype == np.int8
        assert saved.dmsp_electron_conductance_uncertainty_valid.attrs["units"] == "1"
        for name, expected, partial_expected in [("pedersen_conductance_std", 0.2, 0.15), ("hall_conductance_std", 0.3, 0.25), ("pedersen_hall_conductance_covariance", 0.04, 0.03)]:
            variable = saved[f"dmsp_electron_{name}"]
            assert variable.values[0] == pytest.approx(expected)
            assert np.isnan(variable.values[1])
            assert variable.values[2] == pytest.approx(partial_expected)
            assert variable.attrs["units"] == ("S2" if "covariance" in name else "S")
            assert f"dmsp_electron_{name}" in required
        assert np.isfinite(saved.dmsp_electron_hall_conductance.values).all()
        assert np.isfinite(saved.dmsp_electron_pedersen_conductance.values).all()
        assert saved.attrs["dmsp_conductance_uncertainty_method"] == fetch.DMSP_CONDUCTANCE_UNCERTAINTY_METHOD
        assert "NaN channel-error terms are omitted" in saved.attrs["dmsp_conductance_uncertainty_method"]


@pytest.mark.parametrize("module", [fetch, fetch_cs])
@pytest.mark.parametrize("missing", ["electron_pedersen_conductance_std", "electron_hall_conductance_std", "electron_pedersen_hall_conductance_covariance", "electron_conductance_uncertainty_valid", "electron_conductance_uncertainty_channel_count"])
def test_crossing_restart_rejects_missing_conductance_uncertainty_fields(tmp_path, module, missing):
    data = xr.Dataset({name: ("sample", np.ones(2)) for name in module.OUTPUT_FIELDS})
    output = tmp_path / "or_0001.nc"
    data.to_netcdf(output)
    assert fetch.output_is_complete(output, module.OUTPUT_FIELDS)
    data.drop_vars(f"dmsp_{missing}").to_netcdf(output)
    assert not fetch.output_is_complete(output, module.OUTPUT_FIELDS)


@pytest.mark.parametrize("missing", ["electron_raw_counts_valid", "electron_hall_conductance_std", "electron_conductance_uncertainty_valid", "electron_conductance_uncertainty_channel_count"])
def test_load_dmsp_requires_corrected_yearly_fields(tmp_path, missing):
    yearly = tmp_path / "year.nc"
    data = dmsp_sample().drop_vars(missing)
    data.to_netcdf(yearly)
    cache = {}
    try:
        with pytest.raises(ValueError, match="rebuild the yearly DMSP file"):
            fetch.load_dmsp({("f13", 2001): yearly}, cache, "f13", data.time.values[0], data.time.values[0])
    finally:
        for dataset in cache.values():
            dataset.close()


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
        "dmsp_electron_pedersen_conductance": np.full(count, 5.0),
        "dmsp_electron_hall_conductance": np.full(count, 10.0),
        "dmsp_electron_pedersen_conductance_std": np.full(count, 0.2),
        "dmsp_electron_hall_conductance_std": np.full(count, 0.3),
        "dmsp_electron_pedersen_hall_conductance_covariance": np.full(count, 0.04),
        "dmsp_electron_conductance_uncertainty_valid": np.ones(count, dtype=np.int8),
        "dmsp_electron_conductance_uncertainty_channel_count": np.full(count, 19, dtype=np.int16),
        "dmsp_electron_conductance_valid": np.ones(count, dtype=np.int8),
        "dmsp_electron_raw_counts_valid": np.ones(count, dtype=np.int8),
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
