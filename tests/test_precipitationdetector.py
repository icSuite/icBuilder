"""Focused tests for detector-space precipitation Product 2."""

#%% Imports and test inputs

from datetime import datetime, timedelta
import importlib.util
from pathlib import Path

import numpy as np
import pytest
from icphysics import hardy_ion_precipitation
from netCDF4 import Dataset

from icbuilder.fuvdetector import FUVDetector, PREPROCESSING_LABEL
from icbuilder.precipitationdetector import (
    COUNT_UNCERTAINTY_MODE,
    PrecipitationDetector,
    SCHEMA_VERSION,
)


SCRIPT = (
    Path(__file__).parents[1] / "scripts" / "pipeline"
    / "make_precipitation_detector_orbit_files.py"
)
SPEC = importlib.util.spec_from_file_location(
    "make_precipitation_detector_orbit_files", SCRIPT
)
ORBIT_SCRIPT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ORBIT_SCRIPT)

BASE_TIME = datetime(2001, 1, 1, 1)


def regular_camera(sensor, size, spacing, times, value):
    """Make regular geographic detector frames with constant observations."""

    row, column = np.indices((size, size))
    centre = (size - 1) / 2
    latitude = 70 + (row - centre) * spacing
    longitude = (column - centre) * spacing / np.cos(np.deg2rad(70))
    frame_shape = (len(times), size, size)

    camera = {
        "sensor": sensor,
        "source_file": f"{sensor.lower()}.nc",
        "image_field": "dgimg",
        "unsubtracted_image_field": "img",
        "time": np.asarray(times, dtype=object),
        "counts": np.full(frame_shape, value, dtype=float),
        "unsubtracted_counts": np.full(
            frame_shape, value + 100.0, dtype=float
        ),
        "variance": np.full(frame_shape, 4.0, dtype=float),
        "frame_quality": np.full(len(times), 2, dtype=np.int8),
        "quality_weight": np.full(frame_shape, 0.8),
        "glat": np.broadcast_to(latitude, frame_shape).copy(),
        "glon": np.broadcast_to(longitude, frame_shape).copy(),
        "sza": np.full(frame_shape, 80.0),
        "dza": np.full(frame_shape, 20.0),
        "geometry_valid": np.ones(frame_shape, dtype=bool),
    }
    if sensor == "WIC":
        camera.update({
            "mlat": np.broadcast_to(latitude, frame_shape).copy(),
            "mlon": np.broadcast_to(longitude, frame_shape).copy(),
            "mlt": np.mod(
                np.broadcast_to(longitude / 15, frame_shape), 24
            ).copy(),
            "ssalon": np.full(len(times), 10.0),
        })
    return camera


def kp_series():
    return {
        "time": np.array(["2001-01-01T00:00:00"], dtype="datetime64[s]"),
        "kp": np.array([2.0]),
        "status": np.array(["def"]),
        "provenance": {"source": "test", "status": "def"},
    }


def write_fuv_detector(path, two_frames=False):
    wic_times = [BASE_TIME]
    si_times = [BASE_TIME + timedelta(seconds=1)]
    if two_frames:
        wic_times.append(BASE_TIME + timedelta(seconds=10))
        si_times.append(BASE_TIME + timedelta(seconds=11))

    wic = regular_camera("WIC", 8, 0.05, wic_times, 1000.0)
    si12 = regular_camera("SI12", 4, 0.10, si_times, 10.0)
    si13 = regular_camera(
        "SI13", 4, 0.10, [BASE_TIME + timedelta(seconds=1)], 20.0
    )
    product = FUVDetector(wic, si12, si13, software_version="product1-test")
    product.to_nc(path)
    return product


#%% Detector calculation

def test_image_ratio_runs_on_product1_detector_geometry(tmp_path):
    source = tmp_path / "fuv_detector.nc"
    fuv = write_fuv_detector(source)

    precipitation = PrecipitationDetector(
        source,
        kp_series=kp_series(),
        proton_energy_model="constant",
        proton_energy=2.0,
        proton_energy_uncertainty=0.0,
        software_version="product2-test",
    )

    assert precipitation.shape == fuv.shape
    assert precipitation.time.tolist() == fuv.time.tolist()
    np.testing.assert_array_equal(
        precipitation.detector_row, np.arange(fuv.shape[1])
    )
    np.testing.assert_array_equal(
        precipitation.detector_column, np.arange(fuv.shape[2])
    )
    np.testing.assert_allclose(precipitation.Ep, 2.0)
    assert precipitation.method_valid.any()
    assert np.isfinite(precipitation.Fp[precipitation.method_valid]).all()
    assert np.isfinite(precipitation.E0[precipitation.method_valid]).all()
    assert np.isfinite(precipitation.Fe[precipitation.method_valid]).all()
    assert np.isfinite(
        precipitation.dwic_corrected[precipitation.method_valid]
    ).all()
    assert precipitation.count_uncertainty_mode == "measurement"

    expected_weight = (
        precipitation.wic_quality_weight
        * precipitation.si12_quality_weight
        * precipitation.si13_quality_weight
    )
    np.testing.assert_allclose(
        precipitation.method_quality_weight[precipitation.method_valid],
        expected_weight[precipitation.method_valid],
    )


def test_missing_si13_keeps_wic_frame_but_invalidates_ratio(tmp_path):
    source = tmp_path / "fuv_detector.nc"
    write_fuv_detector(source, two_frames=True)

    precipitation = PrecipitationDetector(
        source,
        kp_series=kp_series(),
        proton_energy_model="constant",
    )

    assert precipitation.time.size == 2
    np.testing.assert_array_equal(precipitation.si13_source_index, [0, -1])
    assert precipitation.method_valid[0].any()
    assert not precipitation.method_valid[1].any()
    assert np.isnan(precipitation.E0[1]).all()
    assert np.isnan(precipitation.Fe[1]).all()


def test_product2_rejects_schema1_product1(tmp_path):
    source = tmp_path / "fuv_detector.nc"
    write_fuv_detector(source)
    with Dataset(source, "r+") as nc:
        nc.schema_version = 1

    with pytest.raises(ValueError, match="schema_version 1; expected 2"):
        PrecipitationDetector(source, kp_series=kp_series())


def test_product2_rejects_product1_from_old_time_decoder(tmp_path):
    source = tmp_path / "fuv_detector.nc"
    write_fuv_detector(source)
    with Dataset(source, "r+") as nc:
        nc.delncattr("source_time_decoding")

    with pytest.raises(ValueError, match="missing attributes: source_time_decoding"):
        PrecipitationDetector(source, kp_series=kp_series())


def test_zero_proton_flux_retains_detector_uncertainty(tmp_path):
    source = tmp_path / "fuv_detector.nc"
    write_fuv_detector(source)
    with Dataset(source, "r+") as nc:
        valid = nc["si12_valid"][:].astype(bool)
        counts = nc["si12_counts"][:]
        counts[valid] = -10.0
        nc["si12_counts"][:] = counts

    precipitation = PrecipitationDetector(
        source,
        kp_series=kp_series(),
        proton_energy_model="constant",
    )

    assert precipitation.method_valid.any()
    np.testing.assert_allclose(
        precipitation.Fp[precipitation.method_valid], 0.0
    )
    assert np.isfinite(precipitation.E0[precipitation.method_valid]).all()
    assert np.isfinite(precipitation.Fe[precipitation.method_valid]).all()
    assert np.isfinite(precipitation.dFp[precipitation.method_valid]).all()
    assert np.all(precipitation.dFp[precipitation.method_valid] > 0)


def test_hardy_energy_uses_each_detector_pixels_mlat_and_mlt(tmp_path):
    source = tmp_path / "fuv_detector.nc"
    write_fuv_detector(source)

    precipitation = PrecipitationDetector(
        source, kp_series=kp_series(), proton_energy_model="hardy"
    )
    expected = hardy_ion_precipitation(
        precipitation.kp[:, None, None],
        precipitation.mlt,
        precipitation.mlat,
    )["mean_energy"]

    np.testing.assert_allclose(precipitation.Ep_model, expected)
    assert np.all(
        (precipitation.Ep[np.isfinite(precipitation.Ep)] >= 0.47)
        & (precipitation.Ep[np.isfinite(precipitation.Ep)] <= 46.7)
    )
    assert "approximated" in precipitation.proton_energy_coordinate_note


#%% NetCDF and restart boundary

def test_precipitation_detector_netcdf_is_self_describing(tmp_path):
    source = tmp_path / "fuv_detector.nc"
    write_fuv_detector(source)
    product = PrecipitationDetector(
        source,
        kp_series=kp_series(),
        proton_energy_model="constant",
        proton_energy=5.0,
        proton_energy_uncertainty=0.5,
        software_version="product2-test",
    )
    output = tmp_path / "precipitation_detector.nc"

    ORBIT_SCRIPT.save_precipitation_detector(product, output)

    assert ORBIT_SCRIPT.precipitation_detector_file_status(
        output, source, "constant", 5.0, 0.5
    ) == "complete"
    assert ORBIT_SCRIPT.precipitation_detector_file_status(
        output, source, "constant", 2.0, 0.0
    ) == "mismatch"
    assert not Path(str(output) + ".partial").exists()

    with Dataset(output) as nc:
        assert nc.product_type == "precipitation_detector"
        assert nc.representation == "detector"
        assert nc.schema_version == SCHEMA_VERSION
        assert nc.method == "image_ratio"
        assert nc.proton_flux_source == "SI12"
        assert nc.proton_energy_model == "constant"
        assert nc.source_fuv_detector == str(source)
        assert nc.source_fuv_detector_schema_version == 2
        assert nc.source_preprocessing_label == PREPROCESSING_LABEL
        assert (
            nc.source_fuv_detector_time_decoding
            == "CF date units and calendar"
        )
        assert nc.software_version == "product2-test"
        assert nc.count_uncertainty_mode == COUNT_UNCERTAINTY_MODE
        assert "Product-1 full detector measurement variances" in (
            nc.count_uncertainty_method
        )
        assert "coregistered" in nc.proton_operation_order

        expected = {
            "time", "wic_source_time", "si12_source_time",
            "si13_source_time", "wic_source_index", "si12_source_index",
            "si13_source_index", "wic_frame_quality", "si12_frame_quality",
            "si13_frame_quality", "Kp", "Kp_interval_start", "ssalon",
            "detector_row", "detector_column",
            "glat", "glon", "mlat", "mlon", "mlt", "sza", "dza",
            "wic_quality_weight", "si12_quality_weight",
            "si13_quality_weight", "method_quality_weight",
            "wic_coverage", "si12_coverage", "si13_coverage",
            "wic_valid", "si12_valid", "si13_valid", "method_valid",
            "si12_source_count", "si13_source_count",
            "Ep_model", "Ep", "dEp", "Ep_clipping_flag", "Fp", "dFp",
            "wic_corrected", "dwic_corrected",
            "si13_corrected", "dsi13_corrected",
            "R", "dR", "E0", "dE0", "Fe", "dFe", "varE0Fe",
        }
        assert set(nc.variables) == expected
        for forbidden in ("P", "H", "dP", "dH"):
            assert forbidden not in nc.variables

    broken = tmp_path / "broken.nc"
    product.to_nc(broken)
    with Dataset(broken, "r+") as nc:
        nc.renameVariable("method_valid", "removed_method_valid")
    assert ORBIT_SCRIPT.precipitation_detector_file_status(
        broken, source, "constant", 5.0, 0.5
    ) == "invalid"

    old_schema = tmp_path / "old_schema.nc"
    product.to_nc(old_schema)
    with Dataset(old_schema, "r+") as nc:
        nc.schema_version = 2
    assert ORBIT_SCRIPT.precipitation_detector_file_status(
        old_schema, source, "constant", 5.0, 0.5
    ) == "invalid"

    old_time_decoder = tmp_path / "old_time_decoder.nc"
    product.to_nc(old_time_decoder)
    with Dataset(old_time_decoder, "r+") as nc:
        nc.delncattr("source_fuv_detector_time_decoding")
    assert ORBIT_SCRIPT.precipitation_detector_file_status(
        old_time_decoder, source, "constant", 5.0, 0.5
    ) == "invalid"


def test_orbit_script_writes_and_restarts_from_product_file(tmp_path, monkeypatch):
    input_directory = (
        tmp_path / "fuv_detector" / PREPROCESSING_LABEL
    )
    input_directory.mkdir(parents=True)
    write_fuv_detector(input_directory / "or_0001.nc")
    monkeypatch.setattr(ORBIT_SCRIPT, "load_gfz_kp", kp_series)

    arguments = [
        "--base", str(tmp_path),
        "--orbit", "1",
        "--proton-energy-model", "constant",
        "--proton-energy", "5",
        "--proton-energy-uncertainty", "0.5",
    ]
    assert ORBIT_SCRIPT.main(arguments) == [(1, 1)]

    output = tmp_path / "precipitation_detector" / "IR_constant" / "or_0001.nc"
    assert output.is_file()
    assert ORBIT_SCRIPT.main(arguments) == []


def test_orbit_script_defaults_to_hardy_image_ratio():
    args = ORBIT_SCRIPT.parse_args([])

    assert args.proton_energy_model == "hardy"
    assert args.workers == 1
    assert args.base_output is None
    assert ORBIT_SCRIPT.PRECIPITATION_METHOD == "image_ratio"


def test_orbit_script_supports_separate_bases_and_parallel_runs(
    tmp_path, monkeypatch
):
    base_input = tmp_path / "input"
    base_output = tmp_path / "output"
    input_directory = (
        base_input / "fuv_detector" / PREPROCESSING_LABEL
    )
    input_directory.mkdir(parents=True)
    for orbit in (1, 2):
        (input_directory / f"or_{orbit:04d}.nc").touch()

    calls = []

    def fake_process_orbit(orbit, **settings):
        calls.append((orbit, settings))
        return orbit, 3

    def fake_process_map(function, orbits, **settings):
        assert settings["max_workers"] == 2
        assert settings["chunksize"] == 1
        return [function(orbit) for orbit in orbits]

    monkeypatch.setattr(ORBIT_SCRIPT, "process_orbit", fake_process_orbit)
    monkeypatch.setattr(ORBIT_SCRIPT, "process_map", fake_process_map)
    monkeypatch.setattr(ORBIT_SCRIPT, "load_gfz_kp", kp_series)
    monkeypatch.setattr(ORBIT_SCRIPT, "current_revision", lambda path: "test")

    result = ORBIT_SCRIPT.main([
        "--base-input", str(base_input),
        "--base-output", str(base_output),
        "--workers", "2",
    ])

    assert result == [(1, 3), (2, 3)]
    assert [call[0] for call in calls] == [1, 2]
    for _, settings in calls:
        assert settings["input_directory"] == input_directory
        assert settings["output_directory"] == (
            base_output / "precipitation_detector" / "IR_hardy"
        )
        assert settings["proton_energy_model"] == "hardy"


def test_orbit_script_rejects_invalid_worker_count():
    with pytest.raises(ValueError, match="workers must be at least 1"):
        ORBIT_SCRIPT.main(["--workers", "0"])


def test_orbit_worker_reports_failed_orbit(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise ValueError("broken input")

    monkeypatch.setattr(ORBIT_SCRIPT, "PrecipitationDetector", fail)

    with pytest.raises(
        RuntimeError, match="precipitation_detector orbit 0042 failed"
    ):
        ORBIT_SCRIPT.process_orbit(
            42,
            tmp_path,
            tmp_path,
            kp_series(),
            "hardy",
            2.0,
            0.0,
            "test",
        )
