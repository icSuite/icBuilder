"""Focused tests for the detector-space Product-1 boundary."""

#%% Imports and test helpers

from datetime import datetime, timedelta
import importlib.util
from pathlib import Path
import shutil

import numpy as np
import pytest
from netCDF4 import Dataset, num2date
from scipy.sparse import csr_matrix

from icbuilder.detector_coregistration import map_si_variance
from icbuilder.fuvdetector import (
    FUVDetector,
    PREPROCESSING_LABEL,
    SCHEMA_VERSION,
    load_detector_source,
    match_wic_times,
    source_identity,
)


SCRIPT = (
    Path(__file__).parents[1] / "scripts" / "pipeline"
    / "make_fuv_detector_orbit_files.py"
)
SPEC = importlib.util.spec_from_file_location(
    "make_fuv_detector_orbit_files", SCRIPT
)
ORBIT_SCRIPT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ORBIT_SCRIPT)

BASE_TIME = datetime(2001, 1, 1)


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
        "time": np.asarray(times, dtype=object),
        "counts": np.full(frame_shape, value, dtype=float),
        "variance": np.full(frame_shape, 4.0, dtype=float),
        "frame_quality": np.full(len(times), 2, dtype=np.int8),
        "quality_weight": np.full(frame_shape, 0.8, dtype=float),
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


def write_source_orbit(path, sensor="WIC", missing=None):
    """Write the smallest new-fuvpy source orbit accepted by Product 1."""

    missing = set() if missing is None else set(missing)
    shape = (1, 2, 2)
    with Dataset(path, "w") as nc:
        nc.createDimension("date", 1)
        nc.createDimension("row", 2)
        nc.createDimension("col", 2)
        fields = {
            "dgimg": np.full(shape, 10.0),
            "dgweight": np.full(shape, 0.75),
            "img_variance": np.full(shape, 4.0),
            "glat": np.full(shape, 70.0),
            "glon": np.zeros(shape),
            "sza": np.full(shape, 80.0),
            "dza": np.full(shape, 20.0),
        }
        for name, values in fields.items():
            if name not in missing:
                nc.createVariable(name, "f8", ("date", "row", "col"))[:] = values
        if "frame_quality" not in missing:
            nc.createVariable("frame_quality", "i1", ("date",))[:] = [2]
        model = nc.createVariable("dgmodel", "f8", ("date", "row", "col"))
        model[:] = np.zeros(shape)
        model.spatial_model = "legacy_x_azimuth"
        model.frame_quality_filter = "frame_quality >= 2"
        model.measurement_variance = "img_variance"
        model.variance_weighted = "true"
        model.monotonic_requested = "true"
        model.monotonic_fallback = "true"
        model.spencer_taper = "enabled"
        model.temporal_knot_spacing_target_minutes = 120.0
        model.directional_temporal_knot_spacing_target_minutes = 60.0
        model.tukey_value = 5.0
        model.damping_value = 1e-2
        model.convergence_tolerance = 1e-2
        model.directional_x_knots = np.asarray([-3.5, 0.0, 1.5, 3.5])
        if "date" not in missing:
            date = nc.createVariable("date", "f8", ("date",))
            date[:] = [0.0]
            date.units = "seconds since 2001-01-01 00:00:00"
            date.calendar = "proleptic_gregorian"
        if "t_start" not in missing:
            nc.createVariable("t_start", str, ())[()] = BASE_TIME.isoformat()


def test_source_loader_uses_only_the_new_fuvpy_contract(tmp_path):
    path = tmp_path / "wic.nc"
    write_source_orbit(path)

    source = load_detector_source(path, "WIC")

    assert source["image_field"] == "dgimg"
    np.testing.assert_allclose(source["counts"], 10.0)
    np.testing.assert_allclose(source["quality_weight"], 0.75)
    np.testing.assert_allclose(source["variance"], 4.0)
    np.testing.assert_array_equal(source["frame_quality"], [2])


def test_source_loader_uses_cf_date_units_when_t_start_is_stale(tmp_path):
    path = tmp_path / "wic.nc"
    write_source_orbit(path)
    with Dataset(path, "r+") as nc:
        nc["t_start"][()] = "2000-12-31T21:39:17"

    source = load_detector_source(path, "WIC")

    assert source["time"].tolist() == [BASE_TIME]


@pytest.mark.parametrize(
    "missing", ["dgimg", "dgweight", "img_variance", "frame_quality"]
)
def test_source_loader_rejects_missing_schema2_fields(tmp_path, missing):
    path = tmp_path / f"missing_{missing}.nc"
    write_source_orbit(path, missing=[missing])

    with pytest.raises(ValueError, match=missing):
        load_detector_source(path, "WIC")


def test_source_loader_rejects_wrong_model_and_invalid_uncertainty(tmp_path):
    path = tmp_path / "wrong_model.nc"
    write_source_orbit(path)
    with Dataset(path, "r+") as nc:
        nc["dgmodel"].damping_value = 0.1

    with pytest.raises(ValueError, match="does not match"):
        load_detector_source(path, "WIC")

    path = tmp_path / "invalid_variance.nc"
    write_source_orbit(path)
    with Dataset(path, "r+") as nc:
        nc["img_variance"][0, 0, 0] = -1

    with pytest.raises(ValueError, match="negative"):
        load_detector_source(path, "WIC")

    path = tmp_path / "invalid_frame_quality.nc"
    write_source_orbit(path)
    with Dataset(path, "r+") as nc:
        nc["frame_quality"][:] = [3]

    with pytest.raises(ValueError, match="0, 1, or 2"):
        load_detector_source(path, "WIC")


#%% WIC-led matching

def test_wic_led_matching_is_nearest_deterministic_and_without_reuse():
    wic_times = [BASE_TIME, BASE_TIME + timedelta(seconds=1)]
    sensor_times = [
        BASE_TIME - timedelta(seconds=1),
        BASE_TIME + timedelta(seconds=1),
    ]

    # The first WIC frame has an equal-distance tie and selects the earlier
    # sensor frame. The second frame then uses the remaining exact match.
    np.testing.assert_array_equal(
        match_wic_times(wic_times, sensor_times), [0, 1]
    )

    np.testing.assert_array_equal(
        match_wic_times(wic_times, [BASE_TIME]), [0, -1]
    )


def test_wic_led_matching_keeps_unmatched_wic_frames():
    wic_times = [BASE_TIME, BASE_TIME + timedelta(seconds=120)]
    sensor_times = [BASE_TIME + timedelta(seconds=3)]

    np.testing.assert_array_equal(
        match_wic_times(wic_times, sensor_times), [-1, -1]
    )

    np.testing.assert_array_equal(
        match_wic_times(wic_times[:1], [BASE_TIME + timedelta(seconds=2)]),
        [0],
    )


def test_wic_led_matching_rejects_missing_times_and_declares_order_priority():
    with pytest.raises(ValueError, match="missing values"):
        match_wic_times([BASE_TIME], [np.datetime64("NaT")])

    # Stored WIC order has priority when one SI frame lies within tolerance of
    # two WIC frames. This is explicit experimental behavior, not a global fit.
    np.testing.assert_array_equal(
        match_wic_times(
            [BASE_TIME, BASE_TIME + timedelta(seconds=1.5)],
            [BASE_TIME + timedelta(seconds=1)],
        ),
        [0, -1],
    )


#%% Product construction and serialization

def test_detector_product_coregisters_each_si_channel_independently():
    wic_times = [BASE_TIME, BASE_TIME + timedelta(seconds=120)]
    wic = regular_camera("WIC", 16, 0.05, wic_times, 10.0)
    si12 = regular_camera(
        "SI12", 8, 0.10, [BASE_TIME + timedelta(seconds=1)], 5.0
    )
    si13 = regular_camera(
        "SI13", 8, 0.10,
        [BASE_TIME + timedelta(seconds=121)], 7.0,
    )

    product = FUVDetector(wic, si12, si13, software_version="test")

    assert product.shape == (2, 16, 16)
    np.testing.assert_array_equal(product.si12_source_index, [0, -1])
    np.testing.assert_array_equal(product.si13_source_index, [-1, 0])
    np.testing.assert_allclose(product.si12_counts[0], 5.0)
    np.testing.assert_allclose(product.si13_counts[1], 7.0)
    assert np.isnan(product.si12_counts[1]).all()
    assert np.isnan(product.si13_counts[0]).all()
    assert np.all(product.si12_coverage[0] >= 0.9)
    assert np.all(product.si13_coverage[1] >= 0.9)
    assert np.all(product.si12_source_count[0] >= 1)
    assert np.all(product.si13_source_count[1] >= 1)
    np.testing.assert_allclose(product.wic_variance, 4.0)
    assert np.nanmin(product.si12_variance[0]) > 3.9
    assert np.nanmax(product.si12_variance[0]) <= 4.0
    assert np.nanmin(product.si13_variance[1]) > 3.9
    assert np.nanmax(product.si13_variance[1]) <= 4.0
    np.testing.assert_array_equal(product.wic_frame_quality, [2, 2])
    np.testing.assert_array_equal(product.si12_frame_quality, [2, -1])
    np.testing.assert_array_equal(product.si13_frame_quality, [-1, 2])


def test_si_variance_uses_squared_overlap_weights_and_variance_support():
    mapping = csr_matrix([[0.5, 0.5]])
    variance, coverage = map_si_variance(
        np.array([[4.0, 16.0]]),
        np.array([[True, True]]),
        mapping,
        (1, 1),
    )
    np.testing.assert_allclose(variance, [[5.0]])
    np.testing.assert_allclose(coverage, [[1.0]])

    variance, coverage = map_si_variance(
        np.array([[4.0, np.nan]]),
        np.array([[True, True]]),
        mapping,
        (1, 1),
    )
    assert np.isnan(variance[0, 0])
    np.testing.assert_allclose(coverage, [[0.5]])


def test_detector_product_netcdf_is_self_describing_and_restart_safe(tmp_path):
    wic_file = tmp_path / "wic.nc"
    si12_file = tmp_path / "si12.nc"
    wic_file.write_bytes(b"wic source")
    si12_file.write_bytes(b"si12 source")

    wic_times = [BASE_TIME, BASE_TIME + timedelta(seconds=120)]
    wic = regular_camera("WIC", 16, 0.05, wic_times, 10.0)
    si12 = regular_camera(
        "SI12", 8, 0.10, [BASE_TIME + timedelta(seconds=1)], 5.0
    )
    wic.update(source_identity(wic_file))
    si12.update(source_identity(si12_file))
    product = FUVDetector(wic, si12, software_version="test")
    source_files = {
        "wic": wic_file,
        "si12": si12_file,
        "si13": None,
    }
    output = tmp_path / "or_0001.nc"

    ORBIT_SCRIPT.save_fuv_detector_file(product, output, source_files)

    assert ORBIT_SCRIPT.fuv_detector_file_status(
        output, PREPROCESSING_LABEL, source_files
    ) == "complete"
    assert ORBIT_SCRIPT.fuv_detector_file_status(
        output, "different_preprocessing", source_files
    ) == "mismatch"
    assert not Path(str(output) + ".partial").exists()

    with Dataset(output) as nc:
        assert nc.product_type == "fuv_detector"
        assert nc.representation == "detector"
        assert nc.schema_version == SCHEMA_VERSION
        assert nc.preprocessing_label == PREPROCESSING_LABEL
        assert nc.source_time_decoding == "CF date units and calendar"
        assert nc.wic_image_field == "dgimg"
        assert nc.si12_image_field == "dgimg"
        assert nc.si13_image_field == "dgimg"
        assert "independent detector-pixel" in nc.detector_noise_model
        assert nc.coregistration_overlap_operator_stored == 0
        assert "Kp" not in nc.variables
        for name in ("dE0", "dFe", "dP", "dH", "wic_uncertainty"):
            assert name not in nc.variables

        np.testing.assert_array_equal(nc["wic_source_index"][:], [0, 1])
        np.testing.assert_array_equal(nc["si12_source_index"][:], [0, -1])
        np.testing.assert_array_equal(nc["si13_source_index"][:], [-1, -1])
        np.testing.assert_array_equal(nc["wic_frame_quality"][:], [2, 2])
        np.testing.assert_array_equal(nc["si12_frame_quality"][:], [2, -1])
        assert np.ma.getmaskarray(nc["si13_frame_quality"][:]).all()
        np.testing.assert_allclose(nc["wic_variance"][:], 4.0)
        np.testing.assert_allclose(nc["si12_counts"][0], 5.0)
        assert np.isnan(nc["si12_counts"][1]).all()
        assert np.ma.getmaskarray(nc["si13_source_time"][:]).all()

        decoded = num2date(
            nc["time"][:], nc["time"].units, nc["time"].calendar,
            only_use_cftime_datetimes=False,
        )
        assert decoded.tolist() == wic_times

    documentation_edit = tmp_path / "documentation_edit.nc"
    shutil.copy2(output, documentation_edit)
    with Dataset(documentation_edit, "r+") as nc:
        nc.time_match_rule = "Clearer documentation of the same calculation."
        nc.quality_weight_method = "Clearer quality-weight documentation."
        nc.detector_noise_model = "Clearer detector-noise documentation."
    assert ORBIT_SCRIPT.fuv_detector_file_status(
        documentation_edit, PREPROCESSING_LABEL, source_files
    ) == "complete"

    wrong_units = tmp_path / "wrong_units.nc"
    shutil.copy2(output, wrong_units)
    with Dataset(wrong_units, "r+") as nc:
        nc["wic_counts"].units = "not counts"
    assert ORBIT_SCRIPT.fuv_detector_file_status(
        wrong_units, PREPROCESSING_LABEL, source_files
    ) == "invalid"

    wrong_flags = tmp_path / "wrong_flags.nc"
    shutil.copy2(output, wrong_flags)
    with Dataset(wrong_flags, "r+") as nc:
        nc["wic_frame_quality"].flag_meanings = "wrong"
    assert ORBIT_SCRIPT.fuv_detector_file_status(
        wrong_flags, PREPROCESSING_LABEL, source_files
    ) == "invalid"

    old_time_decoder = tmp_path / "old_time_decoder.nc"
    shutil.copy2(output, old_time_decoder)
    with Dataset(old_time_decoder, "r+") as nc:
        nc.delncattr("source_time_decoding")
    assert ORBIT_SCRIPT.fuv_detector_file_status(
        old_time_decoder, PREPROCESSING_LABEL, source_files
    ) == "invalid"

    missing_diagnostic = tmp_path / "missing_diagnostic.nc"
    shutil.copy2(output, missing_diagnostic)
    with Dataset(missing_diagnostic, "r+") as nc:
        nc.renameVariable(
            "si12_coreg_coverage_maximum", "removed_coverage_maximum"
        )
    assert ORBIT_SCRIPT.fuv_detector_file_status(
        missing_diagnostic, PREPROCESSING_LABEL, source_files
    ) == "invalid"

    changed_sources = source_files | {"si12": tmp_path / "different_si12.nc"}
    assert ORBIT_SCRIPT.fuv_detector_file_status(
        output, PREPROCESSING_LABEL, changed_sources
    ) == "mismatch"


def test_product_status_rejects_missing_required_field(tmp_path):
    path = tmp_path / "broken.nc"
    source_files = {"wic": Path("wic.nc"), "si12": None, "si13": None}
    with Dataset(path, "w") as nc:
        nc.product_type = "fuv_detector"
        nc.representation = "detector"
        nc.schema_version = SCHEMA_VERSION
        nc.preprocessing_label = PREPROCESSING_LABEL
        nc.time_match_tolerance_seconds = 2.0
        for sensor, source in source_files.items():
            nc.setncattr(f"source_{sensor}", "" if source is None else str(source))
            nc.setncattr(
                f"{sensor}_image_field",
                "dgimg",
            )
        nc.createDimension("time", 1)
        nc.createDimension("row", 2)
        nc.createDimension("column", 2)

    assert ORBIT_SCRIPT.fuv_detector_file_status(
        path, PREPROCESSING_LABEL, source_files
    ) == "invalid"


def test_orbit_script_does_not_allow_a_label_to_overstate_preprocessing():
    assert ORBIT_SCRIPT.validate_label(PREPROCESSING_LABEL) == PREPROCESSING_LABEL
    with pytest.raises(ValueError, match="explicit preprocessing branch"):
        ORBIT_SCRIPT.validate_label("historical_fuview")

    wic = regular_camera("WIC", 4, 0.05, [BASE_TIME], 10.0)
    with pytest.raises(ValueError, match="explicit preprocessing branch"):
        FUVDetector(wic, preprocessing_label="historical_fuview")


def test_orbit_script_can_process_orbits_in_parallel(tmp_path, monkeypatch):
    wic_directory = tmp_path / "wic"
    wic_directory.mkdir()
    output_base = tmp_path / "output"
    for orbit in (1, 2):
        (wic_directory / f"wic_or{orbit:04d}.nc").touch()

    calls = []

    def record_orbit(orbit, **settings):
        calls.append((orbit, settings))
        return orbit, 1

    def record_process_map(function, orbits, **settings):
        assert settings["max_workers"] == 2
        assert settings["chunksize"] == 1
        return [function(orbit) for orbit in orbits]

    monkeypatch.setattr(ORBIT_SCRIPT, "process_orbit", record_orbit)
    monkeypatch.setattr(ORBIT_SCRIPT, "process_map", record_process_map)

    result = ORBIT_SCRIPT.main([
        "--base-input", str(tmp_path),
        "--base-output", str(output_base),
        "--workers", "2",
    ])

    assert result == [(1, 1), (2, 1)]
    assert [orbit for orbit, _ in calls] == [1, 2]
    assert all(
        settings["base"] == tmp_path
        and settings["output_directory"]
        == output_base / "fuv_detector" / PREPROCESSING_LABEL
        for _, settings in calls
    )


def test_orbit_script_rejects_nonpositive_worker_count():
    with pytest.raises(ValueError, match="workers must be at least 1"):
        ORBIT_SCRIPT.main(["--workers", "0"])


def test_worker_failure_identifies_the_orbit(tmp_path, monkeypatch):
    wic_directory = tmp_path / "wic"
    wic_directory.mkdir()
    (wic_directory / "wic_or0042.nc").touch()

    def fail(*args, **kwargs):
        raise ValueError("bad geometry")

    monkeypatch.setattr(ORBIT_SCRIPT.FUVDetector, "from_files", fail)

    with pytest.raises(RuntimeError, match="orbit 0042 failed"):
        ORBIT_SCRIPT.process_orbit(
            42,
            tmp_path,
            tmp_path / "output",
            PREPROCESSING_LABEL,
            "test",
        )
