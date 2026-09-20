"""Focused tests for detector-space conductance Product 3."""

#%% Imports and fixtures

from datetime import datetime
import importlib.util
from pathlib import Path

import numpy as np
import pytest
from icphysics import robinson_conductance
from netCDF4 import Dataset, date2num

from icbuilder.conductancedetector import (
    CONDUCTANCE_MODEL,
    SCHEMA_VERSION,
    ConductanceDetector,
)
from icbuilder.fuvdetector import SOURCE_TIME_DECODING
from icbuilder.precipitationdetector import (
    PRECIPITATION_METHOD,
    SCHEMA_VERSION as PRECIPITATION_SCHEMA_VERSION,
)


SCRIPT = (
    Path(__file__).parents[1] / "scripts" / "pipeline"
    / "make_conductance_detector_orbit_files.py"
)
SPEC = importlib.util.spec_from_file_location(
    "make_conductance_detector_orbit_files", SCRIPT
)
ORBIT_SCRIPT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ORBIT_SCRIPT)


def write_precipitation_detector(path, proton_energy_model="hardy"):
    """Write the small schema-2 Product-2 input used by focused tests."""

    shape = (1, 2, 2)
    method_valid = np.array([[[True, True], [True, False]]])
    with Dataset(path, "w") as nc:
        nc.createDimension("time", shape[0])
        nc.createDimension("row", shape[1])
        nc.createDimension("column", shape[2])

        nc.product_type = "precipitation_detector"
        nc.representation = "detector"
        nc.schema_version = PRECIPITATION_SCHEMA_VERSION
        nc.method = PRECIPITATION_METHOD
        nc.proton_flux_source = "SI12"
        nc.proton_energy_model = proton_energy_model
        nc.proton_energy_uncertainty_method = "test uncertainty"
        nc.proton_energy_coordinate_note = "test coordinates"
        nc.proton_response_energy_min = 0.47
        nc.proton_response_energy_max = 46.7
        nc.proton_operation_order = "test operation order"
        nc.count_uncertainty_method = "test count uncertainty"
        nc.source_fuv_detector = "fuv_detector/or_0001.nc"
        nc.source_fuv_detector_sha256 = "source-fuv-sha256"
        nc.source_preprocessing_label = "fuvpy_bs_directional_v1"
        nc.source_fuv_detector_time_decoding = SOURCE_TIME_DECODING
        nc.coordinate_system = "Modified Apex"
        nc.reference_height_km = 130.0
        nc.software_version = "product2-test"
        nc.kp_source = "test"
        if proton_energy_model == "constant":
            nc.proton_energy_constant = 2.0
            nc.proton_energy_uncertainty_constant = 0.0

        time_units = "seconds since 2000-01-01 00:00:00"
        time_value = date2num(
            [datetime(2001, 1, 1, 1)], time_units, calendar="standard"
        )
        for name in (
            "time", "wic_source_time", "si12_source_time",
            "si13_source_time", "Kp_interval_start",
        ):
            variable = nc.createVariable(name, "f8", ("time",))
            variable[:] = time_value
            variable.units = time_units
            variable.calendar = "standard"

        for name in (
            "wic_source_index", "si12_source_index", "si13_source_index"
        ):
            nc.createVariable(name, "i4", ("time",))[:] = [0]
        for name in (
            "wic_frame_quality", "si12_frame_quality", "si13_frame_quality"
        ):
            nc.createVariable(name, "i1", ("time",))[:] = [2]
        nc.createVariable("Kp", "f4", ("time",))[:] = [2.0]
        nc.createVariable("ssalon", "f4", ("time",))[:] = [10.0]
        nc.createVariable("detector_row", "i4", ("row",))[:] = [0, 1]
        nc.createVariable("detector_column", "i4", ("column",))[:] = [0, 1]

        dimensions = ("time", "row", "column")
        geometry = {
            "glat": 70.0,
            "glon": 0.0,
            "mlat": 69.0,
            "mlon": 1.0,
            "mlt": 12.0,
            "sza": 80.0,
            "dza": 20.0,
            "method_quality_weight": 0.7,
            "Ep_model": 1.5,
            "Ep": 1.5,
            "dEp": 0.0,
            "Fp": 0.8,
            "dFp": 0.1,
        }
        for name, value in geometry.items():
            nc.createVariable(name, "f4", dimensions)[:] = np.full(
                shape, value
            )

        precipitation = {
            "E0": np.array([[[2.0, 2.0], [2.0, np.nan]]]),
            "dE0": np.array([[[0.4, 0.4], [np.nan, np.nan]]]),
            "Fe": np.array([[[3.0, 0.0], [4.0, np.nan]]]),
            "dFe": np.array([[[0.5, 0.5], [np.nan, np.nan]]]),
            "varE0Fe": np.array([[[0.0, 0.0], [0.0, np.nan]]]),
        }
        for name, values in precipitation.items():
            nc.createVariable(name, "f4", dimensions)[:] = values

        nc.createVariable("method_valid", "i1", dimensions)[:] = (
            method_valid.astype(np.int8)
        )
        nc.createVariable("Ep_clipping_flag", "i1", dimensions)[:] = (
            np.zeros(shape, dtype=np.int8)
        )

    return method_valid


#%% Calculation and schema

def test_detector_conductance_matches_icphysics_and_keeps_central_values(
    tmp_path
):
    source = tmp_path / "precipitation.nc"
    method_valid = write_precipitation_detector(source)

    product = ConductanceDetector(source, software_version="product3-test")
    expected = robinson_conductance(
        product.E0, product.Fe, product.dE0, product.dFe, product.varE0Fe
    )

    np.testing.assert_allclose(product.P, expected["P"], equal_nan=True)
    np.testing.assert_allclose(product.H, expected["H"], equal_nan=True)
    np.testing.assert_allclose(product.dP, expected["dP"], equal_nan=True)
    np.testing.assert_allclose(product.dH, expected["dH"], equal_nan=True)
    np.testing.assert_array_equal(product.conductance_valid, method_valid)
    np.testing.assert_array_equal(
        product.conductance_uncertainty_valid,
        np.array([[[True, True], [False, False]]]),
    )
    assert np.isfinite(product.P[0, 1, 0])
    assert np.isnan(product.dP[0, 1, 0])
    assert product.dP[0, 0, 1] > 0


def test_detector_conductance_netcdf_is_self_describing(tmp_path):
    source = tmp_path / "precipitation.nc"
    write_precipitation_detector(source)
    product = ConductanceDetector(source, software_version="product3-test")
    output = tmp_path / "conductance.nc"

    ORBIT_SCRIPT.save_conductance_detector(product, output)

    assert ORBIT_SCRIPT.conductance_detector_file_status(
        output, source
    ) == "complete"
    assert not Path(str(output) + ".partial").exists()
    with Dataset(output) as nc:
        assert nc.product_type == "conductance_detector"
        assert nc.representation == "detector"
        assert nc.schema_version == SCHEMA_VERSION
        assert nc.conductance_model == CONDUCTANCE_MODEL
        assert nc.precipitation_method == PRECIPITATION_METHOD
        assert nc.source_precipitation_detector == str(source)
        assert nc.source_precipitation_detector_schema_version == 2
        assert nc.source_precipitation_software_version == "product2-test"
        assert nc.software_version == "product3-test"
        assert "one-sided" in nc.conductance_uncertainty_method
        assert nc.variables["P"].units == "S"
        assert nc.variables["method_quality_weight"].units == "1"
        assert nc.variables["conductance_valid"].shape == (1, 2, 2)

    broken = tmp_path / "broken.nc"
    product.to_nc(broken)
    with Dataset(broken, "r+") as nc:
        nc.renameVariable("conductance_valid", "removed_conductance_valid")
    assert ORBIT_SCRIPT.conductance_detector_file_status(
        broken, source
    ) == "invalid"

    with Dataset(source, "r+") as nc:
        nc.audit_change = "source changed after Product 3 was written"
    assert ORBIT_SCRIPT.conductance_detector_file_status(
        output, source
    ) == "mismatch"


def test_detector_conductance_rejects_wrong_product2_schema(tmp_path):
    source = tmp_path / "precipitation.nc"
    write_precipitation_detector(source)
    with Dataset(source, "r+") as nc:
        nc.schema_version = 1

    with pytest.raises(ValueError, match="supported schema-2"):
        ConductanceDetector(source)


#%% Orbit runner

def test_orbit_script_writes_and_restarts_with_separate_bases(
    tmp_path, monkeypatch
):
    base_input = tmp_path / "input"
    base_output = tmp_path / "output"
    input_directory = base_input / "precipitation_detector" / "IR_hardy"
    input_directory.mkdir(parents=True)
    write_precipitation_detector(input_directory / "or_0001.nc")
    monkeypatch.setattr(ORBIT_SCRIPT, "current_revision", lambda path: "test")

    arguments = [
        "--base-input", str(base_input),
        "--base-output", str(base_output),
        "--orbit", "1",
    ]
    assert ORBIT_SCRIPT.main(arguments) == [(1, 1)]

    output = (
        base_output / "conductance_detector" / "IR_hardy" / "robinson"
        / "or_0001.nc"
    )
    assert output.is_file()
    assert ORBIT_SCRIPT.main(arguments) == []


def test_orbit_script_defaults_and_parallel_routing(tmp_path, monkeypatch):
    defaults = ORBIT_SCRIPT.parse_args([])
    assert defaults.retrieval_label == "IR_hardy"
    assert defaults.conductance_model == "robinson"
    assert defaults.workers == 1
    assert defaults.base_output is None

    base_input = tmp_path / "input"
    base_output = tmp_path / "output"
    input_directory = base_input / "precipitation_detector" / "IR_hardy"
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
            base_output / "conductance_detector" / "IR_hardy" / "robinson"
        )


def test_orbit_script_rejects_invalid_worker_count():
    with pytest.raises(ValueError, match="workers must be at least 1"):
        ORBIT_SCRIPT.main(["--workers", "0"])


def test_orbit_worker_reports_failed_orbit(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise ValueError("broken input")

    monkeypatch.setattr(ORBIT_SCRIPT, "ConductanceDetector", fail)

    with pytest.raises(
        RuntimeError, match="conductance_detector orbit 0042 failed"
    ):
        ORBIT_SCRIPT.process_orbit(42, tmp_path, tmp_path, "test")
