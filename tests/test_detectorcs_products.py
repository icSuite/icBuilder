"""Integration tests for paired detector-to-CS post-processing."""

from collections import Counter
from datetime import datetime, timedelta
import importlib.util
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from icreader.product import NetCDFProduct
from icphysics import robinson_conductance
from netCDF4 import Dataset

from icbuilder.conductancedetector import ConductanceDetector
from icbuilder.detectorcs import build_detector_cs_products
from icbuilder.fuvdetector import FUVDetector, PREPROCESSING_LABEL
from icbuilder.grids import DETECTOR_CS_GRID_ID
from icbuilder.precipitationdetector import PrecipitationDetector


SCRIPT = (
    Path(__file__).parents[1] / "scripts" / "pipeline"
    / "make_detector_cs_orbit_files.py"
)
SPEC = importlib.util.spec_from_file_location(
    "make_detector_cs_orbit_files", SCRIPT
)
ORBIT_SCRIPT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ORBIT_SCRIPT)

BASE_TIME = datetime(2001, 1, 1, 1)


def regular_camera(sensor, size, spacing, times, value):
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


def write_detector_pair(base, orbit=1, nonlinear=False):
    precipitation_directory = base / "precipitation_detector" / "IR_hardy"
    conductance_directory = (
        base / "conductance_detector" / "IR_hardy" / "robinson"
    )
    fuv_directory = base / "fuv_detector" / PREPROCESSING_LABEL
    precipitation_directory.mkdir(parents=True, exist_ok=True)
    conductance_directory.mkdir(parents=True, exist_ok=True)
    fuv_directory.mkdir(parents=True, exist_ok=True)

    wic_times = [BASE_TIME, BASE_TIME + timedelta(seconds=10)]
    si_times = [time + timedelta(seconds=1) for time in wic_times]
    wic = regular_camera("WIC", 8, 0.05, wic_times, 1000.0)
    si12 = regular_camera("SI12", 4, 0.10, si_times, 10.0)
    si13 = regular_camera("SI13", 4, 0.10, si_times, 20.0)

    fuv_filename = fuv_directory / f"or_{orbit:04d}.nc"
    FUVDetector(wic, si12, si13, software_version="product1-test").to_nc(
        fuv_filename
    )
    precipitation_filename = precipitation_directory / f"or_{orbit:04d}.nc"
    PrecipitationDetector(
        fuv_filename,
        kp_series=kp_series(),
        proton_energy_model="constant",
        proton_energy=2.0,
        software_version="product2-test",
    ).to_nc(precipitation_filename)

    if nonlinear:
        with Dataset(precipitation_filename, "r+") as nc:
            valid = nc["method_valid"][:].astype(bool)
            columns = np.indices(valid.shape)[2]
            E0 = np.where(columns < valid.shape[2] / 2, 1.0, 8.0)
            Fe = np.where(columns < valid.shape[2] / 2, 1.0, 9.0)
            nc["E0"][:] = np.where(valid, E0, np.nan)
            nc["Fe"][:] = np.where(valid, Fe, np.nan)

    conductance_filename = conductance_directory / f"or_{orbit:04d}.nc"
    ConductanceDetector(
        precipitation_filename,
        software_version="product3-test",
    ).to_nc(conductance_filename)
    return precipitation_filename, conductance_filename


def test_paired_products_share_mapping_and_round_trip(tmp_path, monkeypatch):
    precipitation_source, conductance_source = write_detector_pair(
        tmp_path / "input", nonlinear=True
    )

    calls = []
    real_mapping = __import__(
        "icbuilder.detectorcs", fromlist=["make_detector_cs_mapping"]
    ).make_detector_cs_mapping

    def count_mapping(*args, **kwargs):
        calls.append(1)
        return real_mapping(*args, **kwargs)

    monkeypatch.setattr(
        "icbuilder.detectorcs.make_detector_cs_mapping", count_mapping
    )
    precipitation, conductance = build_detector_cs_products(
        precipitation_source,
        conductance_source,
        software_version="cs-test",
    )

    assert precipitation.shape == (2, 46, 46)
    assert conductance.shape == (2, 46, 46)
    assert len(calls) == 2
    np.testing.assert_allclose(
        precipitation.E0, conductance.E0, equal_nan=True
    )
    assert precipitation.method_valid.any()
    assert conductance.conductance_valid.any()

    recalculated = robinson_conductance(
        conductance.E0,
        conductance.Fe,
        conductance.dE0,
        conductance.dFe,
        conductance.varE0Fe,
    )["P"]
    finite = np.isfinite(recalculated) & np.isfinite(conductance.P)
    assert np.any(np.abs(recalculated[finite] - conductance.P[finite]) > 1e-5)

    precipitation_output = tmp_path / "precipitation_cs.nc"
    conductance_output = tmp_path / "conductance_cs.nc"
    conductance.companion_precipitation_cs = str(precipitation_output)
    precipitation.to_nc(precipitation_output)
    conductance.to_nc(conductance_output)

    assert ORBIT_SCRIPT.precipitation_cs_file_status(
        precipitation_output, precipitation_source
    ) == "complete"
    assert ORBIT_SCRIPT.conductance_cs_file_status(
        conductance_output, conductance_source, precipitation_output
    ) == "complete"
    with Dataset(conductance_output) as nc:
        assert nc.product_type == "conductance_cs"
        assert nc.representation == "cs"
        assert nc.grid_id == DETECTOR_CS_GRID_ID
        assert nc.variables["P"].shape == (2, 46, 46)
        assert "not recalculated" in nc.nonlinear_ordering
        assert nc.groups["grid"].variables["mlat"].shape == (46, 46)


def test_paired_reduction_reads_each_detector_cube_once(tmp_path, monkeypatch):
    precipitation_source, conductance_source = write_detector_pair(
        tmp_path / "input"
    )
    reads = Counter()
    original_read = NetCDFProduct.read

    def counted_read(product, name, index=None):
        if index is not None:
            raise AssertionError("CS reduction must not decompress frame by frame")
        reads[(product.product_type, name)] += 1
        return original_read(product, name, index)

    monkeypatch.setattr(NetCDFProduct, "read", counted_read)
    build_detector_cs_products(
        precipitation_source,
        conductance_source,
        software_version="cs-test",
    )

    assert reads
    assert max(reads.values()) == 1


def test_combined_orbit_runner_writes_and_restarts(tmp_path, monkeypatch):
    base_input = tmp_path / "input"
    base_output = tmp_path / "output"
    write_detector_pair(base_input)
    monkeypatch.setattr(ORBIT_SCRIPT, "current_revision", lambda path: "test")

    arguments = [
        "--base-input", str(base_input),
        "--base-output", str(base_output),
        "--orbit", "1",
    ]
    assert ORBIT_SCRIPT.main(arguments) == [(1, 2)]

    precipitation_output = (
        base_output / "precipitation_cs" / DETECTOR_CS_GRID_ID
        / "IR_hardy" / "or_0001.nc"
    )
    conductance_output = (
        base_output / "conductance_cs" / DETECTOR_CS_GRID_ID
        / "IR_hardy" / "robinson" / "or_0001.nc"
    )
    assert precipitation_output.is_file()
    assert conductance_output.is_file()
    assert ORBIT_SCRIPT.main(arguments) == []

    with Dataset(precipitation_output, "r+") as nc:
        nc.schema_version = 0
    assert ORBIT_SCRIPT.main(arguments) == [(1, 2)]


def test_combined_orbit_runner_parallel_routing(tmp_path, monkeypatch):
    base_input = tmp_path / "input"
    for directory in (
        base_input / "precipitation_detector" / "IR_hardy",
        base_input / "conductance_detector" / "IR_hardy" / "robinson",
    ):
        directory.mkdir(parents=True)
        for orbit in (1, 2):
            (directory / f"or_{orbit:04d}.nc").touch()

    calls = []

    def fake_process_orbit(task, **settings):
        calls.append((task, settings))
        return task[0], 3

    def fake_process_map(function, tasks, **settings):
        assert settings["max_workers"] == 2
        return [function(task) for task in tasks]

    monkeypatch.setattr(ORBIT_SCRIPT, "process_orbit", fake_process_orbit)
    monkeypatch.setattr(ORBIT_SCRIPT, "process_map", fake_process_map)
    monkeypatch.setattr(ORBIT_SCRIPT, "current_revision", lambda path: "test")

    result = ORBIT_SCRIPT.main([
        "--base-input", str(base_input),
        "--base-output", str(tmp_path / "output"),
        "--workers", "2",
    ])

    assert result == [(1, 3), (2, 3)]
    assert [call[0][0] for call in calls] == [1, 2]


def test_paired_builder_rejects_changed_product2(tmp_path):
    precipitation_source, conductance_source = write_detector_pair(
        tmp_path / "input"
    )
    with Dataset(precipitation_source, "r+") as nc:
        nc.audit_change = "changed after detector Product 3 was written"

    with pytest.raises(ValueError, match="does not identify"):
        build_detector_cs_products(
            precipitation_source,
            conductance_source,
            software_version="cs-test",
        )


def test_combined_orbit_runner_rejects_invalid_worker_count():
    with pytest.raises(ValueError, match="workers must be at least 1"):
        ORBIT_SCRIPT.main(["--workers", "0"])


def test_combined_orbit_runner_executes_two_real_workers(tmp_path):
    base_input = tmp_path / "input"
    base_output = tmp_path / "output"
    write_detector_pair(base_input, orbit=1)
    write_detector_pair(base_input, orbit=2)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--base-input", str(base_input),
            "--base-output", str(base_output),
            "--workers", "2",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "0 complete, 2 pending" in result.stdout
    for orbit in (1, 2):
        assert (
            base_output / "precipitation_cs" / DETECTOR_CS_GRID_ID
            / "IR_hardy" / f"or_{orbit:04d}.nc"
        ).is_file()
        assert (
            base_output / "conductance_cs" / DETECTOR_CS_GRID_ID
            / "IR_hardy" / "robinson" / f"or_{orbit:04d}.nc"
        ).is_file()
