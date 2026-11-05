"""Focused tests for the FAST--IMAGE conjunction dataset builder."""

import importlib.util
from pathlib import Path
import sys
import types

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix


try:
    import icreader  # noqa: F401
except ImportError:
    icreader = types.ModuleType("icreader")
    icreader.open_product = None
    sys.modules["icreader"] = icreader

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts/fast/06_generate_conjunction_dataset.py"
)
SPECIFICATION = importlib.util.spec_from_file_location(
    "fast_conjunction_dataset", SCRIPT
)
BUILDER = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(BUILDER)


class GateGrid:
    """Small deterministic CS grid substitute for support tests."""

    def ingrid(self, longitude, latitude):
        return np.asarray(longitude) < 300

    def bin_index(self, longitude, latitude):
        longitude = np.asarray(longitude)
        return np.zeros(longitude.size, dtype=int), (longitude >= 100).astype(int)


def test_temporal_support_is_closed_at_both_edges():
    centre = 1_000_000_000_000
    times = centre + np.array([-61, -60, 0, 60, 61], dtype=np.int64) * 10**9

    selected = BUILDER.temporal_candidate_indices(times, centre, 60)

    np.testing.assert_array_equal(selected, [1, 2, 3])


def test_cs_gate_requires_grid_coverage_and_finite_binned_data():
    shape = (1, 1, 2)
    wic_coverage = np.ones(shape)
    si13_coverage = np.ones(shape)
    wic = np.ones(shape)
    si13 = np.ones(shape)
    si13_coverage[0, 0, 1] = 0

    selected, values, counts = BUILDER.cs_support_gate(
        GateGrid(), 0,
        np.array([70.0, 71.0, 72.0, np.nan]),
        np.array([1.0, 8.0, 23.0, 1.0]),
        wic_coverage, si13_coverage, wic, si13,
    )

    np.testing.assert_array_equal(selected, [0])
    np.testing.assert_array_equal(values["row"], [0])
    np.testing.assert_array_equal(values["column"], [0])
    assert counts == {
        "nonfinite_footprint_count": 1,
        "outside_cs_grid_count": 1,
        "cs_coverage_data_rejection_count": 1,
    }


class IdentityProjection:
    def geo2cube(self, longitude, latitude, set_points_off_cube_to_nan=True):
        return np.asarray(longitude, dtype=float), np.asarray(latitude, dtype=float)


class DetectorGrid:
    shape = (1, 1)
    projection = IdentityProjection()


def test_detector_match_chooses_nearest_when_footprints_overlap(monkeypatch):
    corners = np.array([
        [[-1, -1], [1, -1], [1, 1], [-1, 1]],
        [[-2, -2], [2, -2], [2, 2], [-2, 2]],
    ], dtype=float)
    centres = np.array([[0.0, 0.0], [0.8, 0.0]])
    mapping = csr_matrix(np.array([[1, 1]], dtype=float))
    matcher = corners, centres, mapping
    monkeypatch.setattr(
        BUILDER, "candidate_footprints", lambda *args: np.array([0, 1])
    )

    row, column, separation, count = BUILDER.match_detector_footprints(
        matcher, DetectorGrid(), np.array([0]), np.array([0]),
        np.array([0.0]), np.array([0.0]),
        np.array([[0.0, 0.0]]), np.array([[0.0, 0.0]]),
    )

    np.testing.assert_array_equal(row, [0])
    np.testing.assert_array_equal(column, [0])
    np.testing.assert_array_equal(count, [2])
    assert separation.item() == pytest.approx(0)


class FakeCS:
    time = np.array(["2000-01-01T00:00:00"], dtype="datetime64[ns]")
    grid = GateGrid()

    def read(self, name):
        return np.ones((1, 1, 2), dtype=float)


class FakePrecipitation:
    def read(self, name):
        if name in ("mlat", "mlt"):
            return np.zeros((1, 1, 2), dtype=float)
        if name == "wic_corrected":
            return np.array([[[10.0, 20.0]]])
        if name == "si13_corrected":
            return np.array([[[5.0, np.nan]]])
        raise KeyError(name)


def test_exact_detector_pixel_without_image_pair_is_not_written(monkeypatch):
    monkeypatch.setattr(BUILDER, "make_detector_footprint_matcher", lambda *args: None)
    monkeypatch.setattr(
        BUILDER,
        "match_detector_footprints",
        lambda *args: (
            np.array([0, 0], dtype=np.int16),
            np.array([0, 1], dtype=np.int16),
            np.array([0.1, 0.2], dtype=np.float32),
            np.array([1, 1], dtype=np.int16),
        ),
    )
    fast = {
        "source_file_id": 0,
        "record": np.array([4, 5], dtype=np.int32),
        "time_ns": BUILDER.datetime64_ns(np.array([
            "2000-01-01T00:00:00", "2000-01-01T00:00:01"
        ], dtype="datetime64[ns]")),
        "mlat": np.array([70.0, 70.0]),
        "mlt": np.array([1.0, 8.0]),
    }
    counts = BUILDER.empty_counts()

    candidates = BUILDER.discover_orbit_candidates(
        85, FakeCS(), FakePrecipitation(), [fast], 60, counts
    )

    assert len(candidates) == 1
    assert candidates[0].fast_source_record.tolist() == [4]
    assert counts["exact_detector_data_miss_count"] == 1
    assert counts["successful_candidate_count"] == 1


def test_unique_match_uses_supported_nearest_frame_and_deterministic_ties():
    frame = pd.DataFrame({
        "fast_source_file_id": [0, 0, 0, 1],
        "fast_source_record": [7, 7, 7, 2],
        "fast_time_ns": [100, 100, 100, 50],
        "image_orbit": [86, 85, 85, 84],
        "image_frame": [3, 4, 2, 1],
        "image_time_ns": [110, 90, 90, 40],
        "time_offset_seconds": [-10.0, 10.0, 10.0, -10.0],
        "detector_row": [1, 1, 0, 1],
        "detector_column": [1, 1, 2, 1],
    })
    for name in BUILDER.CANDIDATE_COLUMNS:
        if name not in frame:
            frame[name] = 1

    selected = BUILDER.select_unique_matches(frame)

    assert len(selected) == 2
    source_zero = selected.loc[selected.fast_source_file_id.eq(0)].iloc[0]
    assert source_zero.image_orbit == 85
    assert source_zero.image_frame == 2
    assert source_zero.successful_frame_candidate_count == 3


class FakeProduct:
    def __init__(self, product_type, filename, attrs, times, shape=(2, 3, 4)):
        self.product_type = product_type
        self.filename = str(filename)
        self.attrs = attrs
        self.time = times
        self.shape = shape
        self.reference_height_km = 130.0
        self.grid_id = BUILDER.GRID_ID
        self.method = attrs.get("method", "image_ratio")


def test_product_source_mismatch_fails_before_matching(tmp_path):
    times = np.array(["2000-01-01"], dtype="datetime64[ns]")
    p2_path = tmp_path / "or_0085.nc"
    cs = FakeProduct(
        "precipitation_cs", tmp_path / "cs.nc",
        {
            "source_precipitation_detector_sha256": "one",
            "source_precipitation_detector": str(p2_path),
        }, times,
    )
    p2 = FakeProduct("precipitation_detector", p2_path, {}, times)
    p3 = FakeProduct(
        "conductance_detector", tmp_path / "p3.nc",
        {
            "source_precipitation_detector_sha256": "two",
            "source_precipitation_detector": str(p2_path),
        }, times,
    )

    with pytest.raises(ValueError, match="different Product 2"):
        BUILDER.validate_orbit_products(85, cs, p2, p3)


def test_atomic_schema_writer_preserves_unique_source_identity(tmp_path):
    output = tmp_path / "conjunctions.nc"
    signatures = {
        "source_fast_collection_signature": "fast",
        "source_precipitation_cs_collection_signature": "cs",
        "source_precipitation_detector_collection_signature": "p2",
        "source_conductance_detector_collection_signature": "p3",
    }
    attrs = {
        "product_type": "fast_image_conjunctions",
        "schema_version": 1,
        **signatures,
    }
    with BUILDER.create_output(
        output, 2, 3, [tmp_path / "a.nc"], {}, attrs
    ) as data:
        data.variables["fast_source_file_id"][:] = [0, 0]
        data.variables["fast_source_record"][:] = [4, 5]

    BUILDER.validate_saved_output(output, 2, 3, signatures)

    with Dataset(output) as data:
        assert data.product_type == "fast_image_conjunctions"
        assert len(data.dimensions["sample"]) == 2
        assert len(data.dimensions["energy_index"]) == 3


def test_duplicate_source_records_fail_saved_output_validation(tmp_path):
    output = tmp_path / "duplicates.nc"
    signatures = {"source_fast_collection_signature": "fast"}
    attrs = {
        "product_type": "fast_image_conjunctions",
        "schema_version": 1,
        **signatures,
    }
    with BUILDER.create_output(
        output, 2, 2, [tmp_path / "a.nc"], {}, attrs
    ) as data:
        data.variables["fast_source_file_id"][:] = [0, 0]
        data.variables["fast_source_record"][:] = [4, 4]

    with pytest.raises(ValueError, match="duplicate FAST records"):
        BUILDER.validate_saved_output(output, 2, 2, signatures)
