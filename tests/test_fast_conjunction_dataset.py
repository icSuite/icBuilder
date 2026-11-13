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
import xarray as xr


try:
    import icreader  # noqa: F401
except ImportError:
    icreader = types.ModuleType("icreader")
    icreader.open_product = None
    sys.modules["icreader"] = icreader

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts/fast/generate_conjunction_dataset.py"
)
SPECIFICATION = importlib.util.spec_from_file_location(
    "fast_conjunction_dataset", SCRIPT
)
BUILDER = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(BUILDER)


def fast_products(paths):
    for path in paths:
        with Dataset(path, "w") as data:
            data.createDimension("time", 1)
            data.createVariable("fast_source_record", "i4", ("time",))[:] = [0]
    return [
        {"path": Path(path), "attrs": {},
         "energy_count": 1, "file_id": index}
        for index, path in enumerate(paths)
    ]


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


def test_cs_gate_requires_finite_wic_not_coverage():
    shape = (1, 1, 2)
    wic_coverage = np.ones(shape)
    wic_coverage[0, 0, 0] = 0
    wic = np.ones(shape)
    wic[0, 0, 1] = np.nan

    selected, values, counts = BUILDER.cs_support_gate(
        GateGrid(), 0,
        np.array([70.0, 71.0, 72.0, np.nan]),
        np.array([1.0, 8.0, 23.0, 1.0]),
        wic_coverage, wic,
    )

    np.testing.assert_array_equal(selected, [0])
    np.testing.assert_array_equal(values["row"], [0])
    np.testing.assert_array_equal(values["column"], [0])
    assert counts == {
        "nonfinite_footprint_count": 1,
        "outside_cs_grid_count": 1,
        "cs_wic_rejection_count": 1,
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
    def read(self, name, index=None):
        if name in ("mlat", "mlt"):
            values = np.zeros((1, 1, 2), dtype=float)
        elif name == "wic_corrected":
            values = np.array([[[10.0, 20.0]]])
        elif name == "si13_corrected":
            values = np.array([[[5.0, np.nan]]])
        else:
            raise KeyError(name)
        if index is None:
            return values
        return values[index]


def test_detector_pixel_with_wic_is_written_even_when_si13_is_missing(monkeypatch):
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
        "product_record": np.array([0, 1], dtype=np.int32),
        "source_record": np.array([4, 5], dtype=np.int32),
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
    assert candidates[0].fast_source_record.tolist() == [4, 5]
    assert counts["detector_wic_rejection_count"] == 0
    assert counts["successful_candidate_count"] == 2


def test_unique_match_uses_supported_nearest_frame_and_deterministic_ties():
    frame = pd.DataFrame({
        "fast_source_file_id": [0, 0, 0, 1],
        "fast_product_record": [7, 7, 7, 2],
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


def test_detector_gather_reads_only_matched_frames():
    class FrameProduct:
        def __init__(self):
            self.frame_runs_read = []

        def read(self, name, frame_slice):
            self.frame_runs_read.append((frame_slice.start, frame_slice.stop))
            frames = np.arange(frame_slice.start, frame_slice.stop)
            return np.stack([
                np.arange(6).reshape(2, 3) + frame * 10
                for frame in frames
            ])

    product = FrameProduct()
    values = BUILDER.gathered_field(
        product,
        "test",
        np.array([1, 1, 2, 3, 5]),
        np.array([0, 1, 0, 0, 1]),
        np.array([2, 1, 0, 0, 2]),
    )

    np.testing.assert_array_equal(values, [12, 14, 20, 30, 55])
    assert product.frame_runs_read == [(1, 4), (5, 6)]


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


def test_product_method_is_irrelevant_when_frames_align(tmp_path):
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

    p2.method = "zhang_paxton"
    BUILDER.validate_aligned_orbit(85, cs, p2, p3)


def test_image_fields_are_selected_only_when_present():
    precipitation = types.SimpleNamespace(
        fields={"wic_corrected", "si12", "E0", "Fe"}
    )
    conductance = types.SimpleNamespace(fields={"P", "H"})

    p2_fields, p3_fields = BUILDER.available_image_fields(
        precipitation, conductance
    )

    assert p2_fields == ("wic_corrected", "si12", "E0", "Fe")
    assert p3_fields == ("P", "H")


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
        output, 2, 3, np.array([70.0, 71.0]),
        fast_products([tmp_path / "a.nc"]),
        (), (), {}, attrs
    ) as data:
        data.variables["fast_source_file_id"][:] = [0, 0]
        data.variables["fast_product_record"][:] = [0, 1]
        data.variables["fast_source_record"][:] = [4, 5]

    BUILDER.validate_saved_output(output, 2, 3)

    with Dataset(output) as data:
        assert data.product_type == "fast_image_conjunctions"
        assert len(data.dimensions["sample"]) == 2
        assert len(data.dimensions["energy_index"]) == 3


def test_mixed_fast_energy_widths_are_padded_and_recorded(tmp_path):
    source_files = []
    products = []
    for file_id, energy_count in enumerate((2, 4)):
        filename = tmp_path / f"fast_{energy_count}.nc"
        source = xr.Dataset(coords={
            "time": np.array(["2000-01-01"], dtype="datetime64[ns]"),
            "energy_index": np.arange(energy_count),
            "msis_altitude_km": np.array([70.0, 71.0]),
        })
        source["fast_source_record"] = ("time", np.array([10 + file_id], dtype=np.int32))
        source["footprint_mlat"] = ("time", np.array([70.0 + file_id]))
        source["footprint_mlt"] = ("time", np.array([2.0 + file_id]))
        source["fast_energy_ev"] = (
            ("time", "energy_index"),
            np.full((1, energy_count), file_id + 1.0),
        )
        source["msis_oxygen_cm3"] = (
            ("time", "msis_altitude_km"),
            np.full((1, 2), 100.0 + file_id),
        )
        source.to_netcdf(filename)
        source_files.append(filename)
        products.append({
            "file_id": file_id,
            "path": filename,
            "energy_count": energy_count,
            "attrs": {},
        })

    matches = pd.DataFrame({
        "fast_source_file_id": [0, 1],
        "fast_product_record": [0, 0],
        "fast_source_record": [10, 11],
    })
    output = tmp_path / "mixed_energy.nc"
    with BUILDER.create_output(
        output,
        sample_count=2,
        energy_count=4,
        msis_altitude=np.array([70.0, 71.0]),
        fast_products=products,
        p2_fields=(),
        p3_fields=(),
        image_metadata={},
        root_attrs={"product_type": "fast_image_conjunctions", "schema_version": 1},
    ) as target:
        BUILDER.write_fast_fields(target, matches, products)

    with Dataset(output) as data:
        np.testing.assert_array_equal(
            data.variables["fast_energy_channel_count"][:], [2, 4]
        )
        spectrum = data.variables["fast_energy_ev"][:]
        np.testing.assert_array_equal(spectrum[0, :2], [1.0, 1.0])
        np.testing.assert_array_equal(np.ma.getmaskarray(spectrum[0]), [0, 0, 1, 1])
        np.testing.assert_array_equal(spectrum[1], [2.0, 2.0, 2.0, 2.0])
        np.testing.assert_array_equal(
            data.variables["msis_oxygen_cm3"][:],
            [[100.0, 100.0], [101.0, 101.0]],
        )


def test_duplicate_source_records_fail_saved_output_validation(tmp_path):
    output = tmp_path / "duplicates.nc"
    signatures = {"source_fast_collection_signature": "fast"}
    attrs = {
        "product_type": "fast_image_conjunctions",
        "schema_version": 1,
        **signatures,
    }
    with BUILDER.create_output(
        output, 2, 2, np.array([70.0, 71.0]),
        fast_products([tmp_path / "a.nc"]),
        (), (), {}, attrs
    ) as data:
        data.variables["fast_source_file_id"][:] = [0, 0]
        data.variables["fast_product_record"][:] = [4, 4]
        data.variables["fast_source_record"][:] = [4, 4]

    with pytest.raises(ValueError, match="duplicate FAST records"):
        BUILDER.validate_saved_output(output, 2, 2)
