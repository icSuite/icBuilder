"""Focused tests for FAST--IMAGE conjunction aggregation."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

from netCDF4 import Dataset
import numpy as np


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts/fast/aggregate_conjunctions.py"
)
SPECIFICATION = importlib.util.spec_from_file_location(
    "aggregate_fast_conjunctions", SCRIPT
)
AGGREGATE = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(AGGREGATE)


def test_conservative_regrid_preserves_integrated_number_flux():
    source_energy = np.array([1.0, 3.0])
    source_width = np.array([2.0, 2.0])
    source_spectrum = np.array([2.0, 4.0])
    target_energy = np.array([0.5, 1.5, 2.5, 3.5])
    target_width = np.ones(4)

    result = AGGREGATE.regrid_spectrum(
        source_energy, source_width, source_spectrum,
        target_energy, target_width,
    )

    np.testing.assert_allclose(result, [2.0, 2.0, 4.0, 4.0])
    assert np.sum(result * target_width) == np.sum(
        source_spectrum * source_width
    )


def test_mixed_energy_tables_are_duration_weighted_on_fine_grid():
    energy = np.array([
        [1.0, 3.0, np.nan, np.nan],
        [0.5, 1.5, 2.5, 3.5],
    ])
    width = np.array([
        [2.0, 2.0, np.nan, np.nan],
        [1.0, 1.0, 1.0, 1.0],
    ])
    earthward = np.array([
        [2.0, 4.0, np.nan, np.nan],
        [6.0, 6.0, 8.0, 8.0],
    ])

    result = AGGREGATE.aggregate_spectra(
        energy, width, earthward, earthward,
        np.array([1.0, 3.0]), 4,
    )

    np.testing.assert_allclose(result[0], energy[1])
    np.testing.assert_allclose(result[1], width[1])
    np.testing.assert_allclose(result[2], [5.0, 5.0, 7.0, 7.0])


def make_conjunction(filename):
    sample_count = 4
    energy_count = 4
    altitude = np.linspace(70.0, 300.0, 226)
    with Dataset(filename, "w") as data:
        data.createDimension("sample", sample_count)
        data.createDimension("energy_index", energy_count)
        data.createDimension("msis_altitude_km", altitude.size)
        data.createDimension("fast_source_file", 1)
        data.processing_complete = "true"
        data.selected_year = 2000
        data.createVariable(
            "fast_source_file", str, ("fast_source_file",)
        )[0] = "fast_source.nc"
        data.createVariable(
            "msis_altitude_km", "f8", ("msis_altitude_km",)
        )[:] = altitude

        scalar = {
            "fast_source_file_id": np.zeros(4, dtype=np.int32),
            "fast_time": np.array([100.0, 101.0, 200.0, 201.0]),
            "image_time": np.array([100.0, 100.0, 200.0, 200.0]),
            "image_orbit": np.array([85, 85, 85, 85], dtype=np.int32),
            "image_frame": np.array([0, 0, 1, 1], dtype=np.int32),
            "detector_row": np.array([2, 2, 3, 3], dtype=np.int16),
            "detector_column": np.array([4, 4, 5, 5], dtype=np.int16),
            "fast_altitude_km": np.array([1000.0, 1100.0, 1200.0, 1300.0]),
            "fast_measurement_valid": np.ones(4, dtype=np.int8),
            "detector_separation_deg": np.array([0.1, 0.2, 0.3, 0.4]),
            "time_offset_seconds": np.array([0.0, 1.0, 0.0, -1.0]),
            "image_sza": np.array([110.0, 110.0, 115.0, 115.0]),
            "cs_wic_coverage": np.ones(4),
            "dombeck_window_record_count": np.array([2, 2, 2, 2], dtype=np.int16),
            "dombeck_window_span_seconds": np.array([2.0, 6.0, 2.0, 2.0]),
            "dombeck_available": np.ones(4, dtype=np.int8),
            "dombeck_intense_precipitation": np.array([0, 1, 1, 1], dtype=np.int8),
            "footprint_geodetic_latitude": np.array([70.0, 72.0, 74.0, 76.0]),
            "footprint_geodetic_longitude": np.array([10.0, 12.0, 14.0, 16.0]),
            "footprint_mlat": np.array([71.0, 73.0, 75.0, 77.0]),
            "footprint_mlt": np.array([23.9, 0.1, 1.0, 1.2]),
            "footprint_magnetic_field_nt": np.full(4, 50_000.0),
        }
        for name, array in scalar.items():
            variable = data.createVariable(name, array.dtype, ("sample",))
            variable[:] = array

        energy = np.array([
            [1000.0, 3000.0, np.nan, np.nan],
            [500.0, 1500.0, 2500.0, 3500.0],
            [500.0, 1500.0, 2500.0, 3500.0],
            [500.0, 1500.0, 2500.0, 3500.0],
        ])
        width = np.array([
            [2000.0, 2000.0, np.nan, np.nan],
            [1000.0, 1000.0, 1000.0, 1000.0],
            [1000.0, 1000.0, 1000.0, 1000.0],
            [1000.0, 1000.0, 1000.0, 1000.0],
        ])
        spectrum = np.array([
            [2e6, 4e6, np.nan, np.nan],
            [6e6, 6e6, 8e6, 8e6],
            [3e6, 3e6, 3e6, 3e6],
            [5e6, 5e6, 5e6, 5e6],
        ])
        for name, array in {
            "fast_energy_ev": energy,
            "fast_energy_width_ev": width,
            "earthward_differential_number_flux_cm2_s_ev": spectrum,
            "anti_earthward_differential_number_flux_cm2_s_ev": spectrum * 0.1,
        }.items():
            data.createVariable(
                name, "f8", ("sample", "energy_index"), fill_value=np.nan
            )[:] = array

        profile = {
            "msis_oxygen_cm3": 1e12 * np.exp(-(altitude - 70) / 30),
            "msis_nitrogen_cm3": 5e11 * np.exp(-(altitude - 70) / 25),
            "msis_molecular_oxygen_cm3": 8e11 * np.exp(-(altitude - 70) / 20),
            "msis_mass_density_g_cm3": 1e-9 * np.exp(-(altitude - 70) / 20),
            "msis_temperature_k": np.full(altitude.size, 700.0),
        }
        for name, array in profile.items():
            data.createVariable(
                name, "f8", ("sample", "msis_altitude_km")
            )[:] = np.broadcast_to(array, (sample_count, altitude.size))


def test_compiler_writes_sample_flags_and_pixel_aggregates(tmp_path):
    source = tmp_path / "conjunction.nc"
    output = tmp_path / "analysis.nc"
    make_conjunction(source)
    args = SimpleNamespace(
        input=[source], annotations=None, output=output,
        maximum_detector_separation_deg=None,
        maximum_altitude_km=None, minimum_sza_deg=None,
        maximum_absolute_time_offset_seconds=None,
        record_batch_size=10, conductance_batch_size=10,
        overwrite=False,
    )

    AGGREGATE.build_dataset(args)

    with Dataset(output) as data:
        assert len(data.dimensions["sample"]) == 4
        assert len(data.dimensions["pixel"]) == 2
        np.testing.assert_array_equal(data.variables["sample_selected"][:], 1)
        np.testing.assert_array_equal(
            data.variables["pixel_fast_sample_count"][:], [2, 2]
        )
        assert np.isclose(data.variables["pixel_footprint_mlt"][0], 0.05, atol=1e-3)
        assert np.isfinite(
            data.variables["pixel_electron_hall_conductance_rees_s"][:]
        ).all()


def test_selection_filters_samples_but_not_dombeck_weak(tmp_path):
    source = tmp_path / "conjunction.nc"
    output = tmp_path / "analysis.nc"
    make_conjunction(source)
    args = SimpleNamespace(
        input=[source], annotations=None, output=output,
        maximum_detector_separation_deg=None,
        maximum_altitude_km=1150.0, minimum_sza_deg=None,
        maximum_absolute_time_offset_seconds=None,
        record_batch_size=10, conductance_batch_size=10,
        overwrite=False,
    )

    AGGREGATE.build_dataset(args)

    with Dataset(output) as data:
        np.testing.assert_array_equal(
            data.variables["sample_selected"][:], [1, 1, 0, 0]
        )
        assert data.variables["sample_dombeck_intense_precipitation"][0] == 0
        assert data.variables["sample_selected"][0] == 1
        assert len(data.dimensions["pixel"]) == 1


def test_annotations_select_only_accepted_tracks(tmp_path):
    source = tmp_path / "conjunction.nc"
    output = tmp_path / "analysis.nc"
    annotations = tmp_path / "annotations.csv"
    make_conjunction(source)
    annotations.write_text(
        "year,image_orbit,image_frame,image_time,fast_source_file,accepted,n_fast_samples\n"
        "2000,85,0,1970-01-01T00:01:40,fast_source.nc,1,2\n"
        "2000,85,1,1970-01-01T00:03:20,fast_source.nc,0,2\n"
    )
    args = SimpleNamespace(
        input=[source], annotations=annotations, output=output,
        maximum_detector_separation_deg=None,
        maximum_altitude_km=None, minimum_sza_deg=None,
        maximum_absolute_time_offset_seconds=None,
        record_batch_size=10, conductance_batch_size=10,
        overwrite=False,
    )

    AGGREGATE.build_dataset(args)

    with Dataset(output) as data:
        np.testing.assert_array_equal(
            data.variables["sample_annotation_status"][:], [1, 1, 0, 0]
        )
        np.testing.assert_array_equal(
            data.variables["sample_selected"][:], [1, 1, 0, 0]
        )
        assert len(data.dimensions["pixel"]) == 1
