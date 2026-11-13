"""Focused tests for plotting the compact FAST--IMAGE dataset."""

import csv
import importlib.util
from pathlib import Path
from types import SimpleNamespace

from netCDF4 import Dataset
import numpy as np


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts/fast/plot_conjunction_analysis.py"
)
SPECIFICATION = importlib.util.spec_from_file_location(
    "plot_fast_conjunction_analysis", SCRIPT
)
PLOT = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(PLOT)


def write_variable(data, name, dimension, array):
    array = np.asarray(array)
    variable = data.createVariable(name, array.dtype, (dimension,))
    variable[:] = array


def make_compact_dataset(filename):
    sample_count = 240
    pixel_count = 120
    sample_year = np.repeat(np.arange(2000, 2004), sample_count // 4)
    pixel_year = np.repeat(np.arange(2000, 2004), pixel_count // 4)
    sample_axis = np.linspace(0, 1, sample_count)
    pixel_axis = np.linspace(0.02, 1, pixel_count)

    with Dataset(filename, "w") as data:
        data.createDimension("sample", sample_count)
        data.createDimension("pixel", pixel_count)
        data.processing_complete = "true"

        sample = {
            "sample_year": sample_year.astype(np.int16),
            "sample_fast_measurement_valid": np.ones(sample_count, dtype=np.int8),
            "sample_annotation_status": np.where(sample_axis > 0.2, 1, 0).astype(np.int8),
            "sample_dombeck_available": np.ones(sample_count, dtype=np.int8),
            "sample_dombeck_intense_precipitation": (sample_axis > 0.5).astype(np.int8),
            "sample_time_offset_seconds": -60 + 120 * sample_axis,
            "sample_detector_separation_deg": 1.2 * sample_axis,
            "sample_fast_altitude_km": 500 + 3000 * sample_axis,
            "sample_footprint_mlt": 24 * sample_axis,
            "sample_footprint_mlat": 60 + 25 * sample_axis,
            "sample_image_sza": 80 + 70 * sample_axis,
            "sample_image_dza": 75 - 60 * sample_axis,
            "sample_earthward_energy_flux_all_valid_mw_m2": 0.01 + 15 * sample_axis,
            "sample_earthward_energy_flux_centre_mw_m2": 0.01 + 12 * sample_axis,
            "sample_robinson_effective_energy_flux_mw_m2": 0.01 + 10 * sample_axis,
            "sample_earthward_mean_energy_all_valid_kev": 0.2 + 8 * sample_axis,
            "sample_earthward_mean_energy_centre_kev": 0.5 + 9 * sample_axis,
            "sample_robinson_effective_mean_energy_kev": 0.7 + 10 * sample_axis,
        }
        for name, array in sample.items():
            write_variable(data, name, "sample", array)

        wic = 50 + 1800 * pixel_axis
        energy = 0.3 + 9 * pixel_axis
        flux = 0.02 + 18 * pixel_axis
        pixel = {
            "pixel_year": pixel_year.astype(np.int16),
            "pixel_time_offset_seconds": -55 + 110 * pixel_axis,
            "pixel_detector_separation_deg": pixel_axis,
            "pixel_fast_altitude_km": 600 + 2800 * pixel_axis,
            "pixel_footprint_mlt": np.mod(np.arange(pixel_count) * 0.8, 24),
            "pixel_footprint_mlat": 62 + 22 * pixel_axis,
            "pixel_image_sza": 85 + 65 * pixel_axis,
            "pixel_image_dza": 75 - 60 * pixel_axis,
            "pixel_annotation_status": np.where(pixel_axis > 0.2, 1, 0).astype(np.int8),
            "pixel_fast_measurement_valid_exposure_fraction": np.ones(pixel_count),
            "pixel_maximum_detector_separation_deg": pixel_axis,
            "pixel_maximum_fast_altitude_km": 600 + 2800 * pixel_axis,
            "pixel_maximum_absolute_time_offset_seconds": np.abs(-55 + 110 * pixel_axis),
            "pixel_dombeck_available_exposure_fraction": np.ones(pixel_count),
            "pixel_earthward_energy_flux_full_mw_m2": flux,
            "pixel_earthward_energy_flux_500ev_mw_m2": 0.85 * flux,
            "pixel_robinson_effective_energy_flux_mw_m2": 0.7 * flux,
            "pixel_earthward_mean_energy_full_kev": energy,
            "pixel_earthward_mean_energy_500ev_kev": 1.1 * energy,
            "pixel_robinson_effective_mean_energy_kev": 1.2 * energy,
            "pixel_electron_hall_conductance_rees_s": 2 + 35 * pixel_axis,
            "pixel_electron_pedersen_conductance_rees_s": 1 + 18 * pixel_axis,
            "pixel_dombeck_intense_exposure_fraction": (pixel_axis > 0.5).astype(float),
            "pixel_image_wic_corrected": wic,
            "pixel_image_R": 5 + 120 * pixel_axis,
            "pixel_image_E0": 0.8 + 8 * pixel_axis,
            "pixel_image_Fe": 0.05 + 17 * pixel_axis,
            "pixel_image_H": 2.5 + 33 * pixel_axis,
            "pixel_image_P": 1.5 + 17 * pixel_axis,
            "pixel_image_method_valid": np.ones(pixel_count, dtype=np.int8),
            "pixel_image_conductance_valid": np.ones(pixel_count, dtype=np.int8),
        }
        for name, array in pixel.items():
            write_variable(data, name, "pixel", array)


def test_plotter_writes_primary_figures_and_correlation_table(tmp_path):
    source = tmp_path / "compact.nc"
    output = tmp_path / "figures"
    make_compact_dataset(source)

    PLOT.build_plots(SimpleNamespace(
        input=source, output=output,
        maximum_distribution_samples=1_000, overwrite=False,
        annotation_status="all", fast_population="all-valid",
        maximum_detector_separation_deg=None, maximum_altitude_km=None,
        minimum_sza_deg=None, maximum_dza_deg=None,
        maximum_absolute_time_offset_seconds=None,
    ))

    expected = (
        "selection_summary.csv", "collection_summary.csv",
        "selection_overview.png", "distribution_before.png",
        "distribution_after.png", "coverage_overview.png",
        "image_fast_comparison.png", "energy_comparison_by_mlt.png",
        "conductance_comparison_by_mlt.png", "ratio_energy_comparison.png",
        "ratio_energy_by_mlt.png", "wic_energy_comparison.png",
        "wic_energy_by_mlt.png", "wic_flux_comparison.png",
        "wic_flux_by_mlt.png", "wic_conductance_comparison.png",
        "wic_conductance_by_mlt.png", "dombeck_diagnostic.png",
        "correlation_summary.csv",
    )
    assert all((output / name).is_file() for name in expected)

    with (output / "correlation_summary.csv").open() as file:
        rows = list(csv.DictReader(file))
    assert {row["group"] for row in rows} == {
        "image_fast", "image_fast_energy_mlt",
        "image_fast_conductance_mlt", "ratio_energy", "wic_energy",
        "wic_flux", "wic_conductance",
    }
    ratio_sectors = {
        row["mlt_sector"] for row in rows
        if row["group"] == "ratio_energy"
    }
    assert ratio_sectors == {
        "all", "22-02", "02-06", "06-10", "10-14", "14-18", "18-22",
    }


def test_plot_filters_are_applied_from_stored_metadata(tmp_path):
    source = tmp_path / "compact.nc"
    make_compact_dataset(source)
    args = SimpleNamespace(
        annotation_status="accepted", fast_population="dombeck-intense",
        maximum_detector_separation_deg=0.9, maximum_altitude_km=3300,
        minimum_sza_deg=100, maximum_dza_deg=50,
        maximum_absolute_time_offset_seconds=45,
    )

    with Dataset(source) as data:
        sample_selected, pixel_selected = PLOT.selection_masks(data, args)
        assert 0 < np.count_nonzero(sample_selected) < sample_selected.size
        assert 0 < np.count_nonzero(pixel_selected) < pixel_selected.size
        assert np.all(
            data.variables["sample_annotation_status"][:][sample_selected] == 1
        )
        assert np.all(
            data.variables["pixel_dombeck_intense_exposure_fraction"][:][pixel_selected] > 0
        )
