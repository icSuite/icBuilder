"""Compile FAST--IMAGE conjunctions into a compact analysis dataset.

The output retains sample-level scalar fields and one physically aggregated
FAST record per IMAGE detector pixel.  It never reopens raw FAST EES
files and never reruns MSIS.
"""

#%% Imports and constants

import argparse
import csv
import json
from pathlib import Path

from icphysics import conductance_from_differential_number_flux
from icphysics import robinson_effective_precipitation
from netCDF4 import Dataset
import numpy as np
from tqdm import tqdm


EV_CM2_S_TO_MW_M2 = 1.602176634e-12
DEFAULT_RECORD_BATCH_SIZE = 50_000
DEFAULT_CONDUCTANCE_BATCH_SIZE = 1_000

PIXEL_KEY_FIELDS = (
    "image_orbit", "image_frame", "detector_row", "detector_column"
)

REQUIRED_SCALARS = (
    "fast_source_file_id", "fast_time", "image_time",
    "image_orbit", "image_frame", "detector_row", "detector_column",
    "fast_altitude_km", "fast_measurement_valid",
    "detector_separation_deg", "time_offset_seconds", "image_sza",
    "cs_wic_coverage",
    "dombeck_window_record_count", "dombeck_window_span_seconds",
    "dombeck_available", "dombeck_intense_precipitation",
    "footprint_geodetic_latitude", "footprint_geodetic_longitude",
    "footprint_mlat", "footprint_mlt", "footprint_magnetic_field_nt",
)

SPECTRUM_FIELDS = (
    "fast_energy_ev", "fast_energy_width_ev",
    "earthward_differential_number_flux_cm2_s_ev",
    "anti_earthward_differential_number_flux_cm2_s_ev",
)

ATMOSPHERE_FIELDS = {
    "altitude_km": "msis_altitude_km",
    "oxygen_cm3": "msis_oxygen_cm3",
    "nitrogen_cm3": "msis_nitrogen_cm3",
    "molecular_oxygen_cm3": "msis_molecular_oxygen_cm3",
    "mass_density_g_cm3": "msis_mass_density_g_cm3",
    "temperature_k": "msis_temperature_k",
}

WEIGHTED_CONTEXT_FIELDS = (
    "fast_time", "fast_altitude_km",
    "footprint_geodetic_latitude", "footprint_mlat",
    "footprint_magnetic_field_nt", "time_offset_seconds",
    "detector_separation_deg", "cs_wic_coverage",
)

CIRCULAR_CONTEXT_FIELDS = (
    ("footprint_geodetic_longitude", 360.0),
    ("footprint_mlt", 24.0),
)


#%% NetCDF and selection utilities

def values(variable):
    """Return a NetCDF variable as an ndarray with masked values made NaN."""

    array = variable[:]
    if np.ma.isMaskedArray(array):
        if np.issubdtype(array.dtype, np.floating):
            return np.asarray(array.filled(np.nan))
        return np.asarray(array.filled(-1))
    return np.asarray(array)


def indexed_values(variable, records):
    """Read records through contiguous runs and restore the requested order."""

    if not variable.dimensions or variable.dimensions[0] != "sample":
        raise ValueError(
            f"{variable.name} is not indexed by the sample dimension"
        )
    records = np.asarray(records, dtype=np.int64)
    if records.size == 0:
        return np.empty((0, *variable.shape[1:]), dtype=variable.dtype)
    order = np.argsort(records, kind="stable")
    sorted_records = records[order]
    breaks = np.flatnonzero(np.diff(sorted_records) > 1) + 1
    starts = np.r_[0, breaks]
    stops = np.r_[breaks, records.size]
    parts = []
    for first, stop in zip(starts, stops):
        source_start = int(sorted_records[first])
        source_stop = int(sorted_records[stop - 1]) + 1
        block = variable[source_start:source_stop]
        local = sorted_records[first:stop] - source_start
        parts.append(block[local])
    array = np.concatenate(parts, axis=0)
    inverse = np.empty_like(order)
    inverse[order] = np.arange(order.size)
    array = array[inverse]
    if np.ma.isMaskedArray(array):
        if np.issubdtype(array.dtype, np.floating):
            array = array.filled(np.nan)
        else:
            array = array.filled(-1)
    return np.asarray(array)


def scalar_fields(data):
    return tuple(
        name for name, variable in data.variables.items()
        if variable.dimensions == ("sample",)
    )


def read_annotations(filename):
    if filename is None:
        return {}
    required = {
        "year", "image_orbit", "image_frame", "image_time",
        "fast_source_file", "accepted",
    }
    with Path(filename).open(newline="") as file:
        reader = csv.DictReader(file)
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{filename} is missing {sorted(missing)}")
        rows = list(reader)

    annotations = {}
    for row in rows:
        accepted = int(row["accepted"])
        if accepted not in (0, 1):
            raise ValueError("FAST annotations must be binary")
        key = (
            int(row["year"]), int(row["image_orbit"]),
            int(np.datetime64(row["image_time"], "s").astype(np.int64)),
            Path(row["fast_source_file"]).name,
        )
        if key in annotations:
            raise ValueError(f"duplicate FAST annotation key {key}")
        annotations[key] = accepted
    return annotations


def annotation_status(data, year, annotations):
    count = len(data.dimensions["sample"])
    status = np.full(count, -1, dtype=np.int8)
    rows = [
        (key, row) for key, row in annotations.items() if key[0] == year
    ]
    if not rows:
        return status

    source_files = [
        Path(name).name
        for name in np.asarray(data.variables["fast_source_file"][:], str)
    ]
    if len(source_files) != len(set(source_files)):
        raise ValueError(f"FAST source basenames are not unique in {year}")
    source_lookup = {name: index for index, name in enumerate(source_files)}

    key_dtype = np.dtype([
        ("orbit", "i4"), ("time", "i8"), ("source", "i4"),
    ])
    annotation_keys = np.empty(len(rows), dtype=key_dtype)
    annotation_value = np.empty(len(rows), dtype=np.int8)
    for index, (key, accepted) in enumerate(rows):
        source_name = key[3]
        if source_name not in source_lookup:
            raise ValueError(
                f"FAST annotation refers to absent source {source_name}"
            )
        annotation_keys[index] = (
            key[1], key[2], source_lookup[source_name],
        )
        annotation_value[index] = accepted

    order = np.argsort(annotation_keys, kind="stable")
    annotation_keys = annotation_keys[order]
    annotation_value = annotation_value[order]

    sample_keys = np.empty(count, dtype=key_dtype)
    sample_keys["orbit"] = values(data.variables["image_orbit"])
    sample_keys["time"] = np.rint(values(data.variables["image_time"])).astype(np.int64)
    sample_keys["source"] = values(data.variables["fast_source_file_id"])

    position = np.searchsorted(annotation_keys, sample_keys)
    within = position < annotation_keys.size
    matched = np.zeros(count, dtype=bool)
    matched[within] = annotation_keys[position[within]] == sample_keys[within]
    status[matched] = annotation_value[position[matched]]
    return status


def validate_input_fields(data):
    required = set(REQUIRED_SCALARS) | set(SPECTRUM_FIELDS)
    required |= set(ATMOSPHERE_FIELDS.values())
    missing = required.difference(data.variables)
    if missing:
        raise ValueError(f"conjunction file is missing {sorted(missing)}")

def pixel_runs(data):
    records = np.arange(len(data.dimensions["sample"]), dtype=np.int64)
    if records.size == 0:
        return records, np.empty(0, int), np.empty(0, int)

    arrays = [
        values(data.variables[name])
        for name in PIXEL_KEY_FIELDS
    ]
    keys = np.rec.fromarrays(arrays, names="orbit,frame,row,column")
    order = np.argsort(keys, order=keys.dtype.names, kind="stable")
    records = records[order]
    keys = keys[order]
    starts = np.r_[0, np.flatnonzero(keys[1:] != keys[:-1]) + 1]
    stops = np.r_[starts[1:], records.size]
    return records, starts, stops


#%% Spectrum aggregation and conductance

def weighted_mean(array, weight):
    array = np.asarray(array, dtype=float)
    weight = np.asarray(weight, dtype=float)
    if array.ndim == 1:
        valid = np.isfinite(array) & np.isfinite(weight) & (weight > 0)
        if not np.any(valid):
            return np.nan
        return np.sum(array[valid] * weight[valid]) / np.sum(weight[valid])
    valid = np.isfinite(array) & np.isfinite(weight[:, None]) & (weight[:, None] > 0)
    numerator = np.sum(np.where(valid, array * weight[:, None], 0), axis=0)
    denominator = np.sum(np.where(valid, weight[:, None], 0), axis=0)
    return np.divide(
        numerator, denominator, out=np.full(array.shape[1], np.nan),
        where=denominator > 0,
    )


def circular_weighted_mean(array, weight, period):
    angle = np.asarray(array, dtype=float) * 2 * np.pi / period
    sine = weighted_mean(np.sin(angle), weight)
    cosine = weighted_mean(np.cos(angle), weight)
    if not np.isfinite(sine) or not np.isfinite(cosine):
        return np.nan
    return np.mod(np.arctan2(sine, cosine) * period / (2 * np.pi), period)


def spectrum_table(energy, width):
    valid = (
        np.isfinite(energy) & (energy > 0)
        & np.isfinite(width) & (width > 0)
    )
    index = np.flatnonzero(valid)
    order = np.argsort(energy[index])
    index = index[order]
    return index, energy[index], width[index]


def regrid_spectrum(source_energy, source_width, source_spectrum,
                    target_energy, target_width):
    """Conservatively map piecewise-constant differential number flux."""

    valid = (
        np.isfinite(source_energy) & (source_energy > 0)
        & np.isfinite(source_width) & (source_width > 0)
        & np.isfinite(source_spectrum) & (source_spectrum >= 0)
    )
    source_energy = source_energy[valid]
    source_width = source_width[valid]
    source_spectrum = source_spectrum[valid]
    if source_energy.size == 0:
        return np.full(target_energy.size, np.nan)
    order = np.argsort(source_energy)
    source_energy = source_energy[order]
    source_width = source_width[order]
    source_spectrum = source_spectrum[order]

    source_lower = source_energy - source_width / 2
    source_upper = source_energy + source_width / 2
    target_lower = target_energy - target_width / 2
    target_upper = target_energy + target_width / 2
    overlap = np.maximum(
        0,
        np.minimum(target_upper[:, None], source_upper[None, :])
        - np.maximum(target_lower[:, None], source_lower[None, :]),
    )
    coverage = np.sum(overlap, axis=1)
    complete = np.isclose(coverage, target_width, rtol=1e-6, atol=1e-8)
    result = np.full(target_energy.size, np.nan)
    result[complete] = (
        overlap[complete] @ source_spectrum / target_width[complete]
    )
    return result


def aggregate_spectra(energy, width, earthward, anti_earthward, weight,
                      output_energy_count):
    tables = [spectrum_table(energy[index], width[index]) for index in range(len(weight))]
    populated = [table for table in tables if table[1].size]
    if not populated:
        return tuple(
            np.full(output_energy_count, np.nan) for _ in range(4)
        )
    lower = min(table[1][0] - table[2][0] / 2 for table in populated)
    upper = max(table[1][-1] + table[2][-1] / 2 for table in populated)
    candidates = []
    for index, table in enumerate(tables):
        if not table[1].size:
            continue
        table_lower = table[1][0] - table[2][0] / 2
        table_upper = table[1][-1] + table[2][-1] / 2
        if table_lower <= lower + 1e-6 and table_upper >= upper - 1e-6:
            candidates.append((table[1].size, index))
    if not candidates:
        raise ValueError("no FAST energy table covers all tables in one IMAGE pixel")
    target_index = max(candidates)[1]
    _, target_energy, target_width = tables[target_index]
    if target_energy.size > output_energy_count:
        raise ValueError("aggregated energy table exceeds output energy dimension")

    earthward_sum = np.zeros(target_energy.size)
    earthward_weight = np.zeros(target_energy.size)
    anti_sum = np.zeros(target_energy.size)
    anti_weight = np.zeros(target_energy.size)
    for index, table in enumerate(tables):
        channel, record_energy, record_width = table
        if not record_energy.size:
            continue
        record_earthward = earthward[index, channel]
        record_anti = anti_earthward[index, channel]
        if (
            np.array_equal(record_energy, target_energy)
            and np.array_equal(record_width, target_width)
        ):
            mapped_earthward = record_earthward
            mapped_anti = record_anti
        else:
            mapped_earthward = regrid_spectrum(
                record_energy, record_width, record_earthward,
                target_energy, target_width,
            )
            mapped_anti = regrid_spectrum(
                record_energy, record_width, record_anti,
                target_energy, target_width,
            )
        valid = np.isfinite(mapped_earthward) & (mapped_earthward >= 0)
        earthward_sum[valid] += mapped_earthward[valid] * weight[index]
        earthward_weight[valid] += weight[index]
        valid = np.isfinite(mapped_anti) & (mapped_anti >= 0)
        anti_sum[valid] += mapped_anti[valid] * weight[index]
        anti_weight[valid] += weight[index]

    aggregate_earthward = np.divide(
        earthward_sum, earthward_weight,
        out=np.full(target_energy.size, np.nan), where=earthward_weight > 0,
    )
    aggregate_anti = np.divide(
        anti_sum, anti_weight,
        out=np.full(target_energy.size, np.nan), where=anti_weight > 0,
    )
    result = [np.full(output_energy_count, np.nan) for _ in range(4)]
    result[0][:target_energy.size] = target_energy
    result[1][:target_energy.size] = target_width
    result[2][:target_energy.size] = aggregate_earthward
    result[3][:target_energy.size] = aggregate_anti
    return tuple(result)


def moments(energy, width, spectrum, minimum_energy_ev=0):
    valid = (
        np.isfinite(energy) & (energy >= minimum_energy_ev)
        & np.isfinite(width) & (width > 0)
        & np.isfinite(spectrum) & (spectrum >= 0)
    )
    number_flux = np.sum(np.where(valid, spectrum * width, 0), axis=1)
    energy_flux = np.sum(
        np.where(valid, spectrum * width * energy, 0), axis=1
    ) * EV_CM2_S_TO_MW_M2
    mean_energy = np.divide(
        energy_flux / EV_CM2_S_TO_MW_M2, number_flux,
        out=np.full(number_flux.shape, np.nan), where=number_flux > 0,
    ) / 1000
    return number_flux, energy_flux, mean_energy


def calculate_conductance(energy, width, spectrum, magnetic_field_nt,
                          atmosphere, maximum_batch_size):
    count = energy.shape[0]
    pedersen = np.full(count, np.nan)
    hall = np.full(count, np.nan)
    channel_count = np.zeros(count, dtype=np.int16)
    groups = {}
    for index in range(count):
        valid = (
            np.isfinite(energy[index]) & (energy[index] > 0)
            & np.isfinite(width[index]) & (width[index] > 0)
            & np.isfinite(spectrum[index]) & (spectrum[index] >= 0)
        )
        if np.count_nonzero(valid) < 2:
            continue
        key = (
            valid.tobytes(), energy[index, valid].tobytes(),
            width[index, valid].tobytes(),
        )
        groups.setdefault(key, []).append(index)

    for indices in groups.values():
        indices = np.asarray(indices, dtype=int)
        channels = (
            np.isfinite(energy[indices[0]]) & (energy[indices[0]] > 0)
            & np.isfinite(width[indices[0]]) & (width[indices[0]] > 0)
            & np.isfinite(spectrum[indices[0]]) & (spectrum[indices[0]] >= 0)
        )
        for start in range(0, indices.size, maximum_batch_size):
            batch = indices[start:start + maximum_batch_size]
            profiles = {
                name: (array if name == "altitude_km" else array[batch])
                for name, array in atmosphere.items()
            }
            result = conductance_from_differential_number_flux(
                energy[batch[0], channels] / 1000,
                spectrum[batch][:, channels], width[batch[0], channels],
                magnetic_field_nt[batch], profiles,
            )
            pedersen[batch] = result["P"]
            hall[batch] = result["H"]
            channel_count[batch] = np.count_nonzero(channels)

    effective = robinson_effective_precipitation(hall, pedersen)
    return {
        "electron_pedersen_conductance_rees_s": pedersen,
        "electron_hall_conductance_rees_s": hall,
        "conductance_spectrum_channel_count": channel_count,
        "robinson_effective_mean_energy_kev": effective["effective_energy"],
        "robinson_effective_energy_flux_mw_m2": effective["effective_energy_flux"],
        "robinson_effective_valid": effective["invertible"].astype(np.int8),
    }


#%% Output schema and pixel batches

def create_variable(data, name, dtype, dimensions, fill=None, attrs=None):
    shape = [len(data.dimensions[dimension]) for dimension in dimensions]
    trailing = int(np.prod(shape[1:])) if len(shape) > 1 else 1
    itemsize = np.dtype(dtype).itemsize
    rows = max(1, min(shape[0], 1_048_576 // max(1, trailing * itemsize)))
    options = {
        "zlib": True, "complevel": 1, "shuffle": True,
        "chunksizes": (rows, *shape[1:]),
    }
    if fill is not None:
        options["fill_value"] = fill
    variable = data.createVariable(name, dtype, dimensions, **options)
    for attribute, value in (attrs or {}).items():
        if attribute != "_FillValue":
            variable.setncattr(attribute, value)
    return variable


def create_output(filename, sample_count, pixel_count, energy_count,
                  msis_altitude, source_files, sample_schema, image_schema,
                  root_attrs):
    output = Dataset(filename, "w", format="NETCDF4")
    output.createDimension("sample", sample_count)
    output.createDimension("pixel", pixel_count)
    output.createDimension("energy_index", energy_count)
    output.createDimension("fast_source_file", len(source_files))
    output.createDimension("msis_altitude_km", len(msis_altitude))

    source = output.createVariable(
        "fast_source_file", str, ("fast_source_file",)
    )
    source[:] = np.asarray(source_files, dtype=object)
    altitude = output.createVariable(
        "msis_altitude_km", "f8", ("msis_altitude_km",)
    )
    altitude[:] = msis_altitude
    altitude.units = "km"

    for name, schema in sample_schema.items():
        create_variable(
            output, f"sample_{name}", schema["dtype"], ("sample",),
            schema["fill"], schema["attrs"],
        )
    create_variable(output, "sample_year", "i2", ("sample",), -1)
    create_variable(output, "sample_annotation_status", "i1", ("sample",), -1)

    integer_pixel_fields = {
        "year": "i2", "image_orbit": "i4", "image_frame": "i4",
        "detector_row": "i2", "detector_column": "i2",
        "annotation_status": "i1",
        "fast_sample_count": "i4", "conductance_spectrum_channel_count": "i2",
        "robinson_effective_valid": "i1",
    }
    float_pixel_fields = (
        "exposure_seconds", "image_time", *WEIGHTED_CONTEXT_FIELDS,
        "footprint_geodetic_longitude", "footprint_mlt",
        "maximum_fast_altitude_km", "maximum_detector_separation_deg",
        "maximum_absolute_time_offset_seconds",
        "fast_measurement_valid_exposure_fraction",
        "dombeck_available_exposure_fraction",
        "dombeck_intense_exposure_fraction",
        "earthward_number_flux_full_cm2_s",
        "earthward_energy_flux_full_mw_m2",
        "earthward_mean_energy_full_kev",
        "earthward_number_flux_500ev_cm2_s",
        "earthward_energy_flux_500ev_mw_m2",
        "earthward_mean_energy_500ev_kev",
        "anti_earthward_number_flux_full_cm2_s",
        "anti_earthward_energy_flux_full_mw_m2",
        "net_earthward_number_flux_full_cm2_s",
        "net_earthward_energy_flux_full_mw_m2",
        "net_earthward_mean_energy_full_kev",
        "electron_pedersen_conductance_rees_s",
        "electron_hall_conductance_rees_s",
        "robinson_effective_mean_energy_kev",
        "robinson_effective_energy_flux_mw_m2",
    )
    for name, dtype in integer_pixel_fields.items():
        create_variable(output, f"pixel_{name}", dtype, ("pixel",), -1)
    for name in dict.fromkeys(float_pixel_fields):
        attrs = sample_schema[name]["attrs"] if name in sample_schema else {}
        create_variable(
            output, f"pixel_{name}", "f8", ("pixel",), np.nan, attrs
        )

    spectrum_units = {
        "fast_energy_ev": "eV",
        "fast_energy_width_ev": "eV",
        "earthward_differential_number_flux_cm2_s_ev": "cm-2 s-1 eV-1",
        "anti_earthward_differential_number_flux_cm2_s_ev": "cm-2 s-1 eV-1",
    }
    for name, unit in spectrum_units.items():
        create_variable(output, f"pixel_{name}", "f8",
                        ("pixel", "energy_index"), np.nan,
                        {"units": unit})

    for name, schema in image_schema.items():
        create_variable(
            output, f"pixel_{name}", schema["dtype"], ("pixel",),
            schema["fill"], schema["attrs"],
        )
    derived_units = {
        "pixel_exposure_seconds": "s",
        "pixel_maximum_fast_altitude_km": "km",
        "pixel_maximum_detector_separation_deg": "degrees",
        "pixel_maximum_absolute_time_offset_seconds": "s",
        "pixel_fast_measurement_valid_exposure_fraction": "1",
        "pixel_dombeck_available_exposure_fraction": "1",
        "pixel_dombeck_intense_exposure_fraction": "1",
        "pixel_earthward_number_flux_full_cm2_s": "cm-2 s-1",
        "pixel_earthward_energy_flux_full_mw_m2": "mW m-2",
        "pixel_earthward_mean_energy_full_kev": "keV",
        "pixel_earthward_number_flux_500ev_cm2_s": "cm-2 s-1",
        "pixel_earthward_energy_flux_500ev_mw_m2": "mW m-2",
        "pixel_earthward_mean_energy_500ev_kev": "keV",
        "pixel_anti_earthward_number_flux_full_cm2_s": "cm-2 s-1",
        "pixel_anti_earthward_energy_flux_full_mw_m2": "mW m-2",
        "pixel_net_earthward_number_flux_full_cm2_s": "cm-2 s-1",
        "pixel_net_earthward_energy_flux_full_mw_m2": "mW m-2",
        "pixel_net_earthward_mean_energy_full_kev": "keV",
        "pixel_electron_pedersen_conductance_rees_s": "S",
        "pixel_electron_hall_conductance_rees_s": "S",
        "pixel_robinson_effective_mean_energy_kev": "keV",
        "pixel_robinson_effective_energy_flux_mw_m2": "mW m-2",
    }
    for name, unit in derived_units.items():
        output.variables[name].units = unit
    for name, value in root_attrs.items():
        output.setncattr(name, value)
    return output


def schema_for_variable(variable):
    fill = getattr(variable, "_FillValue", None)
    if fill is None:
        fill = np.nan if np.issubdtype(variable.dtype, np.floating) else -1
    attrs = {
        name: variable.getncattr(name)
        for name in variable.ncattrs() if name != "_FillValue"
    }
    return {"dtype": variable.dtype, "fill": fill, "attrs": attrs}


def exact_group_value(array, name):
    first = array[0]
    equal = (array == first) | (np.isnan(array) & np.isnan(first)) \
        if np.issubdtype(array.dtype, np.floating) else array == first
    if not np.all(equal):
        raise ValueError(f"IMAGE field {name} differs within one detector pixel")
    return first


def aggregate_pixel_batch(data, records, starts, stops, year, energy_count,
                          msis_altitude, image_fields, conductance_batch_size,
                          annotation_status):
    read_fields = set(SPECTRUM_FIELDS) | {
        source_name for name, source_name in ATMOSPHERE_FIELDS.items()
        if name != "altitude_km"
    }
    read_fields |= set(WEIGHTED_CONTEXT_FIELDS) | {
        name for name, _period in CIRCULAR_CONTEXT_FIELDS
    }
    read_fields |= {
        "image_time", "image_orbit", "image_frame",
        "detector_row", "detector_column",
        "fast_measurement_valid",
        "dombeck_window_record_count", "dombeck_window_span_seconds",
        "dombeck_available", "dombeck_intense_precipitation",
    }
    read_fields |= set(image_fields)
    arrays = {
        name: indexed_values(data.variables[name], records)
        for name in read_fields
    }

    pixel_count = len(starts)
    result = {
        "year": np.full(pixel_count, year, dtype=np.int16),
        "fast_sample_count": (stops - starts).astype(np.int32),
        "annotation_status": np.full(pixel_count, -1, dtype=np.int8),
        "exposure_seconds": np.full(pixel_count, np.nan),
    }
    for name in (*PIXEL_KEY_FIELDS, *WEIGHTED_CONTEXT_FIELDS,
                 *(name for name, _period in CIRCULAR_CONTEXT_FIELDS)):
        result[name] = np.full(pixel_count, np.nan)
    result["image_time"] = np.full(pixel_count, np.nan)
    for name in image_fields:
        result[name] = np.empty(pixel_count, dtype=data.variables[name].dtype)
    for name in (
        "maximum_fast_altitude_km", "maximum_detector_separation_deg",
        "maximum_absolute_time_offset_seconds",
        "fast_measurement_valid_exposure_fraction",
        "dombeck_available_exposure_fraction",
        "dombeck_intense_exposure_fraction",
    ):
        result[name] = np.full(pixel_count, np.nan)

    energy = np.full((pixel_count, energy_count), np.nan)
    width = np.full_like(energy, np.nan)
    earthward = np.full_like(energy, np.nan)
    anti_earthward = np.full_like(energy, np.nan)
    atmosphere = {
        "altitude_km": msis_altitude,
        **{
            name: np.full((pixel_count, len(msis_altitude)), np.nan)
            for name in ATMOSPHERE_FIELDS if name != "altitude_km"
        },
    }

    for output_index, (start, stop) in enumerate(zip(starts, stops)):
        local = slice(start, stop)
        count = arrays["dombeck_window_record_count"][local]
        span = arrays["dombeck_window_span_seconds"][local]
        weight = np.divide(
            span, count, out=np.full(stop - start, np.nan), where=count > 0
        )
        valid_weight = np.isfinite(weight) & (weight > 0)
        if np.any(valid_weight):
            result["exposure_seconds"][output_index] = np.sum(
                weight[valid_weight]
            )

        annotated = annotation_status[local]
        annotated = annotated[annotated >= 0]
        if annotated.size:
            unique = np.unique(annotated)
            if unique.size != 1:
                raise ValueError("annotation flags differ within one IMAGE pixel")
            result["annotation_status"][output_index] = unique[0]

        for name in WEIGHTED_CONTEXT_FIELDS:
            result[name][output_index] = weighted_mean(arrays[name][local], weight)
        for name, period in CIRCULAR_CONTEXT_FIELDS:
            result[name][output_index] = circular_weighted_mean(
                arrays[name][local], weight, period
            )
        for name in PIXEL_KEY_FIELDS:
            result[name][output_index] = exact_group_value(
                arrays[name][local], name
            )
        result["image_time"][output_index] = exact_group_value(
            arrays["image_time"][local], "image_time"
        )
        for name in image_fields:
            result[name][output_index] = exact_group_value(
                arrays[name][local], name
            )

        result["maximum_fast_altitude_km"][output_index] = np.nanmax(
            arrays["fast_altitude_km"][local]
        )
        result["maximum_detector_separation_deg"][output_index] = np.nanmax(
            arrays["detector_separation_deg"][local]
        )
        result["maximum_absolute_time_offset_seconds"][output_index] = np.nanmax(
            np.abs(arrays["time_offset_seconds"][local])
        )
        result["fast_measurement_valid_exposure_fraction"][output_index] = weighted_mean(
            arrays["fast_measurement_valid"][local], weight
        )
        result["dombeck_available_exposure_fraction"][output_index] = weighted_mean(
            arrays["dombeck_available"][local], weight
        )
        result["dombeck_intense_exposure_fraction"][output_index] = weighted_mean(
            arrays["dombeck_intense_precipitation"][local], weight
        )

        if np.any(valid_weight):
            aggregate = aggregate_spectra(
                arrays["fast_energy_ev"][local][valid_weight],
                arrays["fast_energy_width_ev"][local][valid_weight],
                arrays["earthward_differential_number_flux_cm2_s_ev"][local][valid_weight],
                arrays["anti_earthward_differential_number_flux_cm2_s_ev"][local][valid_weight],
                weight[valid_weight], energy_count,
            )
            energy[output_index], width[output_index] = aggregate[:2]
            earthward[output_index], anti_earthward[output_index] = aggregate[2:]

        for name, source_name in ATMOSPHERE_FIELDS.items():
            if name == "altitude_km":
                continue
            atmosphere[name][output_index] = weighted_mean(
                arrays[source_name][local], weight
            )

    result["fast_energy_ev"] = energy
    result["fast_energy_width_ev"] = width
    result["earthward_differential_number_flux_cm2_s_ev"] = earthward
    result["anti_earthward_differential_number_flux_cm2_s_ev"] = anti_earthward

    full = moments(energy, width, earthward)
    threshold = moments(energy, width, earthward, 500)
    anti = moments(energy, width, anti_earthward)
    result["earthward_number_flux_full_cm2_s"] = full[0]
    result["earthward_energy_flux_full_mw_m2"] = full[1]
    result["earthward_mean_energy_full_kev"] = full[2]
    result["earthward_number_flux_500ev_cm2_s"] = threshold[0]
    result["earthward_energy_flux_500ev_mw_m2"] = threshold[1]
    result["earthward_mean_energy_500ev_kev"] = threshold[2]
    result["anti_earthward_number_flux_full_cm2_s"] = anti[0]
    result["anti_earthward_energy_flux_full_mw_m2"] = anti[1]
    net_number = full[0] - anti[0]
    net_energy = full[1] - anti[1]
    result["net_earthward_number_flux_full_cm2_s"] = net_number
    result["net_earthward_energy_flux_full_mw_m2"] = net_energy
    result["net_earthward_mean_energy_full_kev"] = np.divide(
        net_energy / EV_CM2_S_TO_MW_M2, net_number,
        out=np.full(pixel_count, np.nan),
        where=(net_number > 0) & (net_energy > 0),
    ) / 1000

    result.update(calculate_conductance(
        energy, width, earthward,
        result["footprint_magnetic_field_nt"], atmosphere,
        conductance_batch_size,
    ))
    return result


def write_pixel_batch(output, start, result):
    stop = start + len(result["year"])
    for name, array in result.items():
        output.variables[f"pixel_{name}"][start:stop] = array
    return stop


#%% Complete compiler

def inspect_inputs(filenames, annotations):
    information = []
    sample_schema = None
    image_schema = None
    energy_count = None
    msis_altitude = None
    years = set()
    source_files = []

    for filename in filenames:
        with Dataset(filename) as data:
            if data.getncattr("processing_complete") != "true":
                raise ValueError(f"incomplete conjunction file {filename}")
            year = int(data.getncattr("selected_year"))
            if year in years:
                raise ValueError(f"multiple conjunction files supplied for {year}")
            years.add(year)
            fields = scalar_fields(data)
            current_sample_schema = {
                name: schema_for_variable(data.variables[name]) for name in fields
            }
            current_image_schema = {
                name: current_sample_schema[name]
                for name in fields
                if name.startswith("image_")
                and name not in ("image_orbit", "image_frame", "image_time")
            }
            if sample_schema is None:
                sample_schema = current_sample_schema
                image_schema = current_image_schema
                energy_count = len(data.dimensions["energy_index"])
                msis_altitude = values(data.variables["msis_altitude_km"])
            else:
                if set(current_sample_schema) != set(sample_schema):
                    raise ValueError("conjunction files have different scalar fields")
                if set(current_image_schema) != set(image_schema):
                    raise ValueError("conjunction files have different IMAGE fields")
                energy_count = max(
                    energy_count, len(data.dimensions["energy_index"])
                )
                if not np.array_equal(
                    values(data.variables["msis_altitude_km"]), msis_altitude
                ):
                    raise ValueError("conjunction files have different MSIS altitude grids")

            validate_input_fields(data)
            status = annotation_status(data, year, annotations)
            records, run_starts, run_stops = pixel_runs(data)
            local_sources = np.asarray(
                data.variables["fast_source_file"][:], dtype=str
            )
            information.append({
                "path": Path(filename), "year": year,
                "sample_count": len(data.dimensions["sample"]),
                "pixel_count": len(run_starts),
                "annotation_status": status,
                "records": records,
                "run_starts": run_starts, "run_stops": run_stops,
                "source_offset": len(source_files),
                "source_count": len(local_sources),
            })
            source_files.extend(local_sources.tolist())
    return information, sample_schema, image_schema, energy_count, msis_altitude, source_files


def build_dataset(args):
    annotations = read_annotations(args.annotations)
    info, sample_schema, image_schema, energy_count, msis_altitude, source_files = inspect_inputs(
        args.input, annotations
    )
    sample_count = sum(item["sample_count"] for item in info)
    pixel_count = sum(item["pixel_count"] for item in info)
    if pixel_count == 0:
        raise ValueError("conjunction files contain no FAST measurements")

    root_attrs = {
        "title": "Compact FAST--IMAGE conjunction analysis dataset",
        "processing_complete": "true",
        "selected_years": json.dumps(sorted(item["year"] for item in info)),
        "annotations_attached": "true" if annotations else "false",
        "aggregation_selection": "none; every conjunction sample is aggregated",
        "aggregation_unit": "IMAGE orbit/frame/detector row/detector column",
        "aggregation_weight": "dombeck_window_span_seconds divided by dombeck_window_record_count",
        "mixed_energy_grid_rule": "piecewise-constant conservative bin-overlap regridding to the finest covering table in each pixel",
        "dombeck_role": "diagnostic intense/weak stratification only; never a quality selection",
        "sample_count": sample_count,
        "pixel_count": pixel_count,
    }

    output = Path(args.output)
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"{output} exists; pass --overwrite to replace it")
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(output.suffix + ".partial")
    if partial.exists():
        partial.unlink()

    try:
        with create_output(
            partial, sample_count, pixel_count, energy_count, msis_altitude,
            source_files, sample_schema, image_schema, root_attrs,
        ) as target:
            sample_start = 0
            pixel_start = 0
            for item in info:
                with Dataset(item["path"]) as source:
                    sample_stop = sample_start + item["sample_count"]
                    for name in tqdm(
                        sample_schema, desc=f"Copy {item['year']} sample scalars",
                        unit="field",
                    ):
                        array = values(source.variables[name])
                        if name == "fast_source_file_id":
                            array = array + item["source_offset"]
                        target.variables[f"sample_{name}"][sample_start:sample_stop] = array
                    target.variables["sample_year"][sample_start:sample_stop] = item["year"]
                    target.variables["sample_annotation_status"][sample_start:sample_stop] = item["annotation_status"]

                    starts = item["run_starts"]
                    stops = item["run_stops"]
                    group_start = 0
                    progress = tqdm(
                        total=len(starts),
                        desc=f"Aggregate {item['year']} IMAGE pixels", unit="pixel",
                    )
                    while group_start < len(starts):
                        first_record_position = starts[group_start]
                        group_stop = group_start + 1
                        while (
                            group_stop < len(starts)
                            and stops[group_stop - 1] - first_record_position
                            < args.record_batch_size
                        ):
                            group_stop += 1
                        selected_start = starts[group_start]
                        selected_stop = stops[group_stop - 1]
                        records = item["records"][selected_start:selected_stop]
                        local_starts = starts[group_start:group_stop] - selected_start
                        local_stops = stops[group_start:group_stop] - selected_start
                        result = aggregate_pixel_batch(
                            source, records, local_starts, local_stops,
                            item["year"], energy_count, msis_altitude,
                            tuple(image_schema), args.conductance_batch_size,
                            item["annotation_status"][records],
                        )
                        pixel_start = write_pixel_batch(target, pixel_start, result)
                        progress.update(group_stop - group_start)
                        group_start = group_stop
                    progress.close()
                    sample_start = sample_stop

            if sample_start != sample_count or pixel_start != pixel_count:
                raise RuntimeError("compiled output counts do not reconcile")
        partial.replace(output)
    except BaseException:
        if partial.exists():
            partial.unlink()
        raise
    print(f"Saved {sample_count:,} samples and {pixel_count:,} pixels: {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, nargs="+", required=True)
    parser.add_argument("--annotations", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--record-batch-size", type=int, default=DEFAULT_RECORD_BATCH_SIZE)
    parser.add_argument("--conductance-batch-size", type=int, default=DEFAULT_CONDUCTANCE_BATCH_SIZE)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    for filename in args.input:
        if not filename.is_file():
            raise FileNotFoundError(filename)
    if args.annotations is not None and not args.annotations.is_file():
        raise FileNotFoundError(args.annotations)
    if args.record_batch_size < 1 or args.conductance_batch_size < 1:
        raise ValueError("batch sizes must be positive")
    build_dataset(args)


if __name__ == "__main__":
    main()
