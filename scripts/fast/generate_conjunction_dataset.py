"""Build one unannotated FAST--IMAGE conjunction dataset for one year.

The compact IMAGE precipitation CS product is used only as an inexpensive
spatial and temporal gate. Detector Product 2 supplies the matched IMAGE
pixel and precipitation fields; detector Product 3 supplies conductance.
The precipitation method does not affect the matching.
"""

#%% Imports and fields

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from multiprocessing import get_context
from pathlib import Path

from icreader import open_product
from icbuilder.footprints import infer_footprints, overlap_mapping, point_in_quadrilaterals
from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd
from tqdm import tqdm
import xarray as xr


GRID_ID = "image_apex_130km_46x46_v1"
REFERENCE_HEIGHT_KM = 130.0
DEFAULT_FRAME_SUPPORT_SECONDS = 60.0

P2_FIELDS = ("glat", "glon", "mlat", "mlt", "sza", "dza",
             "wic_corrected", "dwic_corrected",
             "si12", "dsi12",
             "si13_corrected", "dsi13_corrected",
             "R", "dR", "E0", "dE0", "Fe", "dFe", "varE0Fe",
             "Ep", "dEp", "Fp", "dFp",
             "wic_valid", "si12_valid", "si13_valid", "method_valid",
             "method_quality_weight", "E0_response_range_flag", "Ep_clipping_flag")

P3_FIELDS = ("P", "H", "dP", "dH", "conductance_valid", "conductance_uncertainty_valid")

INTEGER_IMAGE_FIELDS = {"wic_valid", "si12_valid", "si13_valid", "method_valid",
                        "E0_response_range_flag", "Ep_clipping_flag",
                        "conductance_valid", "conductance_uncertainty_valid"}

CANDIDATE_COLUMNS = ("fast_source_file_id", "fast_product_record", "fast_source_record",
                     "fast_time_ns", "image_orbit", "image_frame", "image_time_ns",
                     "time_offset_seconds", "cs_row", "cs_column", "cs_wic_coverage",
                     "detector_row", "detector_column", "detector_separation_deg",
                     "detector_footprint_match_count")


#%% Small utilities

def orbit_number(filename):
    try:
        return int(Path(filename).stem.split("_")[-1])
    except ValueError as error:
        raise ValueError(f"Cannot read orbit number from {filename}") from error


def index_orbit_files(path, selected_orbits=None):
    selected = None if selected_orbits is None else set(selected_orbits)
    files = {}
    for filename in sorted(Path(path).glob("or_*.nc")):
        orbit = orbit_number(filename)
        if selected is None or orbit in selected:
            files[orbit] = filename
    if not files:
        raise FileNotFoundError(f"No selected IMAGE orbit files found in {path}")
    return files


def datetime64_ns(values):
    return np.asarray(values, dtype="datetime64[ns]").astype(np.int64)


def year_bounds(year):
    start = np.datetime64(f"{year:04d}-01-01", "ns")
    stop = np.datetime64(f"{year + 1:04d}-01-01", "ns")
    return int(start.astype(np.int64)), int(stop.astype(np.int64))


def temporal_candidate_indices(times_ns, image_time_ns, support_seconds):
    support_ns = int(round(float(support_seconds) * 1e9))
    return np.flatnonzero(np.abs(np.asarray(times_ns, dtype=np.int64) - int(image_time_ns)) <= support_ns)


def inspect_image_orbit(filename):
    with Dataset(filename) as data:
        time_variable = data.variables["time"]
        record_count = len(time_variable)
        if record_count == 0:
            raise ValueError(f"IMAGE product has no times: {filename}")
        numeric_time = np.asarray(time_variable[[0, record_count - 1]])
        endpoints = num2date(numeric_time, time_variable.units,
                             getattr(time_variable, "calendar", "standard"),
                             only_use_cftime_datetimes=False)

    endpoints = np.array([
        np.datetime64(value.isoformat(), "ns") for value in endpoints
    ]).astype(np.int64)
    return int(endpoints[0]), int(endpoints[1])


def select_image_orbits_for_year(cs_files, year, support_seconds, workers=1,
                                 return_intervals=False):
    year_start, year_stop = year_bounds(year)
    support_ns = int(round(float(support_seconds) * 1e9))
    items = sorted(cs_files.items())
    intervals = {}

    progress_options = {
        "total": len(items),
        "desc": f"Index IMAGE orbits for {year}",
        "unit": "orbit",
    }
    if workers == 1:
        results = (
            (orbit, *inspect_image_orbit(filename))
            for orbit, filename in items
        )
        for orbit, start_ns, stop_ns in tqdm(results, **progress_options):
            if stop_ns + support_ns >= year_start and start_ns - support_ns < year_stop:
                intervals[orbit] = (start_ns, stop_ns)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(inspect_image_orbit, filename): orbit
                for orbit, filename in items
            }
            for future in tqdm(as_completed(futures), **progress_options):
                orbit = futures[future]
                start_ns, stop_ns = future.result()
                if stop_ns + support_ns >= year_start and start_ns - support_ns < year_stop:
                    intervals[orbit] = (start_ns, stop_ns)

    intervals = dict(sorted(intervals.items()))
    return intervals if return_intervals else list(intervals)


#%% Processed FAST files

def inspect_fast_file(filename, year_start_ns, year_stop_ns):
    with Dataset(filename) as data:
        if data.getncattr("processing_complete") != "true":
            raise ValueError(f"Incomplete FAST product: {filename}")
        if not np.isclose(float(data.getncattr("footprint_height_km")), REFERENCE_HEIGHT_KM):
            raise ValueError(f"FAST product is not mapped to 130 km: {filename}")

        time_variable = data.variables["time"]
        record_count = len(time_variable)
        if record_count == 0:
            raise ValueError(f"FAST times are not nonempty and increasing: {filename}")
        numeric_time = np.asarray(time_variable[[0, record_count - 1]])
        if record_count > 1 and numeric_time[1] <= numeric_time[0]:
            raise ValueError(f"FAST times are not nonempty and increasing: {filename}")
        endpoints = num2date(numeric_time, time_variable.units,
                             getattr(time_variable, "calendar", "standard"),
                             only_use_cftime_datetimes=False)

        endpoints = np.array([np.datetime64(value.isoformat(), "ns") for value in endpoints]).astype(np.int64)
        if endpoints[1] < year_start_ns or endpoints[0] >= year_stop_ns:
            return None

        altitude = np.asarray(data.variables["msis_altitude_km"][:], dtype=float)
        return {"path": Path(filename),
                "start_ns": max(int(endpoints[0]), year_start_ns),
                "stop_ns": min(int(endpoints[1]), year_stop_ns - 1),
                "record_count": record_count,
                "energy_count": len(data.dimensions["energy_index"]),
                "msis_altitude_km": altitude}


def index_fast_files(path, year, selected_files=None, workers=1):
    year_start, year_stop = year_bounds(year)
    filenames = ([Path(filename) for filename in selected_files] if selected_files else sorted(Path(path).rglob("*.nc")))
    if not filenames:
        raise FileNotFoundError(f"No FAST products found in {path}")

    products = []
    progress_options = {
        "total": len(filenames),
        "desc": f"Index FAST products for {year}",
        "unit": "file",
    }
    if workers == 1:
        results = (
            inspect_fast_file(filename, year_start, year_stop)
            for filename in filenames
        )
        for product in tqdm(results, **progress_options):
            if product is not None:
                products.append(product)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(inspect_fast_file, filename, year_start, year_stop)
                for filename in filenames
            ]
            for future in tqdm(as_completed(futures), **progress_options):
                product = future.result()
                if product is not None:
                    products.append(product)
    if not products:
        raise ValueError(f"No FAST products overlap {year}")

    reference_altitude = products[0]["msis_altitude_km"]
    for product in products[1:]:
        if not np.array_equal(product["msis_altitude_km"], reference_altitude):
            raise ValueError(f"MSIS altitude grid differs in {product['path']}")

    products.sort(key=lambda item: (item["start_ns"], str(item["path"])))
    for file_id, product in enumerate(products):
        product["file_id"] = file_id
    return products


def load_fast_interval(product, start_ns, stop_ns):
    with xr.open_dataset(product["path"]) as data:
        times = datetime64_ns(data.time.values)
        records = np.flatnonzero((times >= start_ns) & (times <= stop_ns))
        if not records.size:
            return None
        return {"source_file_id": product["file_id"],
                "product_record": records.astype(np.int32),
                "source_record": np.asarray(data.fast_source_record.values, dtype=np.int32)[records],
                "time_ns": times[records],
                "mlat": np.asarray(data.footprint_mlat.values, dtype=float)[records],
                "mlt": np.asarray(data.footprint_mlt.values, dtype=float)[records]}


#%% CS gate and detector-pixel matching

def cs_support_gate(grid, frame, mlat, mlt, wic_coverage, wic_corrected):
    """Select footprints with an available WIC measurement in the CS product."""

    mlat = np.asarray(mlat, dtype=float)
    mlt = np.asarray(mlt, dtype=float)
    finite = np.isfinite(mlat) & np.isfinite(mlt)
    finite_indices = np.flatnonzero(finite)
    counts = {"nonfinite_footprint_count": int(np.count_nonzero(~finite)),
              "outside_cs_grid_count": 0,
              "cs_wic_rejection_count": 0}

    if not finite_indices.size:
        return np.empty(0, dtype=int), {}, counts

    longitude = mlt[finite] * 15
    latitude = mlat[finite]
    inside = np.asarray(grid.ingrid(longitude, latitude), dtype=bool)
    counts["outside_cs_grid_count"] = int(np.count_nonzero(~inside))
    inside_indices = finite_indices[inside]
    if not inside_indices.size:
        return np.empty(0, dtype=int), {}, counts

    row, column = grid.bin_index(longitude[inside], latitude[inside])
    row = np.asarray(row, dtype=np.int16)
    column = np.asarray(column, dtype=np.int16)
    coverage = np.asarray(wic_coverage[frame, row, column], dtype=float)
    wic = np.asarray(wic_corrected[frame, row, column], dtype=float)
    supported = np.isfinite(wic)
    counts["cs_wic_rejection_count"] = int(np.count_nonzero(~supported))
    return inside_indices[supported], {"row": row[supported], "column": column[supported], "wic_coverage": coverage[supported].astype(np.float32)}, counts


def unit_vectors(mlat, mlt):
    latitude = np.deg2rad(np.asarray(mlat, dtype=float))
    longitude = np.deg2rad(np.asarray(mlt, dtype=float) * 15)
    cosine = np.cos(latitude)
    return np.column_stack((cosine * np.cos(longitude),
                            cosine * np.sin(longitude),
                            np.sin(latitude)))


def make_detector_footprint_matcher(image_mlat, image_mlt, grid):
    corners, _ = infer_footprints(image_mlat, image_mlt, grid)
    mapping, _ = overlap_mapping(image_mlat, image_mlt, grid)
    corners = corners.reshape(-1, 4, 2)
    return corners, np.mean(corners, axis=1), mapping


def candidate_footprints(mapping, grid_shape, row, column):
    ny, nx = grid_shape
    candidates = []
    for nearby_row in range(max(row - 1, 0), min(row + 2, ny)):
        for nearby_column in range(max(column - 1, 0), min(column + 2, nx)):
            cell = nearby_row * nx + nearby_column
            start, stop = mapping.indptr[cell:cell + 2]
            candidates.append(mapping.indices[start:stop])
    return np.unique(np.concatenate(candidates)) if candidates else np.empty(0, dtype=int)


def match_detector_footprints(matcher, grid, cs_row, cs_column, point_mlat, point_mlt, image_mlat, image_mlt):
    corners, centres, mapping = matcher
    point_x, point_y = grid.projection.geo2cube(np.asarray(point_mlt) * 15, point_mlat, set_points_off_cube_to_nan=True)
    point_x = np.asarray(point_x, dtype=float)
    point_y = np.asarray(point_y, dtype=float)

    flat_pixel = np.full(point_x.size, -1, dtype=int)
    match_count = np.zeros(point_x.size, dtype=np.int16)
    for sample in range(point_x.size):
        candidates = candidate_footprints(mapping, grid.shape, int(cs_row[sample]), int(cs_column[sample]))
        inside = point_in_quadrilaterals(point_x[sample], point_y[sample], corners[candidates])
        containing = candidates[inside]
        match_count[sample] = containing.size
        if containing.size == 1:
            flat_pixel[sample] = containing[0]
        elif containing.size > 1:
            distance = ((centres[containing, 0] - point_x[sample]) ** 2 + (centres[containing, 1] - point_y[sample]) ** 2)
            flat_pixel[sample] = containing[np.argmin(distance)]

    matched = flat_pixel >= 0
    row = np.full(point_x.size, -1, dtype=np.int16)
    column = np.full(point_x.size, -1, dtype=np.int16)
    separation = np.full(point_x.size, np.nan, dtype=np.float32)
    if np.any(matched):
        matched_row, matched_column = np.unravel_index(flat_pixel[matched], image_mlat.shape)
        row[matched] = matched_row
        column[matched] = matched_column
        point_vector = unit_vectors(np.asarray(point_mlat)[matched], np.asarray(point_mlt)[matched])
        pixel_vector = unit_vectors(image_mlat[matched_row, matched_column], image_mlt[matched_row, matched_column])
        cosine = np.sum(point_vector * pixel_vector, axis=1)
        separation[matched] = np.rad2deg(np.arccos(np.clip(cosine, -1, 1))).astype(np.float32)
    return row, column, separation, match_count


#%% Candidate discovery

def validate_aligned_orbit(orbit, cs, precipitation, conductance):
    """Check only the geometry and frame alignment needed for matching."""

    if cs.grid_id != GRID_ID:
        raise ValueError(f"Orbit {orbit:04d} uses a different CS grid")
    for product in (cs, precipitation, conductance):
        if not np.isclose(float(product.reference_height_km), REFERENCE_HEIGHT_KM):
            raise ValueError(f"Orbit {orbit:04d} is not mapped to 130 km")
    if precipitation.shape != conductance.shape:
        raise ValueError(f"Orbit {orbit:04d} Product 2/3 shapes differ")
    if not np.array_equal(cs.time, precipitation.time):
        raise ValueError(f"Orbit {orbit:04d} CS/Product 2 times differ")
    if not np.array_equal(precipitation.time, conductance.time):
        raise ValueError(f"Orbit {orbit:04d} Product 2/3 times differ")


def empty_counts():
    return {"temporal_candidate_count": 0,
            "nonfinite_footprint_count": 0,
            "outside_cs_grid_count": 0,
            "cs_wic_rejection_count": 0,
            "detector_footprint_miss_count": 0,
            "detector_wic_rejection_count": 0,
            "ambiguous_detector_footprint_count": 0,
            "successful_candidate_count": 0,
            "duplicate_frame_candidates_removed": 0,
            "final_sample_count": 0}


def discover_orbit_candidates(orbit, cs, precipitation, fast_slices, support_seconds, counts):
    if not fast_slices:
        return []

    wic_coverage = np.asarray(cs.read("wic_coverage"), dtype=float)
    cs_wic = np.asarray(cs.read("wic_corrected"), dtype=float)
    image_times = datetime64_ns(cs.time)
    candidates = []

    for frame, image_time in enumerate(image_times):
        parts = []
        for fast in fast_slices:
            nearby = temporal_candidate_indices(fast["time_ns"], image_time, support_seconds)
            if not nearby.size:
                continue
            counts["temporal_candidate_count"] += nearby.size
            supported, cs_values, gate_counts = cs_support_gate(cs.grid, frame, fast["mlat"][nearby], fast["mlt"][nearby], wic_coverage, cs_wic)
            for name, value in gate_counts.items():
                counts[name] += value
            if not supported.size:
                continue
            chosen = nearby[supported]
            parts.append({"source_file_id": np.full(chosen.size, fast["source_file_id"], dtype=np.int32),
                          "product_record": fast["product_record"][chosen],
                          "source_record": fast["source_record"][chosen],
                          "time_ns": fast["time_ns"][chosen],
                          "mlat": fast["mlat"][chosen],
                          "mlt": fast["mlt"][chosen],
                          "cs_row": cs_values["row"],
                          "cs_column": cs_values["column"],
                          "cs_wic_coverage": cs_values["wic_coverage"]})
        if not parts:
            continue

        frame_data = {name: np.concatenate([part[name] for part in parts]) for name in parts[0]}
        image_mlat = np.asarray(precipitation.read("mlat", frame), dtype=float)
        image_mlt = np.asarray(precipitation.read("mlt", frame), dtype=float)
        matcher = make_detector_footprint_matcher(image_mlat, image_mlt, cs.grid)
        row, column, separation, match_count = match_detector_footprints(matcher, cs.grid, frame_data["cs_row"], frame_data["cs_column"], frame_data["mlat"], frame_data["mlt"], image_mlat, image_mlt)

        contained = match_count > 0
        counts["detector_footprint_miss_count"] += int(np.count_nonzero(~contained))
        counts["ambiguous_detector_footprint_count"] += int(np.count_nonzero(match_count > 1))
        if not np.any(contained):
            continue

        detector_wic = np.asarray(precipitation.read("wic_corrected", frame), dtype=float)
        has_wic = np.zeros(contained.size, dtype=bool)
        has_wic[contained] = np.isfinite(detector_wic[row[contained], column[contained]])
        counts["detector_wic_rejection_count"] += int(np.count_nonzero(contained & ~has_wic))
        keep = contained & has_wic
        if not np.any(keep):
            continue

        count = int(np.count_nonzero(keep))
        counts["successful_candidate_count"] += count
        candidates.append(pd.DataFrame({"fast_source_file_id": frame_data["source_file_id"][keep],
                                        "fast_product_record": frame_data["product_record"][keep],
                                        "fast_source_record": frame_data["source_record"][keep],
                                        "fast_time_ns": frame_data["time_ns"][keep],
                                        "image_orbit": np.full(count, orbit, dtype=np.int32),
                                        "image_frame": np.full(count, frame, dtype=np.int32),
                                        "image_time_ns": np.full(count, image_time, dtype=np.int64),
                                        "time_offset_seconds": ((frame_data["time_ns"][keep] - image_time) / 1e9).astype(np.float32),
                                        "cs_row": frame_data["cs_row"][keep],
                                        "cs_column": frame_data["cs_column"][keep],
                                        "cs_wic_coverage": frame_data["cs_wic_coverage"][keep],
                                        "detector_row": row[keep],
                                        "detector_column": column[keep],
                                        "detector_separation_deg": separation[keep],
                                        "detector_footprint_match_count": match_count[keep]}))
    return candidates


def discover_orbit_file_candidates(orbit, cs_filename, p2_filename, p3_filename,
                                   fast_products, interval_start, interval_stop,
                                   support_seconds):
    counts = empty_counts()
    with open_product(cs_filename) as cs, \
            open_product(p2_filename) as precipitation, \
            open_product(p3_filename) as conductance:
        validate_aligned_orbit(orbit, cs, precipitation, conductance)
        fast_slices = [
            load_fast_interval(product, interval_start, interval_stop)
            for product in fast_products
        ]
        fast_slices = [part for part in fast_slices if part is not None]
        candidates = discover_orbit_candidates(
            orbit, cs, precipitation, fast_slices, support_seconds, counts
        )

    if candidates:
        candidates = pd.concat(candidates, ignore_index=True)
    else:
        candidates = pd.DataFrame({
            name: np.array([], dtype=float) for name in CANDIDATE_COLUMNS
        })
    return orbit, candidates, counts


def select_unique_matches(candidates):
    if candidates.empty:
        result = candidates.copy()
        result["successful_frame_candidate_count"] = np.array([], dtype=np.int16)
        return result

    keys = ["fast_source_file_id", "fast_product_record"]
    count = candidates.groupby(keys)[keys[0]].transform("size")
    ranked = candidates.assign(successful_frame_candidate_count=count.astype(np.int16), absolute_time_offset=np.abs(candidates.time_offset_seconds)).sort_values([*keys, "absolute_time_offset", "image_time_ns", "image_orbit", "image_frame", "detector_row", "detector_column"], kind="stable")
    selected = ranked.drop_duplicates(keys).drop(columns="absolute_time_offset")
    return selected.sort_values(["fast_time_ns", "fast_source_file_id", "fast_product_record"], kind="stable").reset_index(drop=True)


#%% Output file

def create_numeric_variable(output, name, dimensions, dtype, fill=None, attrs=None):
    options = {"zlib": True, "complevel": 4, "shuffle": True}
    if fill is not None:
        options["fill_value"] = fill
    variable = output.createVariable(name, dtype, dimensions, **options)
    for attribute, value in (attrs or {}).items():
        if attribute == "_FillValue":
            continue
        try:
            variable.setncattr(attribute, value)
        except (TypeError, ValueError):
            variable.setncattr(attribute, str(value))
    return variable


def available_image_fields(precipitation, conductance):
    p2 = tuple(name for name in P2_FIELDS if name in precipitation.fields)
    p3 = tuple(name for name in P3_FIELDS if name in conductance.fields)
    return p2, p3


def image_field_metadata(precipitation, conductance, p2_fields, p3_fields):
    metadata = {}
    for product, fields in ((precipitation, p2_fields), (conductance, p3_fields)):
        for name in fields:
            metadata[f"image_{name}"] = dict(product.variable_attrs.get(name, {}))
    return metadata


def create_output(filename, sample_count, energy_count, msis_altitude,
                  fast_products, p2_fields, p3_fields,
                  image_metadata, root_attrs):
    output = Dataset(filename, "w", format="NETCDF4")
    output.createDimension("sample", sample_count if sample_count else None)
    output.createDimension("energy_index", energy_count)
    output.createDimension("msis_altitude_km", len(msis_altitude))
    output.createDimension("fast_source_file", len(fast_products))

    energy = output.createVariable("energy_index", "i2", ("energy_index",))
    energy[:] = np.arange(energy_count)
    altitude = output.createVariable("msis_altitude_km", "f8", ("msis_altitude_km",))
    altitude[:] = msis_altitude
    altitude.units = "km"
    source_files = output.createVariable("fast_source_file", str, ("fast_source_file",))
    source_files[:] = np.asarray([str(product["path"].resolve()) for product in fast_products], dtype=object)

    for name, dtype in {"fast_source_file_id": "i4",
                        "fast_product_record": "i4",
                        "fast_energy_channel_count": "i2",
                        "image_orbit": "i4",
                        "image_frame": "i4",
                        "image_wic_source_index": "i4",
                        "cs_row": "i2",
                        "cs_column": "i2",
                        "detector_row": "i2",
                        "detector_column": "i2",
                        "detector_footprint_match_count": "i2",
                        "successful_frame_candidate_count": "i2",}.items():
        create_numeric_variable(output, name, ("sample",), dtype, fill=-1)

    for name in ("fast_time", "image_time", "image_Kp_interval_start"):
        variable = create_numeric_variable(output, name, ("sample",), "f8", fill=np.nan)
        variable.units = "seconds since 1970-01-01 00:00:00"
        variable.calendar = "standard"
    for name in ("time_offset_seconds", "cs_wic_coverage", "detector_separation_deg", "image_Kp"):
        create_numeric_variable(output, name, ("sample",), "f4", fill=np.nan)

    with Dataset(fast_products[0]["path"]) as source:
        for name, variable in source.variables.items():
            if not variable.dimensions or variable.dimensions[0] != "time" or name == "time":
                continue
            attrs = {key: variable.getncattr(key) for key in variable.ncattrs() if key != "_FillValue"}
            create_numeric_variable(output, name, ("sample", *variable.dimensions[1:]), variable.dtype, getattr(variable, "_FillValue", None), attrs)

    for name in (*p2_fields, *p3_fields):
        output_name = f"image_{name}"
        if name in INTEGER_IMAGE_FIELDS:
            dtype, fill = "i1", -1
        elif name in ("glat", "glon", "mlat", "mlt"):
            dtype, fill = "f8", np.nan
        else:
            dtype, fill = "f4", np.nan
        create_numeric_variable(output, output_name, ("sample",), dtype, fill, image_metadata.get(output_name))

    for name, value in root_attrs.items():
        try:
            output.setncattr(name, value)
        except (TypeError, ValueError):
            output.setncattr(name, json.dumps(value, sort_keys=True, default=str))
    return output


def write_candidate_fields(output, matches):
    output.variables["fast_time"][:] = matches.fast_time_ns.to_numpy(np.int64) / 1e9
    output.variables["image_time"][:] = matches.image_time_ns.to_numpy(np.int64) / 1e9
    for name in ("fast_source_file_id", "fast_product_record", "fast_source_record",
                 "image_orbit", "image_frame", "time_offset_seconds", "cs_row",
                 "cs_column", "cs_wic_coverage", "detector_row", "detector_column",
                 "detector_separation_deg", "detector_footprint_match_count",
                 "successful_frame_candidate_count"):
        output.variables[name][:] = matches[name].to_numpy()


def read_fast_fields(product, positions, records, source_records):
    first_record = int(records.min())
    stop_record = int(records.max()) + 1
    relative_records = records - first_record
    fields = {}
    dimensions = {}

    with Dataset(product["path"]) as source:
        raw_records = np.asarray(
            source.variables["fast_source_record"][first_record:stop_record]
        )[relative_records]
        if not np.array_equal(raw_records, source_records):
            raise ValueError(f"FAST source identity changed in {product['path']}")

        for name, variable in source.variables.items():
            if not variable.dimensions or variable.dimensions[0] != "time" or name == "time":
                continue
            source_key = (
                slice(first_record, stop_record),
                *[slice(None) for _ in variable.dimensions[1:]],
            )
            fields[name] = variable[source_key][relative_records]
            dimensions[name] = variable.dimensions[1:]

    return {"positions": positions,
            "energy_count": product["energy_count"],
            "fields": fields,
            "dimensions": dimensions}


def write_fast_batch(output, results):
    positions = np.concatenate([result["positions"] for result in results])
    order = np.argsort(positions)
    positions = positions[order]
    if positions[-1] - positions[0] + 1 == positions.size:
        output_positions = slice(int(positions[0]), int(positions[-1]) + 1)
    else:
        output_positions = positions
    energy_counts = np.concatenate([
        np.full(result["positions"].size, result["energy_count"], dtype=np.int16)
        for result in results
    ])[order]
    output.variables["fast_energy_channel_count"][output_positions] = energy_counts

    field_names = results[0]["fields"]
    for name in field_names:
        output_variable = output.variables[name]
        values = np.ma.masked_all(
            (len(order), *output_variable.shape[1:]), dtype=output_variable.dtype
        )
        start = 0
        for result in results:
            count = result["positions"].size
            trailing = result["dimensions"][name]
            target_key = [slice(start, start + count)]
            for dimension in trailing:
                if dimension == "energy_index":
                    target_key.append(slice(0, result["energy_count"]))
                else:
                    target_key.append(slice(None))
            values[tuple(target_key)] = result["fields"][name]
            start += count
        output_variable[output_positions] = values[order]


def write_fast_fields(output, matches, fast_products, workers=1, pool=None):
    products = {product["file_id"]: product for product in fast_products}
    groups = []
    for file_id, selected in matches.groupby("fast_source_file_id", sort=True):
        product = products[int(file_id)]
        product = {"path": product["path"],
                   "energy_count": product["energy_count"]}
        groups.append((product,
                       selected.index.to_numpy(dtype=int),
                       selected.fast_product_record.to_numpy(dtype=int),
                       selected.fast_source_record.to_numpy()))

    batch_size = max(8, workers * 2)
    progress = tqdm(total=len(groups), desc="Copy FAST fields", unit="file")
    if workers == 1:
        for start in range(0, len(groups), batch_size):
            batch = groups[start:start + batch_size]
            results = [read_fast_fields(*task) for task in batch]
            write_fast_batch(output, results)
            progress.update(len(batch))
    else:
        if pool is None:
            raise ValueError("Parallel FAST copying requires a worker pool")
        for start in range(0, len(groups), batch_size):
            batch = groups[start:start + batch_size]
            futures = [pool.submit(read_fast_fields, *task) for task in batch]
            results = [future.result() for future in as_completed(futures)]
            write_fast_batch(output, results)
            progress.update(len(batch))
    progress.close()


def gathered_field(product, name, frames, rows, columns):
    gathered = None
    unique_frames = np.unique(frames)
    starts = np.r_[0, np.flatnonzero(np.diff(unique_frames) > 1) + 1]
    stops = np.r_[starts[1:], unique_frames.size]
    for first, stop in zip(starts, stops):
        start_frame = int(unique_frames[first])
        stop_frame = int(unique_frames[stop - 1]) + 1
        selected = (frames >= start_frame) & (frames < stop_frame)
        images = np.asarray(product.read(name, slice(start_frame, stop_frame)))
        if gathered is None:
            gathered = np.empty(frames.size, dtype=images.dtype)
        gathered[selected] = images[
            frames[selected] - start_frame, rows[selected], columns[selected]
        ]
    return gathered


def read_image_fields(orbit, positions, frames, rows, columns,
                      p2_filename, p3_filename, p2_fields, p3_fields):
    fields = {}
    with open_product(p2_filename) as precipitation, \
            open_product(p3_filename) as conductance:
        if hasattr(precipitation, "wic_source_index"):
            fields["image_wic_source_index"] = precipitation.wic_source_index[frames]
        if hasattr(precipitation, "Kp"):
            fields["image_Kp"] = precipitation.Kp[frames]
        if hasattr(precipitation, "Kp_interval_start"):
            kp_time = datetime64_ns(precipitation.Kp_interval_start)
            fields["image_Kp_interval_start"] = kp_time[frames] / 1e9

        for product, names in ((precipitation, p2_fields),
                               (conductance, p3_fields)):
            for name in names:
                if name in product.fields:
                    fields[f"image_{name}"] = gathered_field(
                        product, name, frames, rows, columns
                    )
    return {"orbit": orbit, "positions": positions, "fields": fields}


def write_image_batch(output, results):
    positions = np.concatenate([result["positions"] for result in results])
    order = np.argsort(positions)
    positions = positions[order]
    if positions[-1] - positions[0] + 1 == positions.size:
        output_positions = slice(int(positions[0]), int(positions[-1]) + 1)
    else:
        output_positions = positions

    field_names = set()
    for result in results:
        field_names.update(result["fields"])
    for name in field_names:
        output_variable = output.variables[name]
        values = np.ma.masked_all(len(order), dtype=output_variable.dtype)
        start = 0
        for result in results:
            count = result["positions"].size
            if name in result["fields"]:
                values[start:start + count] = result["fields"][name]
            start += count
        output_variable[output_positions] = values[order]


def write_image_fields(output, matches, p2_files, p3_files,
                       p2_fields, p3_fields, workers=1, pool=None):
    groups = []
    for orbit, selected in matches.groupby("image_orbit", sort=True):
        orbit = int(orbit)
        groups.append((orbit,
                       selected.index.to_numpy(dtype=int),
                       selected.image_frame.to_numpy(dtype=int),
                       selected.detector_row.to_numpy(dtype=int),
                       selected.detector_column.to_numpy(dtype=int),
                       p2_files[orbit], p3_files[orbit],
                       p2_fields, p3_fields))

    batch_size = max(8, workers * 2)
    progress = tqdm(total=len(groups), desc="Copy IMAGE fields", unit="orbit")
    if workers == 1:
        for start in range(0, len(groups), batch_size):
            batch = groups[start:start + batch_size]
            results = [read_image_fields(*task) for task in batch]
            write_image_batch(output, results)
            progress.update(len(batch))
    else:
        if pool is None:
            raise ValueError("Parallel IMAGE copying requires a worker pool")
        for start in range(0, len(groups), batch_size):
            batch = groups[start:start + batch_size]
            futures = [pool.submit(read_image_fields, *task) for task in batch]
            results = [future.result() for future in as_completed(futures)]
            write_image_batch(output, results)
            progress.update(len(batch))
    progress.close()


def validate_saved_output(filename, sample_count, energy_count):
    with Dataset(filename) as data:
        if len(data.dimensions["sample"]) != sample_count:
            raise ValueError("Wrong conjunction sample count")
        if len(data.dimensions["energy_index"]) != energy_count:
            raise ValueError("Wrong conjunction energy dimension")
        if sample_count:
            keys = np.column_stack((data.variables["fast_source_file_id"][:], data.variables["fast_product_record"][:]))
            if np.unique(keys, axis=0).shape[0] != sample_count:
                raise ValueError("Saved output contains duplicate FAST records")


#%% Complete run

def build_dataset(args):
    workers = getattr(args, "workers", 1)
    fast_products = index_fast_files(
        args.fast_path, args.year, args.fast_file, workers
    )
    cs_files = index_orbit_files(args.precipitation_cs_path, args.orbit)
    p2_files = index_orbit_files(args.precipitation_detector_path, args.orbit)
    p3_files = index_orbit_files(args.conductance_detector_path, args.orbit)
    if set(cs_files) != set(p2_files) or set(cs_files) != set(p3_files):
        raise ValueError("Selected CS, Product 2 and Product 3 orbit sets differ")

    image_intervals = select_image_orbits_for_year(
        cs_files, args.year, args.frame_support_seconds, workers,
        return_intervals=True
    )

    if not image_intervals:
        raise ValueError(f"No IMAGE orbits overlap {args.year}")
    cs_files = {orbit: cs_files[orbit] for orbit in image_intervals}
    p2_files = {orbit: p2_files[orbit] for orbit in image_intervals}
    p3_files = {orbit: p3_files[orbit] for orbit in image_intervals}

    output = Path(args.output)
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"{output} exists; pass --overwrite to replace it")

    counts = empty_counts()
    candidates = []
    support_ns = int(round(args.frame_support_seconds * 1e9))
    year_start, year_stop = year_bounds(args.year)
    tasks = []
    for orbit, (image_start, image_stop) in image_intervals.items():
        interval_start = max(image_start - support_ns, year_start)
        interval_stop = min(image_stop + support_ns, year_stop - 1)
        overlapping = [
            {"file_id": product["file_id"],
             "path": product["path"],
             "start_ns": product["start_ns"],
             "stop_ns": product["stop_ns"]}
            for product in fast_products
            if product["stop_ns"] >= interval_start
            and product["start_ns"] <= interval_stop
        ]
        tasks.append((orbit, cs_files[orbit], p2_files[orbit], p3_files[orbit],
                      overlapping, interval_start, interval_stop,
                      args.frame_support_seconds))

    progress_options = {
        "total": len(tasks),
        "desc": "Discover FAST--IMAGE matches",
        "unit": "orbit",
    }
    if workers == 1:
        results = (discover_orbit_file_candidates(*task) for task in tasks)
        progress = tqdm(results, **progress_options)
        for orbit, orbit_candidates, orbit_counts in progress:
            if not orbit_candidates.empty:
                candidates.append(orbit_candidates)
            for name, value in orbit_counts.items():
                counts[name] += value
            progress.set_postfix(
                orbit=f"{orbit:04d}",
                supported=f"{counts['successful_candidate_count']:,}"
            )
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(discover_orbit_file_candidates, *task)
                for task in tasks
            ]
            progress = tqdm(
                as_completed(futures), **progress_options
            )
            for future in progress:
                orbit, orbit_candidates, orbit_counts = future.result()
                if not orbit_candidates.empty:
                    candidates.append(orbit_candidates)
                for name, value in orbit_counts.items():
                    counts[name] += value
                progress.set_postfix(
                    orbit=f"{orbit:04d}",
                    supported=f"{counts['successful_candidate_count']:,}"
                )

    all_candidates = (pd.concat(candidates, ignore_index=True) if candidates else pd.DataFrame({name: np.array([], dtype=float) for name in CANDIDATE_COLUMNS}))
    matches = select_unique_matches(all_candidates)
    counts["duplicate_frame_candidates_removed"] = len(all_candidates) - len(matches)
    counts["final_sample_count"] = len(matches)

    first_orbit = min(cs_files)
    with open_product(p2_files[first_orbit]) as precipitation, open_product(p3_files[first_orbit]) as conductance:
        p2_fields, p3_fields = available_image_fields(precipitation, conductance)
        image_metadata = image_field_metadata(precipitation, conductance, p2_fields, p3_fields)
        precipitation_method = getattr(precipitation, "method", precipitation.attrs.get("method", "unknown"))
        conductance_model = getattr(conductance, "conductance_model", conductance.attrs.get("conductance_model", "unknown"))

    energy_count = max(product["energy_count"] for product in fast_products)
    msis_altitude = fast_products[0]["msis_altitude_km"]
    root_attrs = {"title": "FAST--IMAGE detector-pixel conjunctions",
                  "processing_complete": "true",
                  "selected_year": args.year,
                  "frame_support_seconds": float(args.frame_support_seconds),
                  "reference_height_km": REFERENCE_HEIGHT_KM,
                  "conjunction_gate": "finite corrected WIC at CS and matched detector pixel",
                  "unique_match_rule": "minimum absolute time offset with deterministic ties",
                  "image_precipitation_method": str(precipitation_method),
                  "image_conductance_model": str(conductance_model),
                  "image_product_2_fields": json.dumps(p2_fields),
                  "image_product_3_fields": json.dumps(p3_fields),
                  "fast_input_file_count": len(fast_products),
                  "fast_input_record_count": sum(product["record_count"] for product in fast_products),
                  "fast_output_energy_index_count": energy_count,
                  "image_orbit_count": len(cs_files),
                  **counts}

    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(output.suffix + ".partial")
    if partial.exists():
        partial.unlink()
    copy_pool = (
        ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn"))
        if workers > 1 else None
    )
    try:
        with create_output(partial, len(matches), energy_count, msis_altitude, fast_products, p2_fields, p3_fields, image_metadata, root_attrs) as target:
            write_candidate_fields(target, matches)
            write_fast_fields(
                target, matches, fast_products, workers, copy_pool
            )
            write_image_fields(
                target, matches, p2_files, p3_files, p2_fields, p3_fields,
                workers, copy_pool
            )
        validate_saved_output(partial, len(matches), energy_count)
        partial.replace(output)
    except Exception:
        if partial.exists():
            partial.unlink()
        raise
    finally:
        if copy_pool is not None:
            copy_pool.shutdown(cancel_futures=True)
    print(f"Saved {len(matches):,} FAST--IMAGE conjunctions: {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fast-path", type=Path, required=True)
    parser.add_argument("--precipitation-cs-path", type=Path, required=True)
    parser.add_argument("--precipitation-detector-path", type=Path, required=True)
    parser.add_argument("--conductance-detector-path", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--orbit", type=int, nargs="+")
    parser.add_argument("--fast-file", type=Path, nargs="+")
    parser.add_argument("--frame-support-seconds", type=float, default=DEFAULT_FRAME_SUPPORT_SECONDS)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    for path in (args.fast_path, args.precipitation_cs_path, args.precipitation_detector_path, args.conductance_detector_path):
        if not path.is_dir():
            raise FileNotFoundError(path)
    if not np.isfinite(args.frame_support_seconds) or args.frame_support_seconds < 0:
        raise ValueError("--frame-support-seconds must be finite and nonnegative")
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    if args.fast_file:
        for filename in args.fast_file:
            if not filename.is_file():
                raise FileNotFoundError(filename)
    build_dataset(args)


if __name__ == "__main__":
    main()
