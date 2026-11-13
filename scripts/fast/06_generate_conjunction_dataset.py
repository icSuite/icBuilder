"""Build one unannotated FAST--IMAGE detector-level conjunction dataset."""

#%% Imports

import argparse
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path

from icreader import open_product
from icbuilder.footprints import (
    infer_footprints,
    overlap_mapping,
    point_in_quadrilaterals,
)
from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd
from tqdm import tqdm
import xarray as xr


#%% Product contract

GRID_ID = "image_apex_130km_46x46_v1"
REFERENCE_HEIGHT_KM = 130.0
DEFAULT_FRAME_SUPPORT_SECONDS = 60.0

FAST_SCALAR_FIELDS = (
    "fast_geodetic_latitude",
    "fast_geodetic_longitude",
    "fast_altitude_km",
    "footprint_geodetic_latitude",
    "footprint_geodetic_longitude",
    "footprint_mlat",
    "footprint_mlt",
    "loss_cone_degrees",
    "magnetic_field_ratio_footprint_to_fast",
    "footprint_magnetic_field_nt",
    "spacecraft_potential_v",
    "fast_validity_flag",
    "fast_data_quality_flag",
    "earthward_number_flux_conservative_cm2_s",
    "earthward_energy_flux_conservative_mw_m2",
    "earthward_mean_energy_conservative_kev",
    "earthward_number_flux_centre_cm2_s",
    "earthward_energy_flux_centre_mw_m2",
    "earthward_mean_energy_centre_kev",
    "earthward_number_flux_inclusive_cm2_s",
    "earthward_energy_flux_inclusive_mw_m2",
    "earthward_mean_energy_inclusive_kev",
    "earthward_mean_energy_all_valid_kev",
    "earthward_number_flux_all_valid_cm2_s",
    "earthward_energy_flux_all_valid_mw_m2",
    "anti_earthward_energy_flux_mw_m2",
    "net_earthward_number_flux_cm2_s",
    "net_earthward_energy_flux_mw_m2",
    "net_earthward_mean_energy_kev",
    "earthward_to_anti_earthward_energy_flux_ratio",
    "earthward_angle_bins",
    "electron_pedersen_conductance_rees_s",
    "electron_hall_conductance_rees_s",
    "conductance_spectrum_channel_count",
    "robinson_effective_mean_energy_kev",
    "robinson_effective_energy_flux_mw_m2",
    "robinson_effective_valid",
)

FAST_SPECTRAL_FIELDS = (
    "fast_energy_ev",
    "fast_energy_width_ev",
    "earthward_differential_number_flux_cm2_s_ev",
    "anti_earthward_differential_number_flux_cm2_s_ev",
)

IMAGE_P2_FIELDS = (
    "glat", "glon", "mlat", "mlt", "sza", "dza",
    "method_quality_weight",
    "wic_coverage", "si12_coverage", "si13_coverage",
    "wic_valid", "si12_valid", "si13_valid", "method_valid",
    "Ep_model", "Ep", "dEp", "Fp", "dFp",
    "wic_corrected", "dwic_corrected",
    "si13_corrected", "dsi13_corrected",
    "R", "dR", "E0", "dE0", "Fe", "dFe", "varE0Fe",
    "Ep_clipping_flag",
)

IMAGE_P2_OPTIONAL_FIELDS = (
    "si12", "dsi12", "E0_response_range_flag",
)

IMAGE_P3_FIELDS = (
    "P", "H", "dP", "dH",
    "conductance_valid", "conductance_uncertainty_valid",
)

IMAGE_P3_OPTIONAL_FIELDS = (
    "fang_geographic_latitude_at_reference_height",
    "fang_geographic_longitude_at_reference_height",
    "fang_magnetic_field_150km_microtesla",
    "fang_evaluated_energy",
    "fang_within_model_domain",
    "fang_energy_clipped_low",
    "fang_uncertainty_straddles_low_energy_boundary",
)

BOOLEAN_FIELDS = {
    "robinson_effective_valid",
    "image_wic_valid", "image_si12_valid", "image_si13_valid",
    "image_method_valid", "image_Ep_clipping_flag",
    "image_conductance_valid", "image_conductance_uncertainty_valid",
    "image_fang_within_model_domain",
    "image_fang_energy_clipped_low",
    "image_fang_uncertainty_straddles_low_energy_boundary",
}

INTEGER_FAST_FIELDS = {
    "fast_validity_flag", "fast_data_quality_flag", "earthward_angle_bins",
    "conductance_spectrum_channel_count",
}

P2_COMMON_ATTRIBUTES = (
    "method", "proton_flux_source", "proton_energy_model", "count_source",
    "source_preprocessing_label", "smoothed", "smoothing_method",
    "wic_smoothing_width_pixels", "si12_smoothing_width_pixels",
    "si13_smoothing_width_pixels",
)

P3_COMMON_ATTRIBUTES = ("conductance_model", "conductance_uncertainty_method")

FAST_COMMON_ATTRIBUTES = (
    "particle_dataset", "footprint_height_km", "minimum_energy_ev",
    "conductance_msis_cadence_seconds",
    "conductance_space_weather_indices_sha256", "full_spectrum_definition",
)

CANDIDATE_COLUMNS = (
    "fast_source_file_id", "fast_source_record", "fast_time_ns",
    "image_orbit", "image_frame", "image_time_ns", "time_offset_seconds",
    "cs_row", "cs_column", "cs_wic_coverage", "cs_si13_coverage",
    "detector_row", "detector_column", "detector_separation_deg",
    "detector_footprint_match_count",
)


#%% Small utilities

def collection_signature(files):
    """Identify a file collection without hashing its large data arrays."""

    identities = []
    for filename in sorted(Path(file).resolve() for file in files):
        stat = filename.stat()
        identities.append(f"{filename} {stat.st_size} {stat.st_mtime_ns}")
    return sha256("\n".join(identities).encode()).hexdigest()


def package_version(package):
    """Return an installed package version without making it a requirement."""

    try:
        return version(package)
    except PackageNotFoundError:
        return "unknown"


def orbit_number(filename):
    """Read the four-digit orbit number from an ``or_NNNN.nc`` filename."""

    try:
        return int(Path(filename).stem.split("_")[-1])
    except ValueError as error:
        raise ValueError(f"Cannot read orbit number from {filename}") from error


def index_orbit_files(path, selected_orbits=None):
    """Index one explicitly selected IMAGE product directory by orbit."""

    selected = None if selected_orbits is None else set(selected_orbits)
    files = {}
    for filename in sorted(Path(path).glob("or_*.nc")):
        orbit = orbit_number(filename)
        if selected is None or orbit in selected:
            files[orbit] = filename
    if not files:
        raise FileNotFoundError(f"No selected IMAGE orbit files found in {path}")
    return files


def select_image_orbits_for_year(cs_files, year, support_seconds):
    """Keep IMAGE orbits whose frame support overlaps the selected year."""

    year_start_ns, year_stop_ns = year_bounds(year)
    support_ns = int(round(float(support_seconds) * 1e9))
    selected = []
    files = sorted(cs_files.items())
    for orbit, filename in tqdm(
        files, desc=f"Index IMAGE orbits for {year}", unit="orbit"
    ):
        with open_product(filename) as product:
            times_ns = datetime64_ns(product.time)
        if (
            int(times_ns[-1]) + support_ns >= year_start_ns
            and int(times_ns[0]) - support_ns < year_stop_ns
        ):
            selected.append(orbit)
    return selected


def datetime64_ns(values):
    """Return nanosecond integer timestamps without timezone ambiguity."""

    return np.asarray(values, dtype="datetime64[ns]").astype(np.int64)


def year_bounds(year):
    """Return the half-open UTC nanosecond interval for one year."""

    start = np.datetime64(f"{year:04d}-01-01", "ns")
    stop = np.datetime64(f"{year + 1:04d}-01-01", "ns")
    return int(start.astype(np.int64)), int(stop.astype(np.int64))


def attributes_equal(left, right):
    """Compare scalar or short array product attributes."""

    try:
        return bool(np.array_equal(np.asarray(left), np.asarray(right)))
    except Exception:
        return left == right


def prefixed_image_name(name):
    return f"image_{name}"


def temporal_candidate_indices(times_ns, image_time_ns, support_seconds):
    """Return records inside the closed IMAGE-frame support interval."""

    support_ns = int(round(float(support_seconds) * 1e9))
    times_ns = np.asarray(times_ns, dtype=np.int64)
    return np.flatnonzero(np.abs(times_ns - int(image_time_ns)) <= support_ns)


#%% FAST file index

def inspect_fast_file(filename, year_start_ns, year_stop_ns):
    """Validate one processed FAST product and return its usable time range."""

    with Dataset(filename) as data:
        attrs = {name: data.getncattr(name) for name in data.ncattrs()}
        if attrs.get("processing_complete") != "true":
            raise ValueError(f"{filename} is not a complete FAST archive product")
        if not np.isclose(
            float(attrs.get("footprint_height_km", np.nan)),
            REFERENCE_HEIGHT_KM,
        ):
            raise ValueError(f"{filename} does not use a 130-km FAST footprint")

        required = {"time", *FAST_SCALAR_FIELDS, *FAST_SPECTRAL_FIELDS}
        missing = required.difference(data.variables)
        if missing:
            raise ValueError(f"{filename} is missing {sorted(missing)}")

        time_variable = data.variables["time"]
        numeric_time = np.asarray(time_variable[:])
        if not numeric_time.size or np.any(np.diff(numeric_time) <= 0):
            raise ValueError(f"{filename} FAST times must be nonempty and unique")
        endpoints = num2date(
            numeric_time[[0, -1]],
            time_variable.units,
            getattr(time_variable, "calendar", "standard"),
            only_use_cftime_datetimes=False,
        )
        endpoint_ns = np.array(
            [np.datetime64(value.isoformat(), "ns") for value in endpoints]
        ).astype(np.int64)
        if endpoint_ns[1] < year_start_ns or endpoint_ns[0] >= year_stop_ns:
            return None

        if endpoint_ns[0] >= year_start_ns and endpoint_ns[1] < year_stop_ns:
            selected_count = numeric_time.size
        else:
            decoded = num2date(
                numeric_time,
                time_variable.units,
                getattr(time_variable, "calendar", "standard"),
                only_use_cftime_datetimes=False,
            )
            times_ns = np.array(
                [np.datetime64(value.isoformat(), "ns") for value in decoded]
            ).astype(np.int64)
            selected_count = int(np.count_nonzero(
                (times_ns >= year_start_ns) & (times_ns < year_stop_ns)
            ))
            if not selected_count:
                return None

        energy_count = len(data.dimensions["energy_index"])
        return {
            "path": Path(filename),
            "start_ns": max(int(endpoint_ns[0]), year_start_ns),
            "stop_ns": min(int(endpoint_ns[1]), year_stop_ns - 1),
            "energy_count": energy_count,
            "record_count": int(numeric_time.size),
            "selected_record_count": selected_count,
            "attrs": attrs,
        }


def index_fast_files(path, year, selected_files=None):
    """Index complete FAST products overlapping the requested year."""

    year_start_ns, year_stop_ns = year_bounds(year)
    if selected_files:
        filenames = [Path(file) for file in selected_files]
    else:
        filenames = sorted(Path(path).rglob("*.nc"))
    if not filenames:
        raise FileNotFoundError(f"No FAST NetCDF files found in {path}")

    indexed = []
    for filename in tqdm(
        filenames, desc=f"Index FAST products for {year}", unit="file"
    ):
        product = inspect_fast_file(filename, year_start_ns, year_stop_ns)
        if product is not None:
            indexed.append(product)
    if not indexed:
        raise ValueError(f"No FAST products overlap {year}")

    reference = {
        name: indexed[0]["attrs"].get(name) for name in FAST_COMMON_ATTRIBUTES
    }
    for product in indexed[1:]:
        current = {
            name: product["attrs"].get(name) for name in FAST_COMMON_ATTRIBUTES
        }
        require_common_attributes(
            reference, current, FAST_COMMON_ATTRIBUTES, product["path"]
        )
    indexed.sort(key=lambda item: (item["start_ns"], str(item["path"])))
    for file_id, product in enumerate(indexed):
        product["file_id"] = file_id
    return indexed


def load_fast_interval(product, start_ns, stop_ns):
    """Load the small FAST time/footprint slice needed by one IMAGE orbit."""

    with xr.open_dataset(product["path"]) as data:
        times = datetime64_ns(data.time.values)
        records = np.flatnonzero((times >= start_ns) & (times <= stop_ns))
        if not records.size:
            return None
        return {
            "source_file_id": product["file_id"],
            "record": records.astype(np.int32),
            "time_ns": times[records],
            "mlat": np.asarray(data.footprint_mlat.values, dtype=float)[records],
            "mlt": np.asarray(data.footprint_mlt.values, dtype=float)[records],
        }


#%% CS support gate

def cs_support_gate(grid, frame, mlat, mlt, wic_coverage, si13_coverage,
                    wic_corrected, si13_corrected):
    """Return supported point indices and their CS cells and coverage."""

    mlat = np.asarray(mlat, dtype=float)
    mlt = np.asarray(mlt, dtype=float)
    finite = np.isfinite(mlat) & np.isfinite(mlt)
    finite_indices = np.flatnonzero(finite)
    counts = {
        "nonfinite_footprint_count": int(np.count_nonzero(~finite)),
        "outside_cs_grid_count": 0,
        "cs_coverage_data_rejection_count": 0,
    }
    if not finite_indices.size:
        return np.empty(0, dtype=int), {}, counts

    longitude = mlt[finite] * 15.0
    latitude = mlat[finite]
    inside = np.asarray(grid.ingrid(longitude, latitude), dtype=bool)
    counts["outside_cs_grid_count"] = int(np.count_nonzero(~inside))
    inside_indices = finite_indices[inside]
    if not inside_indices.size:
        return np.empty(0, dtype=int), {}, counts

    row, column = grid.bin_index(longitude[inside], latitude[inside])
    row = np.asarray(row, dtype=np.int16)
    column = np.asarray(column, dtype=np.int16)
    wic_cov = np.asarray(wic_coverage[frame, row, column], dtype=float)
    si13_cov = np.asarray(si13_coverage[frame, row, column], dtype=float)
    wic = np.asarray(wic_corrected[frame, row, column], dtype=float)
    si13 = np.asarray(si13_corrected[frame, row, column], dtype=float)
    supported = (
        np.isfinite(wic_cov) & (wic_cov > 0)
        & np.isfinite(si13_cov) & (si13_cov > 0)
        & np.isfinite(wic) & np.isfinite(si13)
    )
    counts["cs_coverage_data_rejection_count"] = int(
        np.count_nonzero(~supported)
    )
    values = {
        "row": row[supported],
        "column": column[supported],
        "wic_coverage": wic_cov[supported].astype(np.float32),
        "si13_coverage": si13_cov[supported].astype(np.float32),
    }
    return inside_indices[supported], values, counts


#%% Detector-footprint matching

def unit_vectors(mlat, mlt):
    """Convert Modified Apex latitude and MLT to unit-sphere vectors."""

    latitude = np.deg2rad(np.asarray(mlat, dtype=float))
    longitude = np.deg2rad(np.asarray(mlt, dtype=float) * 15.0)
    cos_latitude = np.cos(latitude)
    return np.column_stack((
        cos_latitude * np.cos(longitude),
        cos_latitude * np.sin(longitude),
        np.sin(latitude),
    ))


def make_detector_footprint_matcher(image_mlat, image_mlt, grid):
    """Prepare the same footprint lookup used by the DMSP collector."""

    corners, _ = infer_footprints(image_mlat, image_mlt, grid)
    mapping, _ = overlap_mapping(image_mlat, image_mlt, grid)
    corners = corners.reshape(-1, 4, 2)
    centres = np.mean(corners, axis=1)
    return corners, centres, mapping


def candidate_footprints(mapping, grid_shape, row, column):
    """Return footprints intersecting the CS cell or one neighbour."""

    ny, nx = grid_shape
    candidates = []
    for nearby_row in range(max(row - 1, 0), min(row + 2, ny)):
        for nearby_column in range(max(column - 1, 0), min(column + 2, nx)):
            target = nearby_row * nx + nearby_column
            start, stop = mapping.indptr[target:target + 2]
            candidates.append(mapping.indices[start:stop])
    if not candidates:
        return np.empty(0, dtype=int)
    return np.unique(np.concatenate(candidates))


def match_detector_footprints(matcher, grid, cs_row, cs_column, point_mlat,
                              point_mlt, image_mlat, image_mlt):
    """Find the detector quadrilateral containing each FAST footprint."""

    corners, centres, mapping = matcher
    longitude = np.asarray(point_mlt, dtype=float) * 15.0
    latitude = np.asarray(point_mlat, dtype=float)
    point_x, point_y = grid.projection.geo2cube(
        longitude, latitude, set_points_off_cube_to_nan=True
    )
    point_x = np.asarray(point_x, dtype=float)
    point_y = np.asarray(point_y, dtype=float)

    flat_pixel = np.full(point_x.size, -1, dtype=int)
    match_count = np.zeros(point_x.size, dtype=np.int16)
    for sample in range(point_x.size):
        candidates = candidate_footprints(
            mapping, grid.shape, int(cs_row[sample]), int(cs_column[sample])
        )
        inside = point_in_quadrilaterals(
            point_x[sample], point_y[sample], corners[candidates]
        )
        containing = candidates[inside]
        match_count[sample] = containing.size
        if not containing.size:
            continue
        if containing.size == 1:
            flat_pixel[sample] = containing[0]
        else:
            distance_squared = (
                (centres[containing, 0] - point_x[sample]) ** 2
                + (centres[containing, 1] - point_y[sample]) ** 2
            )
            flat_pixel[sample] = containing[np.argmin(distance_squared)]

    contained = flat_pixel >= 0
    row = np.full(point_x.size, -1, dtype=np.int16)
    column = np.full(point_x.size, -1, dtype=np.int16)
    separation = np.full(point_x.size, np.nan, dtype=np.float32)
    if np.any(contained):
        matched_row, matched_column = np.unravel_index(
            flat_pixel[contained], image_mlat.shape
        )
        row[contained] = matched_row
        column[contained] = matched_column
        point_vectors = unit_vectors(latitude[contained], point_mlt[contained])
        pixel_vectors = unit_vectors(
            image_mlat[matched_row, matched_column],
            image_mlt[matched_row, matched_column],
        )
        cosine = np.sum(point_vectors * pixel_vectors, axis=1)
        separation[contained] = np.rad2deg(
            np.arccos(np.clip(cosine, -1, 1))
        ).astype(np.float32)
    return row, column, separation, match_count


#%% IMAGE collection validation

def require_common_attributes(reference, current, names, filename):
    """Reject a collection that mixes scientific configurations."""

    for name in names:
        if name not in reference and name not in current:
            continue
        if not attributes_equal(reference.get(name), current.get(name)):
            raise ValueError(f"{filename} mixes IMAGE configuration {name!r}")


def validate_orbit_products(orbit, cs, precipitation, conductance,
                            common_identity=None):
    """Validate one aligned CS/Product-2/Product-3 orbit trio."""

    if cs.product_type != "precipitation_cs" or cs.grid_id != GRID_ID:
        raise ValueError(f"Orbit {orbit:04d} does not use the required CS Product 2")
    if precipitation.product_type != "precipitation_detector":
        raise ValueError(f"Orbit {orbit:04d} is not detector Product 2")
    if conductance.product_type != "conductance_detector":
        raise ValueError(f"Orbit {orbit:04d} is not detector Product 3")
    if precipitation.method != "image_ratio":
        raise ValueError(
            f"Orbit {orbit:04d} uses {precipitation.method!r}; the schema-1 "
            "conjunction gate requires image_ratio Product 2"
        )

    for product in (cs, precipitation, conductance):
        if not np.isclose(float(product.reference_height_km), REFERENCE_HEIGHT_KM):
            raise ValueError(f"Orbit {orbit:04d} has a non-130-km IMAGE product")
    if precipitation.shape != conductance.shape:
        raise ValueError(f"Orbit {orbit:04d} Product-2/3 shapes differ")
    if not np.array_equal(precipitation.time, conductance.time):
        raise ValueError(f"Orbit {orbit:04d} Product-2/3 times differ")
    if not np.array_equal(cs.time, precipitation.time):
        raise ValueError(f"Orbit {orbit:04d} CS/detector times differ")

    identity = {
        **{
            name: getattr(precipitation, name, precipitation.attrs.get(name))
            for name in P2_COMMON_ATTRIBUTES
        },
        **{
            name: getattr(conductance, name, conductance.attrs.get(name))
            for name in P3_COMMON_ATTRIBUTES
        },
        "grid_id": cs.grid_id,
        "reference_height_km": float(cs.reference_height_km),
    }
    if common_identity is not None:
        require_common_attributes(
            common_identity, identity, tuple(identity), precipitation.filename
        )
    return identity


#%% Candidate discovery

def empty_counts():
    return {
        "temporal_candidate_count": 0,
        "nonfinite_footprint_count": 0,
        "outside_cs_grid_count": 0,
        "cs_coverage_data_rejection_count": 0,
        "detector_footprint_miss_count": 0,
        "exact_detector_data_miss_count": 0,
        "ambiguous_detector_footprint_count": 0,
        "successful_candidate_count": 0,
        "duplicate_frame_candidates_removed": 0,
        "final_sample_count": 0,
    }


def add_counts(total, addition):
    for name, value in addition.items():
        total[name] = total.get(name, 0) + int(value)


def discover_orbit_candidates(orbit, cs, precipitation, fast_slices,
                              support_seconds, counts):
    """Find fully supported FAST matches for one IMAGE orbit."""

    if not fast_slices:
        return []

    wic_coverage = np.asarray(cs.read("wic_coverage"), dtype=float)
    si13_coverage = np.asarray(cs.read("si13_coverage"), dtype=float)
    cs_wic = np.asarray(cs.read("wic_corrected"), dtype=float)
    cs_si13 = np.asarray(cs.read("si13_corrected"), dtype=float)

    image_times_ns = datetime64_ns(cs.time)
    candidates = []

    for frame, image_time_ns in enumerate(image_times_ns):
        frame_parts = []
        for fast in fast_slices:
            selected = temporal_candidate_indices(
                fast["time_ns"], image_time_ns, support_seconds
            )
            if not selected.size:
                continue
            counts["temporal_candidate_count"] += int(selected.size)
            supported, cs_values, gate_counts = cs_support_gate(
                cs.grid, frame, fast["mlat"][selected], fast["mlt"][selected],
                wic_coverage, si13_coverage, cs_wic, cs_si13,
            )
            add_counts(counts, gate_counts)
            if not supported.size:
                continue
            chosen = selected[supported]
            frame_parts.append({
                "source_file_id": np.full(
                    chosen.size, fast["source_file_id"], dtype=np.int32
                ),
                "record": fast["record"][chosen],
                "time_ns": fast["time_ns"][chosen],
                "mlat": fast["mlat"][chosen],
                "mlt": fast["mlt"][chosen],
                "cs_row": cs_values["row"],
                "cs_column": cs_values["column"],
                "cs_wic_coverage": cs_values["wic_coverage"],
                "cs_si13_coverage": cs_values["si13_coverage"],
            })
        if not frame_parts:
            continue

        frame_data = {
            name: np.concatenate([part[name] for part in frame_parts])
            for name in frame_parts[0]
        }
        image_mlat = np.asarray(precipitation.read("mlat", frame), dtype=float)
        image_mlt = np.asarray(precipitation.read("mlt", frame), dtype=float)
        matcher = make_detector_footprint_matcher(
            image_mlat, image_mlt, cs.grid
        )
        row, column, separation, match_count = match_detector_footprints(
            matcher, cs.grid, frame_data["cs_row"], frame_data["cs_column"],
            frame_data["mlat"], frame_data["mlt"],
            image_mlat, image_mlt,
        )
        contained = match_count > 0
        counts["detector_footprint_miss_count"] += int(
            np.count_nonzero(~contained)
        )
        counts["ambiguous_detector_footprint_count"] += int(
            np.count_nonzero(match_count > 1)
        )
        if not np.any(contained):
            continue

        detector_wic = np.asarray(
            precipitation.read("wic_corrected", frame), dtype=float
        )
        detector_si13 = np.asarray(
            precipitation.read("si13_corrected", frame), dtype=float
        )
        exact_data = np.zeros(contained.size, dtype=bool)
        exact_data[contained] = (
            np.isfinite(detector_wic[row[contained], column[contained]])
            & np.isfinite(detector_si13[row[contained], column[contained]])
        )
        counts["exact_detector_data_miss_count"] += int(
            np.count_nonzero(contained & ~exact_data)
        )
        keep = contained & exact_data
        if not np.any(keep):
            continue

        count = int(np.count_nonzero(keep))
        counts["successful_candidate_count"] += count
        candidate = pd.DataFrame({
            "fast_source_file_id": frame_data["source_file_id"][keep],
            "fast_source_record": frame_data["record"][keep],
            "fast_time_ns": frame_data["time_ns"][keep],
            "image_orbit": np.full(count, orbit, dtype=np.int32),
            "image_frame": np.full(count, frame, dtype=np.int32),
            "image_time_ns": np.full(count, image_time_ns, dtype=np.int64),
            "time_offset_seconds": (
                (frame_data["time_ns"][keep] - image_time_ns) / 1e9
            ).astype(np.float32),
            "cs_row": frame_data["cs_row"][keep],
            "cs_column": frame_data["cs_column"][keep],
            "cs_wic_coverage": frame_data["cs_wic_coverage"][keep],
            "cs_si13_coverage": frame_data["cs_si13_coverage"][keep],
            "detector_row": row[keep],
            "detector_column": column[keep],
            "detector_separation_deg": separation[keep],
            "detector_footprint_match_count": match_count[keep],
        })
        candidates.append(candidate)
    return candidates


def select_unique_matches(candidates):
    """Select one supported IMAGE match for every FAST source record."""

    if candidates.empty:
        result = candidates.copy()
        result["successful_frame_candidate_count"] = np.array([], dtype=np.int16)
        return result

    keys = ["fast_source_file_id", "fast_source_record"]
    candidate_count = candidates.groupby(keys)[keys[0]].transform("size")
    ranked = candidates.assign(
        successful_frame_candidate_count=candidate_count.astype(np.int16),
        absolute_time_offset=np.abs(candidates.time_offset_seconds),
    ).sort_values([
        *keys, "absolute_time_offset", "image_time_ns", "image_orbit",
        "image_frame", "detector_row", "detector_column",
    ], kind="stable")
    selected = ranked.drop_duplicates(keys, keep="first")
    selected = selected.drop(columns="absolute_time_offset")
    return selected.sort_values(
        ["fast_time_ns", "fast_source_file_id", "fast_source_record"],
        kind="stable",
    ).reset_index(drop=True)


#%% NetCDF writer

def source_variable_metadata(fast_products, first_orbit_paths):
    """Collect variable attributes from representative source products."""

    metadata = {}
    with xr.open_dataset(fast_products[0]["path"]) as fast:
        for name in (*FAST_SCALAR_FIELDS, *FAST_SPECTRAL_FIELDS):
            metadata[name] = dict(fast[name].attrs)
    with open_product(first_orbit_paths[0]) as precipitation:
        for name in (*IMAGE_P2_FIELDS, *IMAGE_P2_OPTIONAL_FIELDS):
            if name in precipitation.variable_attrs:
                metadata[prefixed_image_name(name)] = dict(
                    precipitation.variable_attrs[name]
                )
    with open_product(first_orbit_paths[1]) as conductance:
        for name in (*IMAGE_P3_FIELDS, *IMAGE_P3_OPTIONAL_FIELDS):
            if name in conductance.variable_attrs:
                metadata[prefixed_image_name(name)] = dict(
                    conductance.variable_attrs[name]
                )
    return metadata


def create_numeric_variable(output, name, dimensions, dtype="f4", fill=None,
                            attrs=None):
    """Create one compressed numeric variable and copy safe attributes."""

    kwargs = {"zlib": True, "complevel": 4, "shuffle": True}
    if fill is not None:
        kwargs["fill_value"] = fill
    variable = output.createVariable(name, dtype, dimensions, **kwargs)
    for attribute, value in (attrs or {}).items():
        if attribute == "_FillValue":
            continue
        try:
            variable.setncattr(attribute, value)
        except (TypeError, ValueError):
            variable.setncattr(attribute, str(value))
    return variable


def create_output(filename, sample_count, energy_count, source_files, metadata,
                  root_attrs):
    """Create the fixed schema-1 output before orbit-wise filling."""

    output = Dataset(filename, "w", format="NETCDF4")
    output.createDimension("sample", sample_count if sample_count else None)
    output.createDimension("energy_index", energy_count)
    output.createDimension("fast_source_file", len(source_files))

    source_names = output.createVariable("fast_source_file", str, ("fast_source_file",))
    source_names[:] = np.asarray([str(path) for path in source_files], dtype=object)

    integer_variables = {
        "fast_source_file_id": "i4", "fast_source_record": "i4",
        "fast_energy_channel_count": "i2",
        "image_orbit": "i4", "image_frame": "i4", "wic_source_index": "i4",
        "cs_row": "i2", "cs_column": "i2", "detector_row": "i2",
        "detector_column": "i2", "detector_footprint_match_count": "i2",
        "successful_frame_candidate_count": "i2",
    }
    for name, dtype in integer_variables.items():
        create_numeric_variable(output, name, ("sample",), dtype, fill=-1)
    output.variables["fast_energy_channel_count"].units = "1"
    output.variables["fast_energy_channel_count"].description = (
        "Original energy_index size in the selected FAST source product; "
        "spectral values beyond this count are fill padding"
    )

    for name in ("fast_time", "image_time", "image_Kp_interval_start"):
        variable = create_numeric_variable(
            output, name, ("sample",), "f8", fill=np.nan
        )
        variable.units = "seconds since 1970-01-01 00:00:00"
        variable.calendar = "standard"

    scalar_floats = (
        "time_offset_seconds", "cs_wic_coverage", "cs_si13_coverage",
        "detector_separation_deg",
    )
    for name in scalar_floats:
        create_numeric_variable(output, name, ("sample",), "f4", fill=np.nan)

    for name in FAST_SCALAR_FIELDS:
        if name in INTEGER_FAST_FIELDS:
            create_numeric_variable(output, name, ("sample",), "i4", fill=-1,
                                    attrs=metadata.get(name))
        elif name in BOOLEAN_FIELDS:
            create_numeric_variable(output, name, ("sample",), "i1", fill=-1,
                                    attrs=metadata.get(name))
        else:
            create_numeric_variable(output, name, ("sample",), "f8", fill=np.nan,
                                    attrs=metadata.get(name))

    for name in FAST_SPECTRAL_FIELDS:
        create_numeric_variable(
            output, name, ("sample", "energy_index"), "f8", fill=np.nan,
            attrs=metadata.get(name),
        )

    p2_names = (*IMAGE_P2_FIELDS, *IMAGE_P2_OPTIONAL_FIELDS)
    p3_names = (*IMAGE_P3_FIELDS, *IMAGE_P3_OPTIONAL_FIELDS)
    for source_name in (*p2_names, *p3_names):
        output_name = prefixed_image_name(source_name)
        if output_name in output.variables:
            continue
        if output_name in BOOLEAN_FIELDS:
            dtype, fill = "i1", -1
        elif source_name == "E0_response_range_flag":
            dtype, fill = "i1", -1
        elif source_name in ("glat", "glon", "mlat", "mlt"):
            dtype, fill = "f8", np.nan
        else:
            dtype, fill = "f4", np.nan
        create_numeric_variable(
            output, output_name, ("sample",), dtype, fill,
            attrs=metadata.get(output_name),
        )
    create_numeric_variable(output, "image_Kp", ("sample",), "f4", fill=np.nan)

    for name, value in root_attrs.items():
        try:
            output.setncattr(name, value)
        except (TypeError, ValueError):
            output.setncattr(name, json.dumps(value, sort_keys=True))
    return output


def write_candidate_columns(output, matches):
    """Write match identity columns shared by all source groups."""

    mapping = {
        "fast_time": matches.fast_time_ns.to_numpy(dtype=np.int64) / 1e9,
        "image_time": matches.image_time_ns.to_numpy(dtype=np.int64) / 1e9,
    }
    for name in (
        "fast_source_file_id", "fast_source_record", "image_orbit",
        "image_frame", "cs_row", "cs_column", "cs_wic_coverage",
        "cs_si13_coverage", "detector_row", "detector_column",
        "detector_separation_deg", "detector_footprint_match_count",
        "successful_frame_candidate_count", "time_offset_seconds",
    ):
        mapping[name] = matches[name].to_numpy()
    for name, values in mapping.items():
        output.variables[name][:] = values


def write_fast_fields(output, matches, fast_products):
    """Gather each FAST file only for its selected source records."""

    products_by_id = {
        product["file_id"]: product for product in fast_products
    }
    grouped = matches.groupby("fast_source_file_id", sort=True)
    for file_id, selected in tqdm(
        grouped, total=grouped.ngroups, desc="Copy FAST fields", unit="file"
    ):
        product = products_by_id[int(file_id)]
        positions = selected.index.to_numpy(dtype=int)
        records = selected.fast_source_record.to_numpy(dtype=int)
        with xr.open_dataset(product["path"]) as fast:
            energy_count = product["energy_count"]
            output.variables["fast_energy_channel_count"][positions] = energy_count
            for name in FAST_SCALAR_FIELDS:
                output.variables[name][positions] = np.asarray(
                    fast[name].values
                )[records]
            for name in FAST_SPECTRAL_FIELDS:
                output.variables[name][positions, :energy_count] = np.asarray(
                    fast[name].values
                )[records, :]


def gathered_field(product, name, frames, rows, columns):
    """Read contiguous selected frame runs and gather pointwise values."""

    gathered = None
    unique_frames = np.unique(frames)
    run_starts = np.r_[0, np.flatnonzero(np.diff(unique_frames) > 1) + 1]
    run_stops = np.r_[run_starts[1:], unique_frames.size]

    for first, stop in zip(run_starts, run_stops):
        start_frame = int(unique_frames[first])
        stop_frame = int(unique_frames[stop - 1]) + 1
        selected = (frames >= start_frame) & (frames < stop_frame)
        images = np.asarray(product.read(name, slice(start_frame, stop_frame)))
        if gathered is None:
            gathered = np.empty(frames.size, dtype=images.dtype)
        frame_offset = frames[selected] - start_frame
        gathered[selected] = images[
            frame_offset, rows[selected], columns[selected]
        ]
    return gathered


def write_image_fields(output, matches, p2_files, p3_files):
    """Read and write IMAGE data one orbit and one field at a time."""

    grouped = matches.groupby("image_orbit", sort=True)
    for orbit, selected in tqdm(
        grouped,
        total=grouped.ngroups,
        desc="Copy IMAGE detector fields",
        unit="orbit",
    ):
        orbit = int(orbit)
        positions = selected.index.to_numpy(dtype=int)
        frames = selected.image_frame.to_numpy(dtype=int)
        rows = selected.detector_row.to_numpy(dtype=int)
        columns = selected.detector_column.to_numpy(dtype=int)

        with open_product(p2_files[orbit]) as precipitation, open_product(
            p3_files[orbit]
        ) as conductance:
            output.variables["wic_source_index"][positions] = (
                precipitation.wic_source_index[frames]
            )
            output.variables["image_Kp"][positions] = precipitation.Kp[frames]
            kp_start = datetime64_ns(precipitation.Kp_interval_start)
            output.variables["image_Kp_interval_start"][positions] = (
                kp_start[frames] / 1e9
            )

            for name in IMAGE_P2_FIELDS:
                output_name = prefixed_image_name(name)
                output.variables[output_name][positions] = gathered_field(
                    precipitation, name, frames, rows, columns
                )
            for name in IMAGE_P2_OPTIONAL_FIELDS:
                if name in precipitation.fields:
                    output_name = prefixed_image_name(name)
                    output.variables[output_name][positions] = gathered_field(
                        precipitation, name, frames, rows, columns
                    )
            for name in IMAGE_P3_FIELDS:
                output_name = prefixed_image_name(name)
                output.variables[output_name][positions] = gathered_field(
                    conductance, name, frames, rows, columns
                )
            for name in IMAGE_P3_OPTIONAL_FIELDS:
                if name in conductance.fields:
                    output_name = prefixed_image_name(name)
                    output.variables[output_name][positions] = gathered_field(
                        conductance, name, frames, rows, columns
                    )


def validate_saved_output(filename, expected_samples, expected_energy_count,
                          signatures):
    """Reopen the result and validate its descriptor and source identity."""

    with Dataset(filename) as data:
        if data.product_type != "fast_image_conjunctions":
            raise ValueError("saved output has the wrong product_type")
        if int(data.schema_version) != 1:
            raise ValueError("saved output has the wrong schema_version")
        if len(data.dimensions["sample"]) != expected_samples:
            raise ValueError("saved output has the wrong sample count")
        if len(data.dimensions["energy_index"]) != expected_energy_count:
            raise ValueError("saved output has the wrong energy dimension")
        for name, signature in signatures.items():
            if data.getncattr(name) != signature:
                raise ValueError(f"saved output has the wrong {name}")
        if expected_samples:
            file_id = np.asarray(data.variables["fast_source_file_id"][:])
            record = np.asarray(data.variables["fast_source_record"][:])
            keys = np.column_stack((file_id, record))
            if np.unique(keys, axis=0).shape[0] != expected_samples:
                raise ValueError("saved output contains duplicate FAST records")


def output_is_current(filename, year, support_seconds, signatures):
    """Return whether an existing output has the requested source contract."""

    if not Path(filename).is_file():
        return False
    try:
        with Dataset(filename) as data:
            required = {
                "fast_time", "image_time", "fast_source_file_id",
                "fast_source_record", "fast_energy_channel_count",
                "image_orbit", "image_frame",
                "detector_row", "detector_column", "image_E0", "image_Fe",
                "image_P", "image_H",
                "earthward_differential_number_flux_cm2_s_ev",
            }
            return (
                data.product_type == "fast_image_conjunctions"
                and int(data.schema_version) == 1
                and data.getncattr("processing_complete") == "true"
                and int(data.selected_year) == year
                and np.isclose(float(data.frame_support_seconds), support_seconds)
                and required.issubset(data.variables)
                and "sample" in data.dimensions
                and "energy_index" in data.dimensions
                and all(data.getncattr(name) == value for name, value in signatures.items())
            )
    except Exception:
        return False


#%% Complete run

def build_dataset(args):
    """Discover all supported matches and write the consolidated dataset."""

    fast_products = index_fast_files(args.fast_path, args.year, args.fast_file)
    cs_files = index_orbit_files(args.precipitation_cs_path, args.orbit)
    p2_files = index_orbit_files(args.precipitation_detector_path, args.orbit)
    p3_files = index_orbit_files(args.conductance_detector_path, args.orbit)
    if set(cs_files) != set(p2_files) or set(cs_files) != set(p3_files):
        raise ValueError("Selected CS, Product-2 and Product-3 orbit sets differ")
    selected_image_orbits = select_image_orbits_for_year(
        cs_files, args.year, args.frame_support_seconds
    )
    if not selected_image_orbits:
        raise ValueError(f"No selected IMAGE orbits overlap {args.year}")
    cs_files = {orbit: cs_files[orbit] for orbit in selected_image_orbits}
    p2_files = {orbit: p2_files[orbit] for orbit in selected_image_orbits}
    p3_files = {orbit: p3_files[orbit] for orbit in selected_image_orbits}

    signatures = {
        "source_fast_collection_signature": collection_signature(
            [product["path"] for product in fast_products]
        ),
        "source_precipitation_cs_collection_signature": collection_signature(
            cs_files.values()
        ),
        "source_precipitation_detector_collection_signature": collection_signature(
            p2_files.values()
        ),
        "source_conductance_detector_collection_signature": collection_signature(
            p3_files.values()
        ),
    }
    if output_is_current(
        args.output, args.year, args.frame_support_seconds, signatures
    ) and not args.overwrite:
        print(f"Current conjunction dataset already exists: {args.output}")
        return
    if Path(args.output).exists() and not args.overwrite:
        raise FileExistsError(f"{args.output} exists; pass --overwrite to replace it")

    counts = empty_counts()
    candidates = []
    common_identity = None
    support_ns = int(round(args.frame_support_seconds * 1e9))
    year_start_ns, year_stop_ns = year_bounds(args.year)

    orbit_progress = tqdm(
        sorted(cs_files), desc="Discover FAST--IMAGE matches", unit="orbit"
    )
    for orbit in orbit_progress:
        with open_product(cs_files[orbit]) as cs, open_product(
            p2_files[orbit]
        ) as precipitation, open_product(p3_files[orbit]) as conductance:
            identity = validate_orbit_products(
                orbit, cs, precipitation, conductance, common_identity
            )
            if common_identity is None:
                common_identity = identity

            image_times_ns = datetime64_ns(cs.time)
            interval_start = max(int(image_times_ns[0] - support_ns), year_start_ns)
            interval_stop = min(int(image_times_ns[-1] + support_ns), year_stop_ns - 1)
            if interval_stop < interval_start:
                continue
            overlapping = [
                product for product in fast_products
                if product["stop_ns"] >= interval_start
                and product["start_ns"] <= interval_stop
            ]
            fast_slices = [
                load_fast_interval(product, interval_start, interval_stop)
                for product in overlapping
            ]
            fast_slices = [part for part in fast_slices if part is not None]
            candidates.extend(discover_orbit_candidates(
                orbit, cs, precipitation, fast_slices,
                args.frame_support_seconds, counts,
            ))
            orbit_progress.set_postfix(
                orbit=f"{orbit:04d}",
                supported=f"{counts['successful_candidate_count']:,}",
            )

    if candidates:
        all_candidates = pd.concat(candidates, ignore_index=True)
    else:
        all_candidates = pd.DataFrame({
            name: np.array([], dtype=float) for name in CANDIDATE_COLUMNS
        })
    matches = select_unique_matches(all_candidates)
    counts["duplicate_frame_candidates_removed"] = len(all_candidates) - len(matches)
    counts["final_sample_count"] = len(matches)

    source_energy_counts = sorted({
        product["energy_count"] for product in fast_products
    })
    energy_count = max(source_energy_counts)
    first_orbit = min(cs_files)
    metadata = source_variable_metadata(
        fast_products, (p2_files[first_orbit], p3_files[first_orbit])
    )
    root_attrs = {
        "title": "FAST--IMAGE detector-level conjunction samples",
        "product_type": "fast_image_conjunctions",
        "schema_version": 1,
        "representation": "detector_point_samples",
        "workflow": Path(__file__).name,
        "selected_year": args.year,
        "frame_support_seconds": float(args.frame_support_seconds),
        "temporal_support_interval": "closed plus/minus frame_support_seconds",
        "unique_match_rule": (
            "minimum absolute offset; earlier IMAGE time; lower orbit, frame, "
            "detector row and detector column"
        ),
        "cs_support_gate": (
            "inside precipitation_cs grid; finite positive WIC/SI13 coverage; "
            "finite CS wic_corrected/si13_corrected"
        ),
        "detector_support_gate": (
            "FAST 130-km point inside inferred WIC detector footprint with "
            "finite detector wic_corrected/si13_corrected"
        ),
        "annotation_status": "not used; complete unannotated conjunction population",
        "processing_complete": "true",
        "fast_spectral_padding": (
            "Mixed FAST energy_index widths use the maximum output width; "
            "shorter source spectra are right-padded with fill values and "
            "fast_energy_channel_count retains the original width"
        ),
        "reference_height_km": REFERENCE_HEIGHT_KM,
        "source_fast_root": str(Path(args.fast_path).resolve()),
        "source_precipitation_cs_root": str(
            Path(args.precipitation_cs_path).resolve()
        ),
        "source_precipitation_detector_root": str(
            Path(args.precipitation_detector_path).resolve()
        ),
        "source_conductance_detector_root": str(
            Path(args.conductance_detector_path).resolve()
        ),
        "source_signature_basis": (
            "sorted absolute paths, file sizes and modification times; not full hashes"
        ),
        "icBuilder_version": package_version("icbuilder"),
        "icReader_version": package_version("icreader"),
        "icPhysics_version": package_version("icphysics"),
        "NumPy_version": package_version("numpy"),
        "xarray_version": package_version("xarray"),
        "netCDF4_version": package_version("netCDF4"),
        "secsy_version": package_version("secsy"),
        "fast_input_file_count": len(fast_products),
        "fast_source_energy_index_counts": json.dumps(source_energy_counts),
        "fast_output_energy_index_count": energy_count,
        "fast_input_record_count": sum(
            product["selected_record_count"] for product in fast_products
        ),
        "image_orbit_count": len(cs_files),
        "image_configuration": json.dumps(
            common_identity, sort_keys=True, default=str
        ),
        "fast_configuration": json.dumps(
            {
                name: fast_products[0]["attrs"].get(name)
                for name in FAST_COMMON_ATTRIBUTES
            },
            sort_keys=True,
            default=str,
        ),
        **signatures,
        **counts,
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(output.suffix + ".partial")
    if partial.exists():
        partial.unlink()
    try:
        with create_output(
            partial, len(matches), energy_count,
            [product["path"].resolve() for product in fast_products],
            metadata, root_attrs,
        ) as target:
            write_candidate_columns(target, matches)
            write_fast_fields(target, matches, fast_products)
            write_image_fields(target, matches, p2_files, p3_files)
        validate_saved_output(partial, len(matches), energy_count, signatures)
        partial.replace(output)
    except Exception:
        if partial.exists():
            partial.unlink()
        raise

    print(f"Saved {len(matches):,} unique FAST--IMAGE samples: {output}")


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
    parser.add_argument(
        "--frame-support-seconds", type=float,
        default=DEFAULT_FRAME_SUPPORT_SECONDS,
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    for path in (
        args.fast_path, args.precipitation_cs_path,
        args.precipitation_detector_path, args.conductance_detector_path,
    ):
        if not path.is_dir():
            raise FileNotFoundError(path)
    if not np.isfinite(args.frame_support_seconds) or args.frame_support_seconds < 0:
        raise ValueError("--frame-support-seconds must be finite and nonnegative")
    if args.fast_file:
        missing = [file for file in args.fast_file if not file.is_file()]
        if missing:
            raise FileNotFoundError(missing[0])

    build_dataset(args)


if __name__ == "__main__":
    main()
