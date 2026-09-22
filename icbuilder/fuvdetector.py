"""Build coregistered IMAGE-FUV observations on WIC detector pixels.

This module defines the experimental detector-space Product 1 boundary.  WIC
sets the time axis and detector geometry.  SI12 and SI13 are matched
independently and mapped to WIC pixels with detector-footprint overlap.
"""

#%% Imports

from datetime import timedelta
import hashlib
from pathlib import Path

import apexpy
import numpy as np
from apexpy.helpers import subsol
from netCDF4 import Dataset, date2num, num2date

from .detector_coregistration import (
    MAX_COREG_ERROR_KM,
    MIN_SI_COVERAGE,
    make_si_mapping,
    make_wic_transform,
    map_si,
    map_si_variance,
)


#%% Product configuration

SCHEMA_VERSION = 2
TIME_TOLERANCE_SECONDS = 2.0
PREPROCESSING_LABEL = "fuvpy_bs_directional_v1"
IMAGE_FIELDS = {"WIC": "dgimg", "SI12": "dgimg", "SI13": "dgimg"}
UNSUBTRACTED_IMAGE_FIELDS = {"WIC": "img", "SI12": "img", "SI13": "img"}
QUALITY_WEIGHT_METHOD = "fuvpy dgweight on the native detector"
SOURCE_TIME_DECODING = "CF date units and calendar"
BS_MODEL_CONTRACT = {
    "spatial_model": "legacy_x_azimuth",
    "frame_quality_filter": "frame_quality >= 2",
    "measurement_variance": "img_variance",
    "variance_weighted": "true",
    "monotonic_requested": "true",
    "monotonic_fallback": "true",
    "spencer_taper": "enabled",
    "temporal_knot_spacing_target_minutes": 120.0,
    "directional_temporal_knot_spacing_target_minutes": 60.0,
    "tukey_value": 5.0,
    "damping_value": 1e-2,
    "convergence_tolerance": 1e-2,
    "directional_x_knots": np.asarray([-3.5, 0.0, 1.5, 3.5]),
}


def describe_time_match_rule(tolerance_seconds):
    """Describe the complete WIC-led matching rule for product provenance."""

    return (
        "WIC defines the time axis. Each SI channel is matched independently "
        f"to the nearest unused frame within +/- {tolerance_seconds:g} seconds. "
        "WIC frames are processed in stored order. A tie selects the earlier "
        "SI time and then the lower source index. Unmatched SI source index is -1."
    )


TIME_MATCH_RULE = describe_time_match_rule(TIME_TOLERANCE_SECONDS)
DETECTOR_NOISE_MODEL = (
    "img_variance from fuvpy is treated as independent detector-pixel count "
    "variance shared by the background-subtracted and unsubtracted counts. "
    "SI variance is propagated through the area-weighted mean as "
    "sum(a_i^2 variance_i) / sum(a_i)^2. Correlated, background-model, and "
    "coregistration uncertainties are not included."
)
COORDINATE_SYSTEM = "Modified Apex calculated from geodetic coordinates"
REFERENCE_HEIGHT_KM = 130.0


COREGISTRATION_INTEGER_FIELDS = (
    "valid_source_centres",
    "valid_source_footprints",
    "accepted_target_pixels",
)
COREGISTRATION_FLOAT_FIELDS = (
    "roundtrip_median_km",
    "roundtrip_95th_km",
    "coverage_median",
    "coverage_95th",
    "coverage_maximum",
)


#%% Source loading and coordinates

def _as_float(values):
    """Convert a NetCDF variable to a float array with masked values as NaN."""

    values = values[:]
    if np.ma.isMaskedArray(values):
        values = values.filled(np.nan)
    return np.asarray(values, dtype=float)


def file_sha256(filename, chunk_size=8 * 1024 * 1024):
    """Return the SHA-256 digest of one source product without loading it whole."""

    digest = hashlib.sha256()
    with Path(filename).open("rb") as stream:
        for chunk in iter(lambda: stream.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_identity(filename):
    """Return stable content identity plus filesystem provenance."""

    filename = Path(filename)
    stat = filename.stat()
    return {
        "source_file": str(filename),
        "source_sha256": file_sha256(filename),
        "source_size_bytes": int(stat.st_size),
        "source_mtime_ns": int(stat.st_mtime_ns),
    }


def _source_times(nc):
    """Decode source times from the CF metadata attached to ``date``."""

    date = nc.variables["date"]
    encoded = date[:]
    if encoded.ndim != 1:
        raise ValueError("source date must contain one time per frame")
    if np.any(np.ma.getmaskarray(encoded)):
        raise ValueError("source date must not contain missing values")
    if "units" not in date.ncattrs():
        raise ValueError("source date is missing CF units")

    values = np.asarray(encoded, dtype=float)
    if np.any(~np.isfinite(values)):
        raise ValueError("source date must contain finite values")

    calendar = getattr(date, "calendar", "standard")
    try:
        decoded = num2date(
            values,
            date.units,
            calendar,
            only_use_cftime_datetimes=False,
        )
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("source date has invalid CF time metadata") from error

    return np.asarray(decoded, dtype=object)


def _validate_background_model(nc, filename):
    """Require the fuvpy BS configuration named by the product label."""

    if "dgmodel" not in nc.variables:
        raise ValueError(f"{filename} is missing source field: dgmodel")
    model = nc.variables["dgmodel"]
    for name, expected in BS_MODEL_CONTRACT.items():
        if name not in model.ncattrs():
            raise ValueError(f"{filename} dgmodel is missing metadata: {name}")
        actual = model.getncattr(name)
        if isinstance(expected, str):
            matches = str(actual) == expected
        else:
            matches = np.allclose(actual, expected, rtol=0, atol=1e-12)
        if not matches:
            raise ValueError(
                f"{filename} dgmodel metadata {name} does not match "
                f"{PREPROCESSING_LABEL}"
            )


def load_detector_source(filename, sensor):
    """Load corrected and unsubtracted images on native detector geometry."""

    sensor = sensor.upper()
    if sensor not in IMAGE_FIELDS:
        raise ValueError("sensor must be WIC, SI12, or SI13")
    image_field = IMAGE_FIELDS[sensor]
    unsubtracted_image_field = UNSUBTRACTED_IMAGE_FIELDS[sensor]

    filename = Path(filename)
    with Dataset(filename) as nc:
        required = (
            image_field, unsubtracted_image_field,
            "dgweight", "img_variance", "frame_quality",
            "glat", "glon", "sza", "dza", "date", "t_start",
        )
        missing = [name for name in required if name not in nc.variables]
        if missing:
            raise ValueError(
                f"{filename} is missing source fields: {', '.join(missing)}"
            )
        _validate_background_model(nc, filename)

        counts = _as_float(nc.variables[image_field])
        unsubtracted_counts = _as_float(
            nc.variables[unsubtracted_image_field]
        )
        dgweight = _as_float(nc.variables["dgweight"])
        variance = _as_float(nc.variables["img_variance"])
        frame_quality = _as_float(nc.variables["frame_quality"])
        glat = _as_float(nc.variables["glat"])
        glon = _as_float(nc.variables["glon"])
        sza = _as_float(nc.variables["sza"])
        dza = _as_float(nc.variables["dza"])
        time = _source_times(nc)

    shape = counts.shape
    if len(shape) != 3 or time.shape != (shape[0],):
        raise ValueError("source counts must have shape (time, row, column)")
    for name, values in (
        ("unsubtracted_counts", unsubtracted_counts),
        ("dgweight", dgweight), ("img_variance", variance),
        ("glat", glat), ("glon", glon), ("sza", sza), ("dza", dza),
    ):
        if values.shape != shape:
            raise ValueError(f"source {name} does not match count dimensions")
    if frame_quality.shape != (shape[0],):
        raise ValueError("source frame_quality must contain one value per frame")
    if (
        np.any(~np.isfinite(frame_quality))
        or np.any(frame_quality != np.round(frame_quality))
        or np.any(~np.isin(frame_quality, [0, 1, 2]))
    ):
        raise ValueError("source frame_quality values must be 0, 1, or 2")
    if np.any(np.isfinite(variance) & (variance < 0)):
        raise ValueError("source img_variance must not contain negative values")

    source = source_identity(filename)
    source.update({
        "sensor": sensor,
        "image_field": image_field,
        "unsubtracted_image_field": unsubtracted_image_field,
        "time": time,
        "counts": counts,
        "unsubtracted_counts": unsubtracted_counts,
        "variance": variance,
        "frame_quality": frame_quality,
        "quality_weight": dgweight,
        "glat": glat,
        "glon": glon,
        "sza": sza,
        "dza": dza,
        "geometry_valid": np.isfinite(glat) & np.isfinite(glon),
    })
    return source


def add_modified_apex_coordinates(wic, reference_height_km=REFERENCE_HEIGHT_KM):
    """Calculate per-frame Modified Apex coordinates on the WIC detector."""

    shape = wic["counts"].shape
    mlat = np.full(shape, np.nan)
    mlon = np.full(shape, np.nan)
    mlt = np.full(shape, np.nan)
    ssalon = np.full(shape[0], np.nan)

    for frame, frame_time in enumerate(wic["time"]):
        valid = wic["geometry_valid"][frame]
        if not np.any(valid):
            continue

        apex = apexpy.Apex(frame_time, refh=reference_height_km)
        frame_mlat, frame_mlon = apex.convert(
            wic["glat"][frame][valid],
            wic["glon"][frame][valid],
            "geo",
            "apex",
            height=reference_height_km,
        )
        mlat[frame][valid] = frame_mlat
        mlon[frame][valid] = frame_mlon

        subsolar_latitude, subsolar_longitude = subsol(frame_time)
        _, frame_ssalon = apex.geo2apex(
            subsolar_latitude, subsolar_longitude, 318550
        )
        ssalon[frame] = frame_ssalon
        mlt[frame][valid] = (
            180 + np.asarray(frame_mlon, dtype=float) - frame_ssalon
        ) / 15 % 24

    wic = dict(wic)
    wic.update({"mlat": mlat, "mlon": mlon, "mlt": mlt, "ssalon": ssalon})
    return wic


#%% WIC-led time matching

def match_wic_times(wic_times, sensor_times, tolerance_seconds=TIME_TOLERANCE_SECONDS):
    """Match each WIC frame to the nearest unused sensor frame."""

    wic_times = np.asarray(wic_times, dtype=object)
    sensor_times = np.asarray(sensor_times, dtype=object)
    if wic_times.ndim != 1 or sensor_times.ndim != 1:
        raise ValueError("sensor times must be one-dimensional")
    if tolerance_seconds < 0:
        raise ValueError("time-match tolerance must not be negative")
    for name, values in (("WIC", wic_times), ("sensor", sensor_times)):
        converted = np.asarray(values, dtype="datetime64[ns]")
        if np.any(np.isnat(converted)):
            raise ValueError(f"{name} times must not contain missing values")

    tolerance = timedelta(seconds=float(tolerance_seconds))
    used = np.zeros(sensor_times.size, dtype=bool)
    matched = np.full(wic_times.size, -1, dtype=np.int32)

    for wic_index, wic_time in enumerate(wic_times):
        candidates = []
        for sensor_index, sensor_time in enumerate(sensor_times):
            separation = abs(sensor_time - wic_time)
            if not used[sensor_index] and separation <= tolerance:
                candidates.append((separation, sensor_time, sensor_index))

        if candidates:
            _, _, sensor_index = min(candidates)
            matched[wic_index] = sensor_index
            used[sensor_index] = True

    return matched


#%% Detector Product 1

def _empty_channel(shape, time_count):
    """Allocate one SI channel on the WIC detector and time axes."""

    return {
        "counts": np.full(shape, np.nan),
        "unsubtracted_counts": np.full(shape, np.nan),
        "variance": np.full(shape, np.nan),
        "quality_weight": np.full(shape, np.nan),
        "coverage": np.zeros(shape),
        "valid": np.zeros(shape, dtype=bool),
        "unsubtracted_valid": np.zeros(shape, dtype=bool),
        "source_count": np.zeros(shape, dtype=np.int32),
        "source_index": np.full(time_count, -1, dtype=np.int32),
        "source_time": np.full(time_count, None, dtype=object),
        "frame_quality": np.full(time_count, -1, dtype=np.int8),
        "diagnostics": {
            name: np.zeros(time_count, dtype=np.int32)
            for name in COREGISTRATION_INTEGER_FIELDS
        } | {
            name: np.full(time_count, np.nan)
            for name in COREGISTRATION_FLOAT_FIELDS
        },
    }


def _source_count(mapping, valid, shape):
    """Count finite source footprints contributing to each WIC pixel."""

    contributors = mapping.copy()
    contributors.data[:] = 1
    count = contributors @ np.asarray(valid, dtype=np.int32).ravel()
    return np.asarray(count).reshape(shape).astype(np.int32)


class FUVDetector:
    """Coregister subtracted and unsubtracted FUV images on WIC pixels."""

    def __init__(
        self,
        wic,
        si12=None,
        si13=None,
        *,
        preprocessing_label=PREPROCESSING_LABEL,
        time_tolerance_seconds=TIME_TOLERANCE_SECONDS,
        software_version="unrecorded experimental worktree",
    ):

        if wic.get("sensor") != "WIC":
            raise ValueError("wic input must contain WIC detector data")
        for sensor, expected in ((si12, "SI12"), (si13, "SI13")):
            if sensor is not None and sensor.get("sensor") != expected:
                raise ValueError(f"{expected.lower()} input must contain {expected} data")
        for name in ("mlat", "mlon", "mlt", "ssalon"):
            if name not in wic:
                raise ValueError(f"WIC detector data is missing {name}")
        if preprocessing_label != PREPROCESSING_LABEL:
            raise ValueError(
                f"this implementation supports only {PREPROCESSING_LABEL}; "
                "another label requires an explicit preprocessing branch"
            )

        self.product_type = "fuv_detector"
        self.representation = "detector"
        self.schema_version = SCHEMA_VERSION
        self.preprocessing_label = PREPROCESSING_LABEL
        self.time_tolerance_seconds = float(time_tolerance_seconds)
        self.time_match_rule = describe_time_match_rule(
            self.time_tolerance_seconds
        )
        self.software_version = str(software_version)
        self.time = np.asarray(wic["time"], dtype=object).copy()
        self.shape = np.asarray(wic["counts"]).shape
        if len(self.shape) != 3 or self.time.shape != (self.shape[0],):
            raise ValueError("WIC counts must have shape (time, row, column)")

        self.source_products = {"wic": wic.get("source_file", "")}
        self.source_sha256 = {"wic": wic.get("source_sha256", "")}
        self.source_size_bytes = {
            "wic": int(wic.get("source_size_bytes", -1))
        }
        self.source_mtime_ns = {"wic": int(wic.get("source_mtime_ns", -1))}
        self.image_fields = {"wic": wic.get("image_field", "unknown")}
        self.unsubtracted_image_fields = {
            "wic": wic.get("unsubtracted_image_field", "unknown")
        }
        self.wic_source_index = np.arange(self.shape[0], dtype=np.int32)
        self.wic_source_time = self.time.copy()
        self.wic_frame_quality = np.asarray(
            wic["frame_quality"], dtype=np.int8
        ).copy()
        if self.wic_frame_quality.shape != (self.shape[0],):
            raise ValueError("WIC frame_quality must contain one value per frame")

        for name in (
            "counts", "unsubtracted_counts", "variance", "quality_weight",
            "glat", "glon", "mlat", "mlon", "mlt", "sza", "dza",
        ):
            values = np.asarray(wic[name])
            if values.shape != self.shape:
                raise ValueError(f"WIC {name} does not match WIC count dimensions")
            output_name = "wic_counts" if name == "counts" else name
            if name == "unsubtracted_counts":
                output_name = "wic_unsubtracted_counts"
            if name == "variance":
                output_name = "wic_variance"
            if name == "quality_weight":
                output_name = "wic_quality_weight"
            setattr(self, output_name, values.copy())

        self.ssalon = np.asarray(wic["ssalon"], dtype=float).copy()
        if self.ssalon.shape != (self.shape[0],):
            raise ValueError("WIC ssalon must contain one value per frame")

        self.wic_valid = (
            np.isfinite(self.wic_counts)
            & np.asarray(wic["geometry_valid"], dtype=bool)
        )
        self.wic_unsubtracted_valid = (
            np.isfinite(self.wic_unsubtracted_counts)
            & np.asarray(wic["geometry_valid"], dtype=bool)
        )
        self.wic_variance = np.where(
            self.wic_valid
            & np.isfinite(self.wic_variance)
            & (self.wic_variance >= 0),
            self.wic_variance,
            np.nan,
        )
        self.wic_coverage = self.wic_valid.astype(float)

        channels = {
            "si12": _empty_channel(self.shape, self.shape[0]),
            "si13": _empty_channel(self.shape, self.shape[0]),
        }
        sources = {"si12": si12, "si13": si13}

        for name, sensor in sources.items():
            if sensor is None:
                self.source_products[name] = ""
                self.source_sha256[name] = ""
                self.source_size_bytes[name] = -1
                self.source_mtime_ns[name] = -1
                self.image_fields[name] = IMAGE_FIELDS[name.upper()]
                self.unsubtracted_image_fields[name] = (
                    UNSUBTRACTED_IMAGE_FIELDS[name.upper()]
                )
                continue
            self.source_products[name] = sensor.get("source_file", "")
            self.source_sha256[name] = sensor.get("source_sha256", "")
            self.source_size_bytes[name] = int(
                sensor.get("source_size_bytes", -1)
            )
            self.source_mtime_ns[name] = int(sensor.get("source_mtime_ns", -1))
            self.image_fields[name] = sensor.get("image_field", "unknown")
            self.unsubtracted_image_fields[name] = sensor.get(
                "unsubtracted_image_field", "unknown"
            )
            channels[name]["source_index"] = match_wic_times(
                self.time, sensor["time"], self.time_tolerance_seconds
            )
            matched = channels[name]["source_index"] >= 0
            channels[name]["source_time"][matched] = np.asarray(
                sensor["time"], dtype=object
            )[channels[name]["source_index"][matched]]
            channels[name]["frame_quality"][matched] = np.asarray(
                sensor["frame_quality"], dtype=np.int8
            )[channels[name]["source_index"][matched]]

        # Build each WIC-frame transform once and reuse it for both SI cameras.
        for frame in range(self.shape[0]):
            if np.count_nonzero(wic["geometry_valid"][frame]) < 3:
                continue
            transform = make_wic_transform(wic, frame)

            for name, sensor in sources.items():
                source_index = channels[name]["source_index"][frame]
                if sensor is None or source_index < 0:
                    continue

                mapping, _, diagnostics = make_si_mapping(
                    sensor, source_index, transform
                )
                valid_counts = (
                    np.isfinite(sensor["counts"][source_index])
                    & sensor["geometry_valid"][source_index]
                )
                counts, coverage = map_si(
                    sensor["counts"][source_index],
                    valid_counts,
                    mapping,
                    self.shape[1:],
                )
                valid_unsubtracted_counts = (
                    np.isfinite(sensor["unsubtracted_counts"][source_index])
                    & sensor["geometry_valid"][source_index]
                )
                unsubtracted_counts, _ = map_si(
                    sensor["unsubtracted_counts"][source_index],
                    valid_unsubtracted_counts,
                    mapping,
                    self.shape[1:],
                )
                valid_weight = (
                    valid_counts
                    & np.isfinite(sensor["quality_weight"][source_index])
                )
                quality_weight, _ = map_si(
                    sensor["quality_weight"][source_index],
                    valid_weight,
                    mapping,
                    self.shape[1:],
                )
                valid_variance = (
                    valid_counts
                    & np.isfinite(sensor["variance"][source_index])
                    & (sensor["variance"][source_index] >= 0)
                )
                variance, _ = map_si_variance(
                    sensor["variance"][source_index],
                    valid_variance,
                    mapping,
                    self.shape[1:],
                )

                channels[name]["counts"][frame] = counts
                channels[name]["unsubtracted_counts"][frame] = (
                    unsubtracted_counts
                )
                channels[name]["variance"][frame] = variance
                channels[name]["quality_weight"][frame] = quality_weight
                # Retain the raw overlap excess in the frame diagnostics, while
                # storing coverage itself as the requested [0, 1] fraction.
                channels[name]["coverage"][frame] = np.minimum(coverage, 1)
                channels[name]["valid"][frame] = np.isfinite(counts)
                channels[name]["unsubtracted_valid"][frame] = np.isfinite(
                    unsubtracted_counts
                )
                channels[name]["source_count"][frame] = _source_count(
                    mapping, valid_counts, self.shape[1:]
                )

                output_diagnostics = channels[name]["diagnostics"]
                output_diagnostics["valid_source_centres"][frame] = (
                    diagnostics["valid_source_centres"]
                )
                output_diagnostics["valid_source_footprints"][frame] = (
                    diagnostics["valid_source_footprints"]
                )
                output_diagnostics["accepted_target_pixels"][frame] = (
                    diagnostics["target_pixels_coverage_ge_0.9"]
                )
                output_diagnostics["roundtrip_median_km"][frame] = (
                    diagnostics["internal_roundtrip_median_km"]
                )
                output_diagnostics["roundtrip_95th_km"][frame] = (
                    diagnostics["internal_roundtrip_95th_km"]
                )
                for diagnostic_name in (
                    "coverage_median", "coverage_95th", "coverage_maximum"
                ):
                    output_diagnostics[diagnostic_name][frame] = diagnostics[
                        diagnostic_name
                    ]

        for name, channel in channels.items():
            for field, values in channel.items():
                if field == "diagnostics":
                    setattr(self, f"{name}_coregistration", values)
                else:
                    setattr(self, f"{name}_{field}", values)

    @classmethod
    def from_files(
        cls,
        wic_file,
        si12_file=None,
        si13_file=None,
        **configuration,
    ):
        """Load native sensor files and construct an experimental Product 1."""

        wic = load_detector_source(wic_file, "WIC")
        wic = add_modified_apex_coordinates(wic)
        si12 = (
            load_detector_source(si12_file, "SI12")
            if si12_file is not None else None
        )
        si13 = (
            load_detector_source(si13_file, "SI13")
            if si13_file is not None else None
        )
        return cls(wic, si12, si13, **configuration)

    #%% NetCDF output

    @staticmethod
    def _write_time(nc, name, values, units):
        variable = nc.createVariable(name, "f8", ("time",), fill_value=np.nan)
        encoded = np.full(len(values), np.nan)
        matched = np.asarray([value is not None for value in values])
        if np.any(matched):
            encoded[matched] = date2num(
                np.asarray(values, dtype=object)[matched].tolist(),
                units,
                calendar="standard",
            )
        variable[:] = encoded
        variable.units = units
        variable.calendar = "standard"
        variable.time_zone = "UTC"

    def to_nc(self, filename):
        """Write the self-describing detector-space Product 1 file."""

        with Dataset(filename, "w", format="NETCDF4") as nc:
            nc.createDimension("time", self.shape[0])
            nc.createDimension("row", self.shape[1])
            nc.createDimension("column", self.shape[2])

            nc.product_type = self.product_type
            nc.representation = self.representation
            nc.schema_version = self.schema_version
            nc.preprocessing_label = self.preprocessing_label
            nc.software_version = self.software_version
            nc.time_match_tolerance_seconds = self.time_tolerance_seconds
            nc.time_match_rule = self.time_match_rule
            nc.source_time_decoding = SOURCE_TIME_DECODING
            nc.quality_weight_method = QUALITY_WEIGHT_METHOD
            nc.detector_noise_model = DETECTOR_NOISE_MODEL
            nc.coordinate_system = COORDINATE_SYSTEM
            nc.reference_height_km = REFERENCE_HEIGHT_KM
            nc.coregistration_method = (
                "SI footprint overlap area averaged on WIC detector pixels"
            )
            nc.coregistration_footprint_model = "uniform top-hat"
            nc.coregistration_max_roundtrip_error_km = MAX_COREG_ERROR_KM
            nc.coregistration_minimum_coverage = MIN_SI_COVERAGE
            nc.coregistration_overlap_operator_stored = np.int8(0)
            nc.unsubtracted_counts_stored = np.int8(1)

            for sensor in ("wic", "si12", "si13"):
                nc.setncattr(f"source_{sensor}", self.source_products[sensor])
                nc.setncattr(
                    f"source_{sensor}_sha256", self.source_sha256[sensor]
                )
                nc.setncattr(
                    f"source_{sensor}_size_bytes",
                    np.int64(self.source_size_bytes[sensor]),
                )
                nc.setncattr(
                    f"source_{sensor}_mtime_ns",
                    np.int64(self.source_mtime_ns[sensor]),
                )
                nc.setncattr(f"{sensor}_image_field", self.image_fields[sensor])
                nc.setncattr(
                    f"{sensor}_unsubtracted_image_field",
                    self.unsubtracted_image_fields[sensor],
                )

            time_units = "seconds since 2000-01-01 00:00:00"
            self._write_time(nc, "time", self.time, time_units)
            self._write_time(nc, "wic_source_time", self.wic_source_time, time_units)
            self._write_time(nc, "si12_source_time", self.si12_source_time, time_units)
            self._write_time(nc, "si13_source_time", self.si13_source_time, time_units)

            for sensor in ("wic", "si12", "si13"):
                variable = nc.createVariable(
                    f"{sensor}_source_index", "i4", ("time",)
                )
                variable[:] = getattr(self, f"{sensor}_source_index")

                frame_quality = nc.createVariable(
                    f"{sensor}_frame_quality", "i1", ("time",), fill_value=-1
                )
                frame_quality[:] = getattr(self, f"{sensor}_frame_quality")
                frame_quality.long_name = "fuvpy frame-level quality flag"
                frame_quality.flag_values = np.asarray([0, 1, 2], dtype=np.int8)
                frame_quality.flag_meanings = "rejected usable science_ready"

            row = nc.createVariable("detector_row", "i4", ("row",))
            column = nc.createVariable("detector_column", "i4", ("column",))
            row[:] = np.arange(self.shape[1])
            column[:] = np.arange(self.shape[2])

            field_units = {
                "wic_counts": "counts",
                "si12_counts": "counts",
                "si13_counts": "counts",
                "wic_unsubtracted_counts": "counts",
                "si12_unsubtracted_counts": "counts",
                "si13_unsubtracted_counts": "counts",
                "wic_variance": "counts^2",
                "si12_variance": "counts^2",
                "si13_variance": "counts^2",
                "wic_quality_weight": "1",
                "si12_quality_weight": "1",
                "si13_quality_weight": "1",
                "wic_coverage": "1",
                "si12_coverage": "1",
                "si13_coverage": "1",
                "glat": "degrees_north",
                "glon": "degrees_east",
                "mlat": "degrees",
                "mlon": "degrees",
                "mlt": "hours",
                "sza": "degrees",
                "dza": "degrees",
            }
            dimensions = ("time", "row", "column")
            coordinate_fields = {"glat", "glon", "mlat", "mlon", "mlt"}
            for name, units in field_units.items():
                dtype = "f8" if name in coordinate_fields else "f4"
                variable = nc.createVariable(name, dtype, dimensions, zlib=True)
                variable[:] = getattr(self, name)
                variable.units = units

            for sensor in ("wic", "si12", "si13"):
                for name in ("valid", "unsubtracted_valid"):
                    valid = nc.createVariable(
                        f"{sensor}_{name}", "i1", dimensions, zlib=True
                    )
                    valid[:] = getattr(
                        self, f"{sensor}_{name}"
                    ).astype(np.int8)

            for sensor in ("si12", "si13"):
                source_count = nc.createVariable(
                    f"{sensor}_source_count", "i4", dimensions, zlib=True
                )
                source_count[:] = getattr(self, f"{sensor}_source_count")

                diagnostics = getattr(self, f"{sensor}_coregistration")
                for name in COREGISTRATION_INTEGER_FIELDS:
                    variable = nc.createVariable(
                        f"{sensor}_coreg_{name}", "i4", ("time",)
                    )
                    variable[:] = diagnostics[name]
                for name in COREGISTRATION_FLOAT_FIELDS:
                    variable = nc.createVariable(
                        f"{sensor}_coreg_{name}", "f4", ("time",)
                    )
                    variable[:] = diagnostics[name]

            ssalon = nc.createVariable("ssalon", "f4", ("time",))
            ssalon[:] = self.ssalon
            ssalon.units = "degrees"
