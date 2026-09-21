"""Shared WIC-detector to fixed Cubed-Sphere reduction operations."""

#%% Imports

import numpy as np
from icreader import open_product
from netCDF4 import date2num

from .footprints import overlap_mapping


#%% Reduction contract

BINNING_METHOD = "wic_footprint_overlap_area_mean_v1"
UNCERTAINTY_METHOD = (
    "independent detector-pixel variance propagated through the normalized "
    "overlap mean with squared overlap weights; covariance induced by SI-to-"
    "WIC coregistration is unavailable"
)
TIME_FIELDS = (
    "time", "wic_source_time", "si12_source_time", "si13_source_time",
    "Kp_interval_start",
)
INDEX_FIELDS = (
    "wic_source_index", "si12_source_index", "si13_source_index",
)
FRAME_QUALITY_FIELDS = (
    "wic_frame_quality", "si12_frame_quality", "si13_frame_quality",
)
TIME_VALUE_FIELDS = ("Kp", "ssalon")


def make_detector_cs_mapping(mlat, mlt, grid):
    """Build one WIC-footprint mapping for a detector frame."""

    return overlap_mapping(mlat, mlt, grid)


def _flat_valid(values, valid, source_size, *, nonnegative=False):
    """Return flattened values and their accepted source support."""

    values = np.asarray(values, dtype=float).ravel()
    valid = np.asarray(valid, dtype=bool).ravel()
    if values.size != source_size or valid.size != source_size:
        raise ValueError("values and valid must match the mapping source size")

    accepted = valid & np.isfinite(values)
    if nonnegative:
        accepted &= values >= 0
    return values, accepted


def _overlap_for_valid(mapping, valid):
    return np.asarray(mapping @ valid.astype(float)).ravel()


def reduce_area_mean(values, valid, mapping, output_shape):
    """Return the overlap-area-weighted mean and accepted overlap area."""

    values, accepted = _flat_valid(values, valid, mapping.shape[1])
    overlap = _overlap_for_valid(mapping, accepted)
    numerator = np.asarray(mapping @ np.where(accepted, values, 0.0)).ravel()

    mean = np.full(mapping.shape[0], np.nan)
    covered = overlap > 0
    mean[covered] = numerator[covered] / overlap[covered]
    return mean.reshape(output_shape), overlap.reshape(output_shape)


def reduce_measurement_variance(variance, valid, mapping, output_shape):
    """Propagate independent variance through the normalized overlap mean."""

    variance, accepted = _flat_valid(
        variance, valid, mapping.shape[1], nonnegative=True
    )
    overlap = _overlap_for_valid(mapping, accepted)
    squared_mapping = mapping.copy()
    squared_mapping.data **= 2
    numerator = np.asarray(
        squared_mapping @ np.where(accepted, variance, 0.0)
    ).ravel()

    reduced = np.full(mapping.shape[0], np.nan)
    covered = overlap > 0
    reduced[covered] = numerator[covered] / overlap[covered] ** 2
    return reduced.reshape(output_shape), overlap.reshape(output_shape)


def reduce_covariance(covariance, valid, mapping, output_shape):
    """Propagate same-pixel covariance through the normalized overlap mean."""

    covariance, accepted = _flat_valid(
        covariance, valid, mapping.shape[1]
    )
    overlap = _overlap_for_valid(mapping, accepted)
    squared_mapping = mapping.copy()
    squared_mapping.data **= 2
    numerator = np.asarray(
        squared_mapping @ np.where(accepted, covariance, 0.0)
    ).ravel()

    reduced = np.full(mapping.shape[0], np.nan)
    covered = overlap > 0
    reduced[covered] = numerator[covered] / overlap[covered] ** 2
    return reduced.reshape(output_shape), overlap.reshape(output_shape)


def reduce_flag_fraction(flag, valid, mapping, output_shape):
    """Return covered-area flag fraction and whether any contributor is set."""

    flag = np.asarray(flag, dtype=bool).ravel()
    valid = np.asarray(valid, dtype=bool).ravel()
    if flag.size != mapping.shape[1] or valid.size != mapping.shape[1]:
        raise ValueError("flag and valid must match the mapping source size")

    overlap = _overlap_for_valid(mapping, valid)
    flagged_overlap = np.asarray(
        mapping @ (flag & valid).astype(float)
    ).ravel()
    fraction = np.full(mapping.shape[0], np.nan)
    covered = overlap > 0
    fraction[covered] = flagged_overlap[covered] / overlap[covered]
    return (
        fraction.reshape(output_shape),
        (flagged_overlap > 0).reshape(output_shape),
    )


def reduce_support(valid, mapping, cell_area, output_shape):
    """Return contributing footprint count and valid covered-area fraction."""

    valid = np.asarray(valid, dtype=bool).ravel()
    if valid.size != mapping.shape[1]:
        raise ValueError("valid must match the mapping source size")

    overlap = _overlap_for_valid(mapping, valid).reshape(output_shape)
    contributors = mapping.copy()
    contributors.data[:] = 1
    count = np.asarray(
        contributors @ valid.astype(np.int32)
    ).reshape(output_shape)
    coverage = np.minimum(overlap / np.asarray(cell_area), 1.0)
    return count.astype(np.int32), coverage


#%% Shared product I/O

def encode_time(values, units, calendar):
    """Encode decoded reader times while retaining missing source frames."""

    values = np.asarray(values, dtype=object)
    encoded = np.full(values.shape, np.nan)
    valid = np.asarray([value is not None for value in values], dtype=bool)
    if np.any(valid):
        encoded[valid] = date2num(
            values[valid].tolist(), units, calendar=calendar
        )
    return encoded


def read_common_time_fields(source):
    """Copy the one-dimensional frame identity shared by both CS products."""

    result = {}
    for name in TIME_FIELDS:
        encoding = source.time_encoding[name]
        result[name] = encode_time(
            getattr(source, name), encoding["units"], encoding["calendar"]
        )
        result[f"{name}_units"] = encoding["units"]
        result[f"{name}_calendar"] = encoding["calendar"]
    for name in INDEX_FIELDS + FRAME_QUALITY_FIELDS:
        result[name] = np.asarray(getattr(source, name), dtype=int).copy()
    for name in TIME_VALUE_FIELDS:
        result[name] = np.asarray(getattr(source, name), dtype=float).copy()
    return result


def write_common_cs_coordinates(nc, product):
    """Write dimensions, frame identity, and the fixed CS coordinates."""

    time_size, ny, nx = product.shape
    nc.createDimension("time", time_size)
    nc.createDimension("dim1", ny)
    nc.createDimension("dim2", nx)
    nc.createDimension("dim1_edge", ny + 1)
    nc.createDimension("dim2_edge", nx + 1)

    for name in TIME_FIELDS:
        variable = nc.createVariable(name, "f8", ("time",), fill_value=np.nan)
        variable[:] = product.time_fields[name]
        variable.units = product.time_fields[f"{name}_units"]
        variable.calendar = product.time_fields[f"{name}_calendar"]
        variable.time_zone = "UTC"

    for name in INDEX_FIELDS:
        variable = nc.createVariable(name, "i4", ("time",))
        variable[:] = product.time_fields[name]

    for name in FRAME_QUALITY_FIELDS:
        variable = nc.createVariable(name, "i1", ("time",), fill_value=-1)
        variable[:] = product.time_fields[name]
        variable.flag_values = np.asarray([0, 1, 2], dtype=np.int8)
        variable.flag_meanings = "rejected usable science_ready"

    for name, units in (("Kp", "1"), ("ssalon", "degrees")):
        variable = nc.createVariable(name, "f4", ("time",))
        variable[:] = product.time_fields[name]
        variable.units = units

    grid_group = nc.createGroup("grid")
    grid_group.grid_id = product.grid_id
    grid_group.position = np.asarray(product.grid.projection.position, dtype=float)
    grid_group.orientation = np.asarray(
        product.grid.projection.orientation, dtype=float
    )
    grid_group.reference_height_km = product.reference_height_km
    grid_group.radius_metres = product.grid.R

    grid_fields = {
        "xi": (np.asarray(product.grid.xi), "radians", ("dim1", "dim2")),
        "eta": (np.asarray(product.grid.eta), "radians", ("dim1", "dim2")),
        "mlat": (np.asarray(product.grid.lat), "degrees", ("dim1", "dim2")),
        "mlt": (
            np.mod(np.asarray(product.grid.lon) / 15.0, 24.0),
            "hours",
            ("dim1", "dim2"),
        ),
        "xi_edge": (
            np.asarray(product.grid.xi_mesh[0]),
            "radians",
            ("dim2_edge",),
        ),
        "eta_edge": (
            np.asarray(product.grid.eta_mesh[:, 0]),
            "radians",
            ("dim1_edge",),
        ),
    }
    for name, (values, units, dimensions) in grid_fields.items():
        variable = grid_group.createVariable(name, "f8", dimensions, zlib=True)
        variable[:] = values
        variable.units = units


def validate_detector_pair(precipitation, conductance, precipitation_identity):
    """Require Product 2 and Product 3 to describe the same detector orbit."""

    if precipitation.shape != conductance.shape:
        raise ValueError("detector Product 2 and Product 3 shapes do not match")
    if (
        conductance.source_precipitation_detector_sha256
        != precipitation_identity["source_sha256"]
    ):
        raise ValueError("Product 3 does not identify the selected Product 2")

    for name in TIME_FIELDS:
        first = np.asarray(getattr(precipitation, name), dtype=object)
        second = np.asarray(getattr(conductance, name), dtype=object)
        if not np.array_equal(first, second):
            raise ValueError(f"detector Product 2 and Product 3 {name} differ")
        if dict(precipitation.time_encoding[name]) != dict(
            conductance.time_encoding[name]
        ):
            raise ValueError(
                f"detector Product 2 and Product 3 {name} encodings differ"
            )


def build_detector_cs_products(
    precipitation_filename,
    conductance_filename,
    *,
    software_version="unrecorded experimental worktree",
):
    """Reduce paired detector Products 2 and 3 with one mapping per frame."""

    from .conductancecs import ConductanceCS
    from .grids import make_detector_cs_grid
    from .precipitationcs import PrecipitationCS

    grid = make_detector_cs_grid()
    with open_product(precipitation_filename) as precipitation, open_product(
        conductance_filename
    ) as conductance:
        precipitation_product = PrecipitationCS(
            precipitation,
            grid=grid,
            software_version=software_version,
        )
        conductance_product = ConductanceCS(
            conductance,
            grid=grid,
            software_version=software_version,
        )
        validate_detector_pair(
            precipitation,
            conductance,
            precipitation_product.source_identity,
        )
        precipitation_mlat = np.asarray(
            precipitation.read("mlat"), dtype=float
        )
        conductance_mlat = np.asarray(
            conductance.read("mlat"), dtype=float
        )
        if not np.array_equal(
            precipitation_mlat, conductance_mlat, equal_nan=True
        ):
            raise ValueError("detector Product 2 and Product 3 MLAT differ")
        del conductance_mlat

        precipitation_mlt = np.asarray(
            precipitation.read("mlt"), dtype=float
        )
        conductance_mlt = np.asarray(
            conductance.read("mlt"), dtype=float
        )
        if not np.array_equal(
            precipitation_mlt, conductance_mlt, equal_nan=True
        ):
            raise ValueError("detector Product 2 and Product 3 MLT differ")
        del conductance_mlt

        mappings = []
        cell_area = None
        for frame in range(precipitation_product.shape[0]):
            mapping, frame_cell_area = make_detector_cs_mapping(
                precipitation_mlat[frame],
                precipitation_mlt[frame],
                grid,
            )
            mappings.append(mapping)
            if cell_area is None:
                cell_area = frame_cell_area
        del precipitation_mlat, precipitation_mlt

        precipitation_product.reduce(precipitation, mappings, cell_area)
        conductance_product.reduce(conductance, mappings, cell_area)

    return precipitation_product, conductance_product
