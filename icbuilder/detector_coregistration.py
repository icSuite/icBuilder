"""Area-weighted coregistration between native IMAGE detector grids.

The SI pixel centres are first expressed in continuous WIC detector
coordinates. Neighbouring centres then define an approximate quadrilateral
footprint for each SI pixel. Values are averaged into WIC pixels using the
quadrilateral overlap areas.

This is still an approximation because the camera point-spread functions are
unknown. It preserves detector-pixel area much better than treating each SI
pixel centre as an independent point measurement.
"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator
from scipy.sparse import csr_matrix
from scipy.spatial import Delaunay

from .footprints import _overlap_entries_numba, neighbour_step, overlap_mean


IMAGE_RADIUS_KM = 6371.0 + 130.0
MAX_COREG_ERROR_KM = 25.0
MIN_SI_COVERAGE = 0.9


def unit_vectors(latitude, longitude):
    """Convert geographic latitude and longitude to Cartesian unit vectors."""
    latitude = np.deg2rad(np.asarray(latitude))
    longitude = np.deg2rad(np.asarray(longitude))
    return np.stack([
        np.cos(latitude) * np.cos(longitude),
        np.cos(latitude) * np.sin(longitude),
        np.sin(latitude),
    ], axis=-1)


def plane_coordinates(vectors, transform):
    """Project unit vectors onto the local plane used for WIC interpolation."""
    denominator = np.sum(vectors * transform["centre"], axis=-1)
    x = np.divide(
        np.sum(vectors * transform["east"], axis=-1), denominator,
        out=np.full(denominator.shape, np.nan), where=denominator > 0,
    )
    y = np.divide(
        np.sum(vectors * transform["north"], axis=-1), denominator,
        out=np.full(denominator.shape, np.nan), where=denominator > 0,
    )
    return x, y


def make_wic_transform(wic, frame):
    """Map geographic coordinates to continuous WIC detector coordinates."""
    valid = wic["geometry_valid"][frame]
    if np.count_nonzero(valid) < 3:
        raise ValueError("WIC frame needs at least three valid geometry pixels")
    vectors = unit_vectors(wic["glat"][frame], wic["glon"][frame])

    centre = np.mean(vectors[valid], axis=0)
    centre /= np.linalg.norm(centre)
    east = np.cross([0.0, 0.0, 1.0], centre)
    if np.linalg.norm(east) < 1e-8:
        east = np.cross([0.0, 1.0, 0.0], centre)
    east /= np.linalg.norm(east)
    north = np.cross(centre, east)

    transform = {"centre": centre, "east": east, "north": north}
    x, y = plane_coordinates(vectors, transform)
    projected_valid = valid & np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(projected_valid) < 3:
        raise ValueError(
            "WIC frame needs at least three pixels on the local projection"
        )

    row, column = np.indices(valid.shape)
    triangulation = Delaunay(
        np.column_stack([x[projected_valid], y[projected_valid]])
    )
    detector_axes = (np.arange(valid.shape[0]), np.arange(valid.shape[1]))

    transform.update({
        "shape": valid.shape,
        "row": LinearNDInterpolator(
            triangulation, row[projected_valid], fill_value=np.nan
        ),
        "column": LinearNDInterpolator(
            triangulation, column[projected_valid], fill_value=np.nan
        ),
        "vector": [
            RegularGridInterpolator(
                detector_axes, vectors[..., component],
                bounds_error=False, fill_value=np.nan,
            )
            for component in range(3)
        ],
    })
    return transform


def geographic_to_wic(transform, latitude, longitude):
    """Return continuous WIC row/column and the geographic round-trip error."""
    vectors = unit_vectors(latitude, longitude)
    x, y = plane_coordinates(vectors, transform)
    projected_valid = np.isfinite(x.ravel()) & np.isfinite(y.ravel())

    row = np.full(x.size, np.nan)
    column = np.full(x.size, np.nan)
    error = np.full(x.size, np.nan)
    if not np.any(projected_valid):
        return row.reshape(x.shape), column.reshape(x.shape), error.reshape(x.shape)

    projected_points = np.column_stack([
        x.ravel()[projected_valid], y.ravel()[projected_valid]
    ])
    row[projected_valid] = transform["row"](projected_points)
    column[projected_valid] = transform["column"](projected_points)

    detector_valid = projected_valid & np.isfinite(row) & np.isfinite(column)
    if not np.any(detector_valid):
        return row.reshape(x.shape), column.reshape(x.shape), error.reshape(x.shape)

    detector_points = np.column_stack([
        row[detector_valid], column[detector_valid]
    ])
    reconstructed = np.column_stack([
        interpolator(detector_points) for interpolator in transform["vector"]
    ])
    length = np.linalg.norm(reconstructed, axis=1)
    reconstructed = np.divide(
        reconstructed, length[:, None], out=np.full_like(reconstructed, np.nan),
        where=length[:, None] > 0,
    )
    cosine = np.sum(
        reconstructed * vectors.reshape(-1, 3)[detector_valid], axis=1
    )
    error[detector_valid] = (
        IMAGE_RADIUS_KM * np.arccos(np.clip(cosine, -1, 1))
    )
    return row.reshape(x.shape), column.reshape(x.shape), error.reshape(x.shape)


def make_si_mapping(si, frame, transform):
    """Build an SI-footprint to WIC-pixel overlap matrix for one frame."""
    row, column, error = geographic_to_wic(
        transform, si["glat"][frame], si["glon"][frame]
    )
    centre_valid = (
        si["geometry_valid"][frame] & np.isfinite(row) & np.isfinite(column)
        & np.isfinite(error) & (error <= MAX_COREG_ERROR_KM)
    )
    row = np.where(centre_valid, row, np.nan)
    column = np.where(centre_valid, column, np.nan)

    row_step = np.stack([
        neighbour_step(column, 0), neighbour_step(row, 0)
    ], axis=-1)
    column_step = np.stack([
        neighbour_step(column, 1), neighbour_step(row, 1)
    ], axis=-1)
    centre = np.stack([column, row], axis=-1)
    corners = np.stack([
        centre - row_step / 2 - column_step / 2,
        centre - row_step / 2 + column_step / 2,
        centre + row_step / 2 + column_step / 2,
        centre + row_step / 2 - column_step / 2,
    ], axis=-2)

    footprint_valid = np.all(np.isfinite(corners), axis=(-2, -1))
    ny, nx = transform["shape"]
    target, source, area = _overlap_entries_numba(
        np.ascontiguousarray(corners.reshape(-1, 4, 2)),
        np.flatnonzero(footprint_valid).astype(np.int64),
        np.arange(-0.5, nx + 0.5), np.arange(-0.5, ny + 0.5), ny, nx,
    )
    mapping = csr_matrix((area, (target, source)), shape=(ny * nx, row.size))
    coverage = np.asarray(
        mapping @ centre_valid.ravel().astype(float)
    ).reshape(ny, nx)
    accepted_error = error[centre_valid]

    positive_coverage = coverage[coverage > 0]
    diagnostics = {
        "sensor": si["sensor"],
        "time": np.datetime_as_string(
            np.datetime64(si["time"][frame]), unit="ms"
        ),
        "valid_source_centres": int(centre_valid.sum()),
        "valid_source_footprints": int(footprint_valid.sum()),
        "target_pixels_coverage_ge_0.9": int((coverage >= MIN_SI_COVERAGE).sum()),
        "internal_roundtrip_median_km": (
            float(np.median(accepted_error)) if accepted_error.size else np.nan
        ),
        "internal_roundtrip_95th_km": (
            float(np.percentile(accepted_error, 95))
            if accepted_error.size else np.nan
        ),
        "coverage_median": (
            float(np.median(positive_coverage))
            if positive_coverage.size else np.nan
        ),
        "coverage_95th": (
            float(np.percentile(positive_coverage, 95))
            if positive_coverage.size else np.nan
        ),
        "coverage_maximum": (
            float(np.max(positive_coverage))
            if positive_coverage.size else np.nan
        ),
    }
    geometry = {"row": row, "column": column, "coverage": coverage}
    return mapping, geometry, diagnostics


def map_si(values, valid, mapping, shape, minimum_coverage=MIN_SI_COVERAGE):
    """Area-average one SI field into WIC pixels and require adequate coverage."""
    mapped, coverage = overlap_mean(np.where(valid, values, np.nan), mapping, shape)
    mapped[coverage < minimum_coverage] = np.nan
    return mapped, coverage


def map_si_variance(
    variance, valid, mapping, shape, minimum_coverage=MIN_SI_COVERAGE
):
    """Propagate independent SI pixel variances through the area mean."""

    variance = np.asarray(variance, dtype=float)
    valid = np.asarray(valid, dtype=bool) & np.isfinite(variance) & (variance >= 0)
    weights = np.asarray(mapping @ valid.ravel().astype(float)).ravel()
    numerator = np.asarray(
        mapping.power(2) @ np.where(valid, variance, 0).ravel()
    ).ravel()
    mapped = np.divide(
        numerator,
        weights ** 2,
        out=np.full(mapping.shape[0], np.nan),
        where=weights > 0,
    ).reshape(shape)
    coverage = weights.reshape(shape)
    mapped[coverage < minimum_coverage] = np.nan
    return mapped, coverage
