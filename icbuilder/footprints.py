"""Map detector-pixel footprints onto a Cubed-Sphere grid.

Each detector pixel is approximated as a uniform quadrilateral. Its corners
are inferred from the local row and column spacing of the geolocated pixel
centres. The quadrilateral is then split among every CS cell that it overlaps.
"""

#%% Imports

import numpy as np
from numba import njit
from scipy.sparse import csr_matrix


#%% Infer detector-pixel corners

def neighbour_step(values, axis):
    """Estimate one detector-pixel step from adjacent pixel centres."""

    values = np.asarray(values, dtype=float)
    previous = np.full_like(values, np.nan)
    following = np.full_like(values, np.nan)

    if axis == 0:
        previous[1:] = values[:-1]
        following[:-1] = values[1:]
    else:
        previous[:, 1:] = values[:, :-1]
        following[:, :-1] = values[:, 1:]

    step = np.full_like(values, np.nan)
    centre_ok = np.isfinite(values)
    previous_ok = np.isfinite(previous)
    following_ok = np.isfinite(following)

    both = centre_ok & previous_ok & following_ok
    step[both] = (following[both] - previous[both]) / 2

    forward = centre_ok & ~previous_ok & following_ok
    step[forward] = following[forward] - values[forward]

    backward = centre_ok & previous_ok & ~following_ok
    step[backward] = values[backward] - previous[backward]

    return step


def infer_footprints(mlat, mlt, grid):
    """Return four projected corners for each geolocated detector pixel."""

    longitude = np.mod(np.asarray(mlt) * 15, 360)
    xi, eta = grid.projection.geo2cube(
        longitude,
        mlat,
        set_points_off_cube_to_nan=True,
    )

    row = np.stack(
        [neighbour_step(xi, 0), neighbour_step(eta, 0)], axis=-1
    )
    column = np.stack(
        [neighbour_step(xi, 1), neighbour_step(eta, 1)], axis=-1
    )
    centre = np.stack([xi, eta], axis=-1)

    corners = np.stack(
        [
            centre - row / 2 - column / 2,
            centre - row / 2 + column / 2,
            centre + row / 2 + column / 2,
            centre + row / 2 - column / 2,
        ],
        axis=-2,
    )

    valid = np.all(np.isfinite(corners), axis=(-2, -1))
    return corners, valid


def point_in_quadrilaterals(x, y, quadrilaterals):
    """Return which projected quadrilaterals contain one projected point.

    Boundary points count as contained. The quadrilaterals may be clockwise
    or counter-clockwise, but are expected to be convex as produced by
    :func:`infer_footprints`.
    """

    quadrilaterals = np.asarray(quadrilaterals, dtype=float)
    if quadrilaterals.size == 0:
        return np.zeros(0, dtype=bool)
    if quadrilaterals.ndim != 3 or quadrilaterals.shape[1:] != (4, 2):
        raise ValueError("quadrilaterals must have shape (n, 4, 2)")

    following = np.roll(quadrilaterals, -1, axis=1)
    edges = following - quadrilaterals
    relative = np.array([x, y]) - quadrilaterals
    cross = edges[..., 0] * relative[..., 1] - edges[..., 1] * relative[..., 0]

    edge_scale = np.max(np.sum(edges**2, axis=-1), axis=1)
    tolerance = 1e-10 * np.maximum(edge_scale, 1)
    clockwise = np.all(cross <= tolerance[:, None], axis=1)
    counter_clockwise = np.all(cross >= -tolerance[:, None], axis=1)
    finite = np.all(np.isfinite(quadrilaterals), axis=(1, 2))
    return finite & (clockwise | counter_clockwise)


#%% Polygon clipping

def polygon_area(polygon):
    """Return unsigned polygon area in projected xi/eta coordinates."""

    if len(polygon) < 3:
        return 0.0

    polygon = np.asarray(polygon)
    x = polygon[:, 0]
    y = polygon[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def clip_at_boundary(polygon, coordinate, boundary, keep_above):
    """Clip a polygon against one vertical or horizontal boundary."""

    if len(polygon) == 0:
        return []

    clipped = []
    previous = np.asarray(polygon[-1], dtype=float)

    def inside(point):
        if keep_above:
            return point[coordinate] >= boundary
        return point[coordinate] <= boundary

    previous_inside = inside(previous)

    for current in polygon:
        current = np.asarray(current, dtype=float)
        current_inside = inside(current)

        if current_inside != previous_inside:
            difference = current[coordinate] - previous[coordinate]
            if difference != 0:
                fraction = (boundary - previous[coordinate]) / difference
                clipped.append(previous + fraction * (current - previous))

        if current_inside:
            clipped.append(current)

        previous = current
        previous_inside = current_inside

    return clipped


def rectangle_overlap(polygon, xmin, xmax, ymin, ymax):
    """Return polygon overlap area with one rectangular CS cell."""

    clipped = list(np.asarray(polygon, dtype=float))
    clipped = clip_at_boundary(clipped, 0, xmin, keep_above=True)
    clipped = clip_at_boundary(clipped, 0, xmax, keep_above=False)
    clipped = clip_at_boundary(clipped, 1, ymin, keep_above=True)
    clipped = clip_at_boundary(clipped, 1, ymax, keep_above=False)
    return polygon_area(clipped)


#%% Compiled polygon clipping

@njit(cache=True)
def _polygon_area_numba(x, y, size):
    """Shoelace area used inside the compiled overlap loop."""

    if size < 3:
        return 0.0

    area = 0.0
    previous_x = x[size - 1]
    previous_y = y[size - 1]

    for i in range(size):
        current_x = x[i]
        current_y = y[i]
        area += previous_x * current_y - previous_y * current_x
        previous_x = current_x
        previous_y = current_y

    return abs(area) / 2


@njit(cache=True)
def _clip_at_boundary_numba(x, y, size, coordinate, boundary, keep_above,
                            output_x, output_y):
    """Compiled form of ``clip_at_boundary`` for one small polygon."""

    if size == 0:
        return 0

    output_size = 0
    previous_x = x[size - 1]
    previous_y = y[size - 1]
    previous_value = previous_x if coordinate == 0 else previous_y
    if keep_above:
        previous_inside = previous_value >= boundary
    else:
        previous_inside = previous_value <= boundary

    for i in range(size):
        current_x = x[i]
        current_y = y[i]
        current_value = current_x if coordinate == 0 else current_y
        if keep_above:
            current_inside = current_value >= boundary
        else:
            current_inside = current_value <= boundary

        if current_inside != previous_inside:
            difference = current_value - previous_value
            if difference != 0:
                fraction = (boundary - previous_value) / difference
                output_x[output_size] = previous_x + fraction * (current_x - previous_x)
                output_y[output_size] = previous_y + fraction * (current_y - previous_y)
                output_size += 1

        if current_inside:
            output_x[output_size] = current_x
            output_y[output_size] = current_y
            output_size += 1

        previous_x = current_x
        previous_y = current_y
        previous_value = current_value
        previous_inside = current_inside

    return output_size


@njit(cache=True)
def _rectangle_overlap_numba(polygon, xmin, xmax, ymin, ymax):
    """Clip one quadrilateral against one target cell."""

    # A quadrilateral clipped by a rectangle has at most eight vertices.
    # Twelve slots leave a little room while keeping these temporary arrays
    # small inside the heavily repeated calculation.
    x_a = np.empty(12, dtype=np.float64)
    y_a = np.empty(12, dtype=np.float64)
    x_b = np.empty(12, dtype=np.float64)
    y_b = np.empty(12, dtype=np.float64)

    for i in range(4):
        x_a[i] = polygon[i, 0]
        y_a[i] = polygon[i, 1]

    size = _clip_at_boundary_numba(
        x_a, y_a, 4, 0, xmin, True, x_b, y_b
    )
    size = _clip_at_boundary_numba(
        x_b, y_b, size, 0, xmax, False, x_a, y_a
    )
    size = _clip_at_boundary_numba(
        x_a, y_a, size, 1, ymin, True, x_b, y_b
    )
    size = _clip_at_boundary_numba(
        x_b, y_b, size, 1, ymax, False, x_a, y_a
    )

    return _polygon_area_numba(x_a, y_a, size)


@njit(cache=True)
def _overlap_entries_numba(corners, valid_indices, xi_edges, eta_edges,
                           ny, nx):
    """Calculate the nonzero sparse-matrix entries for one frame."""

    # First count the candidate cells so the result arrays can be allocated
    # once. Some candidates will have zero true polygon overlap.
    candidate_count = 0
    for source_index in valid_indices:
        polygon = corners[source_index]
        xmin = polygon[0, 0]
        xmax = polygon[0, 0]
        ymin = polygon[0, 1]
        ymax = polygon[0, 1]

        for corner in range(1, 4):
            xmin = min(xmin, polygon[corner, 0])
            xmax = max(xmax, polygon[corner, 0])
            ymin = min(ymin, polygon[corner, 1])
            ymax = max(ymax, polygon[corner, 1])

        i_start = max(np.searchsorted(xi_edges, xmin, side="right") - 1, 0)
        i_stop = min(np.searchsorted(xi_edges, xmax, side="left"), nx)
        j_start = max(np.searchsorted(eta_edges, ymin, side="right") - 1, 0)
        j_stop = min(np.searchsorted(eta_edges, ymax, side="left"), ny)
        candidate_count += (i_stop - i_start) * (j_stop - j_start)

    target_indices = np.empty(candidate_count, dtype=np.int64)
    source_indices = np.empty(candidate_count, dtype=np.int64)
    overlap_areas = np.empty(candidate_count, dtype=np.float64)
    entry = 0

    for source_index in valid_indices:
        polygon = corners[source_index]
        xmin = polygon[0, 0]
        xmax = polygon[0, 0]
        ymin = polygon[0, 1]
        ymax = polygon[0, 1]

        for corner in range(1, 4):
            xmin = min(xmin, polygon[corner, 0])
            xmax = max(xmax, polygon[corner, 0])
            ymin = min(ymin, polygon[corner, 1])
            ymax = max(ymax, polygon[corner, 1])

        i_start = max(np.searchsorted(xi_edges, xmin, side="right") - 1, 0)
        i_stop = min(np.searchsorted(xi_edges, xmax, side="left"), nx)
        j_start = max(np.searchsorted(eta_edges, ymin, side="right") - 1, 0)
        j_stop = min(np.searchsorted(eta_edges, ymax, side="left"), ny)

        # Most of these pixels still need clipping, but a footprint wholly
        # inside one cell can use its full area directly.
        contained = (
            i_stop - i_start == 1
            and j_stop - j_start == 1
            and xmin >= xi_edges[i_start]
            and xmax <= xi_edges[i_start + 1]
            and ymin >= eta_edges[j_start]
            and ymax <= eta_edges[j_start + 1]
        )

        if contained:
            area = _polygon_area_numba(polygon[:, 0], polygon[:, 1], 4)
            if area > 0:
                target_indices[entry] = j_start * nx + i_start
                source_indices[entry] = source_index
                overlap_areas[entry] = area
                entry += 1
            continue

        for j in range(j_start, j_stop):
            for i in range(i_start, i_stop):
                area = _rectangle_overlap_numba(
                    polygon,
                    xi_edges[i],
                    xi_edges[i + 1],
                    eta_edges[j],
                    eta_edges[j + 1],
                )
                if area > 0:
                    target_indices[entry] = j * nx + i
                    source_indices[entry] = source_index
                    overlap_areas[entry] = area
                    entry += 1

    return (
        target_indices[:entry],
        source_indices[:entry],
        overlap_areas[:entry],
    )


#%% Sparse footprint-to-grid mapping

def overlap_mapping(mlat, mlt, grid):
    """Build one sparse source-pixel to target-cell overlap mapping."""

    corners, valid = infer_footprints(mlat, mlt, grid)
    xi_edges = np.asarray(grid.xi_mesh[0], dtype=float)
    eta_edges = np.asarray(grid.eta_mesh[:, 0], dtype=float)
    ny, nx = grid.shape

    corners = np.ascontiguousarray(corners.reshape(-1, 4, 2))
    valid_indices = np.flatnonzero(valid).astype(np.int64)
    target_indices, source_indices, overlap_areas = _overlap_entries_numba(
        corners,
        valid_indices,
        xi_edges,
        eta_edges,
        ny,
        nx,
    )

    mapping = csr_matrix(
        (overlap_areas, (target_indices, source_indices)),
        shape=(grid.size, np.asarray(mlat).size),
    )
    cell_area = np.diff(eta_edges)[:, None] * np.diff(xi_edges)[None, :]
    return mapping, cell_area


#%% Apply one mapping to source fields

def overlap_mean(values, mapping, output_shape):
    """Return an overlap-weighted mean and contributing overlap area."""

    values = np.asarray(values, dtype=float).ravel()
    finite = np.isfinite(values)
    filled = np.where(finite, values, 0)

    overlap = np.asarray(mapping @ finite.astype(float)).ravel()
    numerator = np.asarray(mapping @ filled).ravel()

    mean = np.full(mapping.shape[0], np.nan)
    covered = overlap > 0
    mean[covered] = numerator[covered] / overlap[covered]
    return mean.reshape(output_shape), overlap.reshape(output_shape)


def overlap_statistics(values, mapping, cell_area, output_shape):
    """Return mean, provisional spread, contributors, and coverage."""

    values = np.asarray(values, dtype=float).ravel()
    finite = np.isfinite(values)
    filled = np.where(finite, values, 0)

    overlap = np.asarray(mapping @ finite.astype(float)).ravel()
    numerator = np.asarray(mapping @ filled).ravel()
    numerator_squared = np.asarray(mapping @ (filled**2)).ravel()

    mean = np.full(mapping.shape[0], np.nan)
    spread = np.full(mapping.shape[0], np.nan)
    covered = overlap > 0
    mean[covered] = numerator[covered] / overlap[covered]
    variance = numerator_squared[covered] / overlap[covered] - mean[covered] ** 2
    spread[covered] = np.sqrt(np.maximum(variance, 0))

    contributors = mapping.copy()
    contributors.data[:] = 1
    count = np.asarray(contributors @ finite.astype(np.int32)).ravel()

    coverage = overlap.reshape(output_shape) / cell_area
    coverage = np.minimum(coverage, 1)

    return (
        mean.reshape(output_shape),
        spread.reshape(output_shape),
        count.reshape(output_shape).astype(np.int32),
        coverage,
    )
