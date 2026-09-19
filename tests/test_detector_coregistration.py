import numpy as np
from scipy.sparse import csr_matrix

from icbuilder.detector_coregistration import (
    geographic_to_wic,
    make_si_mapping,
    make_wic_transform,
    map_si,
)


def regular_camera(sensor, size, spacing):
    """Make a small regular geographic detector for coregistration tests."""
    row, column = np.indices((size, size))
    centre = (size - 1) / 2
    latitude = 70 + (row - centre) * spacing
    longitude = (column - centre) * spacing / np.cos(np.deg2rad(70))
    return {
        "sensor": sensor,
        "time": np.array(["2001-01-01"], dtype="datetime64[ms]"),
        "glat": latitude[None],
        "glon": longitude[None],
        "geometry_valid": np.ones((1, size, size), dtype=bool),
    }


def test_area_mapping_preserves_a_constant_field():
    wic = regular_camera("WIC", size=16, spacing=0.05)
    si13 = regular_camera("SI13", size=8, spacing=0.10)

    transform = make_wic_transform(wic, frame=0)
    mapping, _, diagnostics = make_si_mapping(si13, frame=0, transform=transform)
    values = np.ones((8, 8))
    mapped, _ = map_si(values, np.isfinite(values), mapping, shape=(16, 16))

    assert np.all(np.isfinite(mapped))
    assert np.allclose(mapped, 1)
    assert diagnostics["target_pixels_coverage_ge_0.9"] == 16 * 16


def test_area_mapping_rejects_insufficient_coverage():
    mapping = csr_matrix(np.array([[0.8], [0.9]]))
    values = np.array([[5.0]])
    mapped, coverage = map_si(
        values, np.isfinite(values), mapping, shape=(1, 2)
    )

    assert np.isnan(mapped[0, 0])
    assert mapped[0, 1] == 5
    assert np.allclose(coverage, [[0.8, 0.9]])


def test_wic_transform_excludes_pixels_behind_the_local_projection():
    wic = regular_camera("WIC", size=16, spacing=0.05)

    # A finite but geographically isolated source pixel can lie behind the
    # tangent plane. It is not a usable interpolation point.
    wic["glat"][0, 0, 0] = -70
    wic["glon"][0, 0, 0] = 180

    transform = make_wic_transform(wic, frame=0)

    assert transform["shape"] == (16, 16)


def test_geographic_lookup_skips_missing_detector_coordinates():
    wic = regular_camera("WIC", size=16, spacing=0.05)
    transform = make_wic_transform(wic, frame=0)

    latitude = np.full((16, 16), np.nan)
    longitude = np.full((16, 16), np.nan)
    latitude[4:6, 4:6] = wic["glat"][0, 4:6, 4:6]
    longitude[4:6, 4:6] = wic["glon"][0, 4:6, 4:6]

    queried_points = []
    row_interpolator = transform["row"]

    def record_row_query(points):
        queried_points.append(len(points))
        return row_interpolator(points)

    transform["row"] = record_row_query
    row, column, error = geographic_to_wic(transform, latitude, longitude)

    assert queried_points == [4]
    assert np.isfinite(row[4:6, 4:6]).all()
    assert np.isfinite(column[4:6, 4:6]).all()
    assert np.isfinite(error[4:6, 4:6]).all()
    assert np.isnan(row[:4]).all()
    assert np.isnan(column[:4]).all()
    assert np.isnan(error[:4]).all()

    row, column, error = geographic_to_wic(
        transform,
        np.full((4, 4), np.nan),
        np.full((4, 4), np.nan),
    )
    assert np.isnan(row).all()
    assert np.isnan(column).all()
    assert np.isnan(error).all()
    assert queried_points == [4]
