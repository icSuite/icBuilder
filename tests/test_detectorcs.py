"""Focused tests for shared detector-to-CS geometry and reducers."""

import numpy as np
from scipy.sparse import csr_matrix

from icbuilder.detectorcs import (
    reduce_area_mean,
    reduce_covariance,
    reduce_flag_fraction,
    reduce_measurement_variance,
    reduce_support,
)
from icbuilder.grids import (
    DETECTOR_CS_GRID_ID,
    make_detector_cs_grid,
)


def test_detector_cs_grid_is_frozen_at_46_by_46():
    grid = make_detector_cs_grid()

    assert DETECTOR_CS_GRID_ID == "image_apex_130km_46x46_v1"
    assert grid.shape == (46, 46)


def test_two_footprint_reducers_match_hand_calculation():
    mapping = csr_matrix(np.array([[1.0, 3.0], [0.0, 2.0]]))
    cell_area = np.array([[5.0, 4.0]])
    output_shape = (1, 2)
    valid = np.array([True, True])

    mean, overlap = reduce_area_mean(
        np.array([2.0, 6.0]), valid, mapping, output_shape
    )
    variance, variance_overlap = reduce_measurement_variance(
        np.array([4.0, 9.0]), valid, mapping, output_shape
    )
    covariance, _ = reduce_covariance(
        np.array([2.0, -1.0]), valid, mapping, output_shape
    )
    fraction, any_flag = reduce_flag_fraction(
        np.array([False, True]), valid, mapping, output_shape
    )
    count, coverage = reduce_support(
        valid, mapping, cell_area, output_shape
    )

    np.testing.assert_allclose(mean, [[5.0, 6.0]])
    np.testing.assert_allclose(overlap, [[4.0, 2.0]])
    np.testing.assert_allclose(variance_overlap, overlap)
    np.testing.assert_allclose(variance, [[85 / 16, 9.0]])
    np.testing.assert_allclose(covariance, [[-7 / 16, -1.0]])
    np.testing.assert_allclose(fraction, [[0.75, 1.0]])
    np.testing.assert_array_equal(any_flag, [[True, True]])
    np.testing.assert_array_equal(count, [[2, 1]])
    np.testing.assert_allclose(coverage, [[0.8, 0.5]])


def test_field_specific_missing_values_have_independent_support():
    mapping = csr_matrix(np.array([[1.0, 1.0]]))
    valid = np.array([True, True])

    mean, overlap = reduce_area_mean(
        np.array([2.0, np.nan]), valid, mapping, (1, 1)
    )
    variance, variance_overlap = reduce_measurement_variance(
        np.array([4.0, 9.0]), valid, mapping, (1, 1)
    )

    np.testing.assert_allclose(mean, [[2.0]])
    np.testing.assert_allclose(overlap, [[1.0]])
    np.testing.assert_allclose(variance, [[13 / 4]])
    np.testing.assert_allclose(variance_overlap, [[2.0]])


def test_invalid_or_negative_variance_is_excluded():
    mapping = csr_matrix(np.array([[1.0, 1.0]]))
    variance, overlap = reduce_measurement_variance(
        np.array([-1.0, 9.0]),
        np.array([True, True]),
        mapping,
        (1, 1),
    )

    np.testing.assert_allclose(variance, [[9.0]])
    np.testing.assert_allclose(overlap, [[1.0]])
