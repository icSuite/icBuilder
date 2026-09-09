"""Historical FUVVIEW3 quiet-time background models.

This module is a transparent Python port of the parts of
``image_bckgnd_active_p.pro`` and ``image_bckgnd_active_p2.pro`` needed for
the WIC/SI12/SI13 comparison.  It works on images that are already in the
same flat-fielded count space as the archived quiet-time tables.  F10.7,
clock-angle, detector bias, and raw-image flat-field operations are not part
of this first test.

The IDL tables are logically indexed as ``[SZA, DZA]``.  SciPy restores them
as ``[DZA, SZA]``; all arrays in this module therefore use the Python order
``[DZA, SZA]``.
"""

from pathlib import Path

import numpy as np
from scipy.io import readsav


TABLE_FILES = {
    "WIC": "fit_bin_cts2.idlf",
    "SI12": "fit_s12_bin_cts2.idl",
    "SI13": "fit_s13_bin_cts2.idl",
}


def load_reference_background(sensor, support_directory):
    """Load one FUVVIEW3 quiet-time SZA/DZA background table."""

    sensor = sensor.upper()
    if sensor not in TABLE_FILES:
        raise ValueError(f"unknown IMAGE sensor: {sensor}")

    filename = Path(support_directory) / TABLE_FILES[sensor]
    if not filename.exists():
        raise FileNotFoundError(filename)

    saved = readsav(filename, python_dict=True)
    reference = np.asarray(saved["avg"], dtype=float)
    if reference.shape != (90, 110):
        raise ValueError(
            f"expected a (90, 110) DZA/SZA table, got {reference.shape}"
        )
    return reference


def _idl_round(values):
    """Round halves away from zero, matching IDL's ROUND for these angles."""

    values = np.asarray(values, dtype=float)
    return np.sign(values) * np.floor(np.abs(values) + 0.5)


def _upper_median(values):
    """Return IDL's default median: the upper middle value for an even count."""

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return np.partition(values, values.size // 2)[values.size // 2]


def _idl_median_filter(values, width):
    """Apply IDL MEDIAN's interior-only square filter."""

    values = np.asarray(values, dtype=float)
    result = values.copy()
    radius = width // 2

    for row in range(radius, values.shape[0] - radius):
        for column in range(radius, values.shape[1] - radius):
            window = values[
                row - radius:row + radius + 1,
                column - radius:column + radius + 1,
            ]
            result[row, column] = _upper_median(window)

    return result


def _median_arr(values, radius=5):
    """Port ``median_arr(values, 5, 0)`` from the recovered FUVVIEW3 source.

    Despite the IDL helper's comment, ``radius=5`` executes an 11-by-11
    window.  Zero-valued centres remain zero and zero-valued neighbours are
    excluded from the median.
    """

    values = np.asarray(values, dtype=float)
    result = np.zeros(values.shape, dtype=float)

    for row in range(values.shape[0]):
        row_start = max(0, row - radius)
        row_stop = min(values.shape[0], row + radius + 1)

        for column in range(values.shape[1]):
            if values[row, column] == 0:
                continue

            column_start = max(0, column - radius)
            column_stop = min(values.shape[1], column + radius + 1)
            window = values[row_start:row_stop, column_start:column_stop]
            neighbours = window[(window != 0) & np.isfinite(window)]
            if neighbours.size:
                result[row, column] = _upper_median(neighbours)

    return result


def _fill_nearest_1d(values):
    """Fill invalid entries from the nearest positive entry; ties go lower."""

    values = values.copy()
    good = np.flatnonzero(np.isfinite(values) & (values > 0))
    bad = np.flatnonzero(~np.isfinite(values) | (values <= 0))
    if good.size == 0:
        return values

    for index in bad:
        nearest = good[np.argmin(np.abs(good - index))]
        values[index] = values[nearest]
    return values


def _fill_array(values):
    """Port FUVVIEW3 ``fill_array``: DZA first, then SZA."""

    values = np.asarray(values, dtype=float).copy()

    # IDL array(sza, *) holds DZA at fixed SZA.  In Python that is a column.
    for sza_index in range(values.shape[1]):
        values[:, sza_index] = _fill_nearest_1d(values[:, sza_index])

    # IDL array(*, dza) holds SZA at fixed DZA.  In Python that is a row.
    for dza_index in range(values.shape[0]):
        values[dza_index, :] = _fill_nearest_1d(values[dza_index, :])

    return values


def _bin_detector_image(image, sza, dza, fit_mask, shape):
    """Assign detector pixels to integer SZA/DZA bins in source order."""

    binned = np.zeros(shape, dtype=float)
    image = np.asarray(image, dtype=float)
    sza = np.asarray(sza, dtype=float)
    dza = np.asarray(dza, dtype=float)

    # FUVVIEW3 uses floating subscripts here.  IDL truncates them, and the
    # final detector pixel encountered wins when several pixels share a bin.
    for value, solar_angle, detector_angle in zip(
        image[fit_mask].ravel(), sza[fit_mask].ravel(), dza[fit_mask].ravel()
    ):
        sza_index = int(max(solar_angle, 0.0))
        dza_index = int(max(detector_angle, 0.0))
        if dza_index < shape[0] and sza_index < shape[1]:
            binned[dza_index, sza_index] = value

    return binned


def _map_background(reference, sza, dza, outside_value=0.0):
    """Map a DZA/SZA table back to detector pixels using rounded angles.

    The historical IDL initializes the background to zero before filling the
    part covered by the lookup.  Keep that zero fallback here, but return the
    lookup-support mask separately so unsupported pixels remain identifiable.
    """

    sza = np.asarray(sza, dtype=float)
    dza = np.asarray(dza, dtype=float)
    if sza.shape != dza.shape:
        raise ValueError("SZA and DZA must have the same shape")

    background = np.full(sza.shape, outside_value, dtype=float)
    finite = np.isfinite(sza) & np.isfinite(dza)
    sza_index = np.zeros(sza.shape, dtype=int)
    dza_index = np.zeros(dza.shape, dtype=int)
    sza_index[finite] = _idl_round(sza[finite]).astype(int)
    dza_index[finite] = _idl_round(dza[finite]).astype(int)
    valid = (
        finite
        & (sza_index >= 0) & (sza_index < 110)
        & (dza_index >= 0) & (dza_index < 80)
    )

    # FUVVIEW3 rounds before applying the lookup limits.
    background[valid] = reference[dza_index[valid], sza_index[valid]]
    return background, valid


def _continue_high_sza(background, reference, sza, dza, continuation_sza):
    """Hold one SZA column constant at every larger SZA.

    The continuation preserves the background's DZA dependence.  It is an
    explicit extrapolation, returned as a separate mask, rather than an
    extension of the measured lookup support.
    """

    sza = np.asarray(sza, dtype=float)
    dza = np.asarray(dza, dtype=float)

    finite = np.isfinite(sza) & np.isfinite(dza)
    sza_index = np.zeros(sza.shape, dtype=int)
    dza_index = np.zeros(dza.shape, dtype=int)
    sza_index[finite] = _idl_round(sza[finite]).astype(int)
    dza_index[finite] = _idl_round(dza[finite]).astype(int)

    if continuation_sza < 0 or continuation_sza >= reference.shape[1]:
        raise ValueError("continuation SZA is outside the background table")

    extrapolated = (
        finite
        & (sza_index > continuation_sza)
        & (dza_index >= 0) & (dza_index < 80)
    )
    background[extrapolated] = reference[
        dza_index[extrapolated], continuation_sza
    ]

    return background, extrapolated


def fixed_background(reference, sza, dza):
    """Evaluate the unscaled FUVVIEW3 quiet-time background.

    Outside the table, use the historical zero-background fallback.  The WIC
    wrapper could replace this with a detector-bias surface, but that surface
    is not applied here because its count space is ambiguous after the saved
    image has already been calibrated.
    """

    reference = np.asarray(reference, dtype=float)
    background, support = _map_background(reference, sza, dza)
    return {
        "background": background,
        "support": support,
        "reference": reference.copy(),
        "method": "fuview_fixed",
    }


def active_background_p(image, reference, sza, dza, mlat,
                        continue_high_sza=False):
    """Evaluate the active SZA-rescaled background from the old ``p`` code.

    The recovered routine only fits scales below 105 degrees SZA, and its tail
    is affected by the zero-filled edge of the median-filtered fitting image.
    When ``continue_high_sza=True``, freeze both the scale and the complete
    DZA-dependent background profile at 100 degrees SZA.  This avoids letting
    a falling or zero reference background defeat the scale continuation.  The
    default retains the historical zero fallback so the source-faithful and
    experimental versions remain distinct.

    ``support`` records the original table domain and ``extrapolated`` records
    pixels receiving the new edge continuation.
    """

    image, sza, dza, mlat = _matching_arrays(image, sza, dza, mlat)
    reference = np.asarray(reference, dtype=float).copy()

    fit_mask = (
        np.isfinite(image) & np.isfinite(sza) & np.isfinite(dza)
        & np.isfinite(mlat) & (dza >= 0) & (dza < 70)
        & (sza >= 0) & (sza < 105)
        & ((np.abs(mlat) < 60) | (np.abs(mlat) > 75))
    )
    if not np.any(fit_mask):
        raise ValueError("no detector pixels satisfy the FUVVIEW p fitting mask")

    binned = _bin_detector_image(image, sza, dza, fit_mask, reference.shape)
    filtered = _idl_median_filter(binned, width=9)
    sampled_reference = np.where(filtered != 0, reference, 0.0)

    numerator = np.sum(filtered, axis=0)
    denominator = np.sum(sampled_reference, axis=0)
    scale = np.divide(
        numerator, denominator,
        out=np.zeros(numerator.shape), where=denominator != 0,
    )

    # Retain the complete support description for diagnostics.  The final few
    # fitted bins are not suitable extrapolation anchors because the median
    # window overlaps the zero-filled SZA >= 105 part of the fitting image.
    scale_support = (
        (np.arange(scale.size) < 105)
        & np.isfinite(scale) & (scale > 0)
        & (denominator > 0)
    )
    supported_sza = np.flatnonzero(scale_support)
    if supported_sza.size == 0:
        raise ValueError("active-p scaling produced no positive SZA scale")

    last_scale_sza = int(supported_sza[-1])
    if continue_high_sza:
        continuation_sza = 100
        if not scale_support[continuation_sza]:
            raise ValueError(
                "active-p edge continuation requires a positive scale at "
                "SZA=100"
            )
        scale[continuation_sza + 1:] = scale[continuation_sza]
    else:
        continuation_sza = None

    scaled_reference = reference * scale[None, :]
    if continue_high_sza:
        scaled_reference[:, continuation_sza + 1:] = scaled_reference[
            :, continuation_sza, None
        ]

    background, support = _map_background(scaled_reference, sza, dza)
    extrapolated = np.zeros(support.shape, dtype=bool)
    if continue_high_sza:
        background, extrapolated = _continue_high_sza(
            background, scaled_reference, sza, dza, continuation_sza
        )

    return {
        "background": background,
        "support": support,
        "extrapolated": extrapolated,
        "fit_mask": fit_mask,
        "binned_image": binned,
        "filtered_image": filtered,
        "scale": scale,
        "scale_support": scale_support,
        "last_scale_sza": last_scale_sza,
        "continuation_sza": continuation_sza,
        "reference": scaled_reference,
        "method": (
            "fuview_active_p_edge" if continue_high_sza
            else "fuview_active_p"
        ),
    }


def active_background_p2(image, reference, sza, dza, mlat, sensor,
                         include_polar_cap=False):
    """Evaluate the newer two-dimensional active ``p2`` background.

    The supplied IDL defines ``min_avg`` only in its WIC branch but uses it
    for every camera.  This port explicitly computes the same positive
    SZA>=80, DZA<=70 minimum for SI12 and SI13.  The repair is reported in the
    returned dictionary.
    """

    image, sza, dza, mlat = _matching_arrays(image, sza, dza, mlat)
    sensor = sensor.upper()
    if sensor not in TABLE_FILES:
        raise ValueError(f"unknown IMAGE sensor: {sensor}")

    reference = np.asarray(reference, dtype=float).copy()
    if sensor == "WIC":
        reference[:, 104:] = reference[:, 103, None]

    minimum_region = reference[:71, 80:]
    positive = minimum_region[np.isfinite(minimum_region) & (minimum_region > 0)]
    if positive.size == 0:
        raise ValueError("p2 reference has no positive high-SZA values")
    minimum = float(np.min(positive))

    # The literal SI branches fail because min_avg is undefined.  Apply the
    # complete WIC minimum/floor rule independently to each camera as the
    # smallest coherent repair, and label that repair in the result.
    reference = np.maximum(reference, minimum)

    polar_limit = 80 if include_polar_cap else 90
    fit_mask = (
        np.isfinite(image) & np.isfinite(sza) & np.isfinite(dza)
        & np.isfinite(mlat) & (dza >= 0) & (dza < 70)
        & (sza >= 0) & (sza < 111)
        & ((np.abs(mlat) < 60) | (np.abs(mlat) > polar_limit))
    )
    if not np.any(fit_mask):
        raise ValueError("no detector pixels satisfy the FUVVIEW p2 fitting mask")

    binned = _bin_detector_image(image, sza, dza, fit_mask, reference.shape)
    filtered = _median_arr(binned, radius=5)
    sampled_reference = np.where(filtered != 0, reference, 0.0)
    scale_2d = np.divide(
        filtered, sampled_reference,
        out=np.full(reference.shape, np.nan), where=sampled_reference != 0,
    )
    scale_2d = _fill_array(scale_2d)
    scaled_reference = reference * scale_2d

    scaled_region = scaled_reference[:71, 80:]
    positive = scaled_region[
        np.isfinite(scaled_region) & (scaled_region > 0)
    ]
    if positive.size == 0:
        raise ValueError("p2 scaling produced no positive high-SZA background")
    minimum = float(np.min(positive))

    background, support = _map_background(
        scaled_reference, sza, dza, outside_value=minimum
    )
    background[~np.isfinite(sza) | (sza < 0)] = np.nan

    return {
        "background": background,
        "support": support,
        "fit_mask": fit_mask,
        "binned_image": binned,
        "filtered_image": filtered,
        "scale": scale_2d,
        "reference": scaled_reference,
        "minimum_background": minimum,
        "p2_si_minimum_repair": sensor != "WIC",
        "method": "fuview_active_p2",
    }


def _matching_arrays(image, sza, dza, mlat):
    """Convert detector fields to float arrays and check their shapes."""

    arrays = [np.asarray(value, dtype=float) for value in (image, sza, dza, mlat)]
    if len({value.shape for value in arrays}) != 1:
        raise ValueError("image, SZA, DZA, and MLAT must have matching shapes")
    return arrays
