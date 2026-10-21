"""Estimate WIC--SI13 relative alignment from paired bright Earth limbs.

The calculation uses only raw IMAGE/FUV images and their pointing metadata.
It performs two distinct tests: local-reference fits compare time-dependent
motion between the cameras, while an orbit-wide absolute fit tests whether a
fixed two-angle offset remains identifiable after fitting effective limb
height. The local test cannot measure a constant alignment offset.

The script reads the existing ``wicfiles.h5`` and ``s13files.h5`` catalogs; it
never lists the raw-data directories. Per-orbit files make a stopped run safe
to resume. FUVPY supplies ``frame_quality`` and ``pointing_valid``; only
northern science-ready frames with valid pointing enter either fit. This is a
diagnostic and does not modify IDL or pipeline products.
"""

#%% Imports and fixed detector geometry

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path
import warnings

import fuvpy as fuv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import readsav
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from scipy.spatial.transform import Rotation


EARTH_EQUATOR_KM = 6378.136
EARTH_POLE_KM = EARTH_EQUATOR_KM*(1-1/298.257)
PROFILE_DISTANCE_PX = np.arange(-15, 15.01, .25)


def camera_rotation(record):
    """Reproduce the FUVVIEW3 detector-to-GCI rotation."""
    omega, theta, phi, psi = np.deg2rad([record[key] for key in ['INST_ROLL', 'INST_CO_ELEV', 'INST_AZIMUTH', 'SPIN_PHASE']])
    a, b, c = np.asarray(record['SC_SV'], dtype=float)
    d, e, f = np.asarray(record['A_GCI3'], dtype=float)
    root_a = np.sqrt(1-a*a)
    root_af = np.sqrt(1-a*a-f*f)
    sin_alpha, cos_alpha = b/root_a, c/root_a
    sin_beta, cos_beta = a, root_a
    cos_eta, sin_eta = f/root_a, root_af/root_a
    cos_delta = (a*d+e*root_af)/(1-f*f)
    sin_delta = (a*e-d*root_af)/(1-f*f)
    t1a = np.array([[np.cos(omega), -np.sin(omega), 0], [np.sin(omega), np.cos(omega), 0], [0, 0, 1]])
    t1b = np.array([[1, 0, 0], [0, -np.sin(theta), -np.cos(theta)], [0, np.cos(theta), -np.sin(theta)]])
    t1c = np.array([[np.sin(phi), np.cos(phi), 0], [-np.cos(phi), np.sin(phi), 0], [0, 0, 1]])
    t2a = np.array([[1, 0, 0], [0, cos_alpha, -sin_alpha], [0, sin_alpha, cos_alpha]])
    t2b = np.array([[cos_beta, 0, -sin_beta], [0, 1, 0], [sin_beta, 0, cos_beta]])
    t3a = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])
    t3b = np.array([[cos_beta, 0, sin_beta], [0, 1, 0], [-sin_beta, 0, cos_beta]])
    t3c = np.array([[1, 0, 0], [0, cos_eta, sin_eta], [0, -sin_eta, cos_eta]])
    t3d = np.array([[cos_delta, -sin_delta, 0], [sin_delta, cos_delta, 0], [0, 0, 1]])
    return t3d @ t3c @ t3b @ t3a @ t2b @ t2a @ t1c @ t1b @ t1a


def corrected_rotation(record, x_deg, y_deg):
    delta = Rotation.from_rotvec(np.deg2rad([x_deg, y_deg, 0])).as_matrix()
    return camera_rotation(record) @ delta


def tangent_limb(record, transform=None, height_delta_km=0):
    """Return the projected ellipsoidal limb and inward detector normals."""
    transform = camera_rotation(record) if transform is None else transform
    observer = np.asarray(record['O_GCI'], dtype=float)
    equator = EARTH_EQUATOR_KM + float(record['EMIS_HGT']) + height_delta_km
    pole = EARTH_POLE_KM + float(record['EMIS_HGT']) + height_delta_km
    metric = np.diag([1, 1, (equator/pole)**2])
    w = metric @ observer
    shell = observer @ w-equator**2
    conic = transform.T @ (np.outer(w, w)-shell*metric) @ transform
    centre = np.linalg.solve(conic[:2, :2], conic[:2, 2])
    radius = conic[2, 2]-conic[:2, 2] @ centre
    eigenvalue, eigenvector = np.linalg.eigh(-conic[:2, :2])
    if np.any(eigenvalue <= 0) or radius <= 0:
        raise ValueError('non-elliptical projected limb')
    angle = np.linspace(0, 2*np.pi, 360, endpoint=False)
    tangent = centre[:, None] + eigenvector*np.sqrt(radius/eigenvalue)[None, :] @ np.array([np.cos(angle), np.sin(angle)])
    increment = np.deg2rad([record['ANGL_RES_C'], record['ANGL_RES_R']])
    point = 127 + np.arctan(tangent).T/increment
    derivative = 2*(conic[:2, :2] @ tangent-conic[:2, 2, None]).T
    normal = derivative*increment*(1+tangent.T**2)
    normal /= np.linalg.norm(normal, axis=1, keepdims=True)
    return point, normal, conic


def shifted_limb(record, point, normal, parameters):
    """Intersect adjusted limb with the original fixed sampling normals."""
    x_deg, y_deg, height_100km = parameters
    _, _, conic = tangent_limb(record, corrected_rotation(record, x_deg, y_deg), 100*height_100km)
    increment = np.deg2rad([record['ANGL_RES_C'], record['ANGL_RES_R']])
    shift = np.zeros(len(point))
    for _ in range(10):
        tangent = np.tan((point+shift[:, None]*normal-127)*increment)
        q = np.column_stack([tangent, -np.ones(len(point))])
        value = np.einsum('ni,ij,nj->n', q, conic, q)
        gradient = 2*(q @ conic)[:, :2]*increment*(1+tangent**2)
        slope = np.sum(gradient*normal, axis=1)
        shift -= value/slope
    if not np.all(np.isfinite(shift)):
        raise ValueError('adjusted limb intersection failed')
    return shift


#%% Bright-edge measurement

def edge_samples(record, solar_max_deg):
    point, normal, _ = tangent_limb(record)
    image = np.asarray(record['IMAGE'], dtype=float)
    smooth = gaussian_filter(image, 2)
    gy, gx = np.gradient(smooth)
    scale = max(float(np.nanpercentile(image, 95)), 1)
    rows = []
    for p, n in zip(point, normal):
        inside = p+8*n
        if not np.all((p > 3) & (p < 252)) or not np.all((inside > 2) & (inside < 253)):
            continue
        solar = map_coordinates(np.asarray(record['SZA'], dtype=float), [[inside[1]], [inside[0]]], order=1, mode='nearest')[0]
        if not 0 <= solar < 90:
            continue
        profile = p[None, :] + PROFILE_DISTANCE_PX[:, None]*n[None, :]
        visible = np.all((profile > 2) & (profile < 253), axis=1)
        coordinates = np.array([profile[:, 1], profile[:, 0]])
        derivative = n[0]*map_coordinates(gx, coordinates, order=1, mode='nearest') + n[1]*map_coordinates(gy, coordinates, order=1, mode='nearest')
        candidates = np.flatnonzero(visible)
        if len(candidates) < 20:
            continue
        peak = candidates[np.argmax(derivative[candidates])]
        if peak <= candidates[0]+1 or peak >= candidates[-1]-1 or derivative[peak]/scale <= .005:
            continue
        left, centre, right = derivative[peak-1:peak+2]
        curvature = left-2*centre+right
        correction = .5*(left-right)/curvature if curvature else 0
        if abs(correction) > 1:
            continue
        offset = PROFILE_DISTANCE_PX[peak]+.25*correction
        observed = p+offset*n
        dx = map_coordinates(gx, [[observed[1]], [observed[0]]], order=1, mode='nearest')[0]
        dy = map_coordinates(gy, [[observed[1]], [observed[0]]], order=1, mode='nearest')[0]
        gradient = dx*n[0]+dy*n[1]
        alignment = gradient/max(np.hypot(dx, dy), np.finfo(float).eps)
        sample = lambda distance: map_coordinates(smooth, [[observed[1]+distance*n[1]], [observed[0]+distance*n[0]]], order=1, mode='nearest')[0]
        inside_counts = np.mean([sample(distance) for distance in [3, 5, 7]])
        outside_counts = np.mean([sample(distance) for distance in [-7, -5, -3]])
        normal_angle = np.rad2deg(np.arctan2(-n[0], n[1]))
        rows.append((p[0], p[1], n[0], n[1], normal_angle, offset, solar, gradient, alignment, inside_counts, outside_counts))
    columns = ['predicted_x', 'predicted_y', 'inward_nx', 'inward_ny', 'normal_angle_deg', 'offset_px', 'solar_angle_deg', 'gradient', 'alignment', 'inside_counts', 'outside_counts']
    table = pd.DataFrame(rows, columns=columns).sort_values('normal_angle_deg').reset_index(drop=True)
    if table.empty:
        raise ValueError('no measurable dayside edge')
    dayside = table.solar_angle_deg.between(0, solar_max_deg) & (table.alignment >= .9)
    intensity_limit = .6*table.loc[dayside, 'inside_counts'].quantile(.95)
    gradient_limit = .35*table.loc[dayside, 'gradient'].quantile(.95)
    table['selected'] = (dayside &
                         (table.inside_counts >= intensity_limit) &
                         (table.gradient >= gradient_limit) & (table.inside_counts > table.outside_counts))
    return table


def absolute_arc(table, settings):
    indices = np.flatnonzero(table.selected)
    if len(indices) < settings['minimum_samples']:
        raise ValueError('too_few_selected_samples')
    blocks = np.split(indices, np.flatnonzero(np.diff(table.normal_angle_deg.to_numpy()[indices]) > 3)+1)
    arc = table.iloc[max(blocks, key=len)].copy().reset_index(drop=True)
    span = float(np.ptp(arc.normal_angle_deg))
    if len(arc) < settings['minimum_samples'] or span < settings['minimum_arc_deg']:
        raise ValueError('insufficient_contiguous_arc')
    fraction = (arc.normal_angle_deg-arc.normal_angle_deg.min())/span
    arc['training'] = fraction.between(.25, .75)
    arc['held_out'] = (fraction <= .15) | (fraction >= .85)
    if arc.training.sum() < 6 or arc.held_out.sum() < 4:
        raise ValueError('insufficient_training_or_validation_samples')
    return arc


#%% Absolute orbit fit: fixed camera offset plus orbit-specific limb height

def fit_absolute_sensor(records, arcs, settings):
    result = {'status': 'insufficient_frames', 'frames': len(arcs)}
    if len(arcs) < 3:
        return result
    frame_ids = sorted(arcs)
    bounds = np.array([settings['rotation_bound_deg'], settings['rotation_bound_deg'], settings['height_bound_km']/100])

    def residual(value, held_out=False):
        values = []
        for frame in frame_ids:
            table = arcs[frame]
            point = table[['predicted_x', 'predicted_y']].to_numpy()
            normal = table[['inward_nx', 'inward_ny']].to_numpy()
            mask = table.held_out.to_numpy() if held_out else table.training.to_numpy()
            values.extend((shifted_limb(records[frame], point, normal, value)-table.offset_px.to_numpy())[mask])
        return np.asarray(values)

    fit = least_squares(residual, np.zeros(3), bounds=(-bounds, bounds), loss='soft_l1', f_scale=.5,
                        xtol=1e-9, ftol=1e-9, gtol=1e-9, max_nfev=200)
    held_before = residual(np.zeros(3), True)
    held_after = residual(fit.x, True)
    before_rms = float(np.sqrt(np.mean(held_before**2)))
    after_rms = float(np.sqrt(np.mean(held_after**2)))
    singular = np.linalg.svd(fit.jac*bounds[None, :], compute_uv=False)
    condition = float(singular[0]/singular[-1]) if singular[-1] > 0 else np.inf
    bound_fraction = float(np.max(np.abs(fit.x)/bounds))
    accepted = (fit.success and condition <= settings['maximum_condition'] and bound_fraction < .98 and
                after_rms <= settings['maximum_heldout_rms_px'] and
                after_rms <= before_rms+settings['maximum_heldout_degradation_px'])
    result.update({'status': 'accepted' if accepted else 'validation_failed', 'x_deg': float(fit.x[0]),
                   'y_deg': float(fit.x[1]), 'height_delta_km': float(100*fit.x[2]), 'condition': condition,
                   'bound_fraction': bound_fraction, 'heldout_before_rms_px': before_rms,
                   'heldout_after_rms_px': after_rms})
    return result


#%% Relative time-dependent fit: local limb reference, no absolute alignment

def interpolate_selected(table, angle):
    x = table.normal_angle_deg.to_numpy()
    upper = np.searchsorted(x, angle)
    lower = np.clip(upper-1, 0, len(x)-1)
    upper = np.clip(upper, 0, len(x)-1)
    gap = x[upper]-x[lower]
    valid = table.selected.to_numpy()[lower] & table.selected.to_numpy()[upper]
    valid &= (gap > 0) & (gap < 3) & (angle >= x[0]) & (angle <= x[-1])
    values = np.interp(angle, x, table.offset_px)
    values[~valid] = np.nan
    return values


def local_reference_tables(samples, time, segment, settings):
    result = {}
    cadence = np.median(np.diff(time))
    for frame, source in samples.items():
        own = source.loc[source.selected].copy().reset_index(drop=True)
        if len(own) < settings['minimum_samples']:
            continue
        angle = own.normal_angle_deg.to_numpy()
        neighbours = [other for other in range(max(0, frame-5), min(len(time), frame+6))
                      if other != frame and other in samples and segment[other] == segment[frame]
                      and abs(time[other]-time[frame]) <= 5.5*cadence]
        if len(neighbours) < 3:
            continue
        offsets = np.array([interpolate_selected(samples[other], angle) for other in neighbours])
        count = np.isfinite(offsets).sum(axis=0)
        left = np.isfinite(offsets[np.asarray(neighbours) < frame]).sum(axis=0)
        right = np.isfinite(offsets[np.asarray(neighbours) > frame]).sum(axis=0)
        reference = np.nansum(offsets, axis=0)/np.maximum(count, 1)
        indices = np.flatnonzero((count >= 3) & (left >= 1) & (right >= 1))
        if len(indices) < settings['minimum_samples']:
            continue
        indices = max(np.split(indices, np.flatnonzero(np.diff(angle[indices]) > 3)+1), key=len)
        table = own.iloc[indices].copy().reset_index(drop=True)
        span = float(np.ptp(table.normal_angle_deg))
        if len(table) < settings['minimum_samples'] or span < settings['minimum_arc_deg']:
            continue
        table['observed_difference_px'] = table.offset_px-reference[indices]
        fraction = (table.normal_angle_deg-table.normal_angle_deg.min())/span
        table['training'] = fraction.between(.3, .7)
        table['held_out'] = (fraction <= .15) | (fraction >= .85)
        if table.training.sum() >= 3 and table.held_out.sum() >= 4:
            result[frame] = table
    return result


def relative_prediction(record, table, angles):
    point = table[['predicted_x', 'predicted_y']].to_numpy()
    normal = table[['inward_nx', 'inward_ny']].to_numpy()
    return shifted_limb(record, point, normal, [angles[0], angles[1], 0])


def fit_relative_run(records, tables, time, settings):
    count = len(records)
    step = 1e-5
    weak, conditions = [], []
    for record, table in zip(records, tables):
        training = table.training.to_numpy()
        derivative = np.column_stack([(relative_prediction(record, table, step*np.eye(2)[axis])-
                                       relative_prediction(record, table, -step*np.eye(2)[axis]))/(2*step)
                                      for axis in range(2)])[training]
        _, singular, vectors = np.linalg.svd(derivative, full_matrices=False)
        direction = vectors[-1]
        if weak and direction @ weak[-1] < 0:
            direction = -direction
        weak.append(direction)
        conditions.append(float(singular[0]/singular[-1]))
    weak = np.asarray(weak)
    midpoint = weak[1:]+weak[:-1]
    midpoint /= np.linalg.norm(midpoint, axis=1, keepdims=True)
    cadence = np.median(np.diff(time)) if count > 1 else 120
    sizes = [int(table.training.sum()) for table in tables]
    sparsity = lil_matrix((sum(sizes)+count-1, 2*count), dtype=int)
    lower = 0
    for frame, size in enumerate(sizes):
        sparsity[lower:lower+size, 2*frame:2*frame+2] = 1
        lower += size
    for frame in range(count-1):
        sparsity[lower+frame, 2*frame:2*frame+4] = 1

    def residual(value):
        angles = value.reshape(count, 2)
        terms = []
        for record, table, angle in zip(records, tables, angles):
            training = table.training.to_numpy()
            terms.append((relative_prediction(record, table, angle)-table.observed_difference_px.to_numpy())[training])
        if count > 1:
            terms.append(np.sum(np.diff(angles, axis=0)*midpoint, axis=1)/(.4*np.sqrt(np.diff(time)/cadence)))
        return np.concatenate(terms)

    bound = settings['relative_rotation_bound_deg']
    fit = least_squares(residual, np.zeros(2*count), bounds=(-bound, bound), loss='soft_l1', f_scale=.5,
                        jac_sparsity=sparsity.tocsr(), max_nfev=150)
    return fit.x.reshape(count, 2), np.asarray(conditions), bool(fit.success)


def fit_relative_series(records, samples, time, segment, settings):
    tables = local_reference_tables(samples, time, segment, settings)
    indices = np.array(sorted(tables), dtype=int)
    runs = [] if not len(indices) else np.split(indices, np.flatnonzero((np.diff(indices) != 1) | (np.diff(segment[indices]) != 0))+1)
    result, pending = {}, list(runs)
    while pending:
        run = pending.pop(0)
        if len(run) < 3:
            continue
        values, conditions, success = fit_relative_run([records[frame] for frame in run], [tables[frame] for frame in run], time[run], settings)
        rejected = (np.full(len(run), not success) | (conditions > settings['relative_maximum_condition']) |
                    (np.max(np.abs(values), axis=1) > .98*settings['relative_rotation_bound_deg']))
        if np.any(rejected):
            accepted = run[~rejected]
            if len(accepted):
                pending.extend(np.split(accepted, np.flatnonzero(np.diff(accepted) != 1)+1))
            continue
        for frame, value, condition in zip(run, values, conditions):
            table = tables[frame]
            held = table.held_out.to_numpy()
            observed = table.observed_difference_px.to_numpy()
            after = observed-relative_prediction(records[frame], table, value)
            result[int(frame)] = {'status': 'accepted', 'x_deg': float(value[0]), 'y_deg': float(value[1]),
                                  'condition': float(condition), 'heldout_samples': int(held.sum()),
                                  'heldout_before_rms_px': float(np.sqrt(np.mean(observed[held]**2))),
                                  'heldout_after_rms_px': float(np.sqrt(np.mean(after[held]**2)))}
    return result


#%% Catalog matching, orbit processing and resumable output

def record_time(record):
    year_day, milliseconds = np.asarray(record['TIME']).astype(int)
    return pd.Timestamp(datetime.strptime(str(year_day), '%Y%j')+timedelta(milliseconds=int(milliseconds)))


def read_record(path):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return readsav(path, python_dict=True)['imageinfo'][0]


def load_frame_flags(base, sensor, filenames):
    """Run the production loader and index its two quality flags by exact time."""
    prefix = {'wic': 'wic', 'si13': 's13'}[sensor]
    paths = [str(Path(base)/f'{prefix}_data'/filename) for filename in filenames]
    data = fuv.read_idl(paths, dzalim=75, reflat=(sensor == 'wic'), viewing_geometry=True,
                        pointing_metadata=True, skip_unmapped=True)
    result = {}
    for frame in range(data.sizes['date']):
        timestamp = data.pointing_time.values[frame].astype('datetime64[ns]').astype(np.int64)
        if int(timestamp) in result:
            raise ValueError(f'{sensor} loader returned duplicate pointing_time {data.pointing_time.values[frame]}')
        result[int(timestamp)] = {'hemisphere': str(data.hemisphere.values[frame]),
                                  'frame_quality': int(data.frame_quality.values[frame]),
                                  'pointing_valid': bool(data.pointing_valid.values[frame])}
    audit = {'catalog_files': len(filenames), 'loaded_frames': data.sizes['date'],
             'skipped_unmapped_files': len(data.attrs.get('skipped_unmapped_files', '').splitlines()),
             'skipped_invalid_time_files': len(data.attrs.get('skipped_invalid_time_files', '').splitlines())}
    return result, audit


def flag_rejection(sensor, flags):
    if flags is None:
        return f'{sensor}_loader_skipped'
    if flags['hemisphere'] != 'north':
        return f'{sensor}_not_northern'
    if not flags['pointing_valid']:
        return f'{sensor}_pointing_invalid'
    if flags['frame_quality'] < 2:
        return f'{sensor}_frame_quality_below_2'
    return None


def add_rotation_fields(row, wic_record, si13_record, wic_fit, si13_fit, prefix):
    wic_nominal, si13_nominal = camera_rotation(wic_record), camera_rotation(si13_record)
    wic_adjusted = corrected_rotation(wic_record, wic_fit['x_deg'], wic_fit['y_deg'])
    si13_adjusted = corrected_rotation(si13_record, si13_fit['x_deg'], si13_fit['y_deg'])
    nominal_relative = wic_nominal.T @ si13_nominal
    adjusted_relative = wic_adjusted.T @ si13_adjusted
    relative_adjustment = adjusted_relative @ nominal_relative.T
    matrices = {'wic_gci': wic_adjusted @ wic_nominal.T, 'si13_gci': si13_adjusted @ si13_nominal.T,
                'nominal_relative': nominal_relative, 'relative_adjustment': relative_adjustment}
    for label, matrix in matrices.items():
        vector = np.rad2deg(Rotation.from_matrix(matrix).as_rotvec())
        for component, value in zip('xyz', vector):
            row[f'{prefix}_{label}_{component}_deg'] = float(value)
    row[f'{prefix}_relative_adjustment_magnitude_deg'] = float(np.linalg.norm(
        [row[f'{prefix}_relative_adjustment_{component}_deg'] for component in 'xyz']))


def process_orbit(task):
    orbit, inventory, base, output, settings, overwrite = task
    pairs = inventory['pairs']
    output = Path(output)
    result_file = output/'orbits'/f'or_{orbit:04d}.csv'
    status_file = output/'orbits'/f'or_{orbit:04d}.json'
    if status_file.exists() and not overwrite:
        previous = json.loads(status_file.read_text())
        if previous.get('settings') != settings:
            raise ValueError(f'Orbit {orbit} was produced with different settings; use --overwrite or a new output directory')
        return orbit, 'cached'
    frame_flags, loader_audit, loader_errors = {}, {}, {}
    for sensor in ['wic', 'si13']:
        try:
            frame_flags[sensor], loader_audit[sensor] = load_frame_flags(base, sensor, inventory[f'{sensor}_files'])
        except Exception as error:
            frame_flags[sensor], loader_audit[sensor] = {}, {}
            loader_errors[sensor] = str(error)[:500]
    rows, records, samples, arcs = [], {'wic': {}, 'si13': {}}, {'wic': {}, 'si13': {}}, {'wic': {}, 'si13': {}}
    for frame, pair in enumerate(pairs):
        row = {'orbit': orbit, **pair}
        if loader_errors:
            row['pair_status'] = 'loader_failed'
            rows.append(row)
            continue
        try:
            wic_record = read_record(Path(base)/'wic_data'/pair['wic_filename'])
            si13_record = read_record(Path(base)/'s13_data'/pair['si13_filename'])
            wic_time, si13_time = record_time(wic_record), record_time(si13_record)
            row.update({'wic_time_utc': wic_time.isoformat(), 'si13_time_utc': si13_time.isoformat(),
                        'exact_time_difference_seconds': float((si13_time-wic_time).total_seconds()),
                        'wic_spin': int(wic_record['SPIN']), 'si13_spin': int(si13_record['SPIN'])})
            wic_key = int(wic_time.to_datetime64().astype('datetime64[ns]').astype(np.int64))
            si13_key = int(si13_time.to_datetime64().astype('datetime64[ns]').astype(np.int64))
            wic_flags, si13_flags = frame_flags['wic'].get(wic_key), frame_flags['si13'].get(si13_key)
            for sensor, flags in [('wic', wic_flags), ('si13', si13_flags)]:
                if flags is not None:
                    row.update({f'{sensor}_{name}': value for name, value in flags.items()})
            if abs(row['exact_time_difference_seconds']) > settings['match_seconds']:
                row['pair_status'] = 'exact_time_mismatch'
            else:
                rejection = flag_rejection('wic', wic_flags) or flag_rejection('si13', si13_flags)
                if rejection is not None:
                    row['pair_status'] = rejection
                else:
                    row['pair_status'] = 'edge_candidate'
                    for sensor, record in [('wic', wic_record), ('si13', si13_record)]:
                        records[sensor][frame] = record
                        try:
                            samples[sensor][frame] = edge_samples(record, settings['solar_max_deg'])
                            row[f'{sensor}_edge_status'] = 'measured'
                            row[f'{sensor}_selected_edge_samples'] = int(samples[sensor][frame].selected.sum())
                            arcs[sensor][frame] = absolute_arc(samples[sensor][frame], settings)
                            row[f'{sensor}_absolute_arc_status'] = 'usable'
                        except Exception as error:
                            row[f'{sensor}_edge_status'] = str(error).replace(' ', '_')[:100]
                            row[f'{sensor}_absolute_arc_status'] = 'unusable'
        except Exception as error:
            row['pair_status'] = 'read_or_geometry_error'
            row['error'] = str(error)[:300]
        rows.append(row)
    table = pd.DataFrame(rows)
    if table.empty:
        table = pd.DataFrame(columns=['orbit', 'catalog_time', 'wic_filename', 'si13_filename', 'pair_status'])
        result_file.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(result_file)
        status = {'orbit': orbit, 'candidate_pairs': 0, 'status_counts': {}, 'relative_supported_frames': 0,
                  'absolute_common_candidate_frames': 0,
                  'absolute_fit': {'wic': {'status': 'insufficient_frames', 'frames': 0},
                                   'si13': {'status': 'insufficient_frames', 'frames': 0}},
                  'loader_audit': loader_audit, 'loader_errors': loader_errors, 'settings': settings}
        status_file.write_text(json.dumps(status, indent=2)+'\n')
        return orbit, 0

    # Relative motion uses leave-self-out local references. Its zero point is
    # arbitrary, so these fits test common motion rather than constant offset.
    catalog_time = pd.to_datetime(table.catalog_time)
    time = (catalog_time-catalog_time.iloc[0]).dt.total_seconds().to_numpy()
    cadence = np.median(np.diff(time)) if len(time) > 1 else 120
    relative = {}
    for sensor in ['wic', 'si13']:
        spin_values = table[f'{sensor}_spin'] if f'{sensor}_spin' in table else pd.Series(np.nan, index=table.index)
        spin = pd.to_numeric(spin_values, errors='coerce').to_numpy()
        breaks = np.r_[True, (np.diff(time) > 1.5*cadence) | ~np.isfinite(spin[1:]) |
                       ~np.isfinite(spin[:-1]) | (np.diff(spin) != 1)]
        segment = np.cumsum(breaks)
        relative[sensor] = fit_relative_series(records[sensor], samples[sensor], time, segment, settings)
    common_relative = sorted(set(relative['wic']) & set(relative['si13']))
    for frame in common_relative:
        for sensor in ['wic', 'si13']:
            for name, value in relative[sensor][frame].items():
                table.loc[frame, f'{sensor}_relative_{name}'] = value
        table.loc[frame, 'pair_status'] = 'relative_accepted'
        rotation_values = {}
        add_rotation_fields(rotation_values, records['wic'][frame], records['si13'][frame],
                            relative['wic'][frame], relative['si13'][frame], 'local')
        for name, value in rotation_values.items():
            table.loc[frame, name] = value

    # The absolute fit is deliberately separate. It can identify a constant
    # offset only when one rotation and one effective height explain an orbit.
    common_absolute = sorted(set(arcs['wic']) & set(arcs['si13']))
    absolute = {sensor: fit_absolute_sensor(records[sensor],
                {frame: arcs[sensor][frame] for frame in common_absolute}, settings) for sensor in ['wic', 'si13']}
    if absolute['wic']['status'] == absolute['si13']['status'] == 'accepted':
        for frame in common_absolute:
            rotation_values = {}
            add_rotation_fields(rotation_values, records['wic'][frame], records['si13'][frame],
                                absolute['wic'], absolute['si13'], 'absolute')
            for name, value in rotation_values.items():
                table.loc[frame, name] = value
            for sensor in ['wic', 'si13']:
                for name, value in absolute[sensor].items():
                    if name != 'status':
                        table.loc[frame, f'absolute_{sensor}_{name}'] = value

    result_file.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(result_file.with_suffix('.csv.tmp'), index=False)
    result_file.with_suffix('.csv.tmp').replace(result_file)
    status = {'orbit': orbit, 'candidate_pairs': len(table), 'status_counts': table.pair_status.value_counts().to_dict(),
              'relative_supported_frames': len(common_relative), 'absolute_common_candidate_frames': len(common_absolute),
              'absolute_fit': absolute, 'loader_audit': loader_audit, 'loader_errors': loader_errors,
              'settings': settings}
    status_file.with_suffix('.json.tmp').write_text(json.dumps(status, indent=2)+'\n')
    status_file.with_suffix('.json.tmp').replace(status_file)
    return orbit, len(common_relative)


def catalog_pairs(base, chosen_orbits, match_seconds):
    tables = {}
    for sensor, label in [('wic', 'wic'), ('s13', 'si13')]:
        table = pd.read_hdf(base/f'{sensor}files.h5', key='data').reset_index()
        table = table.rename(columns={table.columns[0]: 'catalog_time', 'filename': f'{label}_filename'})
        table['catalog_time'] = pd.to_datetime(table.catalog_time)
        tables[sensor] = table
    orbits = sorted(set(tables['wic'].orbit) & set(tables['s13'].orbit))
    if chosen_orbits:
        missing = sorted(set(chosen_orbits)-set(orbits))
        if missing:
            print('No paired catalog entries for orbits:', ', '.join(map(str, missing)), flush=True)
        orbits = [orbit for orbit in orbits if orbit in chosen_orbits]
    result = {}
    for orbit in orbits:
        wic = tables['wic'].loc[tables['wic'].orbit == orbit, ['catalog_time', 'wic_filename']].sort_values('catalog_time')
        si13 = tables['s13'].loc[tables['s13'].orbit == orbit, ['catalog_time', 'si13_filename']].sort_values('catalog_time')
        paired = pd.merge_asof(wic, si13, on='catalog_time', direction='nearest', tolerance=pd.Timedelta(seconds=match_seconds))
        paired = paired.dropna().drop_duplicates('si13_filename')
        result[int(orbit)] = {'pairs': paired.to_dict('records'),
                              'wic_files': tables['wic'].loc[tables['wic'].orbit == orbit, 'wic_filename'].tolist(),
                              'si13_files': tables['s13'].loc[tables['s13'].orbit == orbit, 'si13_filename'].tolist()}
    return result


#%% Orbit-level inference and compact diagnostic plots

def bootstrap_orbit_medians(orbit_table, columns, repetitions=5000):
    values = orbit_table[columns].to_numpy(float)
    rng = np.random.default_rng(20261020)
    samples = np.empty((repetitions, len(columns)))
    for i in range(repetitions):
        samples[i] = np.median(values[rng.integers(0, len(values), len(values))], axis=0)
    return np.percentile(samples, [2.5, 97.5], axis=0)


def summarize(output):
    files = sorted((output/'orbits').glob('or_*.csv'))
    tables = [pd.read_csv(filename) for filename in files]
    all_pairs = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    all_pairs.to_csv(output/'paired_frame_results.csv', index=False)
    absolute_fit_rows = []
    for filename in sorted((output/'orbits').glob('or_*.json')):
        status = json.loads(filename.read_text())
        row = {'orbit': status['orbit'], 'common_candidate_frames': status['absolute_common_candidate_frames']}
        for sensor in ['wic', 'si13']:
            row.update({f'{sensor}_{name}': value for name, value in status['absolute_fit'][sensor].items()})
        absolute_fit_rows.append(row)
    absolute_fit_summary = pd.DataFrame(absolute_fit_rows)
    absolute_fit_summary.to_csv(output/'absolute_fit_summary.csv', index=False)
    local = all_pairs.loc[all_pairs.pair_status == 'relative_accepted'].copy() if len(all_pairs) else pd.DataFrame()
    absolute_columns = [f'absolute_relative_adjustment_{component}_deg' for component in 'xyz']
    absolute = all_pairs.dropna(subset=absolute_columns).copy() if len(all_pairs) and set(absolute_columns) <= set(all_pairs) else pd.DataFrame()
    summary = {'processed_orbits': len(files), 'candidate_pairs': len(all_pairs),
               'relative_supported_pairs': len(local), 'absolute_supported_orbits': int(absolute.orbit.nunique()) if len(absolute) else 0}
    if len(absolute_fit_summary):
        summary['wic_absolute_status_counts'] = absolute_fit_summary.wic_status.value_counts().to_dict()
        summary['si13_absolute_status_counts'] = absolute_fit_summary.si13_status.value_counts().to_dict()
    if local.empty and absolute.empty:
        (output/'global_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
        return summary
    rows = []
    local_orbits = set(local.orbit) if len(local) else set()
    absolute_orbits_available = set(absolute.orbit) if len(absolute) else set()
    for orbit in sorted(local_orbits | absolute_orbits_available):
        row = {'orbit': int(orbit)}
        relative_orbit = local.loc[local.orbit == orbit]
        if len(relative_orbit):
            wic = relative_orbit[[f'local_wic_gci_{component}_deg' for component in 'xyz']].to_numpy()
            si13 = relative_orbit[[f'local_si13_gci_{component}_deg' for component in 'xyz']].to_numpy()
            row.update({'relative_supported_pairs': len(relative_orbit),
                        'relative_motion_correlation': float(np.corrcoef((wic-wic.mean(0)).ravel(), (si13-si13.mean(0)).ravel())[0, 1]),
                        'relative_vector_difference_median_deg': float(relative_orbit.local_relative_adjustment_magnitude_deg.median())})
        absolute_orbit = absolute.loc[absolute.orbit == orbit] if len(absolute) else absolute
        if len(absolute_orbit):
            row.update({f'absolute_relative_{component}_median_deg': float(absolute_orbit[column].median())
                        for component, column in zip('xyz', absolute_columns)})
            row['absolute_wic_height_delta_km'] = float(absolute_orbit.absolute_wic_height_delta_km.iloc[0])
            row['absolute_si13_height_delta_km'] = float(absolute_orbit.absolute_si13_height_delta_km.iloc[0])
        rows.append(row)
    orbit_summary = pd.DataFrame(rows).set_index('orbit')
    orbit_summary.to_csv(output/'orbit_summary.csv')
    absolute_summary_columns = [f'absolute_relative_{component}_median_deg' for component in 'xyz']
    absolute_orbits = (orbit_summary.dropna(subset=absolute_summary_columns)
                       if set(absolute_summary_columns) <= set(orbit_summary) else orbit_summary.iloc[:0])
    if len(absolute_orbits):
        columns = [f'absolute_relative_{component}_median_deg' for component in 'xyz']
        ci = bootstrap_orbit_medians(absolute_orbits, columns) if len(absolute_orbits) > 1 else np.full((2, 3), np.nan)
        summary['absolute_relative_adjustment_orbit_median_deg'] = dict(zip('xyz', absolute_orbits[columns].median().to_numpy(float)))
        summary['absolute_relative_adjustment_orbit_bootstrap_95_ci_deg'] = {
            component: [float(ci[0, i]), float(ci[1, i])] for i, component in enumerate('xyz')}
    summary['interpretation'] = ('Local-reference results test whether WIC and SI13 share time-dependent pointing motion; their zero point cannot measure a constant offset. '
                                 'Only orbits passing the separate absolute rotation-plus-height fit contribute to a fixed boresight estimate. The z component is derived, not limb-constrained roll.')
    (output/'global_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=True)+'\n')

    figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    if len(absolute_orbits):
        for component, color in zip('xyz', ['tab:blue', 'tab:orange', 'tab:green']):
            axes[0, 0].plot(absolute_orbits.index, absolute_orbits[f'absolute_relative_{component}_median_deg'], '.', color=color, label=component)
    else:
        axes[0, 0].text(.5, .5, 'No orbit passed the absolute fit', ha='center', va='center', transform=axes[0, 0].transAxes)
    if len(local):
        for component, color in zip('xyz', ['tab:blue', 'tab:orange', 'tab:green']):
            axes[0, 1].scatter(local[f'local_wic_gci_{component}_deg'], local[f'local_si13_gci_{component}_deg'], s=8, alpha=.35, color=color, label=component)
        limits = axes[0, 1].get_xlim()+axes[0, 1].get_ylim()
        axes[0, 1].plot([min(limits), max(limits)], [min(limits), max(limits)], ':', color='black')
        axes[1, 0].hist(local.local_relative_adjustment_magnitude_deg, bins=40)
        axes[1, 1].scatter(local.wic_relative_heldout_after_rms_px, local.si13_relative_heldout_after_rms_px, s=8, alpha=.35)
    axes[0, 0].set(xlabel='Orbit', ylabel='Absolute relative adjustment (deg)', title='Fixed-offset estimates passing absolute validation')
    axes[0, 1].set(xlabel='WIC local GCI correction (deg)', ylabel='SI13 local GCI correction (deg)', title='Common time-dependent motion')
    axes[1, 0].set(xlabel='Local WIC--SI13 vector difference (deg)', ylabel='Frames', title='Relative-motion disagreement')
    axes[1, 1].set(xlabel='WIC held-out RMS (px)', ylabel='SI13 held-out RMS (px)', title='Local-reference fit quality')
    for axis in axes.ravel():
        axis.grid(alpha=.2)
        if axis.get_legend_handles_labels()[0]:
            axis.legend()
    figure.savefig(output/'alignment_summary.png', dpi=180)
    plt.close(figure)
    return summary


#%% Command-line entry point

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=Path('/disk/IMAGE_FUV'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--orbits', type=int, nargs='*', help='Optional bounded orbit list; default is every paired orbit')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--match-seconds', type=float, default=10)
    parser.add_argument('--solar-max-deg', type=float, default=65)
    parser.add_argument('--minimum-samples', type=int, default=12)
    parser.add_argument('--minimum-arc-deg', type=float, default=12)
    parser.add_argument('--rotation-bound-deg', type=float, default=.75)
    parser.add_argument('--height-bound-km', type=float, default=500)
    parser.add_argument('--maximum-condition', type=float, default=200)
    parser.add_argument('--maximum-heldout-rms-px', type=float, default=1.5)
    parser.add_argument('--maximum-heldout-degradation-px', type=float, default=.1)
    parser.add_argument('--relative-rotation-bound-deg', type=float, default=.5)
    parser.add_argument('--relative-maximum-condition', type=float, default=80)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be at least one')
    args.output.mkdir(parents=True, exist_ok=True)
    settings = {name: getattr(args, name) for name in ['match_seconds', 'solar_max_deg', 'minimum_samples',
                'minimum_arc_deg', 'rotation_bound_deg', 'height_bound_km', 'maximum_condition',
                'maximum_heldout_rms_px', 'maximum_heldout_degradation_px', 'relative_rotation_bound_deg',
                'relative_maximum_condition']}
    protocol_file = args.output/'analysis_protocol.json'
    protocol = {'method': 'paired local-reference motion test plus separate orbit-wide absolute rotation-and-height test',
                'settings': settings, 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if protocol_file.exists() and not args.overwrite:
        previous = json.loads(protocol_file.read_text())
        if previous.get('settings') != settings or previous.get('script_sha256') != protocol['script_sha256']:
            parser.error('Output directory was produced by different settings or code; use --overwrite or a new output directory')
    protocol_file.write_text(json.dumps(protocol, indent=2)+'\n')
    pairs = catalog_pairs(args.base, set(args.orbits or []), args.match_seconds)
    tasks = [(orbit, inventory, str(args.base), str(args.output), settings, args.overwrite) for orbit, inventory in pairs.items()]
    print(f"Paired WIC--SI13 limb analysis: {len(tasks)} orbits, {sum(len(inventory['pairs']) for inventory in pairs.values())} catalog pairs", flush=True)
    if args.workers == 1:
        iterator = map(process_orbit, tasks)
    else:
        executor = ProcessPoolExecutor(max_workers=args.workers)
        iterator = executor.map(process_orbit, tasks, chunksize=1)
    try:
        for completed, (orbit, accepted) in enumerate(iterator, 1):
            print(f'[{completed}/{len(tasks)}] orbit {orbit}: {accepted}', flush=True)
    finally:
        if args.workers != 1:
            executor.shutdown()
    summary = summarize(args.output)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
