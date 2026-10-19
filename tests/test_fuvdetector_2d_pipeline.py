"""Integration checks for the opt-in raw-IDL to independent2D Product1 path."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from netCDF4 import Dataset
from icreader import open_product

from icbuilder.fuvdetector import BS2D_MODEL_CONTRACT, BS2D_PREPROCESSING_LABEL, load_detector_source
from test_fuvdetector_product import BASE_TIME, regular_camera, write_source_orbit


SCRIPT = Path(__file__).parents[1]/'scripts/pipeline/make_fuv_detector_orbit_files_2d.py'
SPEC = importlib.util.spec_from_file_location('fuvdetector_2d_pipeline', SCRIPT)
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


def native_camera(sensor='WIC'):
    dates = np.datetime64(BASE_TIME)+np.arange(3)*np.timedelta64(120, 's')
    camera = regular_camera(sensor, 5, .1, dates, 10)
    fields = {name: (('date', 'row', 'col'), camera[name]) for name in ['glat', 'glon', 'sza', 'dza']}
    fields.update({'img': (('date', 'row', 'col'), camera['unsubtracted_counts']), 'img_variance': (('date', 'row', 'col'), camera['variance']), 'frame_quality': ('date', [2, 0, 2]), 'hemisphere': ('date', ['north', 'south', 'north']), 'pointing_time': ('date', dates), 'id': xr.DataArray(sensor)})
    return xr.Dataset(fields, coords={'date': dates, 'row': np.arange(5), 'col': np.arange(5)})


def fitted_camera(data):
    data = data.isel(date=data.frame_quality.values >= 2).copy(deep=True)
    for suffix, model, weight in [('', 100, .75), ('_independent', 90, .25)]:
        data['dgmodel'+suffix] = xr.full_like(data.img, model)
        data['dgimg'+suffix] = data.img-data['dgmodel'+suffix]
        data['dgweight'+suffix] = xr.full_like(data.img, weight)
        data['dgsigma'+suffix] = xr.full_like(data.img, 5)
        data['dgmodel'+suffix].attrs = {**BS2D_MODEL_CONTRACT, 'sigma_model_coefficients': np.array([0., 25.]), 'noise_model': 'supplied_quadratic' if suffix else 'binned', 'shape_prior': 'none', 'temporal_knot_spacing_target_minutes': 120., 'robust_initialization': 'fresh all-ones; updated every fit'}
    data['dgfit_converged'] = True
    data['dgfit_converged_independent'] = ('date', np.ones(data.sizes['date'], dtype=bool))
    data.attrs.update({'pointing_correction_applied': 1, 'pointing_reference_sensor': 'WIC', 'pointing_correction_method': 'clock_phase'})
    return data


def test_old_default_and_new_label_select_their_own_counts_and_weights(tmp_path):
    legacy = tmp_path/'legacy.nc'
    write_source_orbit(legacy)
    old = load_detector_source(legacy, 'WIC')
    np.testing.assert_array_equal(old['counts'], 10)
    np.testing.assert_array_equal(old['quality_weight'], .75)
    new = tmp_path/'two_stage.nc'
    RUNNER.save_native(fitted_camera(native_camera()), new, 'WIC')
    selected = load_detector_source(new, 'WIC', preprocessing_label=BS2D_PREPROCESSING_LABEL)
    assert selected['image_field'] == 'dgimg_independent'
    np.testing.assert_array_equal(selected['counts'], 20)
    np.testing.assert_array_equal(selected['quality_weight'], .25)
    np.testing.assert_array_equal(selected['unsubtracted_counts'], 110)
    with pytest.raises(ValueError, match='does not match'):
        load_detector_source(new, 'WIC')


@pytest.mark.parametrize('field,value', [('sigma', 1), ('pointing', 0), ('temporal', 0), ('independent', 0)])
def test_new_label_rejects_wrong_sigma_pointing_or_failed_fit(tmp_path, field, value):
    path = tmp_path/'invalid.nc'
    RUNNER.save_native(fitted_camera(native_camera()), path, 'WIC')
    with Dataset(path, 'r+') as nc:
        if field == 'sigma':
            nc['dgmodel'].dgsigma_freeze_iteration = value
        elif field == 'pointing':
            nc.pointing_correction_applied = value
        else:
            nc['dgfit_converged'+('_independent' if field == 'independent' else '')][:] = value
    with pytest.raises(ValueError):
        load_detector_source(path, 'WIC', preprocessing_label=BS2D_PREPROCESSING_LABEL)


def test_raw_pipeline_corrects_full_streams_saves_both_and_restarts(tmp_path, monkeypatch):
    events, files = [], {}
    inputs = {sensor: native_camera(sensor) for sensor in ['WIC', 'SI12', 'SI13']}
    originals = {sensor: data.copy(deep=True) for sensor, data in inputs.items()}
    for sensor in inputs:
        path = tmp_path/(sensor+'.idl')
        path.write_bytes(b'synthetic raw input')
        files[sensor] = [path]

    def read(filenames, **settings):
        sensor = Path(filenames[0]).stem
        assert settings['pointing_metadata'] and settings['measurement_variance']
        assert settings['reflat'] == (sensor == 'WIC') and settings['poisson_offset'] == .5
        events.append(('read', sensor))
        return inputs[sensor].copy(deep=True)

    def correct(sensors, method):
        assert method == 'clock_phase' and set(sensors) == set(inputs)
        assert all(data.sizes['date'] == 3 and 'south' in data.hemisphere.values for data in sensors.values())
        events.append(('correct', None))
        output = {sensor: data.copy(deep=True) for sensor, data in sensors.items()}
        for data in output.values():
            data['glat'] = data.glat+.01
            data.attrs['pointing_correction_applied'] = 1
        return output

    def fit(data, **settings):
        assert settings['dgsigma_freeze_iteration'] == 2
        assert data.attrs['pointing_correction_applied'] == 1 and data.sizes['date'] == 2
        assert np.all(data.hemisphere.values == 'north')
        events.append(('fit', str(data.id.item())))
        return fitted_camera(data)

    monkeypatch.setattr(RUNNER.fuv, 'read_idl', read)
    monkeypatch.setattr(RUNNER.fuv, 'correct_pointing', correct)
    monkeypatch.setattr(RUNNER.fuv, 'backgroundmodel_BS2D', fit)
    output = tmp_path/'output'
    assert RUNNER.process_orbit((85, files), output, 'fuv_detector') == (85, 2)
    assert [event[0] for event in events] == ['read']*3+['correct']+['fit']*3
    paths = RUNNER.native_paths(output, 85)
    for sensor, path in paths.items():
        with xr.open_dataset(path) as data:
            np.testing.assert_array_equal(data.dgimg, 10)
            np.testing.assert_array_equal(data.dgimg_independent, 20)
            np.testing.assert_array_equal(data.dgweight_independent, .25)
            np.testing.assert_array_equal(data.img_variance, 4)
            np.testing.assert_allclose(data.glat, originals[sensor].glat.isel(date=[0, 2])+.01)
        xr.testing.assert_identical(inputs[sensor], originals[sensor])
    product = output/'fuv_detector'/BS2D_PREPROCESSING_LABEL/'or_0085.nc'
    with open_product(product) as data:
        assert data.attrs['background_fit_stage'] == 'independent'
        np.testing.assert_array_equal(data.read('wic_counts'), 20)
        np.testing.assert_array_equal(data.read('wic_quality_weight'), .25)
    previous_events = events.copy()
    assert RUNNER.process_orbit((85, files), output, 'fuv_detector') == (85, 'complete')
    assert events == previous_events
    files['WIC'][0].write_bytes(b'changed raw input')
    with pytest.raises(ValueError, match='overwrite'):
        RUNNER.process_orbit((85, files), output, 'fuv_detector')
