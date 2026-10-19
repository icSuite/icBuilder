"""Native2D runner preserves the original CLI, file layout and status arrays."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr


SCRIPT = Path(__file__).parents[1]/'scripts/pipeline/make_orbit_nc_files_2d.py'
SPEC = importlib.util.spec_from_file_location('make_orbit_nc_files_2d', SCRIPT)
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)
SENSORS = {'wic': 'WIC', 's12': 'SI12', 's13': 'SI13'}


@pytest.fixture
def setup_orbit(tmp_path, monkeypatch):
    inputs, events, file_lists = {}, [], {}
    raw, output = tmp_path/'raw', tmp_path/'output'
    for prefix, sensor in SENSORS.items():
        shape = (3, 5, 5)
        dates = np.datetime64('2001-01-01')+np.arange(3)*np.timedelta64(120, 's')
        fields = {name: (('date', 'row', 'col'), np.full(shape, value, dtype=float)) for name, value in [('img', 110), ('img_variance', 4), ('mlat', 70), ('glat', 70), ('sza', 80), ('dza', 20)]}
        fields.update({'frame_quality': ('date', [2, 0, 2]), 'hemisphere': ('date', ['north', 'south', 'north']), 'pointing_time': ('date', dates), 'id': xr.DataArray(sensor)})
        inputs[prefix] = xr.Dataset(fields, coords={'date': dates, 'row': np.arange(5), 'col': np.arange(5)})
        (raw/(prefix+'_data')).mkdir(parents=True)
        name = prefix+'_or0085.idl'
        (raw/(prefix+'_data')/name).write_bytes(b'synthetic raw file')
        file_lists[prefix] = [name]
        (output/prefix).mkdir(parents=True)

    def read(filenames, **settings):
        prefix = Path(filenames[0]).name.split('_')[0]
        assert settings['dzalim'] == 75 and settings['pointing_metadata']
        assert settings['skip_unmapped']
        assert settings['measurement_variance'] and settings['poisson_offset'] == .5
        assert settings['reflat'] == (prefix == 'wic')
        assert Path(filenames[0]) == raw/(prefix+'_data')/file_lists[prefix][0]
        events.append(('read', prefix))
        return inputs[prefix].copy(deep=True)

    def correct(cameras, **settings):
        assert set(cameras) == {'WIC', 'SI12', 'SI13'}
        assert all(data.sizes['date'] == 3 and 'south' in data.hemisphere.values for data in cameras.values())
        events.append(('correct', None))
        result = {sensor: data.copy(deep=True) for sensor, data in cameras.items()}
        for data in result.values():
            data.attrs['pointing_correction_applied'] = 1
            data['glat'] = data.glat+.01
        return result

    def fit(data, **settings):
        assert settings['dgsigma_freeze_iteration'] == 2
        assert data.attrs['pointing_correction_applied'] == 1 and np.all(data.hemisphere == 'north')
        prefix = next(name for name, sensor in SENSORS.items() if sensor == data.id.item())
        events.append(('fit', prefix))
        result = data.isel(date=data.frame_quality.values >= 2).copy(deep=True)
        for suffix, model, weight in [('', 100, .75), ('_independent', 90, .25)]:
            result['dgmodel'+suffix] = xr.full_like(result.img, model)
            result['dgimg'+suffix] = result.img-result['dgmodel'+suffix]
            result['dgweight'+suffix] = xr.full_like(result.img, weight)
            result['dgsigma'+suffix] = xr.full_like(result.img, 5)
        return result

    monkeypatch.setattr(RUNNER.fuv, 'read_idl', read)
    monkeypatch.setattr(RUNNER.fuv, 'correct_pointing', correct)
    monkeypatch.setattr(RUNNER.fuv, 'backgroundmodel_BS2D', fit)
    return raw, output, file_lists, inputs, events


def test_original_cli_defaults_and_flags_are_preserved():
    defaults = RUNNER.parse_args([])
    assert defaults.do_wic and defaults.do_s12 and defaults.do_s13
    assert defaults.parallel is False and defaults.pool_size == 10
    assert defaults.base_output is None
    args = RUNNER.parse_args(['--do_wic', 'false', '--do_s12', 'yes', '--do_s13', '0', '--parallel', 'true', '--pool_size', '3', '--base', '/input', '--base_output', '/output'])
    assert (args.do_wic, args.do_s12, args.do_s13, args.parallel) == (False, True, False, True)
    assert args.pool_size == 3 and args.base_input == '/input' and args.base_output == '/output'


def test_full_stream_pointing_and_both_models_saved_in_original_layout(setup_orbit):
    raw, output, names, inputs, events = setup_orbit
    frozen = {prefix: data.copy(deep=True) for prefix, data in inputs.items()}
    assert RUNNER.process_single_orbit(85, names, raw, output, tuple(SENSORS)) == {prefix: (85, 1) for prefix in SENSORS}
    assert [event[0] for event in events] == ['read']*3+['correct']+['fit']*3
    for prefix in SENSORS:
        with xr.open_dataset(output/prefix/(prefix+'_or0085.nc')) as saved:
            np.testing.assert_array_equal(saved.dgimg, 10)
            np.testing.assert_array_equal(saved.dgimg_independent, 20)
            np.testing.assert_array_equal(saved.dgweight, .75)
            np.testing.assert_array_equal(saved.dgweight_independent, .25)
            np.testing.assert_array_equal(saved.img_variance, 4)
            np.testing.assert_allclose(saved.glat, 70.01)
            assert saved.sizes['date'] == 2 and 't_start' in saved
        xr.testing.assert_identical(inputs[prefix], frozen[prefix])


def test_selected_sensor_still_loads_and_corrects_all_cameras(setup_orbit):
    raw, output, names, inputs, events = setup_orbit
    assert RUNNER.process_single_orbit(85, names, raw, output, ('s12',)) == {'s12': (85, 1)}
    assert [event for event in events if event[0] == 'fit'] == [('fit', 's12')]
    assert {event[1] for event in events if event[0] == 'read'} == set(SENSORS)
    assert not (output/'wic'/'wic_or0085.nc').exists()


@pytest.mark.parametrize('problem', ['quality', 'south', 'mapping'])
def test_wic_no_data_is_zero_and_does_not_block_si(setup_orbit, problem):
    raw, output, names, inputs, events = setup_orbit
    if problem == 'quality':
        inputs['wic'].frame_quality.values[:] = 0
    elif problem == 'south':
        inputs['wic'].hemisphere.values[:] = 'south'
    else:
        inputs['wic'].mlat.values[:] = np.nan
    result = RUNNER.process_single_orbit(83, names, raw, output, tuple(SENSORS))
    assert result == {'wic': (83, 0), 's12': (83, 1), 's13': (83, 1)}
    assert [event for event in events if event[0] == 'fit'] == [('fit', 's12'), ('fit', 's13')]
    assert not (output/'wic'/'wic_or0083.nc').exists()


def test_sensor_fit_failure_is_minus_one_and_other_sensors_continue(setup_orbit, monkeypatch):
    raw, output, names, inputs, events = setup_orbit
    ordinary = RUNNER.fuv.backgroundmodel_BS2D

    def fail_wic(data, **settings):
        if data.id.item() == 'WIC':
            raise ValueError('synthetic inversion failure')
        return ordinary(data, **settings)

    monkeypatch.setattr(RUNNER.fuv, 'backgroundmodel_BS2D', fail_wic)
    result = RUNNER.process_single_orbit(85, names, raw, output, tuple(SENSORS))
    assert result == {'wic': (85, -1), 's12': (85, 1), 's13': (85, 1)}


def test_missing_shared_stream_fails_selected_sensors_without_fitting(setup_orbit):
    raw, output, names, inputs, events = setup_orbit
    names['s13'] = []
    assert RUNNER.process_single_orbit(85, names, raw, output, ('wic', 's12')) == {'wic': (85, -1), 's12': (85, -1)}
    assert not any(event[0] == 'fit' for event in events)


def test_main_reads_hdf_and_writes_original_status_npy_for_selected_sensor(setup_orbit):
    raw, output, names, inputs, events = setup_orbit
    for prefix in SENSORS:
        pd.DataFrame({'orbit': [85], 'filename': names[prefix]}).to_hdf(raw/(prefix+'files.h5'), key='data')
    RUNNER.main(['--base_input', str(raw), '--base_output', str(output), '--do_wic', 'false', '--do_s12', 'true', '--do_s13', 'false', '--parallel', 'false'])
    status = np.load(output/'s12_avail_orbit.npy')
    np.testing.assert_array_equal(status, [[85, 1]])
    assert np.issubdtype(status.dtype, np.integer)
    assert not (output/'wic_avail_orbit.npy').exists() and not (output/'s13_avail_orbit.npy').exists()
    assert [event for event in events if event[0] == 'fit'] == [('fit', 's12')]


def test_orbit_union_and_availability_arrays_only_include_indexed_selected_sensors(setup_orbit, monkeypatch):
    raw, output, names, inputs, events = setup_orbit
    files = {prefix: pd.DataFrame({'orbit': orbits, 'filename': [prefix+'_'+str(orbit)+'.idl' for orbit in orbits]}) for prefix, orbits in [('wic', [83, 85]), ('s12', [85, 86]), ('s13', [85])]}
    calls = []

    def process(orbit, file_lists, inpath, outpath, sensors):
        calls.append((orbit, tuple(sensors)))
        return {prefix: (orbit, 1 if orbit == 85 else -1) for prefix in sensors}

    monkeypatch.setattr(RUNNER, 'process_single_orbit', process)
    result = RUNNER.background_removal(files, raw, output, ('wic', 's12'), parallel=False)
    assert calls == [(83, ('wic',)), (85, ('wic', 's12')), (86, ('s12',))]
    np.testing.assert_array_equal(result['wic'], [[83, -1], [85, 1]])
    np.testing.assert_array_equal(result['s12'], [[85, 1], [86, -1]])
    assert set(result) == {'wic', 's12'}
