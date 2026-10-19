"""Native FUV orbit files with shared metadata pointing and 2D backgrounds."""

#%% Import

from os.path import join as pjoin
from pathlib import Path
import fuvpy as fuv
import pandas as pd
import numpy as np
from multiprocessing import Pool, freeze_support
import argparse

#%% Selected fuvpy preprocessing configuration

BACKGROUND_MODEL_SETTINGS = {
    'dgsigma_freeze_iteration': 2,
    'temporal_knot_spacing': 120,
    'variance_weighted': True,
    'stop': 0.01,
    'tukeyVal': 5,
    'dampingVal': 1e-2
}
SENSOR_IDS = {'wic': 'WIC', 's12': 'SI12', 's13': 'SI13'}

#%% Argument parsing

def str2bool(v):
    if isinstance(v, bool):
        return v
    v = v.lower()
    if v in ('yes', 'true', 't', '1'):
        return True
    elif v in ('no', 'false', 'f', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Process FUV data with 2D backgrounds and shared metadata pointing.")

    parser.add_argument('--do_wic', type=str2bool, default=True, help='Fit and save WIC data (default True)')
    parser.add_argument('--do_s12', type=str2bool, default=True, help='Fit and save SI12 data (default True)')
    parser.add_argument('--do_s13', type=str2bool, default=True, help='Fit and save SI13 data (default True)')
    parser.add_argument('--parallel', type=str2bool, default=False, help='Run in parallel (default False)')
    parser.add_argument('--pool_size', type=int, default=10, help='Size of pool (default 10)')
    parser.add_argument('--base_input', '--base', dest='base_input', type=str,
                        default=str(pjoin(Path(__file__).resolve().parents[2], 'example_data')),
                        help='Base directory containing orbit indices and raw IDL data')
    parser.add_argument('--base_output', type=str, default=None,
                        help='Base directory for orbit products; defaults to base_input')

    return parser.parse_args(argv)

#%% One orbit: shared pointing, then fit each selected sensor

def process_single_orbit(orbit, file_lists, inpath, outpath, sensors):
    results = {prefix: (orbit, -1) for prefix in sensors}
    try:
        print(f'Starting on orbit {orbit}', flush=True)
        images = {}
        for prefix, sensor in SENSOR_IDS.items():
            if not file_lists[prefix]:
                raise ValueError(f'No {sensor} input files for shared pointing in orbit {orbit}')
            filenames = [pjoin(inpath, prefix + '_data', name) for name in file_lists[prefix]]
            images[sensor] = fuv.read_idl(
                filenames, dzalim=75, reflat=(prefix == 'wic'),
                measurement_variance=True, poisson_offset=0.5,
                wic_variance_model='reference', viewing_geometry=True, pointing_metadata=True
            )

        # Correct all sensors together, before selecting hemisphere or quality.
        images = fuv.correct_pointing(images, method='clock_phase')

    except Exception as e:
        print(f'Shared pointing : {orbit} : failed with error {e}', flush=True)
        return results

    for prefix in sensors:
        try:
            s = images[SENSOR_IDS[prefix]]
            s = s.sel(date=s.hemisphere.date[s.hemisphere == 'north'])

            if np.all(np.isnan(s['mlat'].values)) or not np.any(s['frame_quality'].values >= 2):
                print(f'Skipping {prefix} orbit {orbit}, no data', flush=True)
                results[prefix] = (orbit, 0)
                continue

            # FUVPY converges temporal, then fits independent; both are retained.
            s = fuv.backgroundmodel_BS2D(s, **BACKGROUND_MODEL_SETTINGS)
            s = s.assign({'t_start': np.datetime_as_string(s['date'][0], unit='s')})

            outfile = pjoin(outpath, prefix, f"{prefix}_or{str(orbit).zfill(4)}.nc")
            encoding = {var: {"zlib": True, "complevel": 4} for var in s.data_vars}
            s.to_netcdf(outfile, format="NETCDF4", encoding=encoding)
            results[prefix] = (orbit, 1)

        except Exception as e:
            print(f'{prefix} : {orbit} : failed with error {e}', flush=True)

    return results

#%% Existing serial/parallel orbit processing

def background_removal_parallel(args_list, pool_size=10):
    with Pool(pool_size) as pool:
        results = pool.starmap(process_single_orbit, args_list)
    return results

def background_removal_serial(args_list):
    results = []
    for args in args_list:
        results.append(process_single_orbit(*args))
    return results

def background_removal(files, inpath, outpath, sensors, parallel=True, pool_size=10):
    """Process orbits jointly; return each sensor's original (orbit, status) array.

    Status is 0 for no data, 1 for success, and -1 for failure. All three raw
    sensor streams are needed for shared pointing; only selected sensors fit
    backgrounds and write files.
    """
    for prefix in sensors:
        Path(pjoin(outpath, prefix)).mkdir(parents=True, exist_ok=True)

    orbits = np.unique(np.concatenate([files[prefix]['orbit'].unique() for prefix in sensors]))
    args_list = []
    for orbit in orbits:
        file_lists = {prefix: table.loc[table['orbit'] == orbit, 'filename'].tolist() for prefix, table in files.items()}
        orbit_sensors = tuple(prefix for prefix in sensors if file_lists[prefix])
        args_list.append((orbit, file_lists, inpath, outpath, orbit_sensors))

    if parallel:
        results = background_removal_parallel(args_list, pool_size)
    else:
        results = background_removal_serial(args_list)

    avail_orbit = {}
    for prefix in sensors:
        rows = [result[prefix] for result in results if prefix in result]
        avail_orbit[prefix] = np.asarray(rows, dtype=int).reshape(-1, 2)
    return avail_orbit

#%% Run selected sensors

def main(argv=None):
    args = parse_args(argv)
    base_input = args.base_input
    base_output = args.base_output if args.base_output is not None else base_input

    print(f'Data settings:\n WIC: {args.do_wic}\n SI12: {args.do_s12}\n SI13: {args.do_s13}\n')
    print(f'Compute settings:\n Parallel: {args.parallel}\n Pool: {args.pool_size}\n')
    print(f'Input base set to {base_input}')
    print(f'Output base set to {base_output}')

    print('Reading h5 files')
    files = {
        'wic': pd.read_hdf(pjoin(base_input, 'wicfiles.h5'), key='data'),
        's12': pd.read_hdf(pjoin(base_input, 's12files.h5'), key='data'),
        's13': pd.read_hdf(pjoin(base_input, 's13files.h5'), key='data')
    }
    enabled = {'wic': args.do_wic, 's12': args.do_s12, 's13': args.do_s13}
    sensors = tuple(prefix for prefix, selected in enabled.items() if selected)
    if not sensors:
        return

    avail_orbit = background_removal(files, base_input, base_output, sensors,
                                     parallel=args.parallel, pool_size=args.pool_size)
    for prefix in sensors:
        np.save(pjoin(base_output, prefix + '_avail_orbit.npy'), avail_orbit[prefix])


if __name__ == '__main__':
    freeze_support()
    main()
