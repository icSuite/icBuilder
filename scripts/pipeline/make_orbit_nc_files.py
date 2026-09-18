#%% Import

from os.path import join as pjoin
from pathlib import Path
import fuvpy as fuv
import pandas as pd
import numpy as np
from multiprocessing import Pool
import argparse

#%% Selected fuvpy preprocessing configuration

BACKGROUND_MODEL_SETTINGS = {
    'spatial_model': 'legacy_x_azimuth',
    'directional_x_knots': [-3.5, 0, 1.5, 3.5],
    'temporal_knot_spacing': 120,
    'directional_temporal_knot_spacing': 60,
    'spencer_taper': True,
    'variance_weighted': True,
    'monotonic': True,
    'monotonic_fallback': True,
    'stop': 0.01,
    'tukeyVal': 5,
    'dampingVal': 1e-2
}

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

parser = argparse.ArgumentParser(description="Process FUV data.")

parser.add_argument('--do_wic', type=str2bool, default=True, help='Process WIC data (default True)')
parser.add_argument('--do_s12', type=str2bool, default=True, help='Process SI12 data (default True)')
parser.add_argument('--do_s13', type=str2bool, default=True, help='Process SI13 data (default True)')
parser.add_argument('--parallel', type=str2bool, default=False, help='Run in parallel (default False)')
parser.add_argument('--pool_size', type=int, default=10, help='Size of pool (default 10)')
parser.add_argument('--base_input', '--base', dest='base_input', type=str,
                    default=str(pjoin(Path(__file__).resolve().parents[2], 'example_data')),
                    help='Base directory containing orbit indices and raw IDL data')
parser.add_argument('--base_output', type=str, default=None,
                    help='Base directory for orbit products; defaults to base_input')

args = parser.parse_args()

#%% What data to process, how to do it, and where it is

do_wic = args.do_wic
do_s12 = args.do_s12
do_s13 = args.do_s13

parallel = args.parallel
pool_size = args.pool_size

base_input = args.base_input
base_output = args.base_output if args.base_output is not None else base_input

print(f'Data settings:\n WIC: {do_wic}\n SI12: {do_s12}\n SI13: {do_s13}\n')
print(f'Compute settings:\n Parallel: {parallel}\n Pool: {pool_size}\n')
print(f'Input base set to {base_input}')
print(f'Output base set to {base_output}')

#%% Import orbit files file 

print('Reading h5 files')
wicfiles = pd.read_hdf(pjoin(base_input, 'wicfiles.h5'), key='data')
s12files = pd.read_hdf(pjoin(base_input, 's12files.h5'), key='data')
s13files = pd.read_hdf(pjoin(base_input, 's13files.h5'), key='data')

#%%

def process_single_orbit(orbit, files, inpath, outpath, reflat, file_prefix):
    try:
        print(f'Starting on orbit {orbit}')
        file_list = [pjoin(inpath, f) for f in files.loc[files['orbit'] == orbit, 'filename']]
        
        s = fuv.read_idl(file_list, dzalim=75, reflat=reflat,
                         measurement_variance=True, poisson_offset=0.5,
                         wic_variance_model='reference', viewing_geometry=True)
        s = s.sel(date=s.hemisphere.date[s.hemisphere == 'north'])
        s = s.assign({'t_start': np.datetime_as_string(s['date'][0], unit='s')})

        if np.all(np.isnan(s['mlat'].values)):
            print(f'Skipping orbit {orbit}, no data')
            return (orbit, 0)

        s = fuv.backgroundmodel_BS(s, **BACKGROUND_MODEL_SETTINGS)

        outfile = pjoin(outpath, f"{file_prefix}_or{str(orbit).zfill(4)}.nc")
        encoding = {var: 
                    {"zlib": True, "complevel": 4}
                    for var in s.data_vars
                    }
        s.to_netcdf(outfile, format="NETCDF4", encoding=encoding)

        return (orbit, 1)

    except Exception as e:
        print(f'{file_prefix} : {orbit} : failed with error {e}')
        return (orbit, -1)

def background_removal_parallel(files, inpath, outpath, reflat=False):
    file_prefix = files['filename'].iloc[0][:3]
    orbits = files['orbit'].unique()
    args_list = [(orbit, files, inpath, outpath, reflat, file_prefix) for orbit in orbits]

    with Pool(pool_size) as pool: # 500 GB RAM. Max 33 GB per orbit, I think. Max 15 process
        results = pool.starmap(process_single_orbit, args_list)

    return np.array(results)

def background_removal_serial(files, inpath, outpath, reflat=False):
    file_prefix = files['filename'].iloc[0][:3]
    orbits = files['orbit'].unique()
    results = []

    for orbit in orbits:
        result = process_single_orbit(orbit, files, inpath, outpath, reflat, file_prefix)
        results.append(result)

    return np.array(results)

def background_removal(files, inpath, outpath, reflat=False, parallel=True):
    """
    Background removal per orbit

    Parameters
    ----------
    files : DataFrame
        Must contain columns 'orbit' and 'filename'
    inpath : str
        Path to input files
    outpath : str
        Path to save output netCDFs
    reflat : bool
        Whether to reflatten images
    parallel : bool
        If True, use multiprocessing. If False, run serially.

    Returns
    -------
    avail_orbit : np.ndarray
        2D array with orbit numbers and status (0 = no data, 1 = success, -1 = failure)
    """
    Path(outpath).mkdir(parents=True, exist_ok=True)

    if parallel:
        return background_removal_parallel(files, inpath, outpath, reflat)
    else:
        return background_removal_serial(files, inpath, outpath, reflat)

#%% Run WIC
if do_wic:
    inpath  = pjoin(base_input, 'wic_data')
    outpath = pjoin(base_output, 'wic')

    print('Starting work on WIC')
    print('Pulling data from: ' + inpath)
    print('Offlaoding at: ' + outpath)
    avail_orbit = background_removal(wicfiles, inpath, outpath, reflat=True, parallel=parallel)
    np.save(pjoin(base_output, 'wic_avail_orbit.npy'), avail_orbit)

#%% Run s12
if do_s12:
    inpath  = pjoin(base_input, 's12_data')
    outpath = pjoin(base_output, 's12')

    print('Starting work on SI12')
    print('Pulling data from: ' + inpath)
    print('Offlaoding at: ' + outpath)
    avail_orbit = background_removal(s12files, inpath, outpath, parallel=parallel)
    np.save(pjoin(base_output, 's12_avail_orbit.npy'), avail_orbit)

#%% Run s13
if do_s13:
    inpath  = pjoin(base_input, 's13_data')
    outpath = pjoin(base_output, 's13')

    print('Starting work on SI13')
    print('Pulling data from: ' + inpath)
    print('Offlaoding at: ' + outpath)

    avail_orbit = background_removal(s13files, inpath, outpath, parallel=parallel)
    np.save(pjoin(base_output, 's13_avail_orbit.npy'), avail_orbit)
