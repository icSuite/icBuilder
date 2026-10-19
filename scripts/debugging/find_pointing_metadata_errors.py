"""Print IDL filenames behind SPIN and viewing-geometry validation errors.

Read-only diagnostic: no background fitting, pointing correction or file writes.
"""

#%% Imports and the same per-frame DZA comparison used by the FUVPY loader

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import readsav

from fuvpy.src.load import _timestamp_idl
from fuvpy.src.viewing_geometry import calculate_viewing_geometry
from fuvpy.utils.sunlight import subsol


def viewing_status(record):
    coordinates = {}
    for name, fill in [('glat', -1e31), ('glon', -1e31), ('dza', -1)]:
        raw = np.asarray(record[name][0])
        coordinates[name] = np.where(raw == fill, np.nan, raw).astype(float)
    sun_latitude, sun_longitude = subsol([_timestamp_idl(record['time'][0])])
    with np.errstate(invalid='ignore', divide='ignore'):
        geometry = calculate_viewing_geometry(
            coordinates['glat'], coordinates['glon'], np.linalg.norm(np.asarray(record['o_gci'][0], dtype=float)),
            float(record['lat_ctr'][0]), float(record['lon_ctr'][0]), float(record['emis_hgt'][0]),
            sun_latitude[0], sun_longitude[0]
        )
    stored = np.isfinite(coordinates['dza']) & (coordinates['dza'] >= 0)
    reconstructed = np.isfinite(geometry['dza'])
    common = stored & reconstructed
    error = np.max(np.abs(geometry['dza'][common] - coordinates['dza'][common])) if np.any(common) else np.nan
    return int(stored.sum()), int(reconstructed.sum()), int(common.sum()), error


#%% Scan the same sorted filenames as read_idl; retain both sides of spin breaks

def inspect_sensor(base, prefix, orbit, filenames):
    sensor = {'wic': 'WIC', 's12': 'SI12', 's13': 'SI13'}[prefix]
    previous = None
    issues = 0
    for name in sorted(filenames):
        filename = base / (prefix + '_data') / str(name)
        try:
            record = readsav(str(filename))['imageinfo']
            stamp = record['time'][0]
            time = pd.Timestamp(_timestamp_idl(stamp)) + pd.Timedelta(milliseconds=int(stamp[1]) % 1000)
            spin = int(record['spin'][0])
            current = f'{filename} TIME={time.isoformat()} SPIN={spin}'
            if previous is not None and (spin <= previous[0] or time <= previous[1]):
                reasons = []
                if spin <= previous[0]:
                    reasons.append('SPIN_NOT_INCREASING')
                if time <= previous[1]:
                    reasons.append('TIME_NOT_INCREASING')
                print(f'{",".join(reasons)} orbit={orbit} sensor={sensor}\n  previous: {previous[2]}\n  current:  {current}', flush=True)
                issues += 1
            previous = spin, time, current

            stored, reconstructed, common, error = viewing_status(record)
            if common == 0:
                print(f'NO_VALID_DZA orbit={orbit} sensor={sensor}: {current}\n  valid pixels: stored={stored}, reconstructed={reconstructed}, common={common}', flush=True)
                issues += 1
            elif error > 1e-3:
                print(f'DZA_MISMATCH orbit={orbit} sensor={sensor}: {current}\n  maximum error={error:.6g} degrees', flush=True)
                issues += 1
        except Exception as error:
            print(f'READ_OR_GEOMETRY_ERROR orbit={orbit} sensor={sensor}: {filename}\n  {type(error).__name__}: {error}', flush=True)
            issues += 1
    print(f'Orbit {orbit} {sensor}: {len(filenames)} files checked, {issues} issues', flush=True)
    return issues


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, required=True, help='Base containing *files.h5 and *_data folders')
    parser.add_argument('--orbit', type=int, action='append', help='Repeat to select orbits; default: 197 and 365')
    args = parser.parse_args(argv)
    base = args.base.expanduser().resolve()
    for prefix in ['wic', 's12', 's13']:
        files = pd.read_hdf(base / (prefix + 'files.h5'), key='data')
        for orbit in args.orbit or [197, 365]:
            filenames = files.loc[files['orbit'] == orbit, 'filename'].tolist()
            inspect_sensor(base, prefix, orbit, filenames)


if __name__ == '__main__':
    main()
