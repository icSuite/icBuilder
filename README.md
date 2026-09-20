# IMAGE Conductance Builder

`icBuilder` is a tool for processing IMAGE WIC, SI12, and SI13 data to estimate ionospheric conductance with propagated uncertainties.

> ⚠️ `icBuilder` is **not intended** for reading the estimated conductances.  
> For that purpose, use the lightweight companion tool [`icReader`](https://github.com/icSuite/icReader).

## Project Description

This code was developed to robustly estimate ionospheric conductances from IMAGE data. All processed IMAGE data (2000-2001) is available on request (~2.5 TB). Estimated ionospheric conductances with associated uncertainties are available [**here**](https://doi.org/10.5281/zenodo.15579301).

The main purpose of this codebase is to document the data processing procedure. While not primarily designed for external use, the code can be run by others if needed.

## Project Memory

For repository-specific operating guidance and technical continuity, start with
[`AGENTS.md`](AGENTS.md) and
[`vault/START HERE - AI Onboarding.md`](vault/START%20HERE%20-%20AI%20Onboarding.md).

## Dependencies

- [`fuvpy`](https://github.com/aohma/fuvpy) – for FUV image processing **[3]**
- [`secsy`](https://github.com/klaundal/secsy) - for cubed sphere grid generation
- `ZhangPaxton2008` - for generating the bundled electron-energy lookup
- `netCDF4` - for lookup and conductance-product serialization
- `numba` - for compiled detector-footprint overlap calculations
- [`tqdm`](https://github.com/tqdm/tqdm) – for progress bars (optional; can be removed with minor edits)
- [`icReader`](https://github.com/icSuite/icReader) - required for reading
  generated pipeline products and validating detector-first restart state

## Installation

mamba activate your_environment  
git clone https://github.com/BingMM/icBuilder.git  
cd icBuilder  
pip install -e .

The project metadata installs icReader from its current `modular_pipeline`
branch. For joint development, clone `icReader` beside this repository and
install it editable before installing icBuilder. The detector-first pipeline
requires the reader interface that includes lazy detector fields, verified CS
grid reconstruction, and read-only variable metadata.

## Step-by-Step Guide

The code is designed to be run in the following sequence:

### `make_orbit_h5_files.py`

- Scans available WIC, SI12, and SI13 data.
- Assigns each file to an orbit based on `orbitsdates.csv`.

### `make_orbit_nc_files.py`

- Reads WIC, SI12, and SI13 data.
- Applies a background removal algorithm.
- Saves the results into a series of NetCDF files **[3]**.

### `make_background_removal_figures.py` *(optional)*

- Plots the data stored in the NetCDF files for quality inspection.

### `grid_resolution_determination.py` *(optional)*

- Analyzes the NetCDF files to determine the optimal grid resolution.
- Uses the trade-off between standard error and the number of bins with ≥30 measurements (lower limit of the central limit theorem).
- Suggested resolution: 225 km for WIC, 450 km for SI12 and SI13.

### `make_zhang_paxton_lookup.py`

- Builds the collapsed Zhang--Paxton electron-energy lookup on the canonical
  36-by-36 WIC Cubed-Sphere grid.
- Stores `E0` and the latitude-profile spread `dE0` for Kp 0.00 through 9.00
  in steps of 0.01.
- Uses dimensions `(kp, eta, xi)` because MLT is a two-dimensional coordinate
  on the Cubed-Sphere grid.

Generate the bundled table from the repository root:

```bash
python scripts/make_zhang_paxton_lookup.py --workers 1
```

Load one or more Kp levels:

```python
from icbuilder import load_zhang_paxton_lookup

lookup = load_zhang_paxton_lookup(1.519)  # Uses the Kp=1.52 layer.
E0 = lookup["E0"]
dE0 = lookup["dE0"]
```

Kp is rounded to the nearest hundredth; the loader does not interpolate
between table layers. Exact half values round upward. Scalar input returns
36-by-36 arrays, while an array of Kp values returns `(time, 36, 36)` arrays.
The returned dictionary also contains `E0_median`, the selected `kp`, and the
two-dimensional `mlt`, `xi`, and `eta` coordinates.

### `make_conductance_orbit_files.py`

- Ingests the NetCDF files and estimates ionospheric conductance **[1,2,4]**.
- Loads the bundled definitive GFZ Kp series once, matches each final IMAGE
  frame to its enclosing three-hour interval, and selects the corresponding
  Zhang--Paxton lookup layer.
- Uses Zhang--Paxton E0 and the collapsed latitude-profile spread dE0 in place
  of the WIC/SI13 energy inversion. SI13 remains in the stage-1 common-frame
  population and ratio diagnostics, but does not affect E0 or Fe.
- Stores the original GFZ Kp, rounded lookup Kp, Kp interval start, source
  provenance, E0--Fe covariance, and lookup/collapse provenance.

`BinnedImage` groups the source pixels belonging to each populated grid cell
once before calculating medians, standard deviations, weights, and viewing
geometry. Uncertainty multipliers are evaluated once for each distinct sample
count. When coarse SI fields are interpolated to the WIC grid, fields share a
triangulation only when their non-NaN source cells are identical; fields with
different missing-data patterns retain their own valid source cells.

The modular Product-2 builder uses SI12 for the event-specific proton flux and
Hardy et al. (1991) for proton mean energy by default. Hardy is evaluated once
per distinct Kp value on the Product-2 MLT/MLAT grid. The raw model energy is
saved as `Ep_model`; `Ep` is the value clipped to the 0.47--46.7 keV range of
the Frey camera-response tables, and `Ep_clipping_flag` records every changed
cell. `Fp`, `dFp`, and the unmodelled status of Hardy energy uncertainty are
also preserved through Product 3. A fixed-energy comparison remains available
with `--proton-energy-model constant`.

Separate zero- and nonzero-flux masks preserve the one-sided conductance
uncertainty definition without evaluating the singular derivative at zero.
The historical `image_ratio` comparison remains scalar because its piecewise
low-signal rules are not used by the Zhang--Paxton path.

Input and output folder names default to `wic`, `s12`, `s13`, and
`conductance` under `--base`, but can be changed independently. For example,
the Chapman server directories are selected with:

```bash
python scripts/make_conductance_orbit_files.py --base ~/IMAGE_FUV \
    --wic-folder =wic --s12-folder =s12 --s13-folder =s13 \
    --output-folder =conductance --parallel False
```

The bundled `icbuilder/data/gfz_kp_2000_2003.json` contains definitive
(`status=def`) Kp from 2000-01-01 through 2003-07-31, distributed under CC BY
4.0 with DOI `10.5880/Kp.0001`. See `icbuilder/data/README.md` for the exact
query and acquisition record. The loader checks the series and records the
SHA-256 of the file actually used in each output product.

The pipeline interprets the orbit files' timezone-free frame datetimes as UTC.
GFZ timestamps mark interval starts, so a frame at exactly 03:00 uses the
03:00--06:00 Kp value. Time interpolation, nearest matching, gap filling, and
out-of-range clipping are not used. `ConductanceImage` checks again that every
frame is inside the Kp interval serialized with it.

Orbit processing resumes by default. Existing `or_XXXX.nc` files are skipped
only when they open successfully and contain the expected Zhang--Paxton
variables and dimensions. New products are written to a temporary file in the
output directory and moved to their final name only after validation, so a
crash cannot make an incomplete file look finished. Use `--overwrite` to
deliberately recompute all common orbits. Do not run two pipeline instances
against the same output directory at the same time.

If proton subtraction clips electron energy flux to zero, the derivative of
the Robinson conductance with respect to flux is singular. In that case dP
and dH are reported as the one-sided conductance excursions from `Fe=0` to
`Fe=dFe`, rather than using linear uncertainty propagation.

The scientific compatibility of Zhang--Paxton electron mean energy with the
energy quantity assumed by the WIC response and Robinson relations remains an
unresolved publication gate.

### `make_conductance_figures.py`

- *(Not implemented yet)*

### `make_spline_model.py`

- *(Not implemented yet)*

## References

**[1]**. Frey, H. U. et al. (2003). *Summary of Quantitative Interpretation of IMAGE Far Ultraviolet Auroral Data*. In: Burch, J. L. (Ed.), *Magnetospheric Imaging — The IMAGE Prime Mission*. Springer. https://doi.org/10.1007/978-94-010-0027-7_11  
**[2]**. Gasparini, S. et al. (2024). *A quantitative analysis of the uncertainties on reconnection electric field estimates using ionospheric measurements*. *JGR: Space Physics*, 129, e2024JA032599. https://doi.org/10.1029/2024JA032599  
**[3]**. Ohma, A. et al. (2024). *Background removal from global auroral images: Data-driven dayglow modeling*. *Earth Planet. Phys.*, 8(1), 247–257. https://doi.org/10.26464/epp2023051  
**[4]**. Robinson, R. M. et al. (1987). *On calculating ionospheric conductances from the flux and energy of precipitating electrons*. *JGR*, 92(A3), 2565–2569. https://doi.org/10.1029/JA092iA03p02565

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contact

For questions or comments, please contact [michael.madelaire@uib.no].
