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

The optional full spatial **Hardy electron** route evaluates Hardy et al.
(1987), DOI `10.1029/JA092iA11p12275`, Tables 4–5 at each pixel's MLAT,
MLT and frame Kp. It uses the ratio of statistical integral energy flux to
number flux as mean energy `E0`, then derives event electron flux `Fe` from
proton-corrected WIC. SI12 still supplies observed proton flux; the independent
`--proton-energy-model hardy|constant` selects proton energy. SI13 is diagnostic
and missing SI13 does not remove Hardy observations or their uncertainty.
Detector image-ratio and bin-first Zhang–Paxton defaults remain unchanged.

For a detector orbit, use a scratch output base and an explicit label shared
by the subsequent stages:

```bash
python scripts/pipeline/make_precipitation_detector_orbit_files.py \
    --base-input /path/to/input --base-output /tmp/hardy_orbit \
    --precipitation-method hardy --proton-energy-model hardy \
    --retrieval-label HardyElectron_HardyProton --orbit 85 --workers 1
python scripts/pipeline/make_conductance_detector_orbit_files.py \
    --base /tmp/hardy_orbit --retrieval-label HardyElectron_HardyProton \
    --orbit 85 --workers 1
python scripts/pipeline/make_detector_cs_orbit_files.py \
    --base /tmp/hardy_orbit --retrieval-label HardyElectron_HardyProton \
    --orbit 85 --workers 1
```

Without an explicit label, the Hardy Product-2 runner derives
`HardyElectron_HardyProton_smoothed_wic0p8_si120_si131p2` (constant protons
instead produce `HardyElectron_ConstantProton...`). Existing image-ratio labels
retain their `IR_...` form. Restart validation rejects changes in electron
selection, electron model provenance, proton selection and count preprocessing.
The bin-first `scripts/pipeline/make_precipitation_orbit_files.py` also accepts
`--precipitation-method hardy`; select a dedicated `--output-folder`.

Hardy latitude support is northern corrected geomagnetic latitude 50–90°.
The current Modified Apex coordinates at 130 km approximate that convention
and this approximation is recorded. Fractional Kp uses interpolation of
**evaluated log10 fluxes**, clipped to the available 0–6 categories; this is an
implementation policy. A 24-hour Fourier phase includes the normalization
missing from printed equation 5. The mandated latitude-profile log10 floors of 6 equatorward and 7 poleward
apply to both flux quantities and affect background mean energies.

Raw Hardy energy is preserved. The WIC response holds its endpoint outside
0.2–25 keV, with `E0_response_range_flag` codes −1 unavailable, 0 in range,
1 below and 2 above. Detector Product 3 retains the flag; CS products retain
area-weighted `E0_response_range_fraction` and `E0_response_range_any`.
All representations preserve `electron_*` scientific provenance.
Uncertainty is **conditional on fixed Hardy energy**: `dE0=0` and
`varE0Fe=0` exclude model uncertainty, while `dFe` propagates corrected-WIC
measurement error. Missing measurement STD preserves nominal precipitation
and conductance and leaves measurement uncertainty unavailable. Nonpositive
corrected WIC does not define a Hardy retrieval.

The full-spectrum mean energy includes low-energy electrons and extrapolation
to 100 keV. The paper warns that non-Maxwellian low-energy precipitation can
lower this mean relative to the E-region-ionizing component. This WIC/Hardy
mean-energy/Robinson combination is a selectable retrieval requiring scientific
validation; it does not reproduce Hardy's separately fitted conductance maps,
which used a different spectrum-based conductivity calculation.

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

## Corrected WIC and SI13 versus DMSP conductance

`scripts/ratio_relation_validation/conductance_validation.py` analyzes detector
crossings; `conductance_validation_cs.py` analyzes CS crossings. Both run WIC
and SI13 with the same conductance binning and uncertainty/bootstrap analysis,
and additionally estimate paired Hall(WIC, SI13) and Pedersen(WIC, SI13).
The detector entry point defaults to
`~/Dropbox/work/data/IMAGE_FUV/ratio_relation_validation/detector_crossings`.

```bash
python scripts/ratio_relation_validation/conductance_validation.py \
    --matches /path/to/detector_crossings \
    --annotations data/dmsp_frame_annotations.csv
```

Detector outputs use `figures/conductance_relation_validation/`,
`data/wic_conductance_bins.nc` and `data/si13_conductance_bins.nc`.
Detector selection uses finite nonnegative counts, own-camera CS-gate coverage,
valid pixel indices, a containing footprint and detector separation, with the
same DMSP flags/moment/SZA cuts as CS analysis. It adds no ratio, method-valid
or other-camera/pair-valid cut. The extractor's upstream common WIC/SI13 support
and CS coverage gate remain inherited. Error draws are shared by detector
pixel and IMAGE time, so different pixels within one CS cell remain distinct.

`scripts/ratio_relation_validation/conductance_validation_cs.py` reads the
same per-orbit CS crossing files and annotation CSV as the CS ratio analysis.
It produces Hall and Pedersen 2-D histograms for each camera (four relations),
with corrected WIC or SI13 counts on the x axis and inferred electron conductance
in siemens on the y axis. It overlays conductance quantiles within count bins,
with camera-error, joint DMSP conductance-error and block-bootstrap resampling
enabled by default. Each camera uses its own finite nonnegative counts and
positive coverage; WIC selection does not require a usable SI13 measurement.

```bash
python scripts/ratio_relation_validation/conductance_validation_cs.py \
    --matches /path/to/matched_dmsp_image_data_cs \
    --annotations data/dmsp_frame_annotations.csv
```

Outputs are PNG/PDF relation plots and quality diagnostics under
`figures/conductance_relation_validation_cs/`, and counts/conditional quantiles
in `data/wic_conductance_bins_cs.nc` and `data/si13_conductance_bins_cs.nc`.
SI13 figures are `si13_hall_conductance_relation.{png,pdf}` and
`si13_pedersen_conductance_relation.{png,pdf}`, plus SI13 diagnostics.
`--figures`, `--binned-output` (WIC) and `--si13-binned-output` select other
destinations. By default the SI13 NetCDF is saved beside the chosen WIC output,
replacing a leading `wic_` with `si13_` or appending `_si13` to its stem.
`--ignore-annotations` analyzes all available crossings.

Both entry points also produce `joint_{hall,pedersen}_conductance_maps.{png,pdf}`
and `joint_{hall,pedersen}_conductance_slices.{png,pdf}`. These show conditional
medians, 10th–90th percentile spread, sample/orbit occupancy, usable-error
medians, resampling interval widths and realization availability. The third
variable is conductance colour on a fixed WIC–SI13 count plane. Conductance
tails are retained; no Frey curve or ratio cut is added.

Joint defaults are 20 WIC bins × 20 SI13 bins over the camera ranges, with at
least 20 samples (`--minimum-quantile-count`) and 3 distinct original orbits
(`--joint-minimum-orbits`) per estimated cell. Occupancy is retained below
these thresholds. `--joint-wic-bins` and `--joint-si13-bins` control resolution;
`--no-joint` runs only the separate camera analyses. Joint samples are the
intersection of the existing camera selections and inherit their upstream
support. The new output is `data/joint_conductance_bins.nc` for detector or
`data/joint_conductance_bins_cs.nc` for CS. `--joint-binned-output` overrides
it; a custom WIC filename gets `wic_` replaced with `joint_` or `_joint` appended.

Joint realizations perturb both corrected cameras, preserve shared pixel/time
and DMSP-time draws and P/H covariance, and use the chosen block bootstrap.
Cross-camera error covariance is **assumed zero**, since corrected-camera
covariance from shared proton/background/calibration terms is unavailable.
Nominal sample spread and uncertainty in conditional quantiles are separate
outputs. Both camera errors must be usable for the resampled cohort; its
unperturbed quantiles and occupancy are saved separately. If none are usable,
nominal maps remain and the resampling status is explicitly unavailable.
Bootstrap orbit copies never count as additional distinct supporting orbits.

`joint_prediction_comparison.csv` and the joint NetCDF compare WIC-only,
SI13-only and joint fixed-bin median predictions on entirely held-out orbits
(default five folds, `--validation-folds`). All train on the same paired
population and coarse camera bins with the same sample/orbit thresholds.
Scores use identical test rows where every model has training support; no
extrapolation. Shared held-out DMSP records are also removed from training.
`joint_prediction_comparison.{png,pdf}` displays held-out error bars for each
model, with joint-versus-WIC percentage changes written on the panels. Lower
error is better. CSV percentage changes are negative for improvements.
Pooled RMSE/MAE/bias, individual-orbit MAE, equal-orbit mean MAE, prediction
availability and common sample/orbit counts are saved. Fewer than two in-range
orbits gives an explicit unavailable result; supported test coverage should
be reviewed alongside the scores. This nominal comparison does not add
resampled performance intervals or fit a smooth production relation.

To generate the error comparison from an existing joint NetCDF without loading
crossings or repeating resampling (also supported by the CS entry point):

```bash
python scripts/ratio_relation_validation/conductance_validation.py \
    --comparison-only data/joint_conductance_bins.nc \
    --figures /path/to/comparison_figures
```

The default 200 realizations perturb each camera's corrected counts using
`img_wic_std` or `img_si13_std`, draw
P/H jointly using their saved standard deviations and covariance, and
bootstrap whole IMAGE orbits. A shared IMAGE cell (or detector pixel) gets one error draw even
when several DMSP samples match it. Gaussian draws are clipped at zero;
this changes their mean near zero and does not recover the censored signal.
Plots show a 90% resampling interval for the conditional median. NetCDF saves
every realization, all five quantile intervals, and the fraction of usable
realizations per bin; sparse intervals use only populated realizations.
`--realizations`, `--seed`, `--bootstrap-unit none|orbit|crossing`,
`--confidence-level` and `--no-progress` control this calculation.
`crossing` means one orbit/frame/satellite block. `--realizations 0` retains
the descriptive analysis alone. Missing/negative camera errors are excluded from
resampling and counted, while remaining in the descriptive histograms.
Repeated matches to one satellite/DMSP timestamp share one P/H draw, including
when bootstrap repeats an orbit. P/H draws are Gaussian and clipped at zero,
which changes moments near zero. Missing/invalid DMSP uncertainties or invalid
covariance likewise exclude rows only from resampling. Saved cohort counts,
unperturbed eligible-cohort quantiles and plot labels distinguish these rows
from the broader descriptive population. Bands cover that eligible subset.

DMSP errors come from first-order propagation of source per-channel spectral
uncertainties, assuming independent energy-channel errors with fixed MSIS/B.
Source channel errors include their reported calibration contribution;
additional correlated calibration/model errors, low-coverage representativeness,
covariance between different IMAGE cells/frames or DMSP times and common
systematics remain unmodelled. These are not calibrated total-error bounds.

**Missing channel errors:** propagation follows [Redmon et al. (2017)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JA023339):
NaN channel-error terms are omitted from variance/covariance sums. At least
one finite nonnegative error is required; all-missing errors remain unavailable.
`electron_conductance_uncertainty_channel_count` records contributing channel
errors, including defined zeros. Negative/infinite supplied errors and singular
propagation remain invalid. Nominal conductances use the complete measured
spectrum. This convention omits the unknown variance of NaN-error channels.
On the Jan1 F13 check, 19,696 of 20,316 valid spectra have usable uncertainty;
all 870 quality-selected pairs in the CS orbit-85 check can be resampled.

After deploying the updated icPhysics and icBuilder, regenerate yearly DMSP
files with `process_dmsp_yearly.py --overwrite`, then rebuild CS/detector
crossings with their fetch scripts and `--overwrite`. Existing files lack the
new error fields; annotations remain reusable. Rerun this analysis afterward.

Initial display ranges are 0–3000 WIC counts (60 bins), 0–60 SI13 counts
(60 bins), and 0–50 S (100 bins), with 10/25/50/75/90% quantiles requiring at
least 20 samples per camera-count bin.
`--wic-max`, `--si13-max`, `--conductance-max`, `--wic-bins`, `--si13-bins`, `--conductance-bins` and
`--minimum-quantile-count` configure these. Quantiles retain conductance tails
above the display range; excluded histogram counts are reported and saved.

Selection requires finite nonnegative counts/P/H, positive camera coverage, both
DMSP validity flags, and an inside-grid CS match within 1°. Provisional
moment cuts match the ratio analysis: fractional energy/flux uncertainty
≤0.25/0.20, flux 0–1e14, SZA 0–180°. These limits are exposed as command-line
options; see `--help`. WIC selection adds no SI13 or ratio-success requirement;
SI13 selection uses SI13 counts/coverage and the same matched WIC-grid SZA.
**Existing CS WIC, coverage and SZA already inherit ratio-retrieval support
upstream.** Those inputs also inherit common SI13 support from detector proton
correction. The script records this selection limitation; broader WIC support
requires upstream product changes and regeneration.

## DMSP spectra and yearly conductances

`scripts/data/process_dmsp_yearly.py` now writes
`electron_pedersen_conductance`, `electron_hall_conductance` (S),
`electron_conductance_valid`, `electron_raw_counts_valid`,
`reference_magnetic_field` (nT), and the
reference footprint's geodetic latitude/longitude alongside the existing
moments. Install the updated local icPhysics checkout exposing
`conductance_from_spectrum`, plus the optional dependencies with
`python -m pip install '.[dmsp]'`.

The calculation uses full 19-channel electron spectra. It evaluates the
installed ppigrf IGRF model at each Modified-Apex footprint at
`--reference-height` (default 130 km), with coefficients evaluated at the
first source record of each day. The field is held constant through the
conductivity altitude integration, matching the IDL's single-B input.

`electron_raw_counts_valid` is 1 only when all electron count channels are
finite, nonnegative, not equal to the source fill value and within the source
`VALIDMIN`/`VALIDMAX` limits. A failed record retains its moments but receives
NaN Hall/Pedersen and `electron_conductance_valid=0`. This count check is
separate from spectrum completeness and uncertainty selection. Calibrated
flux can remain finite when raw counts are out of range. Existing yearly
files without this flag/policy require rebuilding with `--overwrite` after
ongoing runs finish, or using a separate `--output` directory.

By default, the processor groups valid measurements within each satellite-day
into **2° geodetic latitude × 5° geodetic longitude × 1-hour UTC bins**.
It evaluates NRLMSISE-00 once per occupied bin, at the geographic cell center
and median measurement UTC time. Every spectrum retains its own IGRF field
and conductance calculation; spectra and conductances are not averaged.
Daily indices load once per satellite-year. MSIS evaluates up to 256 occupied
bins at a time, and each bin reuses the shared-profile deposition kernel.

Widths are configurable with `--msis-latitude-bin`, `--msis-longitude-bin`,
and `--msis-time-bin-hours`; `--msis-per-measurement` restores the exact
time/location calculation for sensitivity comparisons. Latitude/longitude
bins begin at -90°/-180°, longitude wraps at the dateline, and time bins restart
at UTC midnight. Widths must divide 180°, 360°, and 24 hours respectively.
`msis_profile_geodetic_lat`, `msis_profile_geodetic_lon`, and `msis_profile_time`
record the actual shared atmosphere used by each valid sample. Invalid spectra
have missing profile coordinates/time and NaN conductance.

On F13 2000-01-01, 868 occupied bins replace 20,316 individual profiles and
initial full daily processing took 6.3 s versus about 49 s per measurement.
For reference conductance ≥0.1 S, the 95th-percentile absolute relative
differences are 0.343% (P) and 0.169% (H); maximum absolute differences over
all valid samples are 0.0705 S and 0.0647 S.

The writer now converts integral CDF_EPOCH milliseconds directly, vectorizes
ApexPy's reference-height algebra, and calls the same MSIS C backend without
`numpy.vectorize` overhead. A serial comparison gave 7.64–7.99 s/day before
and 4.71–4.76 s/day after, with every output scientific array exactly equal.
A complete one-day read/calculate/write took 5.19 s; timings depend on machine
load and are not a measured full-year runtime. The CDF fast path retains
cdflib conversion for fractional milliseconds, TT2000, EPOCH16 and dates
outside datetime64[ns] range.

Two further 2,000-sample comparisons on 2000-07-01 and 2000-07-15 (Ap=7/164)
included every occupied bin, bin edges and strong precipitation. Their
95th-percentile binning differences were below 0.4% in both P/H for exact
conductance ≥0.1 S. These are challenge samples, not unbiased full-day
percentiles or proof of all-season model accuracy.

`icbuilder/data/msis_space_weather_2000_2003.nc` contains CelesTrak historical
indices covering 2000–2003, plus 1999-12-31 for the F10.7 lag. This 242-kB file
contains daily observed F10.7, its centered 81-day average, daily Ap, and eight
three-hour ap values per day. The loaded arrays occupy about 152 kB, plus small
string/metadata overhead. Source qualifiers are retained, including the source's
interpolated F10.7 on 2002-08-30; no new interpolation is performed.

MSIS receives previous-day **observed** F10.7, current-day **observed centered**
81-day F10.7, and current-day daily Ap. Adjusted-to-1-AU flux and trailing averages
are not substituted. The optional three-hour geomagnetic-history mode is not
enabled. The yearly file records `msis_f107`, `msis_f107a`, and `msis_ap`, source
and model provenance. Missing dates, forecasts and invalid values raise errors;
selected source dates are checked before processing writes begin. No downloads
occur during processing.

Run a selected year using the bundled indices:

```bash
python scripts/data/process_dmsp_yearly.py \
    --input /path/to/dmsp_ssj --output /path/to/new_yearly_products \
    --satellite f13 --year 2000
```

For a different interval, download the public CSV once and convert a local slice:

```bash
curl --fail --location https://celestrak.org/SpaceData/SW-All.csv --output /tmp/SW-All.csv
python scripts/data/msis_indices.py --input /tmp/SW-All.csv \
    --output /path/to/indices.nc --start 2000-01-01 --end 2003-12-31
python scripts/data/process_dmsp_yearly.py \
    --satellite f13 --year 2000 --indices /path/to/indices.nc
```

The converter automatically retains the preceding day; source URL/hash,
acquisition-time basis, selected interval and original OBS/INT flags are recorded.
The default bundled interval cannot serve measurements outside 2000–2003.
Public source and conventions:
[CelesTrak](https://celestrak.org/SpaceData/),
[NRLMSISE-00 API](https://pynrlmsise00.readthedocs.io/en/latest/reference/nrlmsise00.html).

An explicit fixed atmosphere remains available with `--atmosphere`, mutually
exclusive with `--indices`. `scripts/data/dmsp_conductance.py` generates it from
supplied conditions. The original `msis120065.xdr` MSIS-90 profile's conditions
remain unknown; using the public successor model does not establish equivalence.

For a **demonstration profile only**, using deliberately specified reference
inputs rather than recovered historical conditions:

```bash
python scripts/data/dmsp_conductance.py --output /tmp/dmsp_reference.nc \
    --date 2000-01-01T12:00:00 --latitude 65 --longitude 0 \
    --f107a 150 --f107 150 --ap 4
```

To use that fixed profile instead of the binned MSIS workflow:

```bash
python scripts/data/process_dmsp_yearly.py \
    --input /path/to/dmsp_ssj --output /path/to/new_yearly_products \
    --satellite f13 --year 2000 --atmosphere /path/to/reference.nc
```

Files record the profile or indices inputs/hash, physics and processing source hashes,
IGRF coefficient identity, and inherited angular/bin/integration assumptions.
Existing moment-only or configuration-mismatched outputs require explicit
`--overwrite` or a new output directory; current files restart without writes.
Missing/negative spectrum channels invalidate the entire conductance estimate;
zero complete spectra give zero conductance. These are inferred electron-only
conductances, without solar or precipitating-ion contributions. Yearly files
also contain `electron_pedersen_conductance_std`, `electron_hall_conductance_std`,
`electron_pedersen_hall_conductance_covariance` and
`electron_conductance_uncertainty_valid` (int8) and
`electron_conductance_uncertainty_channel_count` (int16), using the documented channel-error
propagation above. Missing errors preserve nominal conductances. Original-
atmosphere/IDL-runtime parity is unverified.

## Single-orbit conductance spline test

`scripts/debugging/fit_conductance_cs_splines.py` loads one `conductance_cs`
product through icReader and fits independent Hall/Pedersen tensor-product
B-splines over CS xi, eta and elapsed UTC seconds. Editable defaults use
orbit150, 16 spatial control points per axis, cubic spatial/quadratic temporal
splines, 5-minute interior time-knot spacing and log10 ridge strengths -6.
The 16-control default can suppress a thin auroral oval under uncertainty
weighting: low-STD background dominates a coarse basis. The orbit-150 audit
recovered the tested oval with `--ncp 30`, the legacy spatial default.
Inspect residuals and resolution sensitivity before interpreting the fit or
its conditional posterior STD. `--weighting uniform` is a useful diagnostic.

Selected frames without usable real Hall or Pedersen observations are removed
before spline construction and corner padding. Actual UTC timestamps are
preserved across gaps, and `source_frame_index` records the original orbit
indices. With uncertainty weighting, usable observations require valid positive
errors; uniform weighting only requires nominally valid conductance. Each
component receives corner constraints only on frames with real observations of
that component. At least two usable frames must remain.

```bash
python scripts/debugging/fit_conductance_cs_splines.py \
    --input /path/to/conductance_cs/or_0150.nc \
    --output /path/to/spline_test/or_0150_spline.nc
```

The output stores knots, coefficients, evaluated fits, observations, residuals,
fit-support masks, posterior model variance and STD, and parameters. Separate
Hall and Pedersen PNG sets contain one 2×2 figure every 30 seconds: observed/fit
conductance in the top row and observed/posterior STD in the bottom row.
Observed and fitted conductance share their color scale. Posterior STD has its
own scale, using its finite nonnegative 99th percentile independently of the
observed STD scale. There are no colorbars; all scale limits are written above
the panels. Names end in `_hall_0000.png` or
`_pedersen_0000.png`, with indices on the output-cadence grid. The spline and
posterior STD are evaluated at each output time; observations and their STD
stay at the latest measurement until the next measurement arrives. Titles
show both timestamps. `--plot-cadence` changes the cadence in seconds.

An optional clean PNG sequence evaluates a selected interval at a separate
cadence, using the same fitted coefficients and posterior. Add these options
to the usual command to keep coarse plots and also generate movie frames:

```bash
    --plot-cadence 120 \
    --movie true \
    --movie-start-frame 50 \
    --movie-stop-frame 100 \
    --movie-step 5 \
    --movie-output /path/to/movie_frames
```

Movie bounds use original zero-based input frame numbers (start inclusive,
stop exclusive), including frames subsequently removed from fitting. The
evaluation grid runs from the first to the last selected source timestamp,
anchored at the first, and must remain within fitted time support. Defaults
use the first and last retained fitting frames. `--movie-step` is seconds;
it changes neither the fit nor the coarse output cadence. No MP4 is generated.
PNG files are `hall/hall_0000.png` and `pedersen/pedersen_0000.png` beneath
`--movie-output` (default `<model_stem>_movie` beside the NetCDF).
Each 2×2 figure has only `<type> : <evaluation time>` above the panels,
with titles “Observed conductance”, “Modeled conductance”, “Observed
uncertainty” and “Model variance”. All axes and colorbars are hidden. Colour
limits stay fixed over each component's movie interval. Observed uncertainty
is supplied STD in S; the panel titled “Model variance” keeps the existing
posterior STD display in S, with its own scale. Only the panel title changes. Observations/errors are held at the latest retained
measurement. The NetCDF additionally saves `movie_time`, held observation
indices and H/P movie fits/variance; the existing plot fields are preserved.

`--ncp`, `--k`, `--kt`, `--cpt-step`, `--lH` and
`--lP` adjust the model; `--start-frame`/`--stop-frame` restrict the test and
`--plot-frame` plots only one output-cadence frame; otherwise all output frames
between the first and last selected measurement are plotted. The cadence
starts at the first measurement and does not extrapolate past the last.
`--weighting uncertainty` (default) uses each component's own positive STD
and nominal/uncertainty validity flags. Real observation STD is lower-clipped
to **1 S** by default: `--minimum-std 0.5` selects another floor, and
`--minimum-std 0` restores unfloored uncertainty weighting. This is a hard
minimum (`max(STD, floor)`), not a quadrature model-error term or an upper cap.
Invalid/zero errors remain excluded. Raw `dH`/`dP`, including their figure
panels, remain unchanged; `H/P_fit_std` records the effective fitting errors.
The floor and per-component affected counts are stored in the output metadata.
Artificial corner STD remains controlled separately by `--corner-std`.
`--weighting uniform` ignores the floor and fits all nominally
valid finite nonnegative conductances equally. The legacy median-diagonal-scaled
coefficient ridge is retained. `--wscaling true` is enabled by default, matching
the legacy coefficient-quality scaling formula; `--wscaling false` or
`--no-wscaling` uses identity regularization. Quality scales the coefficient
penalty, rather than the observational error weights. `--lH`/`--lP` set its
log10 strength independently. The current code complements the high-is-auroral
input as `1 - method_quality_weight` before averaging coefficient penalties;
the diagnostic quality variable currently stores this complemented field.
Time knots
have exact requested spacing, rather than the legacy approximate control-count
rule. Fits span gaps and are unconstrained; support masks distinguish data from
extrapolation. This diagnostic output is readable with xarray and is separate
from the legacy icReader spline-file format.

Permanently unsupported corner regions are padded with zero conductance and
0.1 S STD by default (`--corner-std` changes that value; `--no-corner-padding`
disables padding). These are four-connected missing regions touching grid
corners, identified from nominal finite nonnegative H/P support over the whole
input orbit, even when fitting a short interval. Interior holes and temporary
gaps remain missing. Original observations and `dH`/`dP` are preserved;
`H/P_fit_target`, `H/P_fit_std`, `corner_padding_mask` and the observation/fit
masks identify the artificial constraints. Residual RMS excludes padding.
The plots include the zero corners and their assigned STD. Input
`method_quality_weight`, when available, is saved with padded cells set to
0.5; it scales regularization when `wscaling` is enabled. This option requires
the input quality field. Missing quality values use 1, as in the legacy code.
The normalized coefficient penalty is saved as `regularization_diagonal`.

`H_model_variance`/`P_model_variance` (S²) and `H_model_std`/`P_model_std` (S)
describe the fitted field's Gaussian posterior conditional on independent
observation errors, the zero corner constraints and the coefficient ridge prior.
They exclude new-observation noise and model discrepancy. Uniform weighting
assumes a 1 S observation STD, except for the assigned corner STD.
By default, 512 Gaussian coefficient draws estimate variance in batches of 16;
the variance estimate has about 6.25% relative sampling STD at each point.
`--variance-samples` controls that precision and `--variance-seed` fixes the
random draws. `--variance-method exact` uses an exact solve per grid point for
small comparisons; it can be substantially slower for a full orbit.
The NetCDF retains predictions/variance at measurement times and also saves
`H/P_plot_fit`, `H/P_plot_model_variance` and `H/P_plot_model_std` on `plot_time`.
`plot_observation_index` identifies the held measurement at each plot time.

## References

**[1]**. Frey, H. U. et al. (2003). *Summary of Quantitative Interpretation of IMAGE Far Ultraviolet Auroral Data*. In: Burch, J. L. (Ed.), *Magnetospheric Imaging — The IMAGE Prime Mission*. Springer. https://doi.org/10.1007/978-94-010-0027-7_11  
**[2]**. Gasparini, S. et al. (2024). *A quantitative analysis of the uncertainties on reconnection electric field estimates using ionospheric measurements*. *JGR: Space Physics*, 129, e2024JA032599. https://doi.org/10.1029/2024JA032599  
**[3]**. Ohma, A. et al. (2024). *Background removal from global auroral images: Data-driven dayglow modeling*. *Earth Planet. Phys.*, 8(1), 247–257. https://doi.org/10.26464/epp2023051  
**[4]**. Robinson, R. M. et al. (1987). *On calculating ionospheric conductances from the flux and energy of precipitating electrons*. *JGR*, 92(A3), 2565–2569. https://doi.org/10.1029/JA092iA03p02565

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contact

For questions or comments, please contact [michael.madelaire@uib.no].
