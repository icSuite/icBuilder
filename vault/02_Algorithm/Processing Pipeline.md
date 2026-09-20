# Processing Pipeline

Last reviewed: 2026-09-20

> [!NOTE]
> This note describes the live bin-first production path and the implemented
> detector-first Product-1/Product-2 experimental slice. The accepted target
> design is documented in [[Detector-First Product Architecture]].

This note records the durable high-level workflow visible in the README and
live code. It is orientation, not a claim that the full production dataset was
reproduced during the review.

## Detector-first experimental workflow

All reads of generated detector-first and CS products cross the public
icReader boundary. icReader validates the serialized contract, decodes CF
times and fill values, lazily exposes detector cubes, and reconstructs the CS
grid. icBuilder retains the scientific calculations, workflow-specific source
and configuration checks, and every writer. Complete-orbit calculations
materialize each required field once inside the reader context; restart checks
inspect structure and metadata without decompressing the detector cubes.

The new-fuvpy sensor orbits feed
`scripts/pipeline/make_fuv_detector_orbit_files.py`. Schema-2 Product 1 keeps
every WIC frame, matches SI12 and SI13 independently within two seconds, maps
SI footprints onto WIC pixels, copies WIC detector variance, and propagates SI
variance with squared overlap weights. It writes beneath
`fuv_detector/fuvpy_bs_directional_v1` using atomic partial-file publication.

`scripts/pipeline/make_precipitation_detector_orbit_files.py` then consumes
that schema directly. The implemented image-ratio parity branch uses the
square roots of Product-1 variances as count uncertainties, preserves source
frame quality and support, and writes schema-3 `precipitation_detector` files.
The runner defaults to Hardy proton energy and the fixed image-ratio method.
It accepts separate `--base-input` and `--base-output` roots, with the output
defaulting to the input, and processes independent orbits concurrently when
`--workers` is greater than one.
Product 2 explicitly selects the icPhysics `measurement` uncertainty mode
because Product 1 already stores Poisson/compound-Poisson detector variance.
No second image-signal variance is added. Schema 3 records this contract and
invalidates the earlier duplicate-noise outputs. Induced target covariance and
background-model uncertainty are not represented.
This branch is implemented and tested. A complete 187-frame Product-1 orbit
has run in 65.46 seconds with 3.57-GB peak RSS after filtering unavailable SI
coordinates before inverse WIC interpolation. The full-corpus scientific gate
remains pending. The user confirmed that the Product-1 runner subsequently
completed successfully on the UiB `dynamit` server; aggregate product counts
and validation statistics were not yet recorded locally.

`scripts/pipeline/make_conductance_detector_orbit_files.py` consumes one
retrieval directory beneath `precipitation_detector` and writes
`conductance_detector/<retrieval_label>/robinson`. The schema-2 detector
product applies `icphysics.robinson_conductance()` directly to Product-2 E0,
Fe, uncertainty, and covariance without spatial binning or recalculating
precipitation. It preserves the per-frame WIC geometry, precipitation/proton
state, method quality, and source identity. `conductance_valid` describes
finite central Hall/Pedersen values; `conductance_uncertainty_valid` is
separate so unavailable analytic uncertainty does not remove a valid central
estimate. The runner supports separate input/output bases, atomic restart, and
orbit-level `--workers` parallelism. Cubed-Sphere conductance remains a later
product derived by binning this detector conductance.
The schema bump invalidates conductance generated from pre-fix Product 2;
regenerate it after the corrected Product 2. Central P/H remain unchanged.

`scripts/pipeline/make_detector_cs_orbit_files.py` is the shared post-processing
stage for completed detector Products 2 and 3. It freezes the first analysis
grid at 46 by 46 under `image_apex_130km_46x46_v1`, builds one WIC-footprint
mapping per frame, and reuses it to write independent `precipitation_cs` and
`conductance_cs` files. Conductance is averaged from detector P/H and is never
recalculated from binned E0/Fe. Measurement variance and E0--Fe covariance use
squared normalized overlap weights; spatial coverage, contributing footprints,
quality, clipping fraction, and uncertainty validity remain separate fields.
The runner supports separate bases, atomic restart, and orbit-level workers.

The reducer reads each compressed detector variable once per orbit after
building the frame mappings; the rejected frame-wise reader produced identical
science but was about 29 times slower on orbit 0086. Corrected orbits 0085,
0086, 0261, and 0968 completed both 46-by-46 products in 39.60 seconds with two
workers and 1,614,780-KB reported peak RSS. All 1,477,851 central-valid cells
retain propagated uncertainty, and restart skips all four in 0.33 seconds. The
legacy 36-by-36 Zhang--Paxton lookup is not used by this stage and remains
deferred.

## Modular Product-1 binning

`scripts/make_binned_orbit_files.py` defaults to
`--binning-method footprint`. For each frame, `BinnedImage` projects the
geolocated detector-pixel centres into the target CS plane, infers a local
four-corner footprint from row and column neighbour spacing, and clips that
uniform top-hat footprint against every intersected target cell. A sparse
overlap matrix is then reused for signal, weight, SZA, DZA, and LOS geometry.

The stored signal is an overlap-area-weighted mean. `counts` is the integer
number of valid source footprints intersecting a cell, while `coverage` is the
fraction of target-cell area covered by valid footprints. These are separate:
a footprint can observe a cell even when its centre lies elsewhere.

`sigma` is deliberately provisional. It is the overlap-weighted within-cell
spread followed, when requested, by the existing count-based small-sample
inflation. It is not yet a propagated detector measurement uncertainty.

The previous centre-assignment median calculation remains available through
`--binning-method centre`. The chosen method is saved in each Product-1 file
and included in restart validation, making rollback possible without mixing
incompatible files. Common WIC/SI resolution matching is not part of Product
1 and remains a Product-2 task.

## Primary orbit workflow

1. **Orbit indexing**
   `scripts/make_orbit_h5_files.py` assigns raw WIC, SI12, and SI13 filenames
   to orbit intervals and writes sensor-specific HDF5 indices.

2. **Sensor preprocessing**
   `scripts/pipeline/make_orbit_nc_files.py` reads the raw images, selects
   northern observations, reconstructs detector measurement variance and
   viewing geometry, applies the selected directional BS-only `fuvpy`
   background model, and writes compressed per-sensor orbit NetCDF files plus
   availability arrays. The corrected image and quality weight are `dgimg`
   and `dgweight`; the removed SH fields are not produced. `--base_input`
   identifies the orbit-index and raw-data tree. `--base_output` identifies the
   product tree and defaults to the input base when omitted. With neither path
   supplied, both resolve to `<repository>/example_data`.

3. **Optional inspection and grid choice**
   Background-removal plotting and grid-resolution scripts support inspection.
   The README records 225 km for WIC and 450 km for SI12/SI13 as suggested
   resolutions; treat those as documented configuration, not newly validated
   values.

4. **Conductance construction**
   `scripts/make_conductance_orbit_files.py` time-aligns the three instruments,
   removes empty or low-coverage frames, converts coordinates, grids the
   observations, matches each retained frame to definitive GFZ Kp, and builds
   a `ConductanceImage`.

   The fixed Kp series is loaded once before orbit multiprocessing. Matching is
   to the enclosing half-open interval `[start, start + 3 h)` after the final
   fullness filter. Exact interval boundaries select the new value. Missing,
   gapped, or out-of-range data fail; no interpolation or nearest matching is
   used. Orbit-file datetimes have no timezone field and are explicitly
   interpreted as UTC. `ConductanceImage` independently checks that each frame
   lies in the half-open interval stored with it.

   `BinnedImage` can apply the provisional
   pixel-level `cos(DZA)` normalization before calculating image statistics.
   Each binned sensor image retains median SZA, median DZA, median
   `cos(DZA)`, and whether the correction was selected. The geometry uses the
   same valid source pixels as the brightness statistic and follows SI/SI13
   target-grid interpolation and later frame selection. Median
   `cos(DZA)` is diagnostic: because correction precedes the image median, it
   is not generally an exact multiplier between corrected and uncorrected
   binned brightness. The modular Product-1 orbit builder explicitly selects
   `los_correction=False`. Independently, the upstream fuvpy read retains only
   source pixels with `DZA < 75` degrees.

   Binning calculates one flattened cell number per source pixel, stably
   groups the populated cells, and then calculates each cell's statistics in
   original source-pixel order. This replaces repeated full-image scans
   without changing count, NaN, median, standard-deviation, weight, or
   viewing-geometry semantics. Student-t and chi-square multipliers are cached
   by integer sample count. During SI-to-WIC interpolation, fields with the
   same non-NaN source mask share a triangulation; differing masks are kept in
   separate groups so valid information is not discarded merely to gain
   speed.

5. **Serialization**
   `icbuilder/conductanceimage.py` writes sensor statistics, characteristic
   energy, energy flux, Hall and Pedersen conductance, propagated
   uncertainties, weights, times, and grid metadata to per-orbit NetCDF.
   It also writes per-sensor `sza`, `dza`, and `los_factor` arrays with units.
   Root attributes record each sensor's source-image correction (`SH`, `DG`,
   or raw) and whether the pixel-level LOS correction was applied.
   Each stage-1 product also stores original GFZ Kp, rounded lookup Kp, Kp
   interval start, GFZ provenance, Zhang--Paxton lookup/collapse provenance,
   dE0 interpretation, and the induced E0--Fe covariance.

   Existing structurally valid orbit products are skipped by default. A new
   orbit is written to `or_XXXX.nc.partial`, reopened for a lightweight schema
   check, and atomically renamed to `or_XXXX.nc`. Missing or invalid products
   are rerun; `--overwrite` deliberately recomputes all common orbits. The
   output file itself is the completion record, so there is no separate ledger
   that can disagree with the data.

   Product 2 uses SI12 as the event-specific proton-flux source. Proton mean
   energy now defaults to the spatial Hardy et al. (1991) model, evaluated once
   for each distinct Kp value on the Product-2 MLT/MLAT grid. `Ep_model`
   preserves the raw Hardy result; `Ep` is clipped to the 0.47--46.7 keV domain
   tabulated for the Frey camera responses, and `Ep_clipping_flag` records the
   affected cells. The file also stores dEp, Fp, dFp, the energy model, the
   flux source, the Modified-Apex coordinate approximation, and the fact that
   Hardy supplies no formal dEp. A constant-Ep comparison path remains.

   Product 3 carries these proton fields and provenance without recalculating
   them. The production Zhang--Paxton branch then evaluates all supported cells
   with masked float64 arrays. Cells missing any count or count uncertainty
   retain the same NaN support as the scalar calculation. The zero-flux
   Robinson uncertainty is evaluated on a separate mask, so the singular
   nonzero-flux derivative is never evaluated at Fe=0. The old `image_ratio`
   branch deliberately remains a scalar comparison path.

## Secondary spline workflow

`scripts/make_spline_model_files.py` reads conductance products through
`icReader`, constructs `icbuilder.splineimage.SplineImage`, and writes spline
and factor NetCDF products. Historical parameter exploration and resolution
discussion remain in the dated vault notes.

Before trusting a historical parameter or resolution claim, locate its input,
configuration, commit, and generated result. The 2026-07-26 vault migration did
not reproduce those analyses.

## Execution boundary

The processing scripts write or overwrite data products. Use an isolated
copied base directory, ensure output directories exist, start with serial
execution on a representative orbit, and verify dependencies before scaling
up. Do not use the tracked `example_data/` tree as scratch space.

## Exploratory Zhang–Paxton latitude collapse

`scripts/ZhangPaxton2008_collapse.py` reduces the external Zhang–Paxton model
from `ZP(Kp, MLT, MLAT)` to one electron mean energy for each `(Kp, MLT)`.
This is not yet part of the production orbit workflow.

For each fixed Kp and MLT, it evaluates E0 and Q at the centres of 0.01-degree
cells covering 50–90 degrees northern MLAT. It identifies the principal Q
maximum and retains only the contiguous cells around that maximum where Q is
above the selected threshold. The area-weighted mean and median E0 are then
calculated conditionally over those cells with exact spherical-area factors
`sin(latitude_upper) - sin(latitude_lower)`. Consequently, the large
near-zero background outside the selected oval cannot dominate the result,
and poleward cells do not receive the same weight as larger equatorward
cells.

The provisional default is `Q > 0.05 mW m-2`, matching the Figure 8
global-mean inclusion criterion in Zhang and Paxton (2008). The only
sensitivity definition retained is the absolute `Q > 0.25 mW m-2` criterion
used as an auroral-boundary contour in the paper. The user rejected the
relative 10%-of-peak definition because it changes with the strength of each
slice. Empty masks produce NaN rather than a fabricated fallback. Contact
with the 50-degree equatorward limit is flagged as possible truncation;
extension to 90 degrees is separately reported as contact with the physical
pole.

The function also returns weighted spread, a Q-weighted energy mean, latitude
bounds, selected area, and selection diagnostics. The area-weighted mean
remains the primary representative; the area-weighted median describes a
typical unit area and is reported alongside it. The Q-weighted mean emphasizes
the precipitation power and remains a sensitivity result.

The model equations are continuous in MLT because their MLT dependence is a
six-harmonic Fourier series. The collapse therefore evaluates the exact MLT
provided by the caller rather than averaging within a finite MLT sector.
Diagnostic figures use 0.05-hour sampling. Zhang and Paxton fitted the Fourier
parameters from 48 empirical sectors with 0.5-hour width; evaluating the
analytic fit more densely resolves its smooth variation but does not increase
the underlying empirical resolution.

The Epstein latitude profiles are likewise continuous. The 0.01-degree MLAT
grid is deliberately finer than the earlier 0.25-degree diagnostic grid
because a hard threshold on coarse cells made the collapsed curves visibly
jitter as boundary cells entered and left the selection. This is numerical
oversampling, not added physical resolution.

### Fixed-grid lookup

The collapse is now importable from
`icbuilder/zhang_paxton_collapse.py`. It contains only the grid, spherical
weights, contiguous-oval selection, and batched model reduction in execution
order. Results are ordinary dictionaries rather than dataclass result
objects. Diagnostic MLT sampling, plotting, figure saving, and the command
line interface live separately in `scripts/ZhangPaxton2008_collapse.py`.

`icbuilder/data/zhang_paxton_e0_lookup.nc` contains 901 direct Kp layers,
0.00 through 9.00 at 0.01 spacing, on the canonical 36-by-36 WIC
Cubed-Sphere grid. Its dimensions are `(kp, eta, xi)`. Cell-centre MLT is a
two-dimensional coordinate calculated from `(grid.lon / 15) % 24`; neither
Cubed-Sphere axis is a one-dimensional MLT axis.

The stored scientific fields are the area-weighted mean E0, latitude-profile
spread dE0, and area-weighted median E0. The file also records units, the
two-dimensional grid coordinates, threshold, latitude domain and resolution,
and ZhangPaxton2008 package version. The loader checks those coordinates
directly against the canonical IMAGE grid.

Runtime lookup rounds Kp to the nearest hundredth and performs direct indexing
without interpolation; exact half values round upward. Scalar Kp produces
36-by-36 arrays; a vector of frame Kp values produces `(time, 36, 36)` arrays.
`load_zhang_paxton_lookup(kp)` returns a plain dictionary containing the
selected Kp, E0, dE0, median E0, MLT, xi, and eta.

The bundled table uses the still-provisional `Q > 0.05 mW m-2` selection over
50--90 degrees at 0.01-degree numerical resolution. It contains no empty
selections, but 231,123 of 1,167,696 cells touch the 50-degree equatorward
sampling limit. A later change to the threshold or latitude domain requires
regeneration.

The simplified table schema was populated by copying the three scientific
arrays from the previously verified table. All values and coordinates are
exactly equal; this refactor did not repeat all 901 expensive collapses. Fresh
direct-collapse checks at Kp 0.00, 1.52, and 9.00 agree to float32 precision.

### Stage-1 E0 integration

`icbuilder/data/gfz_kp_2000_2003.json` is the official GFZ JSON response from
2000-01-01 through 2003-07-31, one month beyond the last June 2003 IMAGE
orbit. It contains 10,464 definitive three-hour values and is distributed
under CC BY 4.0 with DOI `10.5880/Kp.0001`. The query, acquisition date, and
downloaded-file checksum are recorded in `icbuilder/data/README.md`.
Production processing never contacts GFZ. The loader checks cadence, status,
and physical Kp values, then calculates the checksum of the file actually
used for product provenance.

For each orbit, `ConductanceImage` requests the lookup once with the complete
Kp vector. It validates the resulting `(time, 36, 36)` arrays and grid,
retains original GFZ thirds separately from nearest-hundredth lookup values,
and supplies each pixel's E0/dE0 to `E0_eflux_propagated`.

The override path retains SI12 proton-flux estimation, proton subtraction from
WIC and SI13, WIC-derived Fe/dFe, and R/dR as diagnostics. It does not execute
the WIC/SI13 ratio-to-E0 inversion. Consequently, finite SI13 changes cannot
change E0 or Fe. The old inversion remains accessible only through the
explicit `energy_method="image_ratio"` comparison path.

The collapsed profile spread remains dE0:

`dE0 = sqrt(sum(w * (E0_lat - E0_mean)^2) / sum(w))`.

It represents unresolved MLAT variability, not Zhang--Paxton predictive,
coefficient, or Kp uncertainty. Because `Fe = Wprime / Wm(E0)`, the
first-order induced covariance is

`cov(E0, Fe) = -Wprime * Wm'(E0) / Wm(E0)^2 * dE0^2`.

The implementation passes this covariance to both Robinson uncertainty
functions and stores it. It does not claim to complete the broader uncertainty
model: shared-channel covariance and several upstream uncertainty terms remain
unresolved.

At `Fe=0`, the Robinson conductances are proportional to `sqrt(Fe)`, so their
flux derivatives are singular and the usual linear uncertainty calculation is
undefined. The stored dP and dH instead give the one-sided conductance change
from `Fe=0` to `Fe=dFe`. E0 uncertainty contributes nothing at exactly zero
flux because both conductances are then zero for every E0.

Stage 1 deliberately retains the three-camera common-frame selection, SI13
pixel support, combined weight, and ratio diagnostics. Making SI13 optional is
a separate stage so comparison of old and new E0 does not also change frame
population.

The major scientific gate remains unresolved: Zhang--Paxton reports electron
mean energy, while the IMAGE WIC response and Robinson calculation may assume
a different characteristic-energy definition. The implemented path is ready
for controlled testing, not yet publication-ready production.

### Proposed SI13 role after E0 replacement

Do not use SI13 to invert for E0. Its cleanest remaining use is an independent
radiometric consistency diagnostic for frames where simultaneous SI13 exists:

1. calculate corrected WIC electron counts `Wprime` and corrected SI13
   electron counts `Sprime` with the existing dayglow and SI12-based proton
   corrections;
2. take E0 and dE0 from the collapsed Zhang–Paxton lookup;
3. calculate the primary flux estimate
   `Fe_WIC = Wprime / WIC_response(E0)`;
4. predict the SI13 electron signal
   `Sprime_predicted = SI13_response(E0) * Fe_WIC`;
5. retain `Sprime - Sprime_predicted`, and preferably its normalized residual,
   as a validation variable.

Equivalently, calculate
`Fe_SI13 = Sprime / SI13_response(E0)` and compare it with `Fe_WIC`. This does
not recreate the rejected ratio-to-E0 inversion: E0 remains fixed by
Zhang–Paxton, while the two cameras test whether they recover a consistent
energy flux.

The initial recommendation is to use SI13 for validation only and not require
it for product generation. If the two independently derived flux estimates
prove consistent after proton/background correction, a later joint Fe
estimator could be considered. Such a combination must account for their
shared E0 and SI12-proton uncertainties and should not use naive
inverse-variance weights as if the channels were independent.

For arbitrary Kp, evaluate the published interpolation before performing the
latitude collapse. The six empirical bins are represented at Kp centres 0.75,
2.25, 3.75, 5.25, 7, and 9. Mean energy is interpolated linearly in Kp;
energy flux is interpolated in the paper's nonlinear hemispheric-power
coordinate. Selecting the nearest empirical bin would introduce artificial
discontinuities.
