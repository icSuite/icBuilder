# Handoff - Latest

Last updated: 2026-09-21
Repository snapshot: `modular_pipeline` at `23ea850`
Worktree state: regenerated example products, debugging additions, existing
vault changes, and the implemented icReader migration are uncommitted

## Latest checkpoint: detector DMSP crossing extraction

The new `scripts/ratio_relation_validation/fetch_all_dmsp_crossings.py` has
been migrated off the legacy CS product. It reads schema-3
`precipitation_detector`, uses each frame's native Modified-Apex detector
geometry, and matches F12--F15 SSJ samples within centered +/-60-s temporal
support to their nearest WIC pixel. It writes one compressed NetCDF per orbit
instead of accumulating a multi-gigabyte dataframe. Existing files are skipped
unless `--overwrite` is supplied. Its three data locations can be passed as
`--dmsp-path`, `--image-path`, and `--output-path`.

The saved flat sample contract includes absolute IMAGE and DMSP times,
satellite, source indices, detector indices, angular separation, DMSP energy
and flux with fractional uncertainties, and the relevant stored Product-2
corrected counts, ratio, energy, uncertainties, DZA, method validity, and
quality weight. Extraction deliberately saves all temporal candidates and does
not impose a detector-distance or scientific-quality cut.

The simplified orbit-0364 gate wrote 8,315 rows in about 4.4 s, with exact -60 to +60 s
offsets and no duplicate IMAGE-time/DMSP-time/satellite keys. Median, 95th,
and maximum spatial separations are 0.198, 1.431, and 11.265 degrees. The next
validation stage must therefore establish an explicit spatial-support rule
before fitting the IMAGE-ratio/DMSP-energy relation. Existing CS-era DMSP
annotations should be joined by orbit, absolute image time, and satellite;
legacy `frame_id` is only advisory.

The extractor accepts `--workers N`, with `1` as the serial default.
Multiprocessing uses `spawn`; each process opens its own IMAGE file and keeps
its own independently opened DMSP-file cache. Orbit 0364 and 0968 outputs were
byte-identical between serial and two-worker runs. On this small
unequal-duration gate, wall time fell from 13.7 to 9.7 s; larger pools may be
limited by disk bandwidth.

The initial detector-based analysis is implemented in
`scripts/ratio_relation_validation/ratio_validation.py`. It reads the new
per-orbit files, joins accepted legacy annotations using orbit + absolute
IMAGE time + satellite rather than CS frame index, and writes quality
histograms, selected-distribution histograms, a ratio-energy relation figure,
and `data/ratio_energy_bins.nc`. The latter explicitly stores a 50x40
ratio/energy count grid, conditional DMSP-energy quantiles for each ratio bin,
and conditional IMAGE-ratio quantiles for each DMSP-energy bin. The two
conditional directions are plotted separately over the same 2-D histogram,
with the Frey relation overlaid on both.

The analysis loader treats an empty `sample` dataset as a successfully
processed zero-match orbit and skips it. Required-field validation remains in
place for every non-empty crossing file.

The same validation workflow is now available for the fixed-grid CS Product 2
through `fetch_all_dmsp_crossings_cs.py` and `ratio_validation_cs.py`. CS
crossing files remain per orbit and store the fixed-grid row/column,
centre-to-footprint angular separation, and an explicit inside-grid flag. The
analysis reuses annotations by orbit, absolute IMAGE time, and satellite,
requires the footprint to be inside the CS domain, and applies a provisional
1-degree centre-separation cut. It writes separate CS figures and
`data/ratio_energy_bins_cs.nc`, including conditional quantiles in both axis
directions. Serial and two-worker extraction produced identical results for
orbits 0085 and 0086 (3,847 and 6,938 temporal matches); the annotated
two-orbit analysis completed with 258 selected pairs.

The active 0.25 DMSP-energy fractional-uncertainty, 0.20 flux fractional-
uncertainty, 1-degree detector-separation, 0--150 ratio, and 0--8-keV energy
limits are provisional analysis choices. Inspect the full-data histograms
before treating them as durable quality rules. The two-orbit unannotated gate
loaded 41,275 candidates and retained 5,474 after these cuts.

## Latest checkpoint: detector and CS Coumans diagnostics implemented

`reconstruct_coumans_figure4a_detector.py` now reproduces the Figure-4a track
geometry from detector Product 1. Unlike the older IDL script, the full orbit
product contains the exact Coumans WIC frame at 23:35:29 UTC. The script plots
background-corrected WIC counts, maps the DMSP 130-km geographic footprint to
the time-dependent detector geometry, and writes PNG/PDF output beneath the
paper-reconstruction figure tree. The centered DMSP interval is
23:34:32--23:36:29; median/maximum nearest-pixel angular separations are
0.138/0.233 degrees. The script now also renders the next three WIC frames at
23:37:32, 23:39:34, and 23:41:37 UTC. Each output independently maps the DMSP
track using that frame's detector geometry and selects its own centered
plus-or-minus-60-s support. Green and red points are explicitly labelled as
the support start and end; they are not WIC exposure endpoints.

Orbit-0364 detector precipitation is now available, and one maintained Frey
Figure-16 reconstruction has been revived as
`scripts/paper_reconstructions/reconstruct_frey_figure16_detector.py`. It uses
the 2000-10-28 11:38:24 UTC `precipitation_detector/IR_hardy` frame. Stored
Product-2 `E0` and `Fe` supply the final mean-energy and energy-flux panels.
The added orbit-0364 Product-1 file matches Product 2's recorded SHA-256 and
all frame/source indices. Its exact SI12, WIC, and SI13 counts now supply the
sensor panels. Uncorrected mean energy is calculated from Product-1 WIC/SI13
only where WIC is at least 50 and SI13 at least 3 counts. The earlier attempt
to reverse Product 2's clipped proton correction has been removed. The script
ran successfully and wrote its PNG/PDF pair to
`figures/paper_reconstructions/frey_2003`.

The detector and CS Coumans Figure-4b plots were both truncated to
23:34--23:45 UTC and reran successfully. This changes display extent only.
Both now also overlay the earlier direct-count calculation as a dashed blue
curve: `wic_corrected / si13_corrected`, then the Frey lookup, without the
Product-2 50/3 energy fallback. Detector Product 2 supplies 366 finite stored
ratios versus 519 direct ratios. CS Product 2 supplies 434 stored ratios versus
599 direct ratios; its direct curve is calculated after binning, while its
stored retrieval fields are overlap-reduced from detector Product 2.

The paper-reconstruction folder now contains
`reconstruct_coumans_figure4b_detector.py`, which accepts schema-3
`precipitation_detector` input through icReader. Because detector products do
not have a fixed CS grid, it groups DMSP samples by nearest IMAGE frame and
queries a per-frame spherical KD-tree built from WIC-pixel MLAT/MLT. Stored
Product-2 ratio, energy, and uncertainty fields are sampled at those pixels.

The additional `reconstruct_coumans_figure4b_cs.py` accepts schema-1
`precipitation_cs` and uses its reconstructed 46-by-46 grid. Both scripts now
sample the canonical stored `R`, `dR`, `E0`, and `dE0` fields directly; they
do not recalculate the retrieval from corrected WIC/SI13 counts.

Orbit 0968 completed with 599 centered samples in both scripts: 366 finite
detector ratios and 434 finite CS ratios. Stored E0 remains finite for 580 and
599 samples because the low-signal retrieval branches can define E0 when R is
not finite. Both PNG/PDF pairs rendered beneath
`figures/paper_reconstructions/coumans_2004`.
Both versions shade the native DMSP mean-energy uncertainty using
`ELE_AVG_ENERGY * ELE_AVG_ENERGY_STD`; the legend deliberately calls this an
estimated uncertainty range because the source metadata does not explicitly
identify it as a one-sigma interval.
Both versions shade the stored IMAGE `dR` and `dE0` fields. Large ratio
intervals can dominate the display, so the figure clips only the displayed
ratio band at 1.5 times the Frey table maximum and labels that choice; the
stored numerical uncertainty is not modified.

The plot exposed a Product-2 legacy boundary bug. At the orbit-0968 detector
peak, the file stores `R = 204.74 +/- 79.48` but `E0 = 25 keV` and `dE0 = 0`.
The Frey ratio table stops at 136.49. Its legacy constant upper fill makes the
finite-difference derivative zero outside the table, falsely erasing energy
uncertainty. The 15-keV figure cap is display-only. CS uncertainty reduction
itself is active and uses squared normalized overlap weights, but faithfully
propagates this erroneous zero. Product 2 needs an explicit out-of-domain
ratio policy and downstream regeneration before these uncertainties should be
treated as scientifically valid.

## Latest checkpoint: icReader read migration implemented

Generated products in the active detector-first path now enter icBuilder
through public icReader APIs. Product-1/Product-2 loaders, paired detector-CS
reduction, detector/CS restart validation, and maintained generated-product
diagnostics have been migrated. Builder-specific method, configuration, and
source-identity checks remain in icBuilder; icReader owns schema, dimensions,
required variables, fill handling, CF time decoding, and CS-grid validation.

The change preserves the optimized CS path. A focused regression test rejects
frame-indexed decompression and proves every compressed cube is materialized
at most once. Product-2 and Product-3 orbit 0085 are exactly unchanged for all
variables. All eight CS files for 0085, 0086, 0261, and 0968 are exactly equal
in variables and grids; the two-worker gate took 37.77 seconds and 1,614,272
KB peak RSS versus 39.60 seconds and 1,614,780 KB before migration. Restart
skips all four in 0.42 seconds.

icReader passes 34 tests. icBuilder passes 127 tests with only the same four
known 36-by-36 grid/Zhang--Paxton failures. Legacy binned/precipitation/
conductance restart scans remain raw because the corresponding readers are
eager. Raw fuvpy, HDF, DMSP, lookup, writer, and low-level test access also
remains intentional.

The only closeout dependency is Git identity: icReader's immutable
`variable_attrs` API is currently uncommitted. Commit/push icReader first,
replace icBuilder's temporary `modular_pipeline` dependency with that exact
commit, and then commit icBuilder. See
[[icReader Read Migration Implementation Plan]].

## Latest checkpoint: detector CS post-processing implemented

The user selected 46 by 46 for the first fixed detector-CS representation.
The implementation now validates only the 46-by-46 shape. The former exact
coordinate hash was removed from the writer, restart check, and icReader after
it blocked the server run before orbit 0085 in a different numerical
environment.

One combined post-processing runner creates `precipitation_cs` and
`conductance_cs`, building the WIC detector-footprint mapping once per frame
and applying explicit reducers for means, measurement variance, covariance,
flags, quality, counts, and coverage. Conductance is always reduced from
detector Product 3, never recalculated from binned Product 2. Separate bases,
atomic restart, orbit-level workers, source-pair validation, and schema/grid
shape checks are implemented.

Ten focused tests pass, including a real two-worker run and a nonlinear
ordering guard. The full suite has 126 passes and the same four deferred
36-by-36/Zhang--Paxton failures.

Profiling rejected the first frame-wise compressed-reader implementation:
orbit 0086 spent about 2.2 seconds per frame repeatedly decompressing fields,
while its footprint geometry took about 0.02 seconds. The final reducer builds
all mappings, reads each variable once, and produces files identical to the
original calculation. The complete corrected four-orbit gate finished in
39.60 seconds with two workers and 1,614,780-KB reported peak RSS. All
1,477,851 central-valid cells retain propagated uncertainty. The outputs total
90.84 MB for Product-2 CS and 52.85 MB for Product-3 CS, and all four skip on
restart in 0.33 seconds. See [[Detector CS Post-processing Implementation Plan]].

The legacy 36-by-36 Zhang--Paxton lookup and old bin-first conductance
infrastructure need a later redesign. That work is explicitly deferred and
does not block the image-ratio Product-2 CS and Robinson Product-3 CS stage.
The local implementation, numerical, restart, multiprocessing, and four-orbit
performance gates are complete; server execution is the next operational step.

The actual completed tree at
`/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test` was subsequently
audited. Product 3 reproduces Robinson from detector Product 2 to float32
precision, and the paired CS files have consistent masks, support, provenance,
and spatial structure. Their seven empty frames all have zero detector
Product-2 support because an SI camera is missing; this is not a CS mapping
loss.

Do not yet equate this production success with science acceptance. Hardy mean
energy is clipped below 0.47 keV for 52.1--66.2 percent of detector-valid
pixels in the four orbits, while Hardy uncertainty is unmodelled. `varE0Fe` is
zero everywhere, and 78--97 percent of central-valid CS cells have conductance
uncertainty larger than the central P or H estimate. The next decision is how
clipped-energy and low-signal retrievals should affect validity or downstream
selection.

## Latest checkpoint: Product-2 uncertainty propagation fixed

Measurement variance is carried from Product 1 into Product 2. fuvpy
`img_variance` already contains the SI Poisson or WIC compound-Poisson
counting variance. The pre-fix Product 2 supplied its square root as `dwic`,
`dsi12`, and `dsi13`, after which the extracted legacy
`icphysics.proton_correct_images()` function added WIC/SI signal again as a
Poisson variance term. Since Product-1 `dgimg` is background-subtracted and
can be negative, the duplicate term also created invalid square roots and
unnecessarily missing uncertainties.

On full orbit 0085, the corrected interpretation gives finite `dFp`, corrected
WIC uncertainty, and corrected SI13 uncertainty on every method-valid cell;
the pre-fix files give 96.0%, 96.0%, and 91.9%. For cells finite in both, the
pre-fix uncertainties are inflated by median factors 1.072, 1.001, and 1.075
and 95th-percentile factors 1.336, 1.068, and 1.467. Central Product-2 `Fp`,
`E0`, and `Fe` and central Product-3 `P` and `H` do not depend on this error.
The affected fields are `dFp`, corrected-count uncertainties, `dR`, `dE0`,
`dFe`, `dP`, `dH`, and their validity masks.

The shared function now exposes a backward-compatible `legacy_add_poisson`
default and an explicit `measurement` mode. Product 2 selects `measurement`,
records it as provenance, and no longer adds the duplicate term. Product 2 is
schema 3 and Product 3 is schema 2, so normal restart checks invalidate the
old files. Product 1 remains usable; Products 2 and 3 need regeneration.

The full orbit-0085 regression leaves central Fp/E0/Fe and P/H bit-for-bit
unchanged. All 2,156,871 method-valid cells now have finite dE0/dFe/dP/dH,
versus 1,982,337 before. All 17 icPhysics tests and 19 focused Product-2/3
tests pass. The full icBuilder suite has 116 passes and the same four unrelated
grid/lookup failures.

This does not solve the separately documented omissions of
background-model discrepancy, calibration/response-model uncertainty,
mapped-pixel covariance, or the current zero E0--Fe covariance.

## Latest checkpoint: detector Product 3 implemented

`icbuilder/conductancedetector.py` and
`scripts/pipeline/make_conductance_detector_orbit_files.py` now implement the
new detector-first Product 3 without changing the old fixed-grid conductance
path. Schema-2 `conductance_detector` applies the shared icPhysics Robinson
function to schema-3 `precipitation_detector`, preserves the WIC detector
geometry and precipitation/proton state, and stores separate central and
uncertainty-valid masks. Missing analytic precipitation uncertainty leaves
central Hall/Pedersen conductance intact.

The runner defaults to `IR_hardy/robinson`, accepts separate `--base-input`
and `--base-output`, keeps `--base` as the input alias, supports `--workers`,
adds orbit context to worker failures, and uses source-aware atomic restart.
The old `ConductanceImage` and `make_conductance_orbit_files.py` remain the
distinct bin-first path.

A pre-fix two-worker scratch run generated all four full local examples in
49.89 seconds with 6,947,336-KB peak RSS. Their central-valid mask exactly
equalled all 25,328,544 Product-2 method-valid cells; 23,249,005 cells had
finite propagated uncertainty. The schema bump now makes those Product-3
files stale. Corrected orbit 0085 has finite propagated conductance uncertainty
on all 2,156,871 central-valid cells.

Seven focused Product-3 tests and 37 combined detector Product-1/2/3 tests
pass. The full suite has 116 passes and the same four known grid/lookup
failures. See [[Detector Product 3 Implementation Plan]] for schema details
and deferred `conductance_cs`/icReader work.

## Latest checkpoint: four full-length detector examples

The local pipeline tree at
`/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test` now contains
full-length Product 1 and pre-fix Hardy image-ratio Product 2 for orbits 0085,
0086, 0261, and 0968. Product 1 remains current; the schema-2 Product-2 files
are retained only as comparison evidence and require regeneration.
The Product-2 shapes are `(122, 256, 256)`, `(187, 256, 256)`,
`(272, 256, 256)`, and `(312, 256, 256)`. SI12 and SI13 each match 888 of the
893 WIC frames across the four orbits.

Product 2 contains 25,328,544 method-valid cell-frames. Per-orbit valid support
is 27.0%, 34.1%, 49.9%, and 49.4% of the stored detector cube. All
method-valid cells have finite, non-negative `Ep`, `Fp`, `E0`, and `Fe`.
Non-finite `R` is confined to the intentional legacy branch where corrected
WIC or SI13 is exactly zero and a low-signal central estimate is used.
In those pre-fix files, propagated `dE0` and `dFe` are available for
90.3--92.4% of method-valid cells because of the duplicate-noise defect.

The main scientific warning is the prevalence of proton-energy clipping.
Hardy `Ep_model` is outside the 0.47--46.7-keV Frey response-table domain for
52.1--66.2% of method-valid cells, almost entirely below 0.47 keV. The file
correctly preserves raw `Ep_model`, response-limited `Ep`, and the clipping
flag, so this is not corruption. It is a required interpretation/validation
step before treating the experimental image-ratio products as publication
retrievals or making a full-corpus claim.

## Latest checkpoint: Product-2 operational runner

`scripts/pipeline/make_precipitation_detector_orbit_files.py` now has the same
operational controls needed for the corpus run. `--base-input` locates
`fuv_detector/fuvpy_bs_directional_v1`; `--base-output` selects the Product-2
root and defaults to the input base. The old `--base` spelling remains an
input alias. The default retrieval is the fixed image-ratio method with Hardy
proton energy, producing `precipitation_detector/IR_hardy`.

`--workers N` processes independent orbits with `process_map` when `N > 1`;
one remains the serial default. Worker functions are module-level, execution
remains behind the existing `__main__` guard, each worker owns one atomic
output, and failures include the four-digit orbit number. Invalid worker
counts are rejected before filesystem work.

The focused Product-2 tests pass (`12 passed`), and the combined focused
Product-1/Product-2 tests pass (`30 passed`). A real two-orbit run with
`--workers 2` read one base, wrote another, and completed both Hardy
image-ratio files; the repeat run reported both complete. The full repository
suite has `109 passed, 4 failed`, with the same known 46-by-46 live-grid versus
36-by-36 tests/bundled Zhang--Paxton lookup failures.

## Latest checkpoint: UiB server Product-1 execution

The user confirmed on 2026-09-20 that
`scripts/pipeline/make_fuv_detector_orbit_files.py` completed successfully on
the UiB `dynamit` server against the server IMAGE-FUV products. The editable
installation and required environment dependencies were therefore sufficient
for execution. This closes the server execution blocker for detector Product
1. No exact orbit count, elapsed time, failed-orbit count, or aggregate product
audit was supplied, so those remain unverified rather than inferred.

## Latest checkpoint: multi-orbit detector throughput

The test tree at
`/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test` contains new
fuvpy WIC/SI12/SI13 products for orbits 0085, 0086, 0261, and 0968. During the
serial Product-1 run, orbit 0086 appeared stuck after orbit 0085 completed.
Read-only reproduction confirmed normal forward progress but measured
4.8--9.1 seconds per fully matched frame and 22 minutes 49 seconds for the
complete 235-MB, 187-frame product.

Profiling localized almost all of that cost to `geographic_to_wic()`. It sent
all 65,536 SI detector slots through the inverse WIC Delaunay interpolators even
though only roughly 6,000 had finite projected coordinates. The function now
interpolates and round-trip checks only finite/in-domain points, then places the
results back into full-shaped NaN arrays. No triangulation decimation or
geometry approximation was introduced. A real orbit-0086 frame gives
bit-for-bit identical old/new sparse overlap indices, areas, row, column,
coverage, and diagnostics.

After Numba warm-up, five orbit-0086 frames take 0.12--0.34 seconds for both
cameras. A complete optimized orbit-0086 run wrote the 245,119,212-byte product
in 65.46 seconds, about 21 times faster end to end, with 187 SI12 and 186 SI13
matches and 3.57-GB peak RSS. A repeat run skipped it in 0.56 seconds. At this
measured rate, 1,600 comparable orbits are about 29 serial hours or 2.9 ideal
hours with ten orbit workers; larger orbits and shared I/O will raise the real
parallel wall time.

The clean server environment also exposed that `pyproject.toml` declared no
dependencies and `icbuilder/__init__.py` eagerly imported every workflow.
`icbuilder.fuvdetector` therefore demanded unrelated packages such as `secsy`
before it could run. The package root now lazily resolves its existing public
objects, while the project metadata declares ApexPy, netCDF4, Numba, NumPy,
SciPy, and tqdm. A subprocess regression test confirms that importing the
detector does not import `secsy`, `icreader`, `icphysics`, or `sksparse`.

The initially reported orbit-0085 count of 5 SI12 and 7 SI13 matches was an
icBuilder decoding error, not inter-camera clock drift. `make_orbit_nc_files.py`
set `t_start` before `backgroundmodel_BS()` filtered out low-quality frames.
xarray then wrote each sensor's `date` values with correct CF units relative
to its first retained frame, but the scalar `t_start` remained the earlier,
pre-filter time. The Product-1 loader ignored the CF units and reconstructed
every sensor from that stale scalar, introducing false sensor-dependent offsets
of 11--14 seconds.

The writer now assigns `t_start` after background filtering. The Product-1
loader now decodes `date` from its authoritative CF units and calendar and
rejects missing, masked, non-finite, or invalid time metadata. A regression
test deliberately supplies a stale `t_start` and verifies the CF-decoded time.
Applied directly to the existing test sensor files, the corrected decoder
yields 121/122 independent matches for both SI12 and SI13 at two seconds, all
with 0--1-second separation. Existing detector Product-1 outputs made with the
old decoder are invalid and must be regenerated; the source sensor products do
not need regeneration for the reader fix. New Product-1 files store
`source_time_decoding = CF date units and calendar`, and restart validation
requires it. A normal rerun therefore treats pre-fix Product-1 files as invalid
and regenerates them instead of reporting them complete. Product 2 rejects a
Product-1 input without this marker, propagates it as
`source_fuv_detector_time_decoding`, and requires it in its restart validation,
so pre-fix downstream files are invalidated as well.

The historical old-fuvpy/old-icBuilder conductance file at
`/home/bing/Dropbox/work/data/IMAGE_FUV/conductance/or_0085.nc` confirms the
legacy behavior: it contains 131 timestamps. The old sources contain 235 WIC,
224 SI12, and 232 SI13 frames. The original one-to-one, three-camera rule forms
223 groups whose total span is at most two seconds, and all 131 stored
conductance timestamps belong to those groups; old downstream frame-content
filters remove the other 92. The rule retains the latest timestamp and discards
unmatched frames. Product 1 intentionally retains the WIC time axis and matches
SI12/SI13 independently, but the two-second threshold itself remains supported.

## Latest checkpoint: fuvpy BS-only orbit runner

`scripts/pipeline/make_orbit_nc_files.py` now uses the selected fuvpy
directional BS model directly. `read_idl()` reconstructs measurement variance
and viewing geometry; `backgroundmodel_BS()` uses the selected
`legacy_x_azimuth` model, 120/60-minute temporal spacing, damping `1e-2`,
variance weighting, monotonic fitting with ordinary fallback, and the Spencer
taper. The removed `backgroundmodel_SH()` API and obsolete `n_tKnots` option
are no longer called.

The first remote icProjects run exposed an environment mismatch: its installed
fuvpy `read_idl()` does not expose `wic_variance_model`, while the selected new
API does. Do not weaken the orbit runner to accommodate that older install;
update the icProjects fuvpy checkout/environment before running it.

The script accepts separate `--base_input` and `--base_output` paths. Inputs
always contain the orbit-index HDF files and raw sensor directories; products
and availability arrays are written beneath the output base. Omitting
`--base_output` makes it equal to `base_input`, and `--base` remains a
compatibility alias for `--base_input`. Sensor output directories are created
when absent. The no-argument default correctly resolves to
`<repository>/example_data` using `parents[2]`; the initially implemented
`parents[1]` incorrectly pointed to `scripts/example_data` and was fixed after
reproduction.

The first Python 3.14 parallel run failed during worker bootstrap because the
script executed its parser and pipeline at module import. All executable
workflow code now lives in `main()`, module-level functions remain picklable,
and the `__main__` guard calls `freeze_support()` before `main()`. A forced
two-worker `spawn` test processed WIC orbits 0085, 0086, and 0968 successfully,
verifying the same safe-import requirement used by Python 3.14 `forkserver`.

A serial scratch run completed WIC, SI12, and SI13 for example orbits 0085,
0086, and 0968. All nine products were successful and expose the required
BS-only corrected image, weight, spread, model components, measurement
variance, frame-quality, and viewing-geometry fields. Repository-wide Python
syntax parsing passed. `git diff --check` is clean for the changed orbit
script; the whole-worktree check still reports pre-existing trailing whitespace
in `scripts/dmsp/ratio_validation.py`.

A separate-path scratch run then read WIC orbit indices and IDL inputs from the
tracked example base while writing to a new output base. It created the output
directory, wrote orbits 0085, 0086, and 0968, and stored three successful
availability entries. Omitting `--base_output` was also checked and resolved it
to `base_input` without running a sensor workflow.

The detector integration boundary is now implemented. Product 1 consumes only
the selected BS `dgimg`/`dgweight`/`img_variance`/`frame_quality` contract and
validates the identifying `dgmodel` metadata. The abandoned internal prototype
and legacy-Ohma ingestion are not supported. See
[[New fuvpy Detector Product Implementation Plan]] for the implemented schema
and remaining full-orbit/corpus gates.

## Latest checkpoint: experimental detector Product 2

`icbuilder/precipitationdetector.py` and
`scripts/pipeline/make_precipitation_detector_orbit_files.py` implement the
documented first detector-space Product-2 slice. It is explicitly the
image-ratio parity path, not the candidate publication retrieval. The
calculation remains visible in order: match definitive GFZ Kp to Product-1
WIC times; evaluate Hardy or constant proton energy on the per-frame detector
MLAT/MLT; infer proton flux from already mapped SI12 counts; proton-correct
WIC and SI13; call `icphysics.precipitation_from_ratio`; then store channel and
method validity.

Product 2 preserves every WIC frame and the Product-1 detector geometry. A
pixel without finite WIC, SI12, and SI13 support is method-invalid rather than
removing its frame. Schema 2 stores Kp and interval, raw/response-limited Ep,
clipping, Fp, corrected WIC/SI13, R, E0, Fe, analytic uncertainty fields,
three-channel support/coverage/quality, method support, source times/indices,
and a direct Product-1 reference. Atomic publication and a direct restart
check follow the existing science-script pattern.

Product 2 requires schema-2 Product 1 and supplies the square roots of its
WIC/SI12/SI13 variances to icPhysics. The file records that mapped-SI target
covariance, cross-channel covariance, fitted-background uncertainty, and model
discrepancy remain omitted. The implemented map-counts-then-infer-flux order is
also
experimental and still requires the documented comparison against inferring
proton flux on native SI12 pixels before mapping.

The focused Product-1/Product-2 subset has 30 passing tests. The full suite has
109 passes and the same four known grid/lookup failures.

The fresh schema-2 orbit-0085 Product-2 smoke test completed in 17.97 seconds,
retained 358,925 method-valid cell-frames, carried all three source
frame-quality arrays, used detector variance, and skipped on repeat. This is a
partial-orbit execution/schema check, not publication validation.

## Latest checkpoint: detector Product 1 schema 2

`icbuilder/fuvdetector.py` preserves the WIC detector/time axes and
coregisters new-fuvpy SI12/SI13 observations onto WIC pixels with the
established footprint-overlap mapper.
`scripts/pipeline/make_fuv_detector_orbit_files.py` writes schema-2
`product_type = fuv_detector`, `representation = detector` orbit files with
atomic publication and restart validation. It supports orbit-level parallel
execution through `--workers N`: one remains the serial default, while larger
values use `process_map`, with each worker owning one orbit and output file.

A Halley run reached 1,670 of 1,687 scheduled orbits before a WIC transform
failed because a finite geographic pixel projected to non-finite tangent-plane
coordinates. `make_wic_transform()` now restricts its triangulation to pixels
that are valid in both geographic and projected coordinates. The worker adds
the four-digit orbit number to any future exception, and the progress message
now shows the requested output folder rather than the old hard-coded name.
The atomic files already completed by the interrupted run can be reused by
rerunning without `--overwrite`.

The implemented source boundary is `dgimg`, `dgweight`, `img_variance`, and
`frame_quality` for WIC, SI12, and SI13. The loader also checks the selected
BS-directional `dgmodel` metadata. WIC frames are processed in stored order.
Each SI channel independently selects
the nearest unused frame within two seconds; ties select the earlier SI time
and lower source index. Missing SI retains the WIC frame with NaN
signal/variance/weight, missing frame quality, zero coverage/support, and
source index `-1`. Product 1 stores WIC geographic and
130-km Modified-Apex geometry, per-channel values/weights/validity/coverage,
SI source counts, frame coregistration diagnostics, all source times/indices,
method constants, exact source SHA-256 hashes, and a content-derived Product-1
software identity. WIC variance is preserved directly; SI variance is
propagated with squared overlap weights under an independent-source-pixel
assumption. It explicitly records that induced covariance, fitted-background
uncertainty, and the sparse overlap operator are unavailable.

Restart validation was simplified on 2026-09-03 to match the existing
`binned_file_status()` style. One direct function checks product type/schema,
the requested preprocessing label and time tolerance, selected source paths
and image fields, required units, frame-quality semantics, and array shapes.
It does not parse a configuration
document, construct a processing fingerprint, return detailed diagnostic
objects, or hash source files. Exact source hashes and the content-derived
software identity remain stored as provenance. Descriptive metadata such as
the time-match explanation can be edited without changing restart status.

The fresh isolated orbit-0085 excerpt run first regenerated all three source
orbits with fuvpy `background_model_asymmetry` at `65603e0`, then wrote the
schema-2 `(20, 256, 256)` Product 1 in 142.80 seconds. Both SI cameras matched
all 20 WIC frames. Finite signal and variance support were identical within
each channel: 438,058 WIC, 423,726 SI12, and 410,400 SI13 cell-frames. A repeat
command skipped the completed file in 0.64 seconds. This is an execution/schema
result on a partial tracked excerpt, not complete-orbit scientific validation.

The focused Product-1/Product-2 subset has 30 passing tests. The full repository
suite now has 109 passes and the four already-known 36-by-36/46-by-46 grid and
Zhang--Paxton lookup failures. An independent read-only review found the final
restart identity, source identity, label invariants, matching provenance, and
schema validation adequate for this experimental slice.

The coregistration runtime gate is resolved by the exact finite-coordinate
optimization and complete orbit-0086 benchmark above. Corpus execution remains
a separate decision because scientific gates still include detector radiometry,
detector covariance, SI12 map-versus-infer ordering, and final preprocessing
and frame support.

## Latest checkpoint: detector-first product architecture

The target publication pipeline is documented in
[[Detector-First Product Architecture]]. Experimental detector Products 1,
image-ratio 2, and Robinson 3 now implement the pre-binning chain. The
candidate Zhang--Paxton Product 2 remains bin-first.

WIC detector geometry is the canonical geometry for Products 1--3. SI12 and
SI13 are corrected on their native detector grids and coregistered onto WIC.
Precipitation and conductance are then calculated on that detector geometry.
The fixed Cubed-Sphere products used by the VAE and splines are downstream
representations: Product 2 CS is binned from detector precipitation and
Product 3 CS is independently binned from detector conductance. The latter
must not be recalculated from Product 2 CS because spatial averaging and the
nonlinear conductance forward model do not generally commute.

The candidate Zhang--Paxton Product-2 path remains bin-first. The
36-by-36 Zhang--Paxton table is incompatible with arbitrary
per-frame WIC MLT and will need a grid-independent `(Kp, MLT)` lookup or
callable. Exact schemas, the fixed CS grid, binning uncertainty, and whether a
Product-1 CS representation is needed remain open. Conditional detector-count
variance is now available, while induced covariance is not; the order of SI12
coregistration versus proton-flux inference requires a numerical test.

## Latest checkpoint: historical FUVIEW background comparison

The four-way detector-space comparison is implemented outside the production
pipeline. `icbuilder/fuview_background.py` ports the fixed quiet-time lookup,
the active `p` SZA rescaling, and the active `p2` two-dimensional SZA/DZA
rescaling. `scripts/debugging/compare_fuview_backgrounds.py` compares those
three against current fuvpy on the same stored detector frames, then uses one
common detector coregistration, Hardy proton correction, 50/3 signal mask, and
Frey ratio inversion. It has run on the Frey event (orbit 0364), Coumans/FAST
event (0459), strong orbit 0968 frame, and weak pre-boom orbit 0085 frame.
Figures and exact metrics are under
`figures/debugging/fuview_background/`.

The active-`p` comparison now uses an explicitly labelled SZA-100 continuation.
Its final few fitted bins turn downward because the median window overlaps the
zero-filled SZA >= 105 part of the fitting image. Freezing only the scale was
ineffective where the reference background `B0` had already fallen toward
zero. The revised experiment therefore holds the complete DZA-dependent
product `B0(DZA, 100) * scale(100)` at every larger SZA. The function default
retains the historical zero fallback so the two experiments remain distinct.

On pixels that pass the 50/3 guard under every method, the four-event mean
median ratio / fraction above Frey's `R=136.486` limit is now:

- fixed lookup: 200.4 / 45.3%;
- active `p` with SZA-100 continuation: 85.0 / 21.3%;
- active `p2`: 96.4 / 25.9%;
- current fuvpy: 65.8 / 16.1%.

The revised active-`p` treatment improves substantially over the superseded
scale-only/final-column result (107.9 / 31.6%), but remains worse than current
fuvpy overall. The standard outputs under
`figures/debugging/fuview_background/` have been regenerated with this
treatment. Existing FAST statistics for the older edge treatment are stale
for this comparison; the FAST pass has not yet been rerun with the complete
SZA-100 background continuation.

There is event dependence. `p2` gives a lower shared-support above-limit
fraction than fuvpy for orbit 0968 (14.6% versus 18.7%), but is worse for
orbits 0364, 0459, and 0085. The fixed table is especially poor. The
historical background choice clearly changes the inferred ratio, but this
first result gives no evidence that replacing fuvpy with fixed, `p`, or `p2`
backgrounds generally restores the published ratio regime.

The underlying ports are source-faithful rather than a bit-for-bit FUVVIEW
replay; the active-`p` SZA-100 continuation is an explicitly new experiment.
Saved `img` is already calibrated, WIC was reflattened by fuvpy, and the
ambiguous detector-bias/raw-flat-field wrapper operations are excluded. The
literal `p2` SI code cannot run because `min_avg` is defined only for WIC; the
port transparently applies the WIC positive-minimum/floor rule independently
to SI12 and SI13 and labels this repair. Focused background and detector-
coregistration checks pass (`11 passed`).

The Coumans FAST validation is also implemented in
`compare_fuview_backgrounds_fast.py`. It recalculates the nine IMAGE frames
spanning the pass and samples the WIC-grid results at the mapped 120-km FAST
footprint. Active `p2` retains 155 samples versus 132 for fuvpy. Their own-
support median absolute errors are 1.42 and 1.91 keV, but this is not a clean
`p2` win: on their 124 common samples, `p2` has a slightly smaller central
error but a worse upper tail and places 12.1% at the 25-keV Frey fallback.
Outputs are
`coumans_fast_background_comparison.{png,pdf}`, the complete sample CSV, and
the method-metric CSV in the same figure directory.

The fixed branch uses the IDL zero-background fallback outside the quiet-table
domain. Active `p` with the same fallback retained 226 FAST samples, with
median IMAGE energy 6.83 keV and median absolute error about 4.30 keV. The
earlier FAST edge result (206 samples and 14.71 / 12.21 keV median IMAGE energy
/ median absolute error) used the now-superseded scale-only/final-column
continuation. Do not use it to judge the complete SZA-100 treatment until the
FAST comparison is rerun. Original lookup support and continued pixels remain
separately identified. The optional WIC bias surface remains excluded because
its conversion into the calibrated saved-image count space is ambiguous.

The orbit-0968 `p2` over-subtraction is intrinsic to its recovered 2-D scaling
scheme. Its WIC scale field rises from a median 1.64 to a 95th percentile of
58.2 and maximum 89.8 near SZA 100--110; orbit 0364 reaches only 5.59 and 7.09
at the same quantiles. SI13 behaves similarly. Sparse `binimg/reference`
samples are zero-excluding-median filtered and nearest-filled, so a small quiet
reference can turn residual science counts into an enormous background. The
abrupt switch to the minimum fallback at SZA 110 creates the visible boundary.
Because WIC fails too, this is not caused by the repaired SI `min_avg` bug.

Selected DMSP crossings are the remaining repeatability test. Do not add a
historical-background option to Product 1 unless that comparison shows a
material and repeatable improvement. F10.7 and clock-angle corrections remain
tabled: the source names their routines, but the archive lacks the
implementations and support coefficients, and active `p`/`p2` bypass F10.7 in
their normal path.

The recovered chronology favors, but does not prove, the fixed model for Frey
et al. (2003). Frey's Section 7 describes a quiet-time SZA/DZA model with
sensitivity and F10.7 adjustments and does not mention current-frame scaling;
the recovered active path is marked as new work on 28 June 2002 and explicitly
uses `/nof107`. The change log says the February 2003 `p2`/clock-angle attempt
was left unimplemented because it did not work inside FUVVIEW, and the wrapper
comments out `p2`. Do not describe Frey's exact routine as confirmed without a
versioned processing record, but treat fixed as the best-supported reading.

The Frey Figure-16 script now has a `--image-source fixed` branch. It subtracts
the recovered fixed lookup on each native detector grid, then uses the same
SI-to-WIC coregistration, Hardy proton correction, 50/3 mask, and plotting as
the existing branches. The generated orbit-0364 output is
`figures/debugging/paper_reconstruction/frey_figure16_or_0364_fixed_hardy_mask50_3.{png,pdf}`.
It shows widespread values above the 10-keV display ceiling and does not
recover Frey's final-energy morphology. Interpret this only as a failure of
the available bare lookup: Frey's stated F10.7 adjustment and exact historical
calibration order are still missing.

It also has a source-faithful `--image-source active_p` branch using the
original zero fallback rather than the experimental edge continuation. The
orbit-0364 output is
`frey_figure16_or_0364_active_p_hardy_mask50_3.{png,pdf}`. Active p removes
some scattered high-energy pixels relative to fixed, but the main oval remains
widely saturated above 10 keV and is still unlike Frey's final-energy panel.
The title records that WIC bias and F10.7 are omitted.

The new experimental `--image-source active_p_edge100` branch holds the full
background as well as its scale at SZA 100 degrees. Its orbit-0364 output is
`frey_figure16_or_0364_active_p_edge100_hardy_mask50_3.{png,pdf}`. The
common-support median ratio changes from 103.7 to 72.9 and the fraction above
Frey's ratio ceiling from 32.7% to 20.5%, but the energy map still has extensive
saturation and does not reproduce Frey's final panel.

### Source details retained for the comparison

The recovered FUVIEW3 source identifies a genuinely different background
branch worth testing before abandoning the WIC/SI13 ratio.  The production
wrapper calls `image_bckgnd_active_p.pro`, not the newer `p2` routine.  In
active mode it uses separate 90-by-110 quiet-time SZA/DZA tables for WIC,
SI12, and SI13, estimates an SZA-dependent scale from non-oval pixels in the
same science frame, and subtracts the scaled background.  Its fitting mask is
approximately `DZA < 70`, `SZA < 105`, and `|MLAT| < 60` or `|MLAT| > 75`;
the observed SZA/DZA image is filtered with `median(..., 9)`.  Normal FUVIEW
use disables the available F10.7 correction and leaves the clock-angle call
commented out.

All three lookup tables load directly with `scipy.io.readsav`.  Existing
full-orbit sensor NetCDFs contain native-grid `img`, SZA, DZA, and MLAT, which
is enough for a bounded FUVIEW-style comparison.  Because those `img` values
are already calibrated and WIC has undergone fuvpy reflattening, this is not
automatically an exact replay of FUVIEW's raw-count order.  Native IDL frames
provide the closer test, and the fact that the active model is fitted per
image means the partial-orbit limitation that affects fuvpy refits does not
invalidate this experiment.

The comparison keeps separate functions for loading/mapping the reference,
the `p` active rescaling, and the newer `p2` rescaling. `p2` is not a cosmetic
revision: it replaces the
one-dimensional SZA scaling with a filled SZA/DZA ratio field, calls
`median_arr(..., 5, 0)` rather than `median(..., 9)`, and changes the fitting
limits.  The helper's actual inclusive indexing makes that call an up-to-11 by
11 zero-excluding window despite the misleading comment in its header.  The
supplied `p2` source also defines `min_avg` only for WIC but uses it later for
SI12 and SI13; the implemented repair is exposed explicitly. The retained
outputs include all four backgrounds, corrected detector images, ratios,
energies, fit supports, and scale factors.

F10.7 and clock-angle adjustments are tabled until after this four-way
comparison.  They remain possible downstream sensitivity tests, but not all
combinations reproduce historical processing.  The archive is missing the
F10.7 retrieval/correction routines and coefficient
files named by the source, and it is also missing
`image_bckgnd_clock_correct.pro`.  Daily F10.7 can be sourced externally and
the pixel clock angle can be reconstructed from the saved time and geographic
position, but the exact IMAGE-specific corrections cannot yet be reproduced.
The original `p` and `p2` code also bypasses F10.7 whenever active correction
is selected.  An active-plus-F10.7 result would therefore be an explicitly new
sensitivity experiment.  The Immel (2000) published correction coefficients
belong to DE-1 and must not silently replace the missing IMAGE coefficients.

## Latest checkpoint: Coumans Figure 5 and response feasibility

`scripts/debugging/reconstruct_coumans_figure5b.py` now performs the first
same-event FAST/IMAGE reconstruction for 23 December 2000. The official
`FA_ESA_L2_EES` file contains 529 distributions over the published interval.
An SSCWeb subset supplies the missing spacecraft ephemeris. Mapping to 120 km
gives 60.0--81.7 degrees modified-Apex latitude and 4.58--6.96 MLT, consistent with the
paper, and a local bounce loss cone of 28.9--36.0 degrees.

Coumans does not state its FAST low-energy, weak-precipitation, or smoothing
rules. Using all documented quality/bin flags, the mapped loss cone, energies
at or above 50 eV, and a centred 10-s mean gives a 4.29-keV FAST maximum,
close to the stated 4.5 keV. This is an explicit plausible convention, not yet
proof of exact historical processing. The current Hardy-proton image-ratio
product reaches 11.41 keV along the same track; 137 FAST samples pass its 50/3
count guard. Coumans' published IMAGE curve peaks near 1.5 keV. Outputs are
`figures/debugging/coumans_2004/coumans_2004_figure5b_reconstruction.{png,pdf}`.

The same script now reconstructs Figure 5a from the native orbit-0459 WIC
frame at 21:00:16 UT and overlays the FAST footprint mapped to 120 km. It uses
the fuvpy-corrected WIC `shimg` field and is not used by the quantitative ratio
branch. The output is
`coumans_2004_figure5a_reconstruction.{png,pdf}` in the same directory.
Coumans' blue `<E> model` curves are not fits to FAST or IMAGE: panel b uses
the Kp-dependent Hardy et al. (1985) DMSP electron mean-energy climatology;
panel c uses the Hardy ion climatology (caption: 1991; related 1989 model).
Section 4.1.2 confirms 21:00 UT, making the available 21:00:16 frame the right
choice. The published smooth saturated limb is not reproduced by current
fuvpy `shimg`; it is consistent with the historical quiet-time-airglow
background/display path. Raw `img` retaining a similar crescent is only a
visual coincidence and does not establish that uncorrected counts were used.
The current panel already applies `DZA < 75 degrees`. The brightest 0.5 percent
of the frame has median DZA 70.1 degrees and 95th-percentile 74.5 degrees;
60--65-degree cuts remove the bright crescent. A historical cutoff may shape
the outer boundary, but does not explain the processing mismatch on its own.
The companion full-field DZA diagnostic is
`coumans_2004_figure5a_dza.{png,pdf}` in the same output directory.
The corresponding SZA/terminator diagnostic is
`coumans_2004_figure5a_sza.{png,pdf}`.

The literature search found one independent observed-IMAGE reimplementation,
Gasparini et al. (2024), but no independent in-situ validation of recovered
electron energy. Gasparini assumes the Frey inversion and tests downstream
uncertainty. All located direct validations trace to the original collaboration.
Mende et al. (2003) explicitly limited its February-2002-event ratio retrieval
to SI13 above 10 counts per pixel. That paper states a 2.5-pixel SI12 smoothing,
SI12-based proton subtraction from WIC/SI13, and one-minute smoothing of the
FAST flux, but does not document a WIC/SI13 background algorithm, a WIC/SI13
smoothing width, or whether the threshold preceded proton subtraction. Do not
treat 10 counts as a generally validated IMAGE processing rule.

Composition and viewing-angle corrections are implementable only as forward-
response models. The first bounded test should digitize Meurant's standard
O/N2=0.52 and perturbed 0.37 curves as an uncertainty envelope. A full model
requires separate WIC/SI13 response tables over energy, DZA, and atmosphere,
including line-of-sight O2 absorption. `cos(DZA)` cannot do this. Product 1
retains the needed angles, but the legacy bin-first Product 2 drops them and
the remaining atmospheric, geographic, and spectral-response inputs are not
assembled.

## Latest checkpoint: Hardy proton energy in the modular pipeline

The complete Hardy, image-ratio, footprint-weighted Product-2 corpus is at
`/home/bing/Dropbox/work/data/IMAGE_FUV/=precipitation_IR_hardy_weighted`.
All 1,684 orbit files opened successfully. Together they contain 450,601
frames; every file reports schema 2, `image_ratio`, Hardy proton energy, and
SI12 proton flux, with all required Ep/Fp and electron-precipitation fields.

Product 2 now uses SI12 for event-specific proton flux and Hardy et al. (1991)
for proton mean energy by default. Schema 2 stores raw `Ep_model`, response-
clipped `Ep`, dEp, `Ep_clipping_flag`, Fp, dFp, separate flux-source/energy-
model metadata, coordinate provenance, and the explicit statement that Hardy
dEp is not modelled. Product 3 carries these fields unchanged. The constant-Ep
path remains available with `--proton-energy-model constant`.

Focused verification passed: 15 icPhysics tests, 33 icBuilder tests, and 10
icReader tests. An isolated orbit-0364 image-ratio Product-2/Product-3 run on
the active 46x46 grid wrote and reloaded successfully without touching the
corpus. The existing 36x36 Zhang--Paxton lookup mismatch still blocks that
precipitation branch on the active grid and is not part of this fix.

## Latest checkpoint: high WIC/SI13 ratio diagnosis

**Critical native-data limitation:** the available debugging IDLs are partial
orbit extracts: 20 frames per sensor for 0085/0086 and three for 0968. All
diagnostics that estimate one reflat constant over the loaded stack or refit
the time-dependent fuvpy background on these files are exploratory, not
orbit-level validation. This includes the WIC/SI reflat, histogram-background,
recovered-calibration, and SI-smoothing results below. Repeat them using every
native WIC/SI12/SI13 frame from at least one complete orbit before using them
to rule a mechanism in or out. Complete-corpus Product-2 statistics, source/
table provenance, and synthetic coregistration tests are unaffected.

The large-ratio problem is confirmed across the complete local image-ratio
corpus, not just the Coumans event, and the Hardy rerun has now been tested.
The Hardy and footprint-weighted fixed-2-keV corpora contain the same 1,684
orbits and 450,601 frames. Time and WIC source indices match in every paired
file, while a sampled pair also has identical sensor arrays, weights, grid,
and Kp. On each corpus's own positive support, 36.43% of fixed-2-keV ratios and
38.16% of Hardy ratios exceed the Frey-table maximum `R=136.486`. With the
modern diagnostic `WIC >= 50`, `SI13 >= 3` guard, the corresponding fractions
are 13.55% and 14.37%.

`compare_hardy_fixed_proton_corpus.py` performs the stricter cell-for-cell
test. Among 59,152,481 grid cell-frames that pass the guard under both proton
corrections, the out-of-range fraction increases from 13.47% to 14.36%.
Hardy moves 536,272 cells into the invalid range and 14,485 out of it. It lowers
the ratio in 51.57% of common guarded cells and raises it in 40.80%, but the
upper tail becomes worse. Positive support also contracts: 17,169,634 cells
are valid only under fixed 2 keV and 3,173,793 only under Hardy. The
orbit-balanced sample has median `Ep=4.94 keV`; 29.76% of common guarded cells
are clipped to the 0.47--46.7-keV Frey response-table range. The outputs are in
`figures/debugging/hardy_fixed_proton_comparison/`.

Hardy is the better-documented proton correction and reproduces Frey's stated
model choice, but it does not make the WIC/SI13 electron-energy inversion
reliable. All fractions are descriptive counts of correlated grid cell-frames,
not independent samples, and the thresholds are not historical Coumans rules.

Weak SI13 therefore causes much of the spectacular tail but is not the whole
problem. The earlier fixed-2-keV audit showed orbit 0085 and 0086 to be
especially poor, and the new paired result demonstrates that changing the
proton-energy model cannot resolve the corpus-wide failure.

The orbit-to-orbit failure fraction has a recurring seasonal structure rather
than a monotonic mission-time drift. Kp is only weakly associated with it
(Spearman `rho=-0.052` before the signal guard and `+0.095` after it). Guarded
monthly medians are about 0.25 in June and 0.22 in July, compared with roughly
0.05--0.07 from October through January. Mean positive WIC weight
(`rho=-0.438`) and SI13 positive-weight coverage (`rho=-0.405`) are more
informative. This supports an unresolved seasonally varying
observation/instrument/processing effect with a quality-support contribution;
illumination, background, viewing conditions, and calibration remain
confounded. These are descriptive equal-orbit associations using the
pre-proton ratio. The weight summaries cover the full arrays, neighboring
orbits are temporally correlated, and the nominal correlation p-values should
not be interpreted as causal evidence.

The partial-frame orbit-0968 detector-space stage budget suggests against
several specific explanations but cannot reject them at orbit level. The raw
mapped median ratio is 233 and 71.9% of positive pixels
are beyond the Frey maximum. Current fuvpy background subtraction moves in the
right direction, lowering these to 93 and 23.3% on its smaller remaining
positive support; this is not a fixed-pixel comparison. Common smoothing
lowers the fraction to 21.7%, and the 2-keV proton correction raises it
slightly to 23.8%. Re-estimating the histogram after SI-to-WIC mapping,
following the documented Meurant ordering, makes the final three-frame median
ratio worse (about 116). Halving SI12-derived proton flux, as Coumans did,
changes the median only from 95.0 to 94.1. Neither changes the IMAGE peak
time. The mapped-histogram, half-proton, and current-fuvpy curves peak at 23:38:31;
the 66-s DMSP mean peaks near 23:40:24.

The Frey Figure-16 reconstruction has clarified that its full-polar
preprocessing is not the Meurant nightside histogram method. Frey describes an
instrument-specific, spatial dayglow model based on quiet-time observations,
solar and spacecraft zenith angles, sensitivity, and F10.7. Meurant's scalar
brightness-histogram subtraction is explicitly restricted to nightside data,
where no airglow correction was needed. As expected, applying only scalar
knees to raw `img` (`WIC=1000`, `SI12=15`, `SI13=8` in the latest plot) leaves
strong dayside contamination and is not a valid reproduction of Figure 16.

A complete reference audit found no hidden numerical recipe in Frey for the
remaining Figure-16 operations. Frey specifies flat-field/dayglow and temporal
gain correction, but gives no smoothing kernel, count mask, geomagnetic-grid
resolution, or aggregation rule. Its "temporal changes" are whole-image
sensitivity/gain corrections relative to the 19 June 2000 stellar calibration,
not timing or smoothing. Frey never states that WIC or SI13 was smoothed and
does not mention one-to-two pixels; its only averaging is of fine-resolution
FAST samples for validation. Hubert et al. (2002), the closest cited method
paper, explicitly maps SI13 into WIC image space before a pixelwise ratio but
does not specify interpolation, smoothing, or thresholds. Hubert does use
view-angle-dependent response curves for the actual pixel; forcing all pixels
to nadir changed its test-case hemispheric power by 12--17%. The static nadir
Frey lookup in the present reconstruction does not implement this Hubert step,
but Frey does not state that Figure 16 used Hubert's angle-dependent code.
Coumans et al.'s PSF
smoothing applies to IMAGE--NOAA validation, not the WIC/SI13 inversion.
Meurant et al. explicitly coregisters and smooths the cameras, but Frey does
not cite it, so this remains a related-method sensitivity rather than a
confirmed Figure-16 step. Frey's undocumented disturbed-atmosphere adjustment
is stated only for the active 24 June 2000 FAST validation in Figure 15, where
enhanced oxygen absorption suppressed SI13. No atmosphere specification,
corrected curve, formula, factor, or implementation citation is given, and Frey
does not say that it was applied to Figure 16. The exact Figure-16 implementation
is therefore not fully recoverable from the published papers alone.

The Figure-16 script now generates a separately named 50/3 sensitivity variant
from the same clean fuvpy inputs (`WIC=shimg`, `SI12/SI13=dgimg`). The generated
`frey_figure16_or_0364_mask50_3.{png,pdf}` masks uncorrected energy on the
uncorrected counts and masks final energy/flux on the proton-corrected counts;
the camera panels are unchanged. The 50/3 guard remains a modern diagnostic,
not a documented Frey step.

The weak visual impact of proton correction in the reconstruction is explained
by its fixed 2-keV proton assumption. For a given SI12 signal, the Frey-table
WIC/SI13 subtractions increase from 8.96/0.123 times SI12 at 2 keV to
35.35/0.424 at 25 keV. In the 17--20 MLT sector, the 50/3 support falls from
1,286 to 1,158 pixels at 2 keV but to 373 at 25 keV. Frey says energetic
protons supply most WIC/SI13 signal there, and its general procedure uses Hardy
proton energies, although the Figure-16 proton-energy map is not published.
Thus the dark published final-energy sector is plausibly lost electron support,
not simply a lower-valued map. The higher published uncorrected energy cannot
be caused by proton correction and remains an independent ratio/response issue.

The required 1989/1991 Hardy ion model has now been implemented in icPhysics
from the published coefficient table; Lompe's unrelated 1987 electron routine
is not used. `recreate_frey_figure16.py` now accepts `--proton-model hardy`,
evaluates mean ion energy on the WIC MLT/MLAT grid using definitive GFZ Kp,
and passes that spatial map into the existing SI12-based proton correction.
The old constant-energy path and filenames are preserved. A scratch run on
orbit 0085 completed successfully, followed by the full orbit-0364
reconstruction at definitive `Kp=4.667`. The output is
`frey_figure16_or_0364_hardy.{png,pdf}`, with the 50/3-threshold counterpart
in `frey_figure16_or_0364_hardy_mask50_3.{png,pdf}`. Where mapped SI12 is positive, the
Hardy mean proton energy has a median of 9.06 keV. Relative to the fixed 2-keV
case, the 50/3-supported final map retains 76.5% as many pixels overall and
46.5% as many in the 17--20 MLT sector. This is qualitatively closer to the
strong dusk-sector proton removal in Frey's published figure, but the high
uncorrected WIC/SI13 energy remains unchanged by construction.

The side-by-side diagnostic
`frey_figure16_si12_hardy_comparison.{png,pdf}` now includes the modeled
proton-count terms subtracted from WIC and SI13. Although Hardy mean energy is
finite across a much wider region than the observed SI12 image, the correction
maps follow SI12 because the flux amplitude is
`Fp = clip(SI12, 0) / Tmodel(Ep)`. Hardy selects the energy-dependent camera
response; it does not impose the spatial correction amplitude. There is no
explicit SI12 significance threshold, so small positive SI12 noise still
contributes.

A common multiplicative image-scale sensitivity is now exposed through
`recreate_frey_figure16.py --image-scale`. Hardy+50/3 versions for 0.95, 0.90,
0.85, and 0.80 are saved with `_scale95`, `_scale90`, `_scale85`, and
`_scale80` suffixes. Scaling all three cameras equally leaves the corrected
ratio and retrieved mean energy invariant on common support to numerical
precision, scales electron energy flux by the same factor, and only removes
additional pixels through the fixed 50/3 mask. The retained final-energy
supports are 6,296, 6,166, 5,999, and 5,837 pixels. Common scaling therefore
does not test the unresolved relative WIC/SI13 calibration.

The 0.85 common-scale version is the closest visual match to the published
camera panels so far, but the reconstructed energy flux remains too bright on
the flanks, especially at dusk. This is not a common-gain effect. The WIC
energy-flux response `fWm(E0)` decreases from about 410 at 4 keV to 285 at
8 keV and 223 at 10 keV; consequently, the existing excessive flank energies
inflate inferred flux for the same corrected WIC signal. Stronger spatially
selective proton removal or an undocumented final support threshold remain
additional possibilities.

The Figure-16 script can also footprint-average the three coregistered count
images onto a polar magnetic grid before any ratio or proton calculation via
`--magnetic-grid-km`. At 100 km, the Hardy + 50/3 result is smoother and the
peak WIC/SI12/SI13 counts decrease, but the high-energy discrepancy remains:
the fraction above 10 keV changes from 14.0% to 11.4% before proton correction
and from 18.7% to 17.4% afterward. The energy-flux fraction above 10 mW/m2 is
13.7% versus 13.6%. The output is
`frey_figure16_or_0364_hardy_maggrid100km_mask50_3.{png,pdf}`. Treat these as
descriptive cell fractions because the native and gridded supports differ.

`digitize_frey_figure16_counts.py` extracts the original 992x885 Figure-16
JPEG, classifies confidently colored pixels against the measured 16-color
palette, and now compares against the native-pixel reconstruction rather than
the discarded 50-km version. Omitting the ambiguous black/lowest bin, the mean
color-bin indices are Frey/native = 3.27/4.42 (SI12), 4.70/7.05 (WIC), and
2.92/4.62 (SI13). Published/native white-bin fractions are 0.62%/3.27%,
9.44%/19.28%, and 0.50%/3.47%. All three current count maps are brighter, with
the largest absolute color-bin displacement in WIC. Treat this as binned
raster evidence rather than an exact gain estimate because JPEG compression
and black/missing ambiguity remain.

The script also digitizes the complete published bottom row. Direct
published/native-Hardy+50/3 mean nonblack bins are 6.37/7.00 for uncorrected
energy, 3.36/8.70 for final energy, and 3.10/6.43 for final energy flux.
Published/native white-bin fractions are 31.83%/26.50%, 12.03%/42.22%, and
5.78%/25.70%. These are valid broad distribution comparisons and confirm a
large final-product difference, but do not locate its cause.

Do **not** use the attempted recomputation from the digitized top row as
scientific evidence. The independently rasterized/JPEG-compressed panels do
not preserve a recoverable common numerical grid. Black removal creates an
intensity-selected 7,452-pixel intersection, approximate panel alignment
mispairs sensor counts, and color-bin centres are inadequate inputs to a
nonlinear ratio. The recomputed maps resemble neither Frey's bottom row nor
the native result. The previous inference that this ruled out input-count
differences or demonstrated a proton-correction/response mismatch is
withdrawn. Remove or clearly invalidate the middle row before presenting the
diagnostic; retain the independent panel distributions only.

Hardy source documentation is collected in
`literature/Hardy_ion_model_sources.md`. Both primary papers are now local.
The 1991 functional paper contains the full Fourier coefficient table, Epstein
equations, limiting levels, and reconstruction sequence. It defines mean ion
energy as integral energy flux divided by integral number flux and requires
interpolation of evaluated log-flux values rather than coefficients for
intermediate Kp. No additional document is required before implementation;
the historical FORTRAN would only provide an independent cross-check.

The newer-model check identified OVATION Prime as the leading practical
alternative to Hardy. Its empirical diffuse-ion component fits energy flux and
number flux separately to solar-wind driving, yielding mean ion energy, and is
available in public IDL, Python (`OvationPyme`), and Fortran/GITM code. It is
not species-resolved and requires historical IMF/solar-wind input. The 2015 PGI
ion Auroral Precipitation Model (APMI) also provides average ion energy and
energy flux from AL/Dst, but uses only 1986 F6/F7 data and coarse
precipitation-zone statistics; its advertised API was not verified as usable.
Use Hardy first to reproduce Frey's stated method, then compare OVATION Prime
as the newer pipeline candidate.

The supplied FUVVIEW3 source contains the empirical dayglow calculation in
the live `getudf_var.pro` -> `image_fuview_bckgnd.pro` ->
`image_bckgnd_active_p.pro` path. Its wrapper passes `/nof107` and iteration
zero, disabling the routine's F10.7 branch. Frey's paper says an F10.7
adjustment was applied. Treat this as a verified source/configuration
difference, not yet as proof that the present dataset omitted the correction:
the supplied source may not be the exact historical producer version.

`scripts/debugging/recreate_meurant_figures7_9.py` now reconstructs Meurant
Figures 7 and 9 from 86 native orbit-0364 frames (10:04--12:58 UT). It maps SI
into each WIC detector frame with footprint-overlap weights, applies scalar
histogram backgrounds, smooths WIC and SI13, proton-corrects, and performs the
published MLT/MLAT averages. The WIC keograms reproduce the main published
morphology well. The ratio keograms, however, contain substantially more
values above the paper's 130 upper scale, particularly in Figure 9. This is
evidence that the excessive ratio survives a close approximation of the
historical nightside workflow. The exact knee rule, smoothing kernel, and
event proton energy remain undocumented assumptions exposed as script
arguments.

The white ratio regions are verified over-range values rather than NaNs: NaNs
are black. With the pipeline-consistent 2-keV proton mean energy, 3.4% of
finite Figure-7 bins and 25.4% of finite Figure-9 bins exceed 130. Binning
currently takes the arithmetic mean of finite pixelwise
ratios. Ratios with a corrected SI13 denominator close to zero can therefore
dominate a bin. Test a defensible signal/support mask and check whether the
historical workflow used an undocumented threshold before drawing a stronger
calibration conclusion from this reconstruction.

Do not replace this with a ratio of mean counts when reproducing Meurant:
their method first constructs a pixelwise corrected-ratio image and then says
the pixels are averaged within each keogram bin. The WIC >= 50 and SI13 >= 3
guard is a useful sensitivity test, but it is documented in the modern
Gasparini/legacy-icBuilder handling rather than in Meurant's paper.

The script now saves parallel 50/3-mask figures and variables. At 2-keV proton
mean energy, the guard lowers the fraction of finite bins above 130 from 3.4%
to 1.0% (Figure 7) and 25.4% to 10.3% (Figure 9), and removes the most extreme
near-zero-denominator values, but Figure 9 remains more frequently saturated
than the paper.
Meurant explicitly derives proton flux from SI12 and subtracts the modeled
proton contribution from WIC and SI13, as implemented here. Their histogram
background description contains no knee-selection or fitting rule.

Changing the reconstruction from 8 to 2-keV proton mean energy has only a
modest effect: unmasked Figure-7/9 over-range fractions move from 5.1%/28.1%
to 3.4%/25.4%, and masked fractions from 1.6%/13.4% to 1.0%/10.3%. Frey's
general method obtains proton mean energy from Hardy statistical models; its
June-2000 validation example uses 25 keV, but Figure 16 does not state a value.
Do not pursue the stale `icPhysics` characteristic-energy label as a scientific
issue; it is a documentation error on a mean-energy response axis.

The supplied FUVVIEW3 SI active-background algorithm can be ported to Python.
The saved IDL records and support databases provide the lookup model and
geometry. The original routine is called before calibration, internally
flat-fields raw counts to fit its quiet-time table, unflat-fields the modeled
background, subtracts in raw space, and then calibrates the residual. Since
these operations are linear, a Python port can instead use the saved
calibrated/flat-fielded image and keep both image and background in that space.
Do not feed the saved image to the unmodified IDL wrapper because that would
flat-field it twice. The current records have zero `AIRGLOW_SCALE`, so no prior
FUVVIEW3 airglow removal is indicated.

The numerical mapping checks pass: constant SI fields are preserved exactly,
map/subtract closure is within `8.5e-14`, and the camera morphologies show no
row/column reversal. This does not validate the historical optical
coregistration or pointing. Along the orbit-0968 DMSP track, current fuvpy has
median `R=83.8` and no samples beyond the Frey maximum, so the Coumans mismatch
is a displaced moderate/high-ratio feature rather than the map-wide extreme
tail.

The fuvpy audit found no exposure-time or ratio-unit error. It copies the IDL
`IMAGE` field, selected geometry, instrument ID, and whole-second time; applies
an additional WIC-only reflat; and discards calibration flags, detector high
voltages, source identifiers, WIC `CIMAGE`, and other housekeeping metadata.
It performs no temporal gain, temperature, voltage, or exposure correction
itself. Primary literature says FUVVIEW3 corrected-count products nominally
include flat-field and mission-time/temperature corrections. Ohma et al.
(2024) used the same FUVVIEW3-processed inputs and explicitly state that all
camera images are corrected for lifetime temperature and voltage changes and
that flat-fielding is part of that processing. The simple hypothesis that
SI12/SI13 were never flat-fielded is therefore weakened. The exact FUVVIEW3
version, correction functions, and relative WIC/SI13 gain scale remain
unverified. In particular, fuvpy's `reflat` option is guarded by `id == WIC`;
enabling it for SI12 or SI13 has no effect. This is intentional: it repairs a
WIC-only interaction in which the original flat field also scales a large
constant detector background. Ohma et al. did not observe the corresponding
row artifact in SI12 or SI13, probably because their constant backgrounds are
small. Every local example has `CALIBRATION_FLAG=3` for WIC and `1` for
SI12/SI13. The supplied source defines flag 1 as corrected counts, but a
source-wide search finds no image-calibration meaning for flag 3 and no
`CIMAGE` field. A fair native-image test refitting the same background
model with WIC reflat enabled and disabled changed the median ratio by +2.4%,
-0.9%, and +2.0% for orbits 0085, 0086, and 0968. The fraction above the Frey
limit changed by less than 1.3 percentage points in each case, and the large
failure remains under the original FUVVIEW3 flattening. Reflat can still
amplify individual weak-residual pixels and is not negligible in the tail.

The SI13 denominator is much more sensitive to count-scale systematics. On the
fixed full-corpus guard (`WIC >= 50`, original `SI13 >= 3`), adding 1, 2, or 3
counts to SI13 reduces the out-of-Frey fraction from 14.48% to 10.50%, 7.86%,
and 6.04%. This does not justify adding a constant: it demonstrates that the
current ratio is ill-conditioned near its three-count cutoff. Ohma et al.'s
statement that the SI constant is visually small does not establish that its
uncertainty is negligible in WIC/SI13. An SI13 residual-background estimate and
uncertainty-aware denominator criterion are needed if the ratio is retained as
a diagnostic. SI13 SH background correction raises rather than lowers the
ratio in the tested examples.

The WIC reflat database itself is a legacy calibration artifact rather than a
fuvpy output. Its embedded metadata records creation in 2015 by user `hfrey`;
it contains two `(256,)` WIC row profiles labeled `2000/255` and `2001/205`.
Anders Ohma added the binary unchanged in 2021, and no generation code or
uncertainty is available in fuvpy. The recovered FUVVIEW3 source resolves the
date semantics: labels mark calibration-image dates, while the code selects
the nearest profile after first separating observations before and after the
spacecraft boom loss on day 278 of 2000. fuvpy's October 2000 split represents
that event; the 2001 profile is intentionally reused for later years.

The recovered support directory contains four WIC, twelve SI13, and three SI12
profiles plus mission-sensitivity corrections through 2005. The active code
uses `si_flatfield_dbase.idl` and `si12_flatfield_dbase.idl`; the separately
named `new_si_flatfield_dbase.idl` is not referenced. Calibration flag 1 means
corrected counts, including flat-field and mission correction, which is
consistent with the local SI IDLs. The calibration GUI and calculation path
define only flags 0--2; WIC flag 3 and `CIMAGE` are absent from this code
version, showing that it is not the exact WIC producer tree or that an external
export step altered the record. fuvpy's two WIC curves are a subset of
FUVVIEW3's four curves and its routine omits FUVVIEW3's post-boom
rotation/shift. Consequently `_reflatWIC`
does not exactly reverse the recovered upstream calibration. Before testing an
SI constant-preserving reflat, reproduce the exact FUVVIEW3 date selection and
geometry and establish which active SI13 database generated the dataset.

Across the available examples, all 43 WIC files contain flag 3 and `CIMAGE`,
whereas all 43 SI12 and 43 SI13 files contain flag 1 and no `CIMAGE`.
`CIMAGE` can be negative and is always less than or equal to `IMAGE`; their
difference reaches about 13,390 counts and affects 30--85% of a frame. It is
therefore consistent with an unknown background-subtracted WIC product, but it
does not equal current fuvpy `dgimg` or `shimg`. Do not use it as a defined
calibration product until its producer code is found.

`compare_cimage.py` has tested `CIMAGE` as an alternative WIC input without
changing SI12/SI13, proton correction, footprint mapping, or the comparison
pixels. Fixed-nightside median ratios increase from 199.7 to 216.4 (0085),
146.3 to 164.1 (0086), and 92.0 to 99.1 (0968); above-Frey fractions also
increase in all three cases. `CIMAGE` is therefore not a solution to the high
ratio in the available partial extracts. The diagnostic figures and CSV are
under `figures/debugging/cimage_comparison/`; production remains unchanged.

The newly obtained Frey-event WIC file provides a stronger exact-frame test.
`wic20003021138.idl` matches the orbit-0364 NetCDF geometry and timestamp, and
its `IMAGE` agrees with NetCDF `img` at `r=0.9995`. On the bright polar support
shown in the clean Figure-16 WIC panel, its `CIMAGE` agrees with fuvpy `shimg`
at `r=0.998`; the fit is `CIMAGE = 0.991 shimg + 107 counts` and the median
difference is `+88 counts` (`+4.6%`). The auroral WIC signal is therefore
effectively reproduced by fuvpy for this event, making it unlikely that the
WIC background product explains the excessive WIC/SI13 ratio. See
`scripts/debugging/compare_frey_cimage_event.py` and the associated comparison
figure under `figures/debugging/paper_reconstruction/`.

The same frame now has a direct pre/post-background ratio diagnostic. SI13 is
footprint-mapped onto WIC once, and each paired pixel is compared before
proton correction using WIC `img/shimg` and SI13 `img/dgimg`. Median
`R_after / R_before` is 0.519 over all positive pixels and 0.648 under the
`WIC >= 50`, `SI13 >= 3` guard; the guarded 10--90% range is 0.203--0.887.
The processing-induced scale change is therefore nonuniform and strongly
dependent on SI13 brightness. It generally lowers the ratio for this event,
so it is not an upward gain explanation for the excessive ratios, while the
weak-SI13 tail remains highly sensitive to additive background residuals. See
`scripts/debugging/compare_background_ratio_scale.py` and its figure under
`figures/debugging/paper_reconstruction/`.

An exact sensor decomposition shows why the ratio falls. On the 9,016 guarded
pixels, fuvpy removes median fractions of 50.8% from WIC and 18.8% from SI13;
WIC loses the larger fraction in 99.4% of pixels. This does not by itself imply
WIC over-subtraction because raw WIC contains substantial background/dayglow,
but it rules out SI13 background subtraction as a general upward-rescaling
mechanism in this frame. See
`scripts/debugging/decompose_background_ratio_scale.py` and
`frey_background_fraction_decomposition_or_0364.{png,pdf}`.

The background-code audit found a configuration mismatch. fuvpy's
BS and SH functions use the same equations for WIC, SI12, and SI13; sensor
specificity is supplied mainly by preprocessing and caller parameters. The
fuvpy dayglow-paper workflow uses WIC BS/SH damping `1e-2/1e-4` and SI12/SI13
`1e-1/1e1`. In contrast, icBuilder's current orbit builder passes
`1e-3/1e-4` to every sensor. Product 1 uses WIC `shimg` and SI `dgimg` but
combines DG and SH weights for all sensors.

The SI12 part of that controlled test has been completed over all 247
orbit-0364 frames. Paper-workflow damping changes stored `dgimg` by only
`0.00077 count` median absolute difference (`r=0.9999994` over 2.45 million
pixels). `shimg` changes by `0.360 count` median absolute difference
(`r=0.9955`). Product 1 uses SI12 `dgimg`, so its proton-correction intensity
is effectively unchanged. The corresponding SI13 test gives
`dgimg r=0.9999917`, median absolute difference `0.00233 count`, and a signed
10--90% difference range of `-0.0109--+0.0182 count` over 2.63 million pixels.
Thus the mismatched damping does not explain the high WIC/SI13 ratio through
SI13 intensity on orbit 0364. SI13 `shimg` changes more (`r=0.9718`, median
absolute difference `0.553 count`), so the SH-derived binning-weight effect
remains open. See `scripts/debugging/compare_background_damping.py` and
`si12_background_damping_comparison_or_0364.{png,pdf}` /
`si13_background_damping_comparison_or_0364.{png,pdf}`.

The full-orbit WIC test changes only BS damping from `1e-3` to the paper value
`1e-2`, retaining SH damping `1e-4`. The final WIC `shimg` remains effectively
the same (`r=0.9999971`, median absolute difference `0.156 count`, median
signed difference `+0.0288 count` over 10.44 million pixels). The damping
configuration therefore does not explain the WIC/SI13 ratio through either
intensity on orbit 0364. See
`wic_background_damping_comparison_or_0364.{png,pdf}`.

Code inspection also establishes that fuvpy adds `dampingVal I` directly to
the weighted `G.T @ G`; the parameter is not normalized to the matrix scale.
Both `G` and the data are row-scaled by `w*ws` before forming the normal
equations, so those nominal weights are squared in the least-squares
objective. The effective damping therefore depends on data and matrix scale,
which is a separate numerical issue worth revisiting if the background model
is retained.

The BS scale was quantified on ten common WIC/SI13 orbits spanning 0085--1930.
The median weighted `diag(G.T @ G)` varies `28.4x` for WIC and `12.8x` for
SI13, confirming that fixed damping is not orbit-invariant. It remains tiny on
this typical scale: current damping is `1.0e-8--5.0e-7` of the median diagonal;
paper damping is `1.0e-7--3.0e-6` for WIC and `3.9e-6--5.0e-5` for SI13.
This explains the negligible image changes, while not ruling out effects in
small-eigenvalue directions. See
`scripts/debugging/analyze_background_gtg_scale.py` and
`background_gtg_scale.{csv,png,pdf}`.

Frey does not explicitly say that Figure 16 was smoothed and gives no numerical
count/SNR cutoff. Frey only states that the corrected images were remapped onto
a geomagnetic grid, without documenting its resolution or aggregation.
`recreate_frey_figure16.py` maps SI into native WIC space but then scatters the
native WIC pixels at magnetic coordinates; it does not construct that stated
grid and accepts any positive WIC/SI13 pair. Meurant's separate analysis of the
same event explicitly smooths WIC and SI13, making smoothing a reasonable
sensitivity test but not a confirmed Frey step. Common-grid aggregation,
smoothing, and signal thresholds must therefore be tested and labelled as
hypotheses before attributing the difference to SI13 background or gain.

The legacy propagated retrieval has a `WIC=50` / `SI13=3 count` low-signal
fallback rather than evaluating the ratio when SI13 is below 3 counts. The
Figure 16 reconstruction calls the ratio curve directly and bypasses this rule.
No reviewed paper has yet established those exact numbers as Frey's choices.

That recovered-calibration sensitivity has been carried out on the available
partial extracts of orbits 0085, 0086, and 0968 without changing production
code. It includes the
recovered date selection and post-boom geometry, tests current fuvpy's WIC
mask both as executed and with the apparent intended parentheses, refits the
background independently for each sensor and branch, and uses one fixed
nightside comparison support. SI is mapped onto native WIC detector pixels by
the established Coumans-style footprint-overlap method. Relative to current
fuvpy, applying the recovered WIC+SI13+SI12 treatment changes the median ratio
by `+6.18%`, `+3.17%`, and `-0.91%`; the above-Frey fractions change by
`+2.73`, `+1.17`, and `+0.03` percentage points. It worsens the first two cases
and has little effect on the third partial stack. The available recovered flat-field
conventions therefore do not solve the high-ratio failure.

More specifically, adding recovered SI13 to recovered WIC increases the
90th-percentile ratio and above-Frey fraction in every test orbit; it increases
the median for 0085 and 0086 and changes 0968 only slightly. The subsequent
SI12 reflat is negligible. The ratio-map figures now include full corrected
WIC/SI13 context, a missing/excluded/retained support map, and explicit pixel
counts for the densest fixed-support frame. A zero-based frame can instead be
chosen explicitly with `--plot-frame ORBIT:INDEX`; the regenerated orbit-0968
figure uses `--plot-frame 968:2`, the third of its three available raw example
frames.

The SI-to-WIC mapping defect has been corrected. The Coumans and reflattening
diagnostics now share `icbuilder/detector_coregistration.py`, which represents
SI pixels as quadrilaterals, allocates their values by overlap with WIC pixels,
and applies a 90% coverage threshold. Orbit-0968 frame-2 output is exactly
identical to the original Coumans implementation, and synthetic tests verify
constant-field preservation and coverage rejection. Relative to the former
point interpolation, the corrected mapper lowers the current-fuvpy median by
12.50, 13.25, and 0.67 ratio units for 0085, 0086, and 0968. It changes the
absolute statistics materially but does not reverse the calibration result:
recovered SI13 still raises the 90th-percentile ratio and above-Frey fraction
in all three examples, while SI12 remains negligible.

Separate maps for all three orbit-0968 frames show the same approximately
detector-horizontal ratio band, with recurring row-profile changes around WIC
rows 49--56. The SI flat-field tables are sufficiently jagged to be a credible
cause of weak-signal stripes: orbit-0085 adjacent-row SI13 changes have median,
95th-percentile, and maximum magnitudes of 3.1%, 24.5%, and 69%; SI12 gives
3.7%, 19.7%, and 36%. WIC is much smoother. Orbit 0085 also uses the nearest
available pre-boom profiles from 46 days later for SI13 and 118 days later for
WIC/SI12. This is a real calibration lead, but not yet a demonstrated cause of
the orbit-0968 ratio discontinuity. A direct multiplicative-profile sensitivity
has now divided out the selected recovered SI profile, reapplied either the
original or a median-preserving Gaussian-smoothed profile, and refitted the
background. Smoothing at `sigma=2` and `5` detector rows raises the fixed
nightside median ratio by 4.4%/5.2% for orbit 0085 and 2.2%/6.8% for orbit
0086. Orbit 0968 changes by +3.8%/-0.5%. Above-Frey fractions are unchanged or
worse, except for a negligible 0.3 percentage-point reduction in the strongest
orbit-0968 case. The row plots show altered horizontal structure, not a clean
removal. Simple SI profile smoothing does not correct these partial stacks,
but the result is not orbit-level validation. Production is unchanged.
The profile substitution is exact on positive pixels in these example files:
they use calibration flag 1, have zero `AIRGLOW_SCALE`, and therefore have no
additive term in the recovered FUVVIEW3 calibration equation. Upstream-clipped
zeros cannot be reconstructed and are excluded from the positive fixed
comparison support.

The calibration terminology is now resolved. FUVVIEW3's explicit
`bkg_level` is zero in the example SI files because `AIRGLOW_SCALE=0`, while
fuvpy later estimates the physical constant background `C` from unflattened
WIC dark pixels. That estimate is one scalar per loaded stack, not a prescribed
time series: the partial example stacks give about 411, 417, and 1032 counts for 0085,
0086, and 0968, with a 450-count fallback. The analogous diagnostic SI13
estimates are about 1.88, 1.88, and 5.04 counts. These are inferred constants,
not independent calibration values.

The SI estimate is not yet robust despite the pooled pixel counts. Orbits
0085/0086 have 20 frames but only about 220--320 candidate dark pixels per SI
frame, and 10--24% of the SI13 candidates are upstream-clipped zeros. Orbit
0968 has only three frames but 6,500--8,200 candidates per frame; its SI13
estimate still changes 4.08 -> 5.01 -> 6.42 counts and SI12 changes
5.51 -> 5.98 -> 6.89. This points to systematic mask/emission/clipping
uncertainty rather than insufficient raw pixel count.

The deferred constant test requires a complete native detector-grid orbit but
does not require the original IDLs. SI12/SI13 NetCDF `img` is unchanged from
IDL `IMAGE`. For reflattened WIC `J`, dark pixels obey
`J=C_0+F(C_true-C_0)`; robust regression against the known flat field can thus
recover both the imposed `C_0` and candidate `C_true`. The binned products
cannot support this detector-row test. Estimate `C` per frame from several
independently defined zero-signal masks, prioritizing
off-Earth pixels, and compare the median/mode of `I/F` with robust `I`-versus-
`F` regression. Estimate on held-out detector regions or frames and require the
corrected background to become independent of row and flat-field factor.
Bootstrap whole frames or detector regions. Propagate the resulting interval
through the ratio pipeline only afterward, with DMSP reserved for validation
rather than calibration.

The recovered daily mission-sensitivity correction remains unvalidated and is
distinct from both `C` and the row profile. It derives independent WIC/SI12/
SI13 gains from stellar sensitivity curves and normalizes them to 2000 day
171. The resulting WIC/SI13 multiplier is 0.703 for orbit 0085, 0.802 for 0086,
and 1.456 for 0968. Thus it lowers the ratio in the two weak pre-boom examples
but raises the Coumans-event ratio by about 46%. The exact correction table's
producer compatibility and its consistency with the Frey response calibration
are now the most relevant scalar-calibration questions. The median-preserving
flat-field smoothing did not test this overall gain.

Use the established notation consistently: `B_s(t)` is the multiplicative
stellar/mission correction and `C_s` is the additive detector-background
constant in the Ohma reflat equation. `B_WIC/B_SI13` scales the image ratio.
The complete-orbit dark-pixel experiment estimates `C`; it cannot independently
recover `B_s(t)`.

The calibration history is now plotted in
`flatfield_profiles_over_time.png`, with consecutive-profile differences in
`flatfield_profile_differences.csv`. FUVVIEW3 selects the nearest profile per
sensor on the same side of the boom loss and performs no temporal
interpolation. Exact ties choose the earlier database entry. Orbit 0085 has no
earlier pre-boom calibration and therefore uses future profiles 46 days later
for SI13 and 118 days later for WIC/SI12. Consecutive post-boom WIC profiles
typically differ by about 1% per row, while SI profiles can differ by 10--13%
in the median and much more in individual rows. Frey et al. (2003) explain the
contrast: WIC dayglow-derived deviations were fitted by a smooth parabola,
whereas SI deviations were retained directly in lookup tables. Coumans et al.
(2004) give no evidence of a custom flat field. A full recalibration requires
unflatfielded data and many modeled dayglow exposures; the present IDLs cannot
provide an independent calibration because their `IMAGE` fields are already
flatfielded.

The experiment also confirms a precedence error in fuvpy's WIC background
mask: the missing-SZA clause does not contribute as apparently intended. This
can materially change the fitted detector constant, but independent background
refitting leaves the resulting ratio statistics almost unchanged. It is not a
plausible root cause. The test remains a recovered-FUVVIEW3 sensitivity rather
than an exact producer reconstruction because the recovered source does not
define local WIC flag 3 or `CIMAGE`, WIC and SI examples have different SAVE
provenance, and the post-boom Python geometry has not been checked against an
IDL/GDL golden output.

Diagnostic outputs are under
`figures/debugging/coumans_histogram_background/` and
`figures/debugging/large_ratio_corpus/`, and the new recovered-calibration
outputs are under `figures/debugging/fuview3_reflattening/`. Product 2 and the
production fuvpy processing remain unchanged. Another ordinary reflat patch is
not a promising next step. Exact calibration work would require producer
lineage or a golden output from the actual IDL processing version, especially
for the relative WIC/SI13 gain, flag 3, and `CIMAGE`. Until then, out-of-domain
ratios should be reported as retrieval failures rather than interpreted as
high energies. This reinforces the existing decision to use an alternative E0
source rather than the current IMAGE-ratio retrieval.

The additive `C_s` investigation is now limited to diagnosing detector-row
stripes, not repairing the global ratio scale. WIC has already undergone the
Ohma reflat, and a positive `C_SI13` with `F > 1` would lower much of SI13 and
tend to increase WIC/SI13. The unresolved global discrepancy is therefore best
described as relative WIC/SI13 calibration or compatibility with the Frey
response calibration. `B_WIC/B_SI13` is the only identified mechanism with
direct coherent multiplicative leverage, but the available evidence does not
show that the recovered stellar calibration is wrong.

One final external calibration test is now available. Frey et al. (2003),
Figure 16, provides simultaneous numerically scaled SI12, WIC, and SI13 maps
for orbit 0364 on 28 October 2000 at 11:38 UT; Meurant et al. (2003), Figures 7
and 9, provide WIC-count and corrected-ratio keograms for 10:04--13:00 UT on
the same day. The current `icPhysics` ratio path is only a static interpolation
of the Frey tables and omits the activity-dependent atmospheric correction that
Frey explicitly applied when disturbed composition suppressed SI13 and created
unreasonably high inferred energies. Reproduce these published camera scales
before definitively closing the IMAGE-ratio path.

### Portfolio impact

- Central update needed: Yes.
- The publication-blocking ratio defect is quantified over the complete corpus,
  and the recovered FUVVIEW3 flat-field hypothesis has now been tested without
  resolving it. Exact producer provenance remains the only route to a stricter
  calibration reconstruction; the alternative-E0 path remains justified.
  Project priority and deadlines are unchanged.

## Previous checkpoint: unweighted pointwise Bayesian ratio inversion

The exploratory Bayesian calculation was simplified to one readable script:

```bash
python scripts/dmsp/bayesian_ratio.py
```

The script loads manually accepted F12--F15 frames from the 3.5-GB
`matches.nc`; the +/-60-s IMAGE--DMSP association was already applied when
that file was generated. It keeps finite positive DMSP energy and flux,
requires `dE/E <= 0.25`, `dQ/Q <= 0.20`, positive WIC, corrected `SI13 > 3`,
and restricts the diagnostic to `E <= 5 keV` and `R <= 150`. This leaves
44,861 measurements.

The scientific calculation is direct and uses linear bins with no smoothing:

1. count the empirical prior `p(E)`;
2. row-normalize the joint `E--R` histogram to obtain `p(R | E)`;
3. calculate `p(E | R) proportional to p(R | E) p(E)` and normalize each
   ratio bin;
4. plot the raw modal-bin maxima as dots through `p(R | E)` and `p(E | R)`,
   masking empty bins, and overlay a centred five-bin rolling mean;
5. independently obtain `p(E | R)` by column-normalizing the joint histogram;
6. plot the four distributions in
   `figures/debugging/dmsp_bayesian_ratio.png`.

The Bayesian and direct-histogram posteriors agree to `5.55e-17`. Their modal
trace remains near 0.25 keV over much of the ratio range because it is the
maximum bin of a skewed conditional distribution. Earlier binned relations
used medians or other quantiles and therefore need not follow the mode. Binwise
medians are now shown as dashed cyan lines on every conditional panel; the
posterior median is generally about 0.7--1.5 keV and varies with ratio.

The prior itself is a selected and truncated distribution, not a general
auroral climatology. Its modal 0.25-keV bin contains 9.7% of rows; median and
mean energy are 1.12 and 1.44 keV. Equal-frame and equal-orbit medians remain
near 1 keV, while DMSP flux-quartile medians increase from 0.36 to 2.16 keV.
The modal bin moves to about 1.69 keV if equal-width log-energy bins are used,
so the mode must not be promoted to a characteristic energy.

Every measurement has equal weight and is treated as independent, as requested.
The visible isolated high-energy likelihood cells are supported by very few
measurements and are not a stable response. No cross-validation, score
sensitivity, frame/orbit weighting, or log-energy calculation remains in this
first diagnostic. The earlier three-script exploratory framework was removed;
its generated files may still exist locally but are superseded.

### Portfolio impact

- Central update needed: Yes.
- This is a descriptive view of the pointwise conditional distributions, not
  evidence of predictive skill. Do not change Product 2 on this basis.

## Historical checkpoint: crossing-averaged Bayesian ratio inversion

This older orbit-held-out result is retained here only to explain the change
to pointwise measurements. Its overextended Bayesian scripts were removed when
the diagnostic was consolidated into `scripts/dmsp/bayesian_ratio.py`.

The analysis uses the outcome-independent `image_ratio` paired with
`response_dmsp_energy`; it deliberately avoids the target-conditioned
`frey_image_ratio` subset. Complete IMAGE orbits are assigned to five
year-balanced folds, each training orbit has equal total weight, and 66-s and
120-s support are evaluated independently. The primary standard histogram is
checked against coarse, fine, wide-ratio, and fixed six-energy-by-five-ratio
models shrunk toward the empirical energy prior with strengths 0.5, 2, and 5.
Bayesian normalization is checked numerically in every fold.

The 66-s sample contains 251 rows from 113 orbits. Relative to the empirical
prior, the standard ratio posterior changes held-out log score by -0.1158
nats (95% conditional orbit-resampling interval -0.2097 to -0.0333), CRPS by
-0.0052 dex (-0.0148 to 0.0039), and posterior-median absolute log error by
only +0.0024 dex. Positive values mean improvement. The 120-s sample contains
233 rows from 102 orbits and likewise shows no established improvement. The
sparse/shrunk sensitivities do not reverse the decision. Standard-posterior
50%, 80%, and 90% coverage also falls relative to the prior at both supports.

The calculation therefore answers the mathematical question but not in the
desired scientific direction: the observed WIC/SI13 ratio does not improve
held-out mean-energy prediction beyond the empirical DMSP prior in this
selected sample. This is stronger than noting that two conditional plots are
not inverses, but it is not a mission-wide result. The prior is learned from
manually accepted F13/F15 crossings in 2000--2001, measurement-error
uncertainty is not propagated, and the confidence intervals condition on the
fixed cross-validated predictions rather than refitting every bootstrap draw.

Important scope correction: the 251-row 66-s analysis is a crossing-level
test. The compact input starts with 700 accepted frame--satellite pairs, only
392 of which have any finite common response, and the `response_coverage >=
0.5` rule leaves 251. Those are track-averaged WIC, SI13, ratio, and DMSP
energy values. This was a defensible independence unit for comparing one
energy model per crossing, but it discards most of the spatially matched
one-second samples and can wash out a pointwise `R--E` relationship. Do not
interpret the current negative result as rejecting a pixel-level Bayesian
retrieval. Before that decision, repeat the likelihood analysis with valid
pointwise matches, unique nearest-frame assignment, orbit-held-out splitting,
and frame/orbit weights so correlated seconds add spatial information without
being counted as independent experiments.

Reusable tables are under `data/dmsp_energy_validation/`. The explanatory
figures are `bayesian_ratio_components.{png,pdf}` (prior, likelihood, and
posterior), `bayesian_ratio_posterior.{png,pdf}` (posterior structure and Frey
comparison), and `bayesian_ratio_validation.{png,pdf}` (proper-score gains and
coverage), all under
`figures/debugging/dmsp_energy_validation/`.

### Portfolio impact

- Central update needed: Yes.
- The Bayesian inversion is feasible, but the IMAGE ratio fails the held-out
  predictive gate. Do not add it to Product 2. The remaining research decision
  is whether to use an explicit empirical energy baseline, continue with
  Zhang--Paxton for physical/coverage reasons despite weak pooled validation,
  or first investigate the unresolved cross-camera preprocessing/calibration
  failure.

## Previous checkpoint: DMSP energy decision analysis

The requested constant/Frey/Zhang--Paxton plan is implemented and verified as
standalone analysis only. No file under `icbuilder/`, no production orbit
script, and no icPhysics routine was changed. Run in this order:

```bash
python scripts/dmsp/build_energy_validation.py
python scripts/dmsp/validate_energy_models.py
python scripts/dmsp/validate_frey_response.py
python scripts/dmsp/plot_energy_validation.py
```

The builder streams only manually accepted F13/F15 rows out of the 3.5-GB
`data/matches.nc`, then makes one record per IMAGE frame, satellite, and
temporal support. The 66-s and 120-s supports are half-open centered windows
and are analyzed separately. DMSP energy is calculated from resolution-matched
energy and number-flux proxies as `sum(Q) / sum(Q / E)`. Camera-specific DMSP
moments use exactly those seconds with valid WIC and corrected SI13 support.
The forward predictors integrate `Q * response(E)` per DMSP second and enforce
the 0.2--25-keV Frey table domain before averaging; they do not apply the
response nonlinearity to one already-averaged energy.
The product stores definitive Kp matched at the IMAGE time and a
grid-independent collapsed Zhang--Paxton value interpolated from a 0.05-h MLT
grid; it deliberately does not use the obsolete 36x36 lookup.

Primary 66-s conclusions:

- 531 energy records from 173 orbits enter the energy comparison; 227 records
  from 107 orbits enter the response-supported Frey test.
- A physical constant gives 1.36 keV. Its orbit-weighted 50% range is
  0.68--2.13 keV and its 90% range is 0.23--4.27 keV.
- The Frey shape gives small WIC improvement and clear SI13 improvement over
  flux alone, but no WIC/SI13 improvement. The ratio gain is fitted directly;
  its improvement is -0.0080 dex with 95% conditional orbit-resampling
  interval -0.0208--0.0050.
- Raw collapsed Zhang--Paxton is worse than the log-median constant. Scaling
  its amplitude by about 0.61 gives only 0.0056 dex held-out improvement, with
  95% conditional orbit-resampling interval -0.0096--0.0207. A free log-space
  slope is also not established.
- The pooled result is not universal across precipitation strength. In the
  highest DMSP energy-flux quartile, raw, scaled, and calibrated Zhang--Paxton
  improve by 0.1185, 0.0344, and 0.0220 dex at 66-s support, with conditional
  intervals above zero. The lower three quartiles show no fitted-model gain.
  At 120 s only the raw high-flux improvement remains established. The
  physical constant also improves in the highest quartile, so this does not
  isolate Zhang--Paxton's Kp/MLT structure from a general
  strong-precipitation energy shift. This is a diagnostic stratification by
  held-out DMSP flux, not yet a usable production gate.
- Year and satellite transfer scores are saved, but this annotated sample
  covers only 2000--2001 and is not a production climatology. No extra DZA/MLT
  cut is applied after manual annotation. The predictive ranges describe
  empirical residual scatter; per-second DMSP fractional uncertainties are
  quality cuts rather than a propagated measurement-error model.

Reusable outputs are under `data/dmsp_energy_validation/`. The key figures are
`figures/debugging/dmsp_energy_validation/frey_forward_closure.{png,pdf}` and
`energy_model_comparison.{png,pdf}`. The new
`energy_flux_strata.{png,pdf}` figure and corresponding CSV tables expose the
strong-precipitation exception. The production decision remains open: test an
IMAGE-observable strong-precipitation gate, enlarge the accepted DMSP support,
or accept the constant as the honest pooled benchmark. Do not introduce any
of these models into Product 2 without an explicit user decision.

### Portfolio impact

- Central update needed: Yes.
- The user authorized and completed the DMSP-calibrated Zhang--Paxton analysis,
  but not its production adoption. It does not beat a constant across the
  pooled sample. The highest DMSP-flux quartile has a repeatable energy-level
  shift, but the present test does not show that Zhang--Paxton structure rather
  than precipitation strength explains it.

## Project state

The experimental 100-km orbit-0085 and orbit-0968 products reveal that viewing
coverage and source sampling differ far more than the legacy downstream
quality fields indicate. For frame 000, WIC occupied-grid coverage is 73.5%
versus 93.5%, SI12 is 40.1% versus 73.6%, and final three-sensor image-ratio
support is 22.5% versus 70.7%. Orbit 0085's occupied SI cells are almost all
single-pixel bins, whereas orbit 0968 contains many two-pixel cells. The
background-correction weight is nearly unchanged between these cases and does
not include sample density; the bin-first Product 2 drops all three native
count fields. In the current experimental single-pixel branch, `sigma = 0` also
means missing within-bin replication is represented as perfect precision.
Do not solve this by changing the stored grid per frame: the VAE/covariance
workflow needs a fixed spatial basis. The next scientific design decision is
how to retain per-sensor support and estimate single-cell uncertainty, and
whether a fixed-grid adaptive smoothing scale is justified.

Do not attempt to repair the sparse-cell problem by substituting another
small-sample statistic. At one sample there is no within-cell variance
estimate, and at two samples the median is merely the midpoint. Moreover, the
current within-cell spread conflates count noise with true unresolved auroral
gradients. Treat the current 100-km result as an oversampled remapping. The
next comparison should test a coarser common effective resolution and/or a
footprint- or kernel-based remap that propagates the known uncertainty of the
individual detector measurements.

The standard multi-resolution ratio workflow is now the leading proposed
replacement: geometrically register the channels, degrade each to a common
PSF/projected footprint no sharper than the worst required channel,
area-resample the matched surface-brightness fields to the fixed CS grid, and
form the ratio only on adequate common support. The present SI-to-WIC
interpolation performs registration/resampling but not resolution matching.
The next investigation should establish or approximate the projected spatial
responses of WIC, SI12, and SI13, including their viewing-angle dependence.

Raw-IDL inspection confirms that this investigation can begin without a new
ephemeris source. The records retain spacecraft position, attitude/look
vectors, camera orientation, FOV, angular sampling, per-pixel GLAT/GLON and
DZA, and `EMIS_HGT = 130 km`; fuvpy currently keeps only the already-derived
pixel centres and angles. Frame 000 of orbit 0085 is sampled from about 7.88
Earth radii, giving approximate nadir footprints of 51 km (WIC) and 98--104 km
(SI), whereas orbit 0968 is at about 5.83 Earth radii and gives 36 km and
69--73 km respectively. This provides a geometric explanation for the denser
orbit-0968 bins. Before implementation, verify the original metadata field
definitions and whether the 128-by-128 SI arrays are 2-by-2 binned from the
reported 256-pixel angular sampling. A centre-derived polygon remap is a
viable independent approximation if the exact camera model cannot be found.

The likely practical route no longer depends on recovering every raw IDL.
Before binning, derive each footprint from the local spacing and orientation
of the per-pixel `glat`/`glon` centre lattice already stored in the fuvpy orbit
NetCDF. Construct corners from half row/column differences on the 130-km
shell, use one-sided differences at valid edges, and compare the inferred area
against DZA. This captures the range difference that a DZA-only scaling misses.
Use an assumed DZA-scaled top-hat footprint only where the centre lattice lacks
enough neighbours.

`scripts/debugging/test_pixel_footprint_regridding.py` implements and verifies
that experiment without modifying products. It reads processed fuvpy frames,
infers four corners from local row/column centre differences, clips the
polygons against every intersected CS cell, and writes PNG/PDF diagnostics to
`figures/debugging/pixel_footprints`. The footprint lattice is internally
consistent: interior coverage is approximately one with less than 5% overlap
error everywhere tested. Orbit 0085 has median square-equivalent WIC/SI13
footprint widths of 58.6/114.2 km, versus 39.9/76.2 km in orbit 0968. On the
current 200-km target grid, footprint overlap increases orbit-0085 SI13 support
from 71.4% centre occupancy to 77.3% any coverage. This verifies the geometry
concept but not yet the radiometric response, uncertainty propagation, or
common-resolution ratio.

The inferred uniform footprint response is now implemented in `BinnedImage`.
`icbuilder/footprints.py` contains the corner inference, polygon clipping, and
sparse overlap application. Product 1 stores overlap-weighted values,
intersecting-footprint counts, and valid fractional coverage. The former
centre-bin median remains available as `binning_method="centre"`; the orbit
CLI, NetCDF descriptor, atomic-save check, and restart check all preserve that
choice. `sigma` remains explicitly provisional so its uncertainty semantics
can be redesigned separately. Common-resolution WIC/SI treatment remains a
Product-2 task.

The current Halley Product-1 run takes about 225 s per orbit. No profiler has
yet been run, but `overlap_mapping` is the strong suspected bottleneck: it
loops over valid detector footprints and candidate grid cells and performs
four-boundary polygon clipping for every overlap. Sparse matrix application
to the image, geometry, and weight fields occurs only after that mapping is
built and is unlikely to dominate.

`scripts/make_orbit_h5_files.py` no longer contains the obsolete DTU-local
absolute path. It accepts `--base`, defaults to the repository `example_data`,
discovers the three IDL folders with `pathlib`, assigns files using
`orbitdates.csv`, and rewrites the three historical `*files.h5` indices while
reporting file and orbit counts.

The first modular product is now implemented. `PreImage` carries a validated
canonical sensor name and `BinnedImage` requires an explicit timestamp for
each frame. Its dedicated NetCDF writer stores sensor-native binned signal,
uncertainty, weights, source counts, geometry, corrections, time, subsolar
longitude, and grid metadata. It contains no Kp, precipitation, proton-energy,
or conductance variables. `BinnedImage` is now native-grid only: interpolation
and `target_grid` were removed, counts are integer-valued, and NetCDF includes
the actual xi, eta, MLAT, and MLT arrays.

`scripts/make_binned_orbit_files.py` now processes WIC, SI12, and SI13
independently instead of copying the three-camera conductance workflow. It
writes native-grid files under `binned/wic`, `binned/si12`, and
`binned/si13`, retains sensor-specific time support, skips files that already
exist, and publishes each new file atomically through a temporary `.partial`
path. Orbit discovery now scans the actual NetCDF files in each sensor
directory and no longer depends on the stale `*_avail_orbit.npy` manifests.
The legacy conductance script was updated for the explicit sensor/time API.

`PrecipitationImage` now defines the initial Product-2 boundary. Image ratio
matches WIC/SI12/SI13 one-to-one within two seconds. Zhang--Paxton matches WIC
and SI12, then attaches SI13 without reducing that time support; missing SI13
uses source index `-1` and NaN arrays. Two explicit preparation functions keep
these method rules separate. The verified regular xi/eta grids use
explicit four-cell bilinear interpolation instead of Delaunay triangulation.
All four surrounding SI cells must be finite. Targets in the outer half of a
physical SI boundary cell inherit the nearest SI value; targets beyond the
physical grid and internal missing-data gaps remain NaN. Variance uses squared
bilinear weights internally and the source uncertainty at the boundary.
Zhang--Paxton
requires WIC/SI12 only, while the ratio method requires all three sensors.
The fixed grid mapping is calculated once per SI sensor and reused for its
signal/uncertainty and quality-weight interpolation; SI12 and SI13 mappings
remain separate. The interpolation function accepts one array at a time;
uncertainty propagation is selected explicitly and signal/uncertainty support
is not forced to match at this stage. Native sample counts are never
interpolated or saved in Product 2.

Both methods now carry regridded SI13, uncertainty, weight, and the
SI12-proton-corrected SI13 quantities. Zhang--Paxton accepts an absent SI13
orbit; SI13 remains diagnostic and cannot change its E0 or WIC-derived Fe.
Orbit-0085 smoke tests produced 20 frames for both methods with SI13 present.
A Zhang--Paxton run without an SI13 file retained all 20 frames with index
`-1` and NaN raw/corrected SI13.

Proton correction is now independent of method-specific precipitation
inference. `PrecipitationImage` calls the shared SI12 correction once, stores
the corrected WIC/SI13 fields, and passes corrected counts to either the ratio
or Zhang--Paxton routine. The corrected-count ratio entry point reproduces the
legacy combined calculation in a frozen reference case; no zero-SI12 shortcut
is used. Shared SI12-induced WIC/SI13 covariance remains deferred.

The combined sensor weight is now finalized in Product 2. Image ratio
multiplies WIC, SI12, and SI13 weights; Zhang--Paxton multiplies WIC and SI12
because SI13 is diagnostic. Product 3 carries this weight unchanged.

The user will choose explicit dataset directories such as `ZP_P2` to compare
precipitation methods and proton energies. Do not force method-specific
subdirectories. Instead, restart checks must validate that every existing file
matches the requested method, proton settings, schema, and source stage before
skipping it; a mismatch should fail clearly rather than silently mix datasets.

The short modular fix list is complete. Product 3 carries `ssalon` and the
reader reconstructs time-dependent magnetic longitude. Product-1 and Product-2
restart checks validate structure and requested configuration; Product 3 also
checks agreement with its precipitation input. The stale weight and output
directory tests were updated. The focused icBuilder suite passes 62 tests and
the icReader loader suite passes 10 tests. Existing Product-3 files made before
this change lack `ssalon` and must be regenerated.

Product-2 processing exposes `--proton-method`, `--proton-energy`, and
`--proton-energy-uncertainty`. The startup summary reports the selected `Ep`
and `dEp`; Product 2 serializes them, Product 3 carries them unchanged, and
both modular readers expose them.

The SI-grid boundary loss is fixed. Bilinear interpolation remains the general
interior rule, while target centres between the outer SI centres and physical
SI edges use the nearest edge cell. On the canonical grids this fills exactly
the missing 140-cell WIC boundary ring and maps a complete 18-by-18 input to
all 1,296 WIC cells. It does not fill internal source-data gaps or targets
outside the physical SI domain. A live orbit-0085 calculation recovered 1,900
finite SI12 boundary values that were all NaN in the saved centre-only file;
remaining boundary NaNs correspond to missing SI source cells. Restart
validation now includes the regridding rule, so old Product-2 files fail as
configuration mismatches instead of being silently skipped. Near-zero
Robinson uncertainty, possible count-noise double counting, ratio covariance,
and detailed publication provenance are also tabled for later review.

**Confirmed by the user:** near-zero Robinson uncertainty does not block the
current regeneration. Product-3 `dP` and `dH` are provisional and Product 3
may be regenerated after the analytic uncertainty method is corrected. The
icAnalyzer uncertainty workflow will instead use Monte Carlo sampling of
upstream inputs and repeated nonlinear conductance conversion.

`ConductanceImage` and `make_conductance_orbit_files.py` now form a compact
Product-3 pipeline. They load precipitation products, call the array-based
icPhysics Robinson model, and save the precipitation state, Hall/Pedersen
conductance and uncertainty, Kp, weight, grid, precipitation/proton choices,
model provenance, and source filename. The old preprocessing calculations
were removed from this stage. A temporary orbit-0085 run produced a
20-by-36-by-36 conductance file successfully.

`scripts/make_conductance_figures.py` now uses the public `icReader.load()`
interface for modular products. It follows each conductance file's recorded
precipitation source and adds proton-corrected WIC, proton-corrected SI13, and
their ratio to the Product-3 precipitation state, conductance, uncertainty,
covariance, weight, and Hall/Pedersen panels. It also follows the precipitation
product's recorded Product-1 WIC source and plots WIC detector zenith angle on
valid WIC/SI13 support on a fixed 0--75 degree scale. This avoids displaying
high-DZA geometry in cells where no ratio exists. Input/output directories and
orbit/frame selection are command-line options. A scratch render of orbit
0085, frame 0 verified the 4-by-4 layout before this display-only mask change.

The regenerated corpus is mounted locally from Halley. Product 1 contains
1,687 WIC, 1,694 SI12, and 1,693 SI13 orbit files. The Zhang--Paxton,
2-keV-proton run is complete with 1,686 precipitation and 1,686 conductance
files. The image-ratio precipitation run remains active. A fixed 188-orbit
snapshot, spanning orbits 0085--0291, was used for the DZA sensitivity test;
no image-ratio conductance directory existed at inspection time.

The class contains no proton, energy, flux, covariance, Hall, or Pedersen
equations. It selects the crude icPhysics ratio or Zhang--Paxton function by
default, while explicit function injection remains available. Its NetCDF
writer stores the prepared fields, returned precipitation quantities, target
grid, Kp, source indices, and provenance. The new and old Zhang--Paxton
calculations match exactly for all saved precipitation fields over complete
example orbit 0085.

`PrecipitationImage` now also accepts binned-product filenames and loads them
through `icreader.load()`. Omitted Kp defaults to the bundled local definitive
GFZ series; explicit loaded images and a preloaded `kp_series` remain supported
for efficient bulk processing. Filename/default-Kp construction and explicit
object/Kp construction produced exactly identical Product-2 fields for example
orbit 0085. Ten focused Product-2 tests pass.

All three modular writers now declare a root `product_type` and
`schema_version = 1`: `binned_fuv`, `precipitation`, and `conductance`. The
focused writer/schema tests pass (19 tests). The actual Product-1 writer also
passed a temporary end-to-end round trip through `icReader.load()` on the
nested SI grid. The precipitation and conductance reader branches remain
explicit placeholders until those builder products are complete.

The legacy conductance orbit script is temporarily not runnable because its
SI12 and SI13 constructors still pass the removed `target_grid` argument. It
must later consume the explicit precipitation product. A Product-2 orbit
reader also remains to be implemented.

`scripts/make_precipitation_image_orbit_files.py` is now the Product-2 orbit
builder. It discovers the sensor-specific Product-1 NetCDF files directly,
uses the correct method-specific sensor intersection, loads Kp once before
orbit processing, and writes to a user-selected output folder. Existing files
are skipped only after their method, proton settings, schema, arrays, and grid
have been validated. Mismatched configurations stop with a clear error;
`--overwrite` deliberately replaces them. New files are atomically published
only after the partial file passes the same validation.

`icBuilder` now contains an importable Zhang–Paxton latitude collapse and a
complete fixed-grid lookup. The 8.4-MB bundled NetCDF has dimensions
`(901 Kp, 36 eta, 36 xi)` and stores direct layers from Kp 0.00 through 9.00
at 0.01 spacing. The loader rounds Kp to the nearest hundredth and directly
indexes E0, profile-derived dE0, and the area-weighted median without runtime
interpolation. It directly compares the stored two-dimensional coordinates
with the active Cubed-Sphere grid.

The lookup implementation has been reduced from production-style
infrastructure to transparent research code. The loader is a 79-line module
with one public function returning a dictionary. The 150-line generation
script shows the serial Kp loop, optional `process_map`, collapse, and NetCDF
writing in execution order. The table stores only the three scientific
fields, coordinates, units, collapse settings, and ZhangPaxton2008 package
version.

The latitude-collapse implementation is now separated in the same style. The
reusable module is 306 lines of direct NumPy calculation (down from a
1,065-line calculation/plotting/CLI module), returns ordinary dictionaries,
and contains no dataclasses or array type machinery. The 520-line diagnostic
script retains the necessarily verbose construction of four explanatory
figures and the CLI. This is deliberate plotting content rather than
scientific-framework complexity.

The orbit pipeline now loads the bundled definitive GFZ Kp series once before
multiprocessing and matches each final retained frame to the enclosing
half-open three-hour interval. `ConductanceImage` requests all lookup layers
once per orbit, checks their shape and grid, verifies that every frame remains
inside its serialized interval, and supplies E0/dE0 to the count conversion.
Missing, gapped, or out-of-range Kp fails explicitly. The local GFZ series now
contains 10,464 definitive values through 2003-07-31, one month beyond the
last June 2003 IMAGE orbit. Its actual SHA-256 is recorded in provenance
without a hard-coded checksum gate.

The orbit pipeline now resumes by default. It skips structurally valid current
products, reruns missing or invalid outputs, and provides `--overwrite` for a
deliberate full rebuild. Workers write and validate `or_XXXX.nc.partial` before
atomically publishing `or_XXXX.nc`; the final file is the completion record.

The paired E0/dE0 override bypasses the WIC/SI13 energy inversion while keeping
SI12 proton correction, WIC-derived Fe/dFe, and R/dR diagnostics. The induced
first-order E0--Fe covariance is passed to both Robinson uncertainty functions
and stored. Original GFZ thirds, rounded lookup Kp, interval starts, source
metadata, lookup/collapse provenance, and dE0 interpretation are serialized.
At `Fe=0`, dP and dH are the one-sided conductance excursions to `Fe=dFe`;
this replaces an invalid minimum-uncertainty calculation that failed at and
above 4 keV.

The modular Zhang--Paxton product now uses WIC/SI12 frame support and carries
SI13 only where available. SI13 changes do not affect E0 or Fe under the
override.

The Zhang–Paxton implementation did not read or modify IMAGE data or tracked
`example_data/` products. The user subsequently fixed the independent
SI-grid construction defect and regenerated orbit 0085/0086 products and
figures.

`BinnedImage` has been optimized without changing its scientific calculation.
It stably groups source pixels once by populated grid cell, evaluates the
Student-t and chi-square uncertainty multipliers once per distinct sample
count, and shares target-grid triangulations only among fields with identical
non-NaN source masks. Fields with different missing-data support remain
independent. Comments in the calculation explain the historical count, NaN,
weight, and geometry semantics that must be preserved.

`ConductanceImage` now calculates the combined weight and the six
Ep/dEp-dependent proton response values once per orbit, and the immutable
camera response tables use persistent interpolators. The Zhang--Paxton
production path applies the same count correction, flux, covariance, and
Robinson equations to masked float64 arrays. Its zero-flux mask preserves the
one-sided dP/dH definition without evaluating singular derivatives. The
historical `image_ratio` comparison remains scalar and reuses the same cached
proton response.

A subsequent full read-only audit did read the five tracked example
conductance products and the representative raw/intermediate metadata. It did
not regenerate or modify them. The audit establishes that the current
pipeline and products are not publication-ready. The Frey table transcription
is correct, but most stored E0 values are hard fallbacks or saturation rather
than interior ratio retrievals. See
`vault/02_Algorithm/Audit - 2026-07-29.md`.

## Implemented definition

For each latitude slice:

1. evaluate Zhang–Paxton E0 (keV) and Q (mW m-2) in cell centres;
2. find the principal Q maximum;
3. retain the contiguous cells around that maximum where Q exceeds the chosen
   threshold;
4. calculate the area-weighted mean and median E0 over those cells using exact
   spherical weights `sin(latitude_upper) - sin(latitude_lower)`.

The reusable API supports broadcastable Kp/MLT inputs. It returns the
area-weighted mean and median, weighted spread, Q-weighted sensitivity mean,
selected latitude bounds and area, threshold cutoff, peak diagnostics, empty
masks, possible sampling-limit contact, and physical-pole contact.

The provisional default is `Q > 0.05 mW m-2`, based on the Zhang–Paxton Figure
8 mean-energy criterion. The implementation and figures also compare the
absolute `Q > 0.25 mW m-2` boundary criterion. **Confirmed by the user:** the
relative 10%-of-peak threshold is removed, the area mean remains primary, and
the area median is calculated and displayed alongside it.

The collapse evaluates the exact MLT supplied by the caller and does not
average within an MLT bin. Diagnostic maps use 0.05-hour sampling. The
published model is continuous in MLT through its Fourier representation,
although it was fitted from empirical 0.5-hour sectors.

MLAT cells now use 0.01-degree spacing over 50–90 degrees. This fine grid
numerically resolves the continuous Epstein profiles and prevents the hard
Q-threshold boundary from producing the visible jitter present at
0.25-degree spacing. It does not imply 0.01-degree empirical accuracy.

## Verification

- Footprint geometry, NaN support, Product-1 schema, rollback centre binning,
  and reader tests pass: 17 focused icBuilder tests and 10 icReader tests.
- One real orbit-0085 WIC frame completed on the current 46-by-46 grid with
  1,655 finite cells; coverage stayed in `[0, 1]`. A temporary Product-1
  footprint file round-tripped through `icreader.load()` with matching signal
  and coverage.
- The full icBuilder suite has 61 passing and 5 failing tests. The five
  failures are the known current-grid mismatch: tests and the bundled
  Zhang--Paxton lookup expect 36-by-36, while the uncommitted grid setting is
  46-by-46. They are not footprint-binning failures.

- The conductance-orbit CLI now accepts independent WIC, SI12, SI13, and
  output folder names beneath `--base`, while retaining the original defaults.
  AST parsing, `--help`, and direct assertions for Chapman's `=wic`, `=s12`,
  `=s13`, and `=conductance` paths pass without running the data pipeline.
- The lookup contains 901 direct Kp layers and 1,167,696 grid values per
  field. It has no empty selections. E0 spans 0.211--4.383 keV and dE0 spans
  0.212--2.233 keV.
- The simplified NetCDF was populated from the previous verified table rather
  than recomputing all 901 layers. The Kp, xi, eta, MLT, E0, dE0, and median
  E0 arrays are exactly equal before and after the schema rewrite.
- One 36-by-36 Kp layer took 10.67 seconds and 259 MB RSS. A full run was
  estimated at 2.7 hours serial or roughly 40--60 minutes with four workers,
  which did not justify recomputing scientifically unchanged values.
- Fresh direct collapses checked nine cells spanning Kp 0.00, 1.52, and 9.00.
  The maximum difference across E0, dE0, and median E0 was `2.4e-7 keV`, as
  expected from float32 table storage.
- All 45 focused tests pass. They cover canonical grid shape/nesting,
  two-dimensional MLT, nearest-hundredth Kp quantization, lookup shape,
  scalar/vector access, direct-collapse agreement, definitive-Kp integrity and
  boundary matching and checksum, paired E0 override, SI13 invariance, induced
  covariance, zero-flux Robinson uncertainty, binned geometry, and NetCDF
  provenance.
- The refactored collapse is numerically identical to the previous
  implementation for synthetic disconnected and empty intervals, strict
  threshold boundaries, non-finite inputs, equatorward/non-polar/polar domain
  flags, scalar inputs, broad broadcast Kp/MLT grids, and both absolute
  thresholds. The comparison used exact tolerances (`rtol=0`, `atol=0`) with
  aligned NaNs.
- All four collapse diagnostics regenerate as valid one-page PNG/PDF pairs
  and were visually inspected after the refactor. Their 4,800-slice summary
  is unchanged: 0 empty, 1,001 equatorward-limit contacts, 0 non-polar upper
  contacts, and 1,378 physical-pole contacts.
- `icbuilder/data/zhang_paxton_e0_lookup.nc` remains byte-for-byte unchanged
  (SHA-256 `53c06e2c2ebb5ac8c605185aab7bba65575630be5339bb0d7e1c789262d07e65`);
  no lookup regeneration is required.
- A serial scratch run of orbit 0085 completed in 21 seconds without touching
  tracked examples. The new and old products both contain 20 frames. Original
  Kp 2.667/4.333 map to lookup Kp 2.67/4.33; all 18,958 finite E0 pixels
  exactly match the selected lookup, and their covariance values are finite.
- `figures/zhang_paxton_lookup_kp1_52.png` and PDF show the MLT coordinate,
  E0, and dE0 in native `(xi, eta)` coordinates and were visually inspected.
- The table uses `Q > 0.05 mW m-2` over 50--90 degrees with 0.01-degree
  numerical MLAT cells. 231,123 table cells touch the 50-degree equatorward
  sampling limit; changing this provisional domain requires regeneration.

- All Python files under `icbuilder/`, `scripts/`, and `tests/` parse, and
  `git diff --check` passes.
- Frozen-reference tests show exact equality for grouped binning, cached
  uncertainty inflation, and grouped interpolation, including NaNs and masks.
- On tracked orbit-0085 intermediate inputs, complete `BinnedImage`
  construction improved from 8.687 to 3.782 seconds for WIC, 1.850 to 0.815
  seconds for SI12, and 1.824 to 0.818 seconds for SI13. All seven output
  fields were exactly equal. A full isolated orbit run improved from 25.83 to
  18.94 seconds, and the resulting 32-variable NetCDF file was byte-for-byte
  identical. Tracked example products were not modified.
- After the ConductanceImage optimization, another isolated orbit-0085 run
  improved from 18.84 to 7.52 seconds (2.51x). The complete 32-variable
  NetCDF remained byte-for-byte identical. Maximum RSS changed from 349.2 to
  349.5 MB. Under cProfile, `_compute_conductance` fell from 20.41 seconds and
  208,538 SciPy interpolator constructions to 0.004 seconds and one cached
  proton-response evaluation. Focused vector/scalar tests require exact E0,
  dE0, Fe, R, covariance, P, H, weights, and NaN support; propagated
  uncertainties permit only float64 roundoff. Across all orbit-0085 cells,
  maximum absolute differences were `4.4e-16` for dFe, `2.9e-11` for dR, and
  `1.8e-15` for both dP and dH; dR is large in magnitude near a small SI13
  denominator, so this absolute difference is still last-place arithmetic.
- The earlier four collapse figure sets remain valid; their polar orientation,
  threshold curves, and dE0 companion were visually checked previously.

For Kp integers 0–9 and 480 MLT samples at 0.05-hour spacing (4,800 slices),
using 0.01-degree MLAT cells over 50–90 degrees and the default threshold:

- 0 selections were empty;
- 1,001 reached the 50-degree equatorward sampling limit, indicating possible
  domain truncation;
- 0 reached a non-polar upper sampling limit;
- 1,378 reached the physical 90-degree pole, which is not sampling truncation.

Against a 0.01-hour reference grid, interpolation of the 0.05-hour area-mean
product had a 0.021 keV 99th-percentile error and a 0.037 keV maximum error.
For the four Kp=2/Kp=5 threshold curves in the result figure, refining MLAT
from 0.25 to 0.01 degrees reduced RMS second differences from
0.016–0.023 keV to 0.00095–0.00124 keV and removed the visible jitter.

## Remaining scientific boundary

Preserve these issues while testing the implemented stage-1 path:

- `ConductanceImage.Ep` is proton mean energy. The Frey response tables are
  explicitly tabulated against `<E>`; `icphysics.image` currently mislabels
  this table axis as characteristic energy in comments and docstrings.
- **Confirmed by the user:** DMSP tests did not reproduce the Frey et al.
  (2003) WIC/SI13-to-E0 relationship. Stage 1 now replaces IMAGE-derived E0
  and dE0 entirely with Zhang–Paxton rather than using a fallback or blend.
- Zhang–Paxton electron mean energy versus icBuilder electron
  characteristic-energy compatibility requires scientific confirmation.
- **Confirmed by the user:** use the 0.01-degree collapse to generate a lookup
  once for each supported Kp state on the fixed 36-by-36 IMAGE Cubed-Sphere
  grid. Production frames should index this table rather than recompute
  latitude profiles.
- Since the collapsed model is `E0(Kp, MLT)`, grid-cell MLAT does not enter the
  table value; it would matter only for a future return to the uncollapsed
  Zhang–Paxton model.
- **Confirmed by the user:** do not add a new precipitation-support mask; no
  robust automatic construction is available. Corrected WIC brightness
  supplies the spatial amplitude through Fe, while existing IMAGE coverage
  and background handling remain in place.
- Do not choose a single empirical Kp bin. The published method and current
  package interpolate E0 linearly in Kp and Q in hemispheric-power space
  between fitted-bin centres before the latitude collapse.
- Zhang–Paxton publishes no coefficient covariance or predictive-error model.
  **Confirmed by the user:** do not introduce DMSP calibration. Use the
  spherical-area-weighted E0 spread within the same selected latitude profile
  as dE0. This is uncertainty from discarding MLAT in the collapse, not formal
  Zhang–Paxton predictive or coefficient uncertainty.
- The stage-1 path now propagates the first-order E0–Fe covariance. This does
  not resolve shared-channel covariance or missing model-prediction
  uncertainty.
- Removing the WIC/SI13 ratio from E0 inference removes SI13 as a required
  Zhang--Paxton input. **Confirmed implementation:** retain simultaneous SI13
  as an optional proton-corrected diagnostic channel. At fixed Zhang–Paxton E0, compare
  WIC-derived and SI13-derived Fe, or predict corrected SI13 counts from
  WIC-derived Fe. Do not initially require SI13 for product generation or
  combine it into Fe.

## Audit findings that constrain integration

`scripts/debugging/debug_image_ratio.py` now reproduces the modular ratio path
in memory without overwriting any data product. It writes detector/background
stage plots, common-grid signal traces, response-domain and retrieval-category
summaries, and the numerical values behind the plots to `figures/debugging/`.
The figures were generated and visually inspected for example orbits 0085 and
0086.

- The IDL-to-fuvpy check confirms that saved SI13 `img` reproduces the IDL
  `IMAGE`, while WIC changes because the upstream read uses `reflat=True`.
- `scripts/debugging/debug_fuvpy_background.py` reruns the same input orbits
  into scratch space using both current and published fuvpy regularization,
  then compares DG/SH choices with identical grids and time support. It does
  not overwrite scientific products.
- WIC SH with SI12/SI13 DG is the best tested processing choice. After proton
  correction, current versus publication settings give median WIC/SI13 ratios
  of 242.9 versus 244.0 for orbit 0085 and 193.5 versus 194.2 for orbit 0086;
  the corresponding table-interior fractions are 8.9% versus 8.9% and 26.8%
  versus 26.6%.
- The saved SI SH products were generated with `1e-4` damping, whereas the
  publication uses 10. Correctly damped SI SH is less aggressive, but still
  worsens the WIC/SI13 ratio relative to SI DG. The user's decision not to use
  SH correction on SI12 and SI13 is therefore supported.
- Gasparini et al. (2024) provides a positive-control case using all three
  cameras. The paper first grids them together at 250 km, averages three raw
  images and their dayglow models over six minutes, clips after dayglow and
  proton removal, derives E0, and finally smooths the conductances through a
  regularized SECS inversion.
- The corresponding local orbit 0478 interval is much healthier than the May
  examples. The regenerated archive has median strong-pixel uncorrected ratio
  40.8 with 56.9% inside the Frey table. The historical corrected `R` has
  median 41.5 with 56.3% inside. The failure is therefore not universal.
- A three-frame mean applied approximately to the already binned May examples
  does not repair them: their median ratios remain about 207/244 for orbit
  0085 and 171/189 for orbit 0086 before/after proton correction. Temporal
  averaging alone cannot explain the positive-control difference.
- The paper's Zenodo record `10.5281/zenodo.10203397` contains only a processed
  WIC subset, not the SI12/SI13 products required for direct reproduction.
- Before proton correction, strong-pixel median WIC/SI13 is 208.2 in orbit
  0085 and 177.5 in orbit 0086; 14.6% and 30.4% fall inside the 34.84--136.49
  response-table interval.
- After proton correction, the medians rise to 242.9 and 193.5; table-interior
  fractions fall to 8.9% and 26.8%.
- Considering every finite corrected pixel, only 1.6% and 5.4% are interior
  retrievals. Roughly half use the weak-SI13 `E0=1` fallback and about 14--16%
  use the high-ratio `E0=25` saturation.
- Coumans et al. (2004) capped the result at 15 keV, while the active code
  saturates at the 25-keV Frey table endpoint. This changes the reported upper
  energy but cannot explain why the supplied ratio is already too high.
- The current code has no explicit Coumans factor-of-two adjustment of the
  SI12-derived proton flux. Since the high ratio precedes proton correction,
  this is a secondary discrepancy rather than the first failure boundary.
- Orbit 0968 is the exact 2001-10-21 DMSP F15 interval in Coumans et al. (2004)
  Figure 4 and is now the best direct benchmark. Their processing mapped SI
  into WIC image space, co-registered common fields of view, smoothed WIC and
  SI13 before taking the ratio, used spatial proton estimates, and interpreted
  the ratio with modeled camera responses. Their forward model treats viewing
  geometry and absorption, but the electron table is a nadir calculation and
  a per-pixel angle correction of the Figure-4 ratio is not documented. The current
  pipeline does not yet match the WIC/SI13 effective point-spread function and
  uses a fixed lookup response.
- Coumans Figure 4 is not a uniformly successful validation: the paper reports
  a greater-than-50% underestimate over part of the pass and artificial later
  energy peaks caused by WIC background subtraction. A controlled orbit-0968
  reconstruction should compare each processing stage rather than only the
  final capped energy curve.
- The historical background was independently estimated per camera, not by a
  joint WIC/SI13 fit. Frey and Coumans used camera-specific quiet-time dayglow
  responses with viewing/sensitivity corrections; Meurant instead used a
  per-image histogram on the nightside. fuvpy uses the same BS/SH model form
  for all three cameras, apart from WIC re-flattening and camera masks. The
  plausible historical difference is therefore the old camera-specific
  reference/calibration, not missing cross-camera statistical coupling.
- Historical safeguards and caveats are now explicit. Frey's response is a
  fixed-atmosphere nadir calculation and one active FAST case required an
  undocumented disturbed-atmosphere correction. Coumans caps electron energy
  at 15 keV and concludes that WIC/SI13 is useful for morphology while its
  quantitative values are often overestimated because of background
  sensitivity. Meurant's published atmosphere perturbation changes a ratio-120
  retrieval by about 16 percent, too little by itself to explain the corpus tail.
- The direct IMAGE inverse validations are one connected processing lineage.
  Gasparini et al. (2024) independently reimplemented the observed-IMAGE
  inversion but did not validate retrieved energy against particles. No
  peer-reviewed independent in-situ validation was found. Independent
  Polar-UVI work and Galand--Lummerzheim forward modelling support the
  two-colour physics, not IMAGE preprocessing or inversion accuracy.
- The first Coumans Figure-5b reconstruction is complete. Its FAST track and
  plausible loss-cone moment agree with the published geometry and peak scale,
  while the current IMAGE product reaches 11.41 keV instead of the roughly
  1.5-keV published IMAGE peak. Finish the undocumented FAST selection
  sensitivity and a detector-space IMAGE reconstruction before treating the
  numerical difference as final.
- Meurant's GLOW WIC/SI13 curve is numerically close to the active Frey
  response over their common range (roughly ratio 100 at 10 keV and 140--150
  at 25 keV). Its extension to 50 keV and Coumans' separate 15-keV output cap
  cannot explain why the measured ratio is already too high.
- fuvpy handles camera-specific background removal, WIC re-flattening, masks,
  and upstream FUVIEW3 geolocation, but not SI-to-WIC detector-space
  co-registration or common WIC/SI13 PSF matching. Its publication demonstrates
  WIC, SI12, and SI13 background removal but does not compare against Meurant's
  histogram method or validate the derived electron-energy ratio.
- DMSP F12--F15 SSJ4 one-second CDF products cover the IMAGE interval and
  provide calibrated electron spectra, total energy flux, average energy,
  uncertainty, and ephemeris. A compact comparison archive can retain only
  northern magnetic latitudes above 40 degrees and samples within an expanded
  window around simultaneous IMAGE frames. Preserve at least two minutes of
  padding so the published 66--120 s DMSP smoothing can be reproduced.
- WIC has a nominal 10-s exposure and SI13 a 5-s exposure. Those exposure
  differences are already incorporated in Frey's tabulated WIC and SI13 count
  responses and therefore in the active `E1/E3` ratio lookup. No extra
  exposure normalization should be applied; this issue is closed.
- A small common-grid smoothing test does not repair orbit 0968. Gaussian
  smoothing of both corrected channels at 0.5 current grid cells (about
  100 km) changes the fraction above the Frey table from 16.6% to 13.3%.
  Even an extreme 400-km sigma leaves 8.4% above the table, while WIC-only
  smoothing eventually makes the tail worse. Retain detector-space
  registration and common-PSF treatment as a deferred refinement if the
  ratio-scale problem is otherwise solved.
- `scripts/debugging/download_dmsp_ssj.py` now accepts an explicit UTC
  interval and downloads and combines all daily SSJ files touched by it. The
  complete F15 pass crosses midnight: 2001-10-21 23:33:15 through 2001-10-22
  00:01:58. The two daily CDFs produce a compact NetCDF with 1,724 one-second
  records above 40 degrees AACGM latitude. Its supplied
  average electron energies span approximately 0.034--5.23 keV; this moment
  has not yet been independently reconstructed from the energy channels. The
  reducer converts the spacecraft position to geodetic coordinates, maps it
  with ApexPy to the 130-km IMAGE shell, and stores geodetic, quasi-dipole,
  and MLT footprint coordinates for the same-event comparison.
- `scripts/debugging/compare_dmsp_image_ratio_pass.py` compares an arbitrary
  reduced pass with a processed image-ratio precipitation file. For orbit
  0968, 1,310 of 1,724 samples lie inside the IMAGE grid and within 75 s of a
  frame. The complete-pass figure displays both auroral crossings and the
  polar cap. The first crossing retains recognizable but displaced IMAGE
  structure; the return crossing has ratios of several hundred to several
  thousand, far beyond the Frey response table.
- `scripts/debugging/compare_dmsp_orbit.py` is the new orbit-number interface.
  It derives the IMAGE interval and stored grid, discovers every locally
  available DMSP satellite/day, maps each track to the grid height, separates
  continuous on-grid passes, and writes one comparison figure per pass plus
  `pass_summary.csv`. It recalculates WIC/SI13 and ratio-response energy from
  `wic_corrected` and `si13_corrected`, never stored `R` or `E0`, so
  Zhang--Paxton products remain valid inputs. The orbit-0968 smoke test found
  seven F15 passes and gave identical inventories for the IR and ZP products.
  `scripts/download_dmsp.py` retains only the local
  `ssj/precipitating-electrons-ions` archive.
- The bulk DMSP download is safely resumable. Existing daily CDFs are checked
  before being skipped; new files are streamed to `.partial`, checked against
  the HTTP byte count and expected 86,400-record SSJ variables, and atomically
  renamed. Each transfer is attempted up to five times and a persistent
  failure is reported without terminating the remaining downloads.
  Validation after the first server disconnect found 746 valid CDFs and one
  zero-byte file, which the next run will replace automatically.
- Current archive access is unstable during directory discovery: connections
  are accepted but frequently closed without an HTTP response or time out
  after 60 seconds. Cache the discovered URL manifest before optimizing file
  transfers. A later optional threaded downloader should be tested at two and
  then four workers only; 100 public-server connections are explicitly ruled
  out.
- `scripts/process_dmsp_yearly.py` converts the raw archive into compressed,
  atomic satellite-year files containing only northern records with
  recalculated modified-Apex latitude at or above 40 degrees. It uses a
  130-km reference height, keeps original geocentric and AACGM coordinates,
  and retains electron/ion energy and flux moments with fractional
  uncertainties. Counts and differential spectra remain in the immutable raw
  CDFs. A one-day F13 smoke test produced 24,802 sorted unique records and
  verified the northern latitude threshold and units.
- A direct source-to-reduced audit found no DMSP import defect in the earlier
  606-record focused subset. Its times are monotonic and exactly match the
  source CDF, and the retained spacecraft
  coordinates, AACGM latitude, electron mean energy and uncertainty, and total
  energy flux are numerically unchanged. The Figure-4 interval contains 251
  finite one-second energies in eV and the reconstruction correctly converts
  them to keV. Geocentric radius is explicitly in kilometres, the derived
  spacecraft height is 855--857 km, and ApexPy's reported mapping error is
  below 0.000007 degrees. The source contains no omitted quality-flag field.
- The coordinate path in the three-frame plots is correct and does not use
  spacecraft-altitude AACGM as an ionospheric footprint. DMSP is mapped along
  the field to 130 km and then expressed as QD/MLT; this agrees with modified
  Apex at the spacecraft using `refh=130` to within 1.2e-5 degrees. The CDF
  AACGM track is only 0.31 degrees poleward in the Figure-4b interval and
  moves the sampled ratio peak five seconds earlier. Direct `geo2qd` at the
  approximately 855-km spacecraft altitude is the misleading alternative: it
  is about 1.80 degrees equatorward and shifts the sampled peak by 30--50 s.
  Coumans explicitly mapped the track to 120 km; changing the present mapping
  from 130 to 110 km moves the peak by only one second. The figure reads the
  36-by-36 grid stored in the mounted product (`L=20,000 km`, 225-km cells,
  older 110-km radius); it does not use the current 46-by-46 `L=50,000 km`,
  200-km grid definition. Reconstructing the saved image on that newer grid is
  not a like-for-like coordinate comparison.
- `scripts/debugging/reconstruct_coumans_figure4b.py` performs that first
  comparison. The supplied DMSP mean-energy curve agrees qualitatively with
  Coumans Figure 4b, validating the event, time interval, and particle
  product. The current IMAGE reconstruction does not: it peaks near 12.4 keV
  around 23:38:30 and falls to the 0.2-keV floor after 23:40, with correlation
  -0.12 against the nine-second-smoothed DMSP series. Ratios sampled along the
  mapped pass remain below 111, so this particular track does not encounter
  the whole-image high-ratio saturation. Source AACGM versus ApexPy-mapped
  coordinates change energy by only about 0.3 keV in the median. The next
  comparison should therefore focus on current versus Coumans preprocessing,
  not the DMSP ephemeris or Frey lookup orientation. The comparison now uses
  nearest-cell values rather than bilinear spatial interpolation; a companion
  WIC/ratio figure marks all 16 grid cells that supply the plotted IMAGE series.
- The reconstructed IMAGE peak precedes the visually similar Coumans peak by
  about two minutes, but a simple one-frame indexing error does not reproduce
  the shift: previous, nearest, and next-frame selection all retain the main
  peak at 23:38:22. The corresponding Coumans time lies about 7.3 degrees
  farther poleward along the mapped track (58.8 versus 66.1 degrees QD), much
  larger than ordinary AACGM/QD or 120/130-km mapping differences. Test the
  three ratio frames with time-marked track positions before deciding between
  serious geolocation/co-registration error, upstream image/time association,
  and genuinely different morphology from Coumans' preprocessing.
- The three-frame alignment figure rules out selecting the wrong IMAGE frame,
  but a new along-track experiment reveals an unexpectedly precise two-minute
  position/time offset. Keeping DMSP energies at their original times while
  sampling IMAGE at the footprint from `time - 2 min` moves the IMAGE peak to
  23:40:26 with ratio 109 and retrieved energy 11.6 keV, close to Coumans'
  23:40:20 and 12.5--13 keV peak. The three fixed-frame profiles peak at
  23:40:26, 23:40:26, and 23:40:43, and the later shifted structure also
  resembles the published two-peak curve. This is strong diagnostic evidence
  for a position/time-association difference, not justification to hard-code
  a two-minute correction. Identify whether the discrepancy lies in the DMSP
  ephemeris association, the historical Coumans method, or the present
  reconstruction before changing production data.
- `scripts/debugging/plot_coumans_figure4a_idl.py` now provides the raw
  geolocation test. It bypasses fuvpy and all gridded products, plots the native
  WIC IDL `IMAGE`, and maps the 130-km DMSP footprint directly through the IDL
  `GLAT`/`GLON` arrays. The resulting track crosses the left auroral oval in
  the same detector-space location and orientation as Coumans Figure 4a. Its
  median/maximum nearest-pixel separations are 0.136/0.229 degrees, far below
  the roughly 7-degree shift under consideration. The exact Coumans detector
  frame at 23:35:29 is absent from the local IDL subset, but the mounted Halley
  fuvpy orbit files contain simultaneous WIC, SI12, and SI13 products at that
  time. Its processed image-ratio profile still peaks at 23:38:51; at the
  published late peak near 23:40:20 the ratio is only 23.9 and maps to the
  0.2-keV floor. Nearby 23:33--23:41 frames also retain the earlier peak.
  Before proton correction, the exact-frame early peak is already about 137
  (`WIC/SI13 = 7322/53.4`), while the late position is only about 30. The
  discrepancy therefore enters no later than the background-corrected channel
  images; proton correction and the energy lookup do not create it.
  **Confirmed by the user:** an independent reconstruction by a colleague
  using an older IMAGE dataset shows the same discrepancy. Audit initial
  WIC/SI13 handling and historical response compatibility next; DMSP import,
  gross geolocation, missing-frame selection, and downstream frame assignment
  are now low-priority suspects.
- The user has corrected the separate Product-1 height inconsistency in the
  live worktree: the ApexPy input height, modified-apex reference height, and
  CS-grid shell are now all 130 km. Numerical comparison with native FUVIEW3
  `MLAT` strongly supports treating IDL `GLAT` as geodetic. Direct use at
  130 km gives a signed median latitude residual of +0.0005 degrees and RMSE
  0.271 degrees; converting it from geocentric to WGS84 geodetic first gives
  +0.178 degrees and RMSE 0.283 degrees. `GLON` requires no such conversion.
  This remains an evidence-based inference because the original FUVIEW3
  source/manual has not been recovered.
- In the earlier mixed-height DMSP comparison, re-expressing the track in the
  Product-1 convention changed 64/251 nearest-cell assignments but shifted the
  ratio peak only from 23:38:22 to 23:38:17 and left its value unchanged. The
  former height inconsistency is not the explanation for the two-minute
  Coumans disagreement.
- The DMSP comparison correctly treats the CS-grid longitude-like coordinate
  as `MLT * 15`: Product 1 bins IMAGE pixels that way, DMSP footprint MLT is
  calculated with ApexPy, and the comparison passes `footprint_mlt * 15` to
  `geo2cube`. The selected cells differ from the track by at most 1.10 degrees
  latitude and 0.20 MLT hours over the Figure-4 interval, ruling out an
  accidental magnetic-longitude/MLT substitution.

- Orbit 0085 has a median uncorrected binned WIC/SI13 ratio near 208 for SI13
  above three counts, already beyond the table maximum near 136.5. Fixed
  proton correction raises the comparable median to about 244.
- Across five tracked conductance products, only about 1.7–11.2% of valid
  pixels are interior-table E0 retrievals. Most are exact 0.2-, 1-, or
  25-keV fallback/saturation values.
- The former coarse-SI-grid defect is fixed in the current uncommitted work.
  Numerical inspection confirms exact 18-by-18 nesting: every SI boundary is
  every second 36-by-36 WIC edge, and corresponding SI centres are the
  pairwise WIC-centre midpoints.
- Even after that bug is fixed, 450-km SI13 is upsampled to 225-km WIC without
  matching effective spatial resolution.
- Fixed-2-keV proton correction, independent clipping, a three-count SI13
  threshold, and hard fallbacks amplify the failure. They are not sufficient
  explanations for the already-high raw/uncorrected ratio.
- The user is investigating the SI-to-WIC mapping. Published IMAGE-FUV work
  did map SI into WIC image space, so interpolation itself is not an error.
  The confirmed nested-grid construction defect appears modest in aggregate
  and is not the leading E0 explanation.
- The upstream Laundal/Østgaard FUVIEW3 product uses corrected counts,
  including detector flat-field and mission-time/temperature correction.
  The generated products still do not retain enough calibration,
  detector-state, source, unit, processing, sample-count, or retrieval-status
  metadata to prove the exact recipe or compatibility with the Frey response.
- A fixed standard-atmosphere/nadir response is applied to frames admitted out
  to 75 degrees detector zenith angle. Published work identifies atmospheric
  oxygen variability as a direct cause of high WIC/SI13 and unreasonable
  inferred energy.
- A focused literature review found no validated IMAGE WIC/SI13 retrieval
  cutoff. Ohma et al. (2024) says an upper DZA limit is required because
  `1/cos(DZA)` diverges in the optically thin dayglow model and large-angle
  geolocation is uncertain; 70--80 degrees is only an example range. This supports the existing
  near-limb exclusion but does not prove the ratio is quantitative up to 75
  degrees. Frey et al. (2003) and Meurant et al. (2003) use nadir/vertical
  response curves and require viewing-dependent line-of-sight and absorption
  treatment. A common cosine multiplier cancels in WIC/SI13 and is not a
  solution to the ratio problem.
- Follow-up inspection confirms that fuvpy uses DZA in its dayglow model but
  does not normalize the remaining auroral counts to a nadir-equivalent
  response. In an exploratory same-225-km-grid orbit-0085 test, the median
  uncorrected strong-pixel ratio increased from about 56 for DZA below 30
  degrees to 93, 130, and 156 for cutoffs of 45, 60, and 75 degrees. Changing
  geographic coverage prevents a causal interpretation, but viewing angle is
  now a concrete diagnostic rather than only a literature concern.
- A mission-scale post-hoc test now applies cumulative WIC DZA masks to all
  1,686 regenerated Zhang--Paxton precipitation orbits. A 70-degree mask
  retains 97.27% of product pixels, 94.65% of summed energy flux, 96.23% of
  summed Pedersen conductance, and 96.31% of summed Hall conductance relative
  to the existing 75-degree support. A 60-degree mask retains 86.88% of
  product pixels but only 76.73% of summed energy flux.
- The same test on the fixed, incomplete 188-orbit image-ratio snapshot does
  not support high DZA as the cause of the ratio problem. Changing 75 to 70
  degrees raises the median orbit ratio from 89.81 to 92.54, while the fraction
  inside the Frey table changes only from 49.01% to 49.38%. It retains 97.45%
  of product pixels and 94.22% of summed energy flux. These image-ratio values
  are provisional until the full run finishes.
- Median orbit ratios hid the pathological upper tail. At 75 degrees, 25.18%
  of supported ratio pixels exceed 150, 14.96% exceed 200, and 5.10% exceed
  300; the 90th percentile is 237.84. Every sampled orbit has at least one
  value above 150. At 30 degrees, 26.86% still exceed 150. The DZA sweep
  therefore does contain the previously observed failure, but the first plot
  summarized it with the wrong statistic.
- `ratio_vs_dza.png` directly displays the 4,339,186 supported pixels as a
  log-density field with conditional quantiles and tail fractions. It shows no
  positive ratio-versus-DZA trend: the median is nearly flat through most of
  the range and falls near 75 degrees, while ratios above 150--300 occur at all
  DZA. This pooled result is still confounded by auroral location and event.
- The strongest immediate lead is SI13 denominator sensitivity. Of supported
  ratio pixels, 69.7% have corrected SI13 between 3 and 7.5 counts. Ratios
  above 300 have median SI13 4.44 counts and median WIC about 1,780 counts.
  However, 20.9% of pixels with SI13 between 20 and 30 counts still exceed
  150, so low SI13 alone is not the full cause. Test matched effective spatial
  resolution and aggregation statistics before modifying fuvpy background
  subtraction.
- A 150-km common-grid experiment exposed two distinct support losses. Product
  1 requires two finite detector pixels per bin; in orbit 0085 only 25.1% of
  SI12 and 30.6% of SI13 cells satisfy this, with 43.5% of SI12 cells holding
  exactly one pixel. Product 2 then unnecessarily applies four-neighbor
  bilinear interpolation between identical grids, reducing SI12 support to
  3.7%, SI13 to 7.9%, and final image-ratio support to 1.6%. Add an identical-
  grid direct-copy path before using common grids. For the scientific
  resolution test, coarsen WIC to 450-km SI resolution rather than refining SI.
- The identical-grid direct-copy path is implemented and tested. It compares
  shape plus xi, eta, latitude, and longitude, then copies signal, uncertainty,
  and weight without changing support. On current 100-km orbit-0085 files it
  preserves all 65,373 finite SI12 and 69,322 finite SI13 cells. Remaining SI
  gaps come from empty detector sampling rather than Product-2 interpolation.
- This screening test uses the median WIC DZA already stored in Product 1 and
  leaves retained values unchanged. A bin below the threshold can contain raw
  source pixels above it, and the SI-camera viewing angles are not jointly
  masked. A final cutoff decision therefore requires Product 1 to be rerun
  from source pixels at the candidate limit.
- The uncertainty chain is not publication-grade: saturation can return
  dE0=0, corrected channels share omitted uncertainty, E0–Fe covariance is
  forced to zero, and binned-median uncertainty is not estimated coherently.
  The user asked that uncertainty redesign be recorded and deferred.
- `SplineImage.solverP` passes Hall uncertainty into the Pedersen solver, a
  confirmed secondary-product bug. The user asked that it be retained for
  later correction.
- With Zhang–Paxton replacing ratio-derived E0, the weak-SI13 fallback problem
  no longer affects the primary E0 product. Revisit it only if SI13 remains a
  validation channel.
- Østgaard et al. (2018) used solar- and satellite-zenith-angle-dependent
  dayglow subtraction but deliberately avoided absolute auroral intensity
  comparisons; it warns about oblique views and supplies no auroral LOS
  normalization. The correct quantitative treatment is an angle-dependent
  atmosphere/instrument forward response. Multiplying each
  dayglow-subtracted pixel by `cos(DZA)` is only the plane-parallel,
  optically-thin first-order diagnostic and must be validated over a chosen
  moderate-DZA range. LOS remains relevant after Zhang–Paxton because WIC
  intensity still determines Fe.
- The repository contains an SZA/DZA and cosine-correction path.
  Inspection confirms that all three example input channels contain matching
  degree-valued geometry arrays and that the operation correctly multiplies
  the background-subtracted pixels by `cos(DZA)` before binning.
  `BinnedImage` now stores median SZA, median DZA, median `cos(DZA)`, and the
  correction flag, and carries them through interpolation and discard.
  `PreImage.discard` now also keeps raw `img` aligned. `ConductanceImage`
  copies and serializes the geometry separately for WIC, SI12, and SI13 and
  records each channel's correction mode and LOS-applied flag. The modular
  Product-1 builder explicitly sets `los_correction=False`; the optional
  correction remains unsuitable as a general three-camera quantitative
  treatment, especially for SI12. A one-frame orbit-0085 WIC test gave
  corrected/uncorrected binned ratios of 0.296, 0.623, and 0.932 at the 5th,
  50th, and 95th percentiles.

## Next action

Review and then execute [[icReader Read Migration Implementation Plan]],
starting with the dependency/baseline gate and the Product-1-to-Product-2
loader. Keep each slice numerically identical before migrating the next
boundary. The separate scientific decision about Hardy clipping and
low-signal validity remains unresolved; this interface cleanup must not alter
those products.

The remaining text in this section records the older ratio/background
debugging backlog. It is not the implementation sequence for the new product
architecture.

Finish Figure 5b before adding another empirical fit. Establish how weak FAST
precipitation and the lowest EESA channels were excluded historically, then
repeat the IMAGE branch directly from the unbinned orbit-0459 WIC/SI12/SI13
fields using detector-space area coregistration. Keep the present 50-eV result
as an explicit sensitivity, not a silently selected match to the paper.

Visual inspection of the new F13 and F15 crossing panels finds poor agreement
for most or all inspected crossings, not merely the Coumans event. Stop fitting
a time offset to that single event. First use orbit 0968 as a positive-control
comparison between `data/matches.nc` and the focused Figure-4 script. Then
restrict the corpus comparison to finite, non-default, adequately illuminated
WIC/SI13 samples and compare observed ratio with the Frey-table ratio predicted
from DMSP energy. Coumans-style smoothing can be tested separately, but it
cannot reasonably account for missing broad correspondence across nearly all
crossings.

The selected DMSP variable and units are verified: `electron_mean_energy` is
raw SSJ `ELE_AVG_ENERGY` (total energy flux divided by total number flux),
scaled from eV to keV. Values above 10--15 keV are genuinely rare in the
northern F13/F15 products. The regenerated match product now also carries
total electron energy flux and the fractional uncertainties of both energy
and flux.

Do not use `good_orbit_segments.csv` as the primary DMSP validation filter.
The current `<503` plus accepted-frame selection reduces 23.43 million matches
to 472,863 and nightside matches from 4.76 million to 43,577; finite IMAGE
ratio support leaves 10,990, dominated by F13 and limited to 04--06 and 18--22
MLT. The strict filter at line 185 leaves only 2,426 nightside rows. The
annotations describe globally usable VAE frames and can discard frames whose
local DMSP-track pixels are still usable. Build the baseline from local
signal/geometry/SSJ quality and apply the annotation selection only as a
sensitivity check.

Use `scripts/dmsp/annotate_dmsp_frames.py` for the new local quality filter.
It shows corrected WIC/SI13 with all matched DMSP tracks for one frame and
annotates one satellite at a time. The active track is thick magenta and any
simultaneous tracks are thin gray. Figure-window `0`/`1` keypresses are saved
immediately, without Enter, in a resumable CSV keyed by orbit, IMAGE timestamp,
and satellite. `q` or closing the figure exits without annotating the active
satellite. The existing 22,397 single-satellite rows were verified to map
uniquely to the new resume key. The no-data `0` and dayside-only `2` tests now
operate per satellite. An isolated two-satellite test produced independent
`0` and `2` rows, and orbit 0466 frame 034 verified the magenta/gray display.

The uncertainty histograms show clear transitions near 0.25 for DMSP mean
energy and 0.20 for total energy flux. Treat these as candidate quality cuts
and verify the conclusions against looser and stricter cuts. Coumans et al.
(2004) compared IMAGE with the same moment mean stored in the SSJ product.
Frey's response tables are also tabulated by mean energy `<E>`, although the
current `icphysics.image` comment calls the table axis characteristic energy.
A peak-channel or fitted-Maxwellian characteristic energy is a separate,
model-dependent diagnostic rather than a replacement for this validation
quantity.

With finite positive IMAGE ratio, those uncertainty cuts leave 69,184 samples
from accepted annotations. Do not impose a lower flux-amplitude cutoff. The
uncertainty cuts already reject almost the entire lowest flux decile, while the
retracted `3e11` limit came from overgeneralizing Meurant et al.'s
event-specific 1.5-mW/m2 comment and retained only 27.7% of the
quality-filtered sample.

The upper tail contains a distinct data defect. Seven surviving samples have
electron energy flux at or above `1e14 eV cm-2 sr-1 s-1`; all are F12, and raw
CDF inspection identifies isolated one-second spectral spikes, including
counts beyond `ELE_COUNTS_OBS` `VALIDMAX`. The next-largest normal-looking
accepted value is `7.06e13`. Use `flux < 1e14` as a temporary gross-artifact
guard for the current reduced archive, not as a physical upper limit. The
proper fix is to preserve raw count-validity and record-quality flags when
building the yearly files. DMSP instrument notes independently warn that all
F12 instruments were noisy from 1999 onward. Use F13/F15 for the primary
validation, show F12 separately, and stratify F14 by year because its electron
detector degraded progressively after 2001.

Restrict to `0.2 <= E <= 25 keV` only for a quantitative comparison with the
Frey table; this is model support, not a DMSP quality selection. Keep the
current matches as raw 1-s samples until resolution matching. Before final
regression, smooth total energy flux and inferred number flux to IMAGE support,
then form their ratio to obtain the matched mean energy.

The preliminary ratio plot does not apply the complete DMSP quality/satellite
mask or corrected `SI13 >= 3` support rule. A clean F13/F15 repetition with
those selections leaves 18,906 unique timestamps; median DMSP energy across
eight Frey-domain ratio bins is only 1.17--1.96 keV and the 90th percentile is
3.59--6.04 keV. Quality filtering does not recover the Frey inverse. Reverse
conditioning still shows weak physical sensitivity: median ratio rises from
about 66 near 0.35 keV to about 112 near 6.8 keV, but pooled Spearman
correlation is only 0.148 and crossing-level slopes are mixed.

An orbit-grouped monotonic ratio model improves log-energy prediction over a
constant by only 0.1% on Frey support, with an orbit-balanced confidence
interval spanning zero. Treat this as evidence that a universal deterministic
ratio inverse is currently non-predictive, not proof that the camera signals
contain no energy sensitivity. The fuvpy robust background weight is also not
a simple quality probability: low values often identify auroral foreground but
can equally identify model errors or artifacts.

**Proposed next decision:** run one final resolution-matched forward-closure
test of `WIC / Q` and `SI13 / Q` separately against their Frey response shapes.
If it fails, stop trying to repair the inverse. Compare a probabilistic constant,
collapsed Zhang--Paxton, and DMSP-calibrated Zhang--Paxton under blocked
orbit/year/satellite validation. Test an IMAGE-ratio residual only if it adds
repeatable held-out skill and has production-domain support. This is proposed,
not yet a confirmed pipeline choice.

Accepted-frame distributions now demonstrate that IMAGE's high-energy tail is
mostly a retrieval artifact. DMSP exceeds 10/15 keV in only 0.21%/0.025% of
111,818 finite samples, whereas IMAGE exceeds those levels in 13.0%/11.3% of
its finite matches and returns exactly 25 keV in 8.7%. Exact-25-keV IMAGE
matches have DMSP median 1.02 keV. The Frey table maps only ratios 34.8--136.5;
29.8% of accepted finite ratios exceed 136.5 and are converted to 25 keV by
the interpolator's upper `fill_value`. Treat 25 keV as out-of-domain status,
not a physical estimate. The two products both report mean energy, so this is
not a characteristic-versus-mean factor-of-two error.

The low DMSP histogram is not evidence of a unit/conversion bug. Raw CDF
metadata gives `ELE_AVG_ENERGY` in eV as total energy flux divided by total
number flux, the reducer applies the correct single `1e-3` conversion to keV,
and matched values equal the yearly product. The unconditional accepted-frame
sample includes weak polar-cap/background precipitation. Conditioning on
electron energy flux moves median DMSP mean energy from 0.53 keV to 2.16,
2.67, and 3.24 keV for thresholds of `2e11`, `5e11`, and
`1e12 eV cm-2 sr-1 s-1`; their 95th percentiles are 6.91, 7.96, and 8.69 keV.
This supports 2--3 keV as a plausible simple auroral baseline for these times,
but not yet as a universal constant for all IMAGE orbits. The 15--25-keV part
of the Frey response curve covers rare physical/model inputs rather than the
typical observed distribution.

Use that energy flux to interpret the possible WIC-count/DMSP-energy trend in
the accepted annotations. Plot WIC against DMSP electron energy flux, then
WIC divided by that flux against DMSP mean energy and the Frey response table.
Do not fit mean energy directly from WIC counts until this test shows that the
trend survives conditioning on flux; single-band brightness and mean energy
are not uniquely related.

Annotation quality has a sharp May--June 2001 hole. The CSV currently reaches
orbit 0838; acceptance among manually judged `0`/`1` rows is about 0.18% in
May and 0.11% in June, recovering to 7.7% in July and 12.8% in early August.
The failure affects both F12 and F13. From orbit 0600 to 0725 the valid WIC
grid changes from 50% sunlit/36% dark to 82% sunlit/5% dark, while median DZA
stays near 36--39 degrees. Frey et al. (2003) also place an IMAGE thermal
maximum around May 15 and document temperature-dependent SI gain. The leading
interpretation is therefore summer dayglow/background difficulty compounded
by SI thermal/gain behavior, not changing DZA or one DMSP satellite. Because
the upstream FUVIEW3 product is nominally mission-time/temperature corrected,
this is not proof that gain correction is absent. A focused raw/BS/SH and
calibration-metadata comparison across orbits 0675--0800 is the next causal
test.

Investigate orbit 0341 frame 087 as the first WIC/SI13 registration case. Its
apparent displacement exists before proton correction, and the active
three-sensor matcher limits source-time spread to two seconds. The pipeline
trusts each sensor's supplied geolocation and has no residual feature-based
coregistration step. A preliminary integer-shift test favored about five
200-km WIC cells but still correlated poorly, so do not apply a global shift.
First compare the original per-sensor geographic coordinate maps and degrade
WIC to SI resolution; then test whether displacement follows detector
position/DZA (mapping) or auroral morphology (physical spectral response).
Also compare raw, BS-only, and SH-corrected fields on that common grid. The
published methods estimated background independently by instrument, but the
active product uses WIC SH and SI13 DG, so correction-stage asymmetry is a more
specific concern than independent fitting alone.

The full DMSP SSJ archive and yearly northern reductions are now available,
and `scripts/dmsp/ratio_validation_data.py` has produced the first full-corpus
IMAGE-ratio/DMSP match file. Before drawing scientific conclusions from
`scripts/dmsp/ratio_validation.py`, replace the disproven `[time, time + 120
s]` assignment. Match each DMSP sample to its nearest IMAGE central snapshot
within `+/- 60 s`, retain the signed separation, include the final frame's
centred support, and regenerate `data/matches_2min.nc` before treating its
statistics as final.

The instrument literature now resolves exposure, cadence, and timestamp
meaning. IMAGE spins once per 120 s, while WIC views an Earth location for
about 10 s and SI for about 5 s during each spin. Frey et al.'s calibration
values already include those nominal exposure durations. NASA's IMAGE
data-management plan assigns the WIC image time to the centre of integration;
Coumans et al. (2002) also use the central snapshot time and extract the
satellite track over `time +/- 1 min`. fuvpy passes the IDL `TIME` field
through unchanged. Replace the asymmetric `[time, time + 120 s]` matcher with
a centred nearest-frame association and retain the signed time separation.

Do not interpret the apparent success of the current -1-min Coumans shift as
a timestamp result. `reconstruct_coumans_footprint_minus_1min.py` shifts only the
DMSP footprint used to sample IMAGE while leaving DMSP energy and time
unchanged. That diagnostic tests along-track spatial association, not IMAGE
clock time.

The focused Figure-4 scripts now apply the published convention directly.
The Figure-4b reconstruction was reduced from more than 600 lines to a focused
251-line calculation that keeps nearest-frame matches within `+/- 60 s`; the
native-IDL Figure-4a map highlights the same centred track segment while
retaining the full pass as context. A run against the tracked orbit-0968
products retained 248 centred time matches, including 243 finite IMAGE ratios,
but left the IMAGE-ratio maximum at 23:38:26. The separate footprint-shift
script is now a 29-line spatial experiment: using the DMSP position from 60 s
earlier moves the IMAGE maximum to 23:39:26 while leaving DMSP energy and time
unchanged. Do not describe that result as an IMAGE timestamp correction.

Resolve the Zhang--Paxton lookup/grid mismatch before rerunning Product 2. The
active binned WIC grid is 46x46, whereas the lookup used by icPhysics and its
loader validation are fixed at 36x36. Regenerate the table from the active
grid and remove the split ownership in which the generator/default output live
in icBuilder but the production loader and bundled table live in icPhysics.
Changing only the hardcoded shape is insufficient because the stored MLT,
xi, eta, E0, and dE0 arrays are themselves 36x36.

Deploy the verified Numba Product-1 footprint kernel before resuming the
roughly 300-hour Halley rerun. The change is local to polygon clipping and
matches the readable reference overlap matrix to floating-point precision.
Locally it reduced a 20-frame WIC `BinnedImage` from 24.9 s to 0.40 s and a
complete WIC orbit to 1.71 s. Confirm Numba is installed on Halley, pull the
change, benchmark one orbit there, and then resume without `--overwrite` so
already completed matching products are skipped. Apex conversion and NetCDF
I/O are now likely to dominate.

Let the image-ratio precipitation run finish, then rerun
`scripts/debugging/dza_threshold_sensitivity.py` so that its image-ratio
statistics cover the full corpus. The current Zhang--Paxton result is already
complete. The sensitivity does not justify lowering the operational limit
solely to repair WIC/SI13. Ohma et al. does not prescribe 70 degrees; choosing
that lower end of its example range would remove about 5% of summed energy
flux while leaving the ratio problem essentially unchanged.

Before treating the result as publication-ready, resolve whether
Zhang--Paxton electron mean energy is compatible with the energy quantity
assumed by the WIC response and Robinson relations. In parallel, continue the
channel-specific LOS investigation because WIC brightness still determines
Fe.

For the IMAGE-ratio investigation, do not change the selected WIC SH and SI
DG background products or tune the fuvpy regularization. The controlled rerun
shows that the published settings leave the preferred-path ratio essentially
unchanged. Orbit 0086 has a clear increase in median ratio with WIC DZA on a
fixed common-pixel comparison, but orbit 0085 is not monotonic and has too few
low-DZA pixels for the same test. Viewing angle remains a concrete diagnostic,
not a demonstrated correction.

The next defensible boundary is the provenance and compatibility of the
calibrated FUVIEW3 counts, including WIC reflat/calibration, cross-camera
response, and the original Frey/Coumans processing recipe. An original
reference event or processed count product is needed for absolute validation;
tuning the E0 interpolation cannot fix incompatible input counts.

`scripts/debugging/recreate_frey_figure16.py` implements the first published
benchmark from the unbinned sensor NetCDFs. It performs the established
footprint-overlap SI-to-WIC coregistration, applies the SI12 proton correction,
and reproduces the Figure-16 count/energy/flux layout without reading binned
products or refitting fuvpy. The orbit-0364 run selected the simultaneous
11:38:24 frames and completed successfully. The exact published proton-energy
assumption remains unresolved; the diagnostic currently exposes a fixed value.

Use orbit 0478 as the next positive-control experiment. If its Product-1
WIC/SI12/SI13 files can be recovered from the remote corpus, run the modular
pipeline for the paper interval and compare each processing stage with the
archived `conductance/or_0478.nc`. This separates an event-dependent physical
response from a processing difference without requiring the unavailable
fuvpy demonstration event.

Only after the controlled comparison should stage 2 decide whether SI13
becomes optional and which validation diagnostic, if any, is scientifically
defensible. Do not combine that support change with the initial E0 comparison.

## Portfolio impact

- Central update needed: No
- Changes: a repository-local implementation plan now defines where icReader
  should replace duplicate reads of generated icBuilder products and where raw
  access must remain.
- No scientific result, publication priority, or deadline changed.
- Next technical action: implement the migration in parity-checked slices,
  preserving detector-to-CS read-once performance and restart behavior.

## Entry points

- `vault/02_Algorithm/icReader Read Migration Implementation Plan.md`
- `icbuilder/conductancedetector.py`
- `icbuilder/detectorcs.py`
- `icbuilder/precipitationcs.py`
- `icbuilder/conductancecs.py`
- `scripts/pipeline/make_conductance_detector_orbit_files.py`
- `scripts/pipeline/make_detector_cs_orbit_files.py`
- `scripts/ZhangPaxton2008_collapse.py`
- `icbuilder/zhang_paxton_collapse.py`
- `icbuilder/zhang_paxton_lookup.py`
- `icbuilder/kp.py`
- `icbuilder/grids.py`
- `icbuilder/conductanceimage.py`
- `icbuilder/precipitationimage.py`
- `icbuilder/fuvdetector.py`
- `icbuilder/precipitationdetector.py`
- `icbuilder/detector_coregistration.py`
- `scripts/pipeline/make_fuv_detector_orbit_files.py`
- `scripts/pipeline/make_precipitation_detector_orbit_files.py`
- `icbuilder/imagesat_e0_eflux_estimates.py`
- `scripts/download_gfz_kp.py`
- `scripts/make_conductance_orbit_files.py`
- `scripts/make_zhang_paxton_lookup.py`
- `scripts/plot_zhang_paxton_lookup.py`
- `scripts/debugging/debug_image_ratio.py`
- `scripts/debugging/recreate_frey_figure16.py`
- `scripts/debugging/digitize_frey_figure16_counts.py`
- `scripts/debugging/dza_threshold_sensitivity.py`
- `scripts/debugging/plot_ratio_vs_dza.py`
- `scripts/debugging/download_dmsp_ssj.py`
- `scripts/debugging/compare_dmsp_orbit.py`
- `scripts/debugging/compare_dmsp_image_ratio_pass.py`
- `figures/debugging/image_ratio_response_domain.png`
- `figures/debugging/image_ratio_retrieval_categories.png`
- `figures/debugging/image_ratio_summary.csv`
- `figures/debugging/paper_reconstruction/frey_figure16_digitized_count_distributions.png`
- `figures/debugging/paper_reconstruction/frey_figure16_digitized_bottom_comparison.png`
- `figures/debugging/dza_sensitivity/dza_threshold_summary.csv`
- `figures/debugging/dza_sensitivity/dza_band_summary.csv`
- `figures/debugging/dza_sensitivity/ratio_vs_dza.png`
- `literature/README.md`
- `icbuilder/data/zhang_paxton_e0_lookup.nc`
- `icbuilder/data/gfz_kp_2000_2003.json`
- `icbuilder/data/README.md`
- `tests/test_kp.py`
- `tests/test_e0_override.py`
- `tests/test_conductanceimage_zhang_paxton.py`
- `tests/test_zhang_paxton_lookup.py`
- `tests/test_grids.py`
- `tests/test_zhang_paxton_collapse.py`
- `figures/README.md`
- `figures/zhang_paxton_collapse_process_map.png`
- `figures/zhang_paxton_collapse_latitude_slice.png`
- `figures/zhang_paxton_collapse_result.png`
- `figures/zhang_paxton_collapse_dE0_result.png`
- `figures/zhang_paxton_lookup_kp1_52.png`
- `vault/02_Algorithm/Processing Pipeline.md`
- `vault/02_Algorithm/Detector-First Product Architecture.md`
- `vault/02_Algorithm/Proposed Modular Pipeline Redesign.md`
- `vault/02_Algorithm/Audit - 2026-07-29.md`
