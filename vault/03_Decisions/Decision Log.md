# Decision Log

Last reviewed: 2026-10-11

## 2026-10-11 — Optional full spatial Hardy electron energy

**Confirmed by the user:** implement the optional Hardy electron mean-energy
model using the supplied functional paper and each pixel's MLAT, MLT and Kp.
Electron/proton choices remain independent; no default replacement requested.

**Implementation choices:** correct the Fourier phase to 24 h periodicity;
interpolate evaluated log10 flux for fractional Kp; expose the Northern 50–90°
source domain and Modified Apex approximation. Retain raw E0 and flag endpoint
WIC response evaluation. Propagate measurement errors conditional on fixed
Hardy energy, explicitly excluding prediction error. These conventions are
documented implementation choices, not newly confirmed scientific policies.
Broader particle/coordinate validation precedes production model selection.


## 2026-10-09 — Plan dynamic orbit-specific IMAGE response curves

**Confirmed by the user:** add a dynamic MSIS/GLOW and sensor-response pipeline
to the TODO list. Generate orbit-specific WIC/SI13 ratio relations conditioned
on F10.7, F10.7a, Ap, geographic latitude/longitude, DZA and UT. Reduce cost
through occupied model-evaluation bins, with daily solar and daily/3-hourly
geomagnetic drivers. Accept GLOW's supplied excitation physics.

The user hopes Benoit Hubert will supply DZA-dependent sensor response curves;
receipt and exact content remain unconfirmed. Bin sizes, calibration inputs,
runtime/accuracy and production adoption are unresolved. This is planned work,
not an implementation or an instruction to launch a full-corpus regeneration.
[[TODO]] owns the checklist; [[Audit - Published IMAGE Ratio Limitations]] owns
the source/provenance evidence.

## 2026-10-06 — Default real-observation STD floor in the CS spline prototype

**Confirmed by the user:** lower-clip observed uncertainty to1 S by default
in fit_conductance_cs_splines.py and expose the value. Implemented
--minimum-std;0 disables. Use max(STD,floor) for valid real observations,
preserving raw dH/dP, invalid masks, quality convention and separate artificial
corner STD. Uniform weighting is unchanged. This local fitting policy is not
a correction of shared physics or approval of a calibrated error model.
67 focused tests and a seven-frame real execution gate pass; intended40-control
fit-quality evaluation remains pending. [[Processing Pipeline]] owns details.

## 2026-10-06 — Add paired WIC/SI13 conductance diagnostics

**Confirmed by the user:** implement the joint Hall/Pedersen analysis discussed
in the conversation. Retain independent camera analyses and use existing paired
crossings. Begin with fixed count-plane conditional quantiles and maps/slices,
paired measurement-error/orbit resampling, and compare single/joint median
predictors on identical wholly held-out orbit data. No Frey curve or ratio
cut. Sparse support and unavailable scores must remain explicit.

**Implementation choices:** default20×20 count cells,20 samples/3 distinct
original orbits,5 held-out folds; configurable. Assume zero cross-camera error
covariance because it is unavailable, and record that limitation. This is a
diagnostic empirical relation, not approval of a new production conductance
model or proof of improved real-data prediction. Method/evidence owner:
[[DMSP Spectrum Conductance Validation]].

## 2026-10-06 — Provide separate detector and CS conductance entry points

User identified the absent detector conductance analysis. Add
conductance_validation.py for detector crossings while retaining the existing
_cs entry point. Share count binning and analytical-error/block-bootstrap code,
with explicit representation-specific geometry and error keys. Detector camera
draws share a pixel/frame, not a CS cell. Keep independent camera selection,
document inherited common support and retain separate detector/CS output paths.

## 2026-10-06 — Restore legacy weight scaling and evaluate plots every30s

**Confirmed by the user:** expose wscaling and make figures default to30s
cadence, holding observations unchanged while the spline updates. The previous
identity-only ridge did not implement the legacy batch default.

**Implementation:** enable normalized coefficient-average method_quality_weight
regularization by default, matching the actual legacy LTL_diag calculation;
allow --wscaling false/--no-wscaling. Use padded quality0.5 and legacy NaN->1.
Apply the same precision matrix to posterior variance, including predictions
between measurement timestamps. Plot observations/STD use latest-at-or-before
holding, without interpolating them; cadence is anchored to the first selected
measurement and stays within the fitted interval. Save output times, held
indices and continuous fit/variance/STD. --plot-cadence exposes seconds.

## 2026-10-06 — Constrain unsupported spline corners to zero

**Confirmed by the user:** pad unsupported corners with zero conductance and
low uncertainty; use0.5 for padded method_quality_weight if used. Separate
Hall/Pedersen figures have2×2 panels (corrected from the original4×4 request),
observed/fit conductance above observed/posterior STD, without colorbars.

**Implementation choice:** identify corner-connected cells without any valid
finite nonnegative H/P over the full input orbit, not just selected frames.
Use configurable0.1S STD; retain raw observations/errors separately and label
padding in masks and plots. Exclude artificial constraints from residual RMS.
Stored quality is0.5 on padding; quality remains unused in fitting. The posterior
is conditional on this imposed boundary assumption. Sensitivity comparison can
disable padding. This does not fill transient or isolated interior gaps.

## 2026-10-06 — Add spline posterior model variance

**Confirmed by the user:** add model variance to the single-orbit conductance_cs
spline prototype. Reuse the fitted sparse precision matrix and save separate
Hall/Pedersen fitted-field variance and STD alongside the unchanged fits.

**Implementation choice:** estimate diag(G A^-1 G.T) with512 batched Gaussian
coefficient draws by default, with exact evaluation available for small checks.
Record sampling precision and seed. Interpret this as the Gaussian posterior
conditional on independent observation errors and the coefficient ridge prior;
it excludes new-observation noise and model discrepancy. Uniform fits assume
1S measurement STD. This remains a diagnostic prototype, not a production
spline-format migration or full-orbit scientific validation.

## 2026-10-06 — Analyze corrected SI13 alongside WIC

**Confirmed by the user:** add SI13 to the conductance-validation script.
Run Hall/Pedersen versus corrected SI13 counts using the same analytical DMSP
errors and orbit bootstrap as WIC. Each camera uses its own counts, coverage
and count-error cohort, preserving WIC selection. Reuse existing crossing
fields; retain separate per-camera outputs and explicit camera provenance.

## 2026-10-06 — Follow published DMSP NaN-error omission

**Confirmed by the user:** follow the published method and ignore NaN channel-
error terms. This supersedes our strict complete-channel-error policy and the
unselected reversion/low-count/moment alternatives. Retain the full nominal
spectrum, propagate available finite nonnegative errors with P/H covariance,
and save channel availability count plus separate uncertainty validity.
All-missing errors remain unavailable; supplied invalid errors and singular
propagation remain invalid. Missing terms are omitted, not a recovered noise
model. Source/method evidence and end-to-end checks are owned by icBuilder's
DMSP Spectrum Conductance Validation note.

## 2026-10-06 — Propagate DMSP spectral errors analytically

**Confirmed by the user:** Implement the proposed analytical propagation,
with independent numerical/Monte Carlo checks, rather than per-record production
Monte Carlo. Reuse the spectrum/MSIS calculation and retain nominal P/H.

**Implementation:** Propagate supplied source per-channel fractional errors
through ionization and equilibrium electron density to P/H standard deviations
and joint covariance. Assume independent energy-channel errors explicitly;
source STD contains its reported calibration contribution, but add no extra
common calibration, model, atmosphere, B or energy errors. Missing/negative
source errors or singular zero-ionization derivatives invalidate only errors.
Keep source low-count NaNs; do not substitute zero or invent a count-error model.
Separate validity flags and uncertainty-cohort diagnostics accompany saved
yearly/crossing outputs and joint Gaussian analysis draws.

The source intentionally leaves errors undefined below one corrected count,
so this is incomplete coverage: Jan1 F13 has889 uncertainty-valid of20,316
nominal-valid spectra. Extending to censored channels requires a separate
scientific low-count model. No such model or production corpus rerun was adopted.

## 2026-10-05 — Reuse MSIS profiles within spatial and hourly bins

**Confirmed by the user:** Use 2° geodetic latitude, 5° geodetic longitude,
and 1-hour UTC bins to reduce repeated MSIS evaluation of 1-second DMSP data.

**Implementation:** For each occupied bin within a satellite-day, evaluate
the neutral profile at its geographic cell center and median valid-sample
time. Keep each spectrum and IGRF magnitude individual, share the neutral
profile/deposition kernel, and record the representative coordinates/time.
Retain per-measurement mode for sensitivity comparisons; record bin widths
in yearly provenance and restart identity. A single-day benchmark reduces
20,316 profiles to 868 and full processing from about 49 s to 6.3 s. This
does not establish all-season approximation accuracy.

**Equivalent speedups:** Keep those bin dimensions, and reduce conversion/
wrapper overhead using exact integral CDF_EPOCH arithmetic, vectorized ApexPy
reference-height algebra and the same C gtd7 backend. Serial old/new comparisons
match every scientific array exactly; fifty focused index/integration tests
pass. Two additional seasonal/disturbed challenge samples support the selected
resolution without establishing a universal binning-error bound.

## 2026-10-05 — Add spectrum-inferred DMSP conductance to yearly files

**Confirmed by the user:** Translate the provided IDL spectrum calculation
into icPhysics and integrate it into `process_dmsp_yearly.py`, with a
reference magnetic-field model such as IGRF. Investigate the missing MSIS
table before interpreting or regenerating the analysis dataset.

**Implementation:** Keep physics numerical and array-based in icPhysics;
keep CDF handling, footprint mapping, IGRF, atmosphere generation and yearly
I/O in icBuilder. Retain model and driving-input provenance.

**Follow-up implementation:** The user requested local historical index
series for measurement-specific MSIS. Use daily observed prior-day F10.7,
current-day centered observed 81-day F10.7 and daily Ap from CelesTrak, stored
locally with source OBS/INT qualifiers and hashes. NRLMSISE-00 profiles vary
at each measurement time/footprint; explicit fixed profiles remain available.
This is a successor-model workflow, not recovery of the original MSIS-90
table. No IMAGE conductance estimator changed.
See [[DMSP Spectrum Conductance Validation]].

## 2026-10-01 — Apply optional fixed Gaussian smoothing inside Product 2

**Decision:** Keep Product 1 as the unsmoothed coregistered observation
boundary. Apply optional per-sensor Gaussian smoothing in detector Product 2,
enabled by default with WIC sigma 0.8, SI12 disabled, and SI13 sigma 1.2 in WIC
detector pixels. Do not create a separate smoothed observation product. Remove
the experimental boxcar and adaptive processing paths.

Convert each sensor's fuvpy background-fit weight to an auroral residual score
as `1 - dgweight`. Propagate that score for each smoothed sensor as a
normalized Gaussian-weighted geometric product over the same contributors,
using a recorded `1e-6` floor, and then multiply the working WIC, SI12, and
SI13 scores. Store only the combined `method_quality_weight` in Product 2;
retain the separate source weights in Product 1. Keep unsubtracted-product
quality uniform because the background-fit Tukey weight does not describe
those counts.

**Rationale:** Fixed smoothing is cheap enough to repeat in Product 2 and
avoids another stored product. Sensor multiplication preserves the intended
strict all-sensors corroboration requirement, while the geometric spatial
propagation prevents a weak auroral score from being spread as readily as it
would by an arithmetic mean and does not depend on the number of kernel pixels
as a raw product would. The schema-4 boundary and
explicit provenance prevent old unsmoothed/schema-3 products from being
silently treated as the new default.

**Downstream amendment:** Carry the complete count-source, smoothing,
variance-propagation, and Tukey-weight recipe into detector Product 3 and
Product-3 CS. Bump those contracts to schemas 3 and 2, respectively, while
keeping their readers compatible with legacy schemas 2 and 1. The combined
quality field is physically carried downstream, so its construction must be
understandable without reopening Product 2.

**Support amendment:** Define `method_quality_weight` wherever all three
working sensor weights have finite common detector support, independently of
ratio, E0, Fe, or conductance success. Reduce it into Product-2 CS and
Product-3 CS on that same independent finite support. Retrieval and
conductance validity remain separate masks.

**Semantic amendment:** `method_quality_weight` is high-is-auroral, not the
original background-fit inlier weight. icAnalyzer must use it in the stored
orientation and must not apply another `1 - weight` transformation. Because
artifacts can also receive low fuvpy weights, interpret it as a corroborated
three-sensor residual score rather than a calibrated probability.

## 2026-09-20 — Freeze detector CS post-processing at 46 by 46

**Decision:** Produce `precipitation_cs` and `conductance_cs` as a shared
post-processing stage on one explicitly frozen 46-by-46 grid. Reuse the WIC
detector-footprint mapping and field reducers, but reduce conductance directly
from `conductance_detector` rather than calculating it from binned
precipitation.

Defer the legacy 36-by-36 Zhang--Paxton lookup and bin-first conductance
infrastructure. Their redesign is required later but does not block the
image-ratio and Robinson CS products.

**Rationale:** Current binned IMAGE and conductance products are 46 by 46, and
Products 2 and 3 share identical WIC detector geometry. Freezing the grid
prevents secsy-version drift; sharing the mapping avoids duplicate geometry
work; binning after the nonlinear Robinson calculation preserves the accepted
detector-first architecture.

**Amendment, 2026-09-21:** Validate only that the CS grid is 46 by 46. Do not
hash or require exact byte-level equality of derived coordinates. The exact
coordinate hash blocked production across otherwise valid numerical
environments before the server could process orbit 0085.

Implementation details and verification gates are maintained in
[[Detector CS Post-processing Implementation Plan]].

## 2026-09-19 — Replace the abandoned detector preprocessing prototype

**Decision:** Treat only Ohma's legacy fuvpy implementation and the new fuvpy
implementation currently on `background_model_asymmetry` as meaningful
preprocessing states. Remove the internal `current_fuvpy_v1` prototype from
live icBuilder code, tests, and current documentation. Do not preserve it as a
runtime compatibility mode or create parallel Product-1 branches for it.

Only the new fuvpy implementation is a detector Product-1 ingestion target.
Legacy Ohma fuvpy remains available for comparison and historical validation;
do not implement a legacy detector-product loader or conversion path.

Implement one new BS-directional detector Product-1 contract using `dgimg`,
`dgweight`, detector measurement variance, and source frame quality. The
detailed implementation and verification sequence is maintained in
[[New fuvpy Detector Product Implementation Plan]].

**Rationale:** The internal prototype no longer exists as a supported
implementation and has no meaning to users outside its development history.
Retaining its labels and field assumptions would add dead complexity and make
the current preprocessing boundary harder to understand. Versioning should
describe real product contracts, not preserve abandoned scaffolding.

## 2026-09-03 — Keep Product-1 restart status direct

**Decision:** Follow the existing `binned_file_status()` pattern for detector
Product 1. In one direct function, check product type/schema, the requested
preprocessing label and time tolerance, selected source paths and image fields,
and required array shapes. Keep the four statuses `missing`, `invalid`,
`mismatch`, and `complete`.

Do not add a processing-configuration document, fingerprint layer, detailed
status API, or source-file hashing to restart discovery. Source SHA-256 values
and software identity remain stored in the product as provenance, but do not
turn the orbit script into an integrity-validation framework. Explanatory
attributes do not decide compatibility.

**Rationale:** This is a scientific processing script used by a small research
group. The earlier validation code obscured the processing workflow and
duplicated information already recorded in the NetCDF file and source-control
history. Restart status only needs enough evidence to avoid mistaking an
obviously incomplete or differently configured orbit file for finished work.

## 2026-09-01 — Make the publication pipeline detector-first

**Decision:** Treat WIC detector geometry as the canonical geometry for
Products 1--3. Correct each sensor on its native detector grid, coregister
SI12 and SI13 onto WIC detector pixels, then calculate precipitation and
conductance before mapping completed products to a fixed Cubed-Sphere grid.

The regular-grid products are derived representations. In particular,
`conductance_cs` is the binned form of `conductance_detector`; it is not
conductance recalculated from `precipitation_cs`. A separate SI Cubed-Sphere
stage is removed from the target architecture. Storage is not a limiting
design constraint.

**Rationale:** Coregistration is an instrument-geometry problem, whereas the
Cubed-Sphere grid is an analysis representation needed by the VAE, covariance,
and spline workflows. Separating them prevents an SI interpolation grid from
becoming part of the physical retrieval, preserves the best available sensor
geometry, and allows detector products to be rebinned without repeating the
physics. Binning after the nonlinear conductance model also avoids assuming
that spatial averaging and the forward model commute.

The fixed grid, file names, schemas, uncertainty propagation, and arbitrary-
MLT Zhang--Paxton lookup remain open implementation decisions. See
[[Detector-First Product Architecture]].

## 2026-08-05 — Use orbit products as restart records

**Decision:** Resume conductance processing from structurally valid orbit
NetCDF files. Write each new product to a same-directory partial file,
validate it after closing, and atomically rename it to the final name. Do not
maintain a separate completed-orbit ledger.

**Rationale:** The product is the authoritative evidence that an orbit
finished. A concurrent ledger needs locking and can disagree with the file if
a process stops between the two writes. Atomic publication leaves a final
name only after a valid save, while missing or invalid outputs naturally enter
the next run. `--overwrite` remains available for deliberate recomputation.

## 2026-08-05 — Keep Kp validation scientific and replaceable

**Decision:** Bundle definitive GFZ Kp through 2003-07-31, validate its
cadence, status, dimensions, and physical range, and calculate the loaded
file's checksum for provenance. Do not reject a valid replacement because it
does not match a checksum or record count hard-coded in the source.

**Rationale:** The original 2000--2001 bundle caused the full IMAGE run to
fail on its first 2002 frame. Git versions the downloaded JSON, while a
dynamically recorded checksum preserves exact product provenance without
turning an intentional time-range extension into coordinated code changes.

## 2026-07-31 — Keep the collapse calculation separate from diagnostics

**Decision:** Keep the reusable Zhang--Paxton latitude reduction as direct
NumPy functions returning ordinary dictionaries. Keep MLT diagnostic
sampling, plots, figure output, and the command line interface in the
diagnostic script.

**Rationale:** The previous 1,065-line package module combined the scientific
calculation with three dataclasses, extensive type machinery, four figure
products, internal checks, and CLI handling. That obscured the short
scientific sequence a student needs to inspect. The separated code retains
meaningful domain and shape checks and batching of the expensive model
evaluation, while old/new numerical comparisons confirm that the scientific
outputs and bundled lookup did not change.

## 2026-07-30 — Use a one-sided conductance uncertainty at zero flux

**Decision:** When proton subtraction clips electron energy flux to zero,
report dP and dH as the conductance reached at the one-sigma upper flux
`Fe=dFe`.

**Rationale:** Robinson conductance is proportional to `sqrt(Fe)`, making its
flux derivative singular at `Fe=0`. Linear propagation is therefore undefined
at that physical boundary. The one-sided excursion is finite for every E0 and
states directly what is being reported. E0 uncertainty contributes nothing at
exactly zero flux because conductance is zero for every E0 there.

## 2026-07-30 — Integrate definitive GFZ Kp and replace IMAGE E0

**Decision:** Use the fixed definitive GFZ Kp series to select collapsed
Zhang--Paxton E0/dE0 for every retained IMAGE frame. Keep the current
three-camera population and SI13 diagnostics during the first comparison
stage, but do not allow SI13 to affect E0 or Fe.

**Rationale:** This implements the user's decision to replace the unvalidated
WIC/SI13 energy inversion without simultaneously changing data support. Kp is
matched to enclosing half-open three-hour intervals with no temporal
interpolation or gap filling. The original GFZ thirds and rounded lookup layer
are both preserved.

Use the collapsed MLAT-profile spread as dE0 and describe it only as unresolved
latitude variability. Propagate its induced first-order covariance with
WIC-derived Fe through the Robinson relations. Do not invent Zhang--Paxton
coefficient, Kp, or predictive uncertainties.

The old ratio inversion remains only as an explicit regression/comparison
path. Making SI13 optional is a later, separate decision. Electron
mean-energy versus characteristic-energy compatibility remains a scientific
gate before publication use.

## 2026-07-30 — Use direct nearest-hundredth Kp lookup layers

**Decision:** Store one directly evaluated collapse for every Kp value from
0.00 through 9.00 in steps of 0.01. At runtime, round Kp to the nearest
hundredth and select that layer without interpolation. Exact half values round
upward.

**Rationale:** The fixed 36-by-36 Cubed-Sphere grid makes the complete table
small enough to bundle with the package. Direct hundredth-Kp layers remove the
need for a second approximation after Zhang--Paxton's published Kp
interpolation and the nonlinear Q-defined latitude selection. The two spatial
axes are the grid's `(eta, xi)` indices; MLT remains a two-dimensional
coordinate associated with those cells.

The lookup stores only the representative area-weighted E0, profile-derived
dE0, area-weighted median, direct grid coordinates, units, and scientific
provenance. The loader checks the stored coordinates directly against the
active grid. Generation remains parallelizable because a measured layer takes
about 10.7 seconds, but the implementation uses the repository's ordinary
`process_map` pattern rather than custom task, checkpoint, or version
machinery.

## 2026-07-29 — Replace IMAGE-derived electron E0 and dE0

**Decision:** Replace the WIC/SI13-ratio-derived electron E0 and dE0 entirely
with a Zhang–Paxton-based estimate. Do not retain the IMAGE estimate as a
fallback, blend component, or calibration target.

**Rationale:** The purpose of introducing Zhang–Paxton is not merely to fill
low-signal gaps. The user reports that comparison with DMSP did not reproduce
the empirical relationship between the WIC/SI13 ratio and E0 proposed by Frey
et al. (2003). Without validation of that relationship, the current IMAGE E0
and its propagated dE0 do not have an adequate basis.

This decision sets the scientific direction but does not authorize production
integration yet. The replacement still requires:

- confirmation that Zhang–Paxton mean energy is compatible with the energy
  quantity used by the WIC response and Robinson relations;
- use of the spherical-area-weighted E0 spread within the selected
  Zhang–Paxton latitude profile as dE0, explicitly interpreted as uncertainty
  from collapsing away MLAT rather than model-fit error;
- nonzero E0–Fe covariance or ensemble propagation, because Fe is calculated
  using E0; and
- a decision on whether SI13 remains for quality control, validity checks, or
  validation after it leaves the E0 inference.

**Clarification later on 2026-07-29:** Do not add a new precipitation-support
mask; the user does not consider one robustly automatable. Corrected WIC
brightness supplies the spatial amplitude. Do not introduce DMSP into the new
estimator or uncertainty calculation. The DMSP result is the reason to reject
the Frey relation, not an input to its replacement.

## 2026-07-29 — Precompute the collapse on the fixed IMAGE grid

**Decision:** If the Zhang–Paxton collapse is integrated, calculate it once
for each supported Kp state on the fixed 36-by-36 Cubed-Sphere grid and use a
lookup during IMAGE processing.

**Rationale:** Every IMAGE frame uses the same grid-cell MLT and MLAT
coordinates. Repeating the expensive 0.01-degree latitude calculation per
frame would reproduce identical results. A small Kp-indexed table turns the
production calculation into array indexing while preserving the fine-grid
offline result.

**Superseded in part on 2026-07-29:** The lookup architecture did not initially
settle the integration role. The user subsequently chose full replacement of
IMAGE-derived E0 and dE0. Compatibility with the downstream energy definition
remains unresolved.

## 2026-07-29 — Use absolute oval thresholds and report area mean and median

**Decision:** Implement a provisional `ZP(Kp, MLT, MLAT)` collapse by selecting
the contiguous Q-above-threshold interval around the principal Q maximum and
computing the exact spherical-area-weighted mean electron energy inside it.
Retain that mean as the primary representative and report the area-weighted
median alongside it.

**Rationale:** Averaging the full latitude slice would let the extensive
near-zero region outside the auroral oval dominate. Exact latitude-cell areas
also avoid overrepresenting equal-width cells near the pole. The median gives
a robust typical-area diagnostic without replacing the spatial mean.

The default `Q > 0.05 mW m-2` follows the Zhang–Paxton Figure 8 mean-energy
criterion. `Q > 0.25 mW m-2` remains an absolute sensitivity check. The user
rejected and removed the relative 10%-of-peak threshold because its meaning
changes with every slice. This decision defines the exploratory reduction
only; it does not approve conductance-pipeline integration.

Diagnostic MLT sampling is 0.05 hour. This is numerical sampling of a
continuous Fourier model, not a claim of 0.05-hour empirical resolution; the
published fit used 0.5-hour observational sectors.

Default MLAT sampling is 0.01 degree. This numerical oversampling of the
continuous Epstein profiles removes the threshold-selection jitter seen at
0.25-degree spacing; it is not a claim of 0.01-degree empirical resolution.

## 2026-07-26 — Treat the project as complete for now

**Decision:** Record `icBuilder` as completed for the current scope while
preserving an explicit need to revisit it later.

**Rationale:** The user directly confirmed both points. No revisit date or
trigger was supplied, so none is inferred.

## 2026-07-26 — Preserve and modernize the established vault

**Decision:** Keep `log/icBuilder/` as the project vault, add the shared live
note architecture within it, and retain the older dated notes and attachments
as historical evidence.

**Rationale:** The repository already used this Obsidian-backed location. A
second renamed vault would split project memory and break continuity.

The architecture migration is documentation-only and does not imply
scientific, dataset, or implementation progress.

**Superseded 2026-07-26:** The user subsequently standardized all project
memory at the repository-root `vault/` path. The complete existing vault was
moved intact; only its location and operational references changed.
