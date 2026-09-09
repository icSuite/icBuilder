# Decision Log

Last reviewed: 2026-09-03

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
