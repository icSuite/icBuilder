# Hardy Electron Mean Energy - Implementation Plan

Last reviewed: 2026-10-11
Status: Implemented optional full spatial route; scientific qualification remains.
Repository evidence: icBuilder `modular_pipeline` at `eecefeb`; icPhysics
`feat/spectrum-conductance` at `eecc263`, with preserved uncommitted work.

## Accepted scope and implemented calculation

The user requested the option and selected full spatial MLAT/MLT/Kp evaluation
on 2026-10-11. Electron and proton choices are independent. Existing defaults
remain unchanged. The pipeline retains its selected camera counts, smoothing
and SI12 proton correction, then computes:

1. Hardy number and energy flux at each pixel's coordinates and frame Kp;
2. raw mean electron energy E0 = energy flux / number flux, in keV;
3. event electron flux Fe = proton-corrected WIC / WIC response(E0);
4. existing Robinson electron conductance from E0 and Fe.

Hardy's statistical flux supplies the energy ratio only. Observed WIC still
supplies event amplitude. WIC/SI12 define Hardy support; SI13 is diagnostic
and cannot invalidate nominal values or their measurement uncertainty.

## Source and numerical ownership

The user supplied Hardy 1985 and the complete Hardy 1987 functional paper.
Copies are indexed in `literature/README.md`. icPhysics owns the pure numerical
implementation and its source note
`../icPhysics/vault/02_Algorithm/Hardy Electron Model.md` (path from the
repository root).

Hardy 1987, printed pp. 12284–12285, equations 4–5 and tables 4–5 define the
four-parameter Epstein profile and sixth-order Fourier expansions. All 728
coefficients were independently visually audited. This electron profile is
separate from the existing six-parameter Hardy ion model.

`hardy_electron_precipitation(kp, mlt, mlat)` returns number flux, energy flux
and mean energy with NumPy broadcasting. MLT wraps 24 h; northern corrected
geomagnetic latitude 50–90° is supported, with unavailable outputs elsewhere
and no implicit hemisphere reflection. The original spectra span 50 eV–20 keV
and their integrals extrapolate to 100 keV.

The Fourier phase uses `n*pi*MLT/12`: the printed equation omits pi, so this
24-hour normalization is explicitly inferred. The paper's background floors are
log10 flux 6 equatorward and 7 poleward for both flux quantities. Fractional Kp
linearly interpolates evaluated log10 flux after floors; finite Kp clips 0–6.
This interpolation is an implementation convention because the electron paper
does not specify it. The highest published activity group is Kp>=6-.

icBuilder supplies Modified Apex coordinates at 130 km as a documented
approximation to the source's corrected geomagnetic coordinates.

## Selection, products and compatibility

- Detector: `PrecipitationDetector(method="hardy")`; runner option
  `--precipitation-method hardy`. Proton selection remains
  `--proton-energy-model hardy|constant`.
- Bin-first: `PrecipitationImage(method="hardy")` and the equivalent runner
  option. Optional SI13 remains attached for diagnostics.
- Dedicated default retrieval labels begin `HardyElectron_HardyProton` or
  `HardyElectron_ConstantProton` and retain count/smoothing suffixes.
  `README.md` owns exact stage commands and label conventions.
- Product 2, Product 3 and paired CS outputs carry `electron_*` model, source,
  coefficient, Kp, coordinate, uncertainty and response-policy provenance.
  Existing schemas are retained because these additions are optional.
- Detector/bin-first `E0_response_range_flag` is signed int8:
  -1 unavailable, 0 within 0.2–25 keV, 1 below, 2 above. Raw E0 is retained and
  the existing WIC interpolator holds endpoint responses outside the table.
- CS products carry `E0_response_range_fraction` and
  `E0_response_range_any` from valid overlap support.
- icReader exposes provenance and additive fields, preserving older files.
  Restart checks reject missing fields/metadata, malformed shapes and
  electron provenance inconsistent with upstream products.

## Uncertainty interpretation

First version propagates measurement errors conditional on fixed Hardy E0.
Computational dE0=0 and varE0Fe=0 exclude prediction uncertainty; they do not
claim physically exact energy. Corrected WIC STD gives dFe through the same
response. Missing measurement errors preserve nominal precipitation and
conductance; only uncertainty is unavailable. Nonpositive corrected WIC
is unavailable. No DMSP calibration or predictive Hardy error model is added.

The 1987 paper fits conductivities using moments above 661 eV. Full-spectrum
Hardy E0 combined with WIC and Robinson does not reproduce those fitted
conductance maps.

## Verification checkpoint

- icPhysics: 72 tests pass. Independent scalar evaluation of 10,175 points
  agrees within 4.33e-14 relative for mean energy. All 728 coefficients audited;
  an independently found table-transcription digit was corrected before gate.
- Existing image-ratio detector outputs match HEAD exactly across 8 combinations
  of proton choice, count source and smoothing: 800 arrays/attributes compared.
- Builder: 75 focused integration/restart tests plus 18 CS/particle-diagnostic
  tests pass; reader full suite 46 passes. Restart shape checks inspect lazy
  metadata without decompressing detector/CS flag arrays.
- Real orbit 0085, 122 × 256 × 256, completed serially in isolated scratch:
  Product 2 ~40.8 s, Product 3 ~8.2 s, paired CS ~7.9 s. All 756,975 nominal-valid
  detector pixels also have conditional conductance errors. E0 spans
  0.129–9.542 keV; 9,071 valid pixels flag below the response table. Both CS
  products are 122 × 46 × 46; 2,541 flagged CS cells retain overlap fractions.
  Standalone E0/Fe checks agree within 1e-6 after float32 serialization.
  All three runners accept the complete products on restart.
- Particle tools now evaluate uncalibrated pointwise Hardy moments over the
  same complete DMSP temporal windows, with common support and orbit-held-out
  comparisons. Legacy files without Hardy remain supported. Synthetic
  old/new workflows and real two-orbit build/evaluation/plot execution pass.
  The two-orbit smoke test is not a scientific model-selection result.

Scratch evidence: `/tmp/icbuilder_hardy_gate/gate_summary.json` and
`/tmp/icbuilder_hardy_energy_diagnostics/`. Scratch artifacts are disposable.
No tracked example products were regenerated for this task, no full corpus
was processed, and no commit/default-model change was made.

## Next scientific action

Inspect spatial morphology and coordinate sensitivity, then run a broader
accepted DMSP/FAST comparison on matched support before selecting Hardy for
production. Quantify the Modified Apex approximation and investigate the
out-of-table WIC response policy. Predictive energy uncertainty remains a
separate scientific decision. Existing response/calibration and particle
integration-limit differences still apply.
