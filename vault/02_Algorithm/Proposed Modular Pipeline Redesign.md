# Proposed Modular Pipeline Redesign

Last reviewed: 2026-09-01
Status: Implemented bin-first redesign retained as a description of the live
`modular_pipeline` branch. Its future architecture is superseded by
[[Detector-First Product Architecture]], which coregisters SI onto WIC before
precipitation and delays Cubed-Sphere binning until after the detector-level
physical products.

## Motivation

The current orbit pipeline combines sensor binning, proton correction,
electron energy and flux inference, Robinson conductance conversion, and
serialization in one processing path. This makes a change to one scientific
model look like a reason to regenerate every upstream result.

The user proposes a clean restart that reuses the verified calculations but
separates stable observations from replaceable scientific models. The goal is
not a wholesale rewrite: retain useful binning, response, Zhang--Paxton, and
Robinson calculations while giving each stage a clear input and output.

## Proposed product stages

### 1. Binned FUV observations

Write a model-independent per-orbit product containing the binned WIC, SI12,
and SI13 observations, count uncertainties, masks and support, timestamps,
grid coordinates, viewing geometry, and preprocessing provenance.

This product should contain observations and instrument preprocessing, not a
choice of electron-energy or conductance model. In particular, proton-corrected
WIC and SI13 depend on an assumed proton energy and response model; storing
them in the next stage prevents a changed proton assumption from invalidating
or disguising the reusable binned observations.

The three channels do not need to be forced onto identical frame support in
this product. Later methods can select the channels they require.

The existing `BinnedImage` represents one sensor, not a combined WIC/SI12/SI13
orbit. Preserving that meaning suggests one binned file per sensor and orbit,
or otherwise distinct sensor time dimensions. A single three-camera time axis
would preserve the current SI13 restriction and prevent the Zhang--Paxton path
from using otherwise valid WIC/SI12 frames.

The implemented Product-1 layout is:

```text
binned/wic/or_XXXX.nc
binned/si12/or_XXXX.nc
binned/si13/or_XXXX.nc
```

WIC remains on its native 36-by-36 grid and SI12/SI13 on their native
18-by-18 grid. Each sensor orbit is selected and processed independently.
Files contain binned signal and uncertainty, weights, native-grid source
counts, sensor time, SZA/DZA/LOS diagnostics, subsolar longitude, correction
provenance, and reconstructable grid metadata. Kp is intentionally excluded;
it belongs to precipitation inference. `BinnedImage` is strictly native-grid
and contains no interpolation path. Its NetCDF stores the actual xi, eta,
MLAT, and MLT arrays alongside integer source-pixel counts.

### 2. Precipitation inference

Read the binned product, perform the SI12-based proton correction, and store
the corrected WIC and SI13 signals together with electron energy and flux.
Support at least two explicit methods:

- `image_ratio`: infer E0 from corrected WIC/SI13 and then infer Fe;
- `zhang_paxton`: obtain E0 from Kp and MLT and infer Fe from corrected WIC,
  retaining corrected SI13 as an optional diagnostic rather than a required
  input to E0.

The methods should produce separate method-labelled products or directories
rather than silently replacing variables in one ambiguous file. This also
allows their frame support to differ: the ratio method needs WIC, SI12, and
SI13, whereas the Zhang--Paxton method principally needs WIC and SI12.

Each product should preserve method, coefficients or lookup provenance,
assumed proton energy, masks, units, uncertainties, covariance terms, and the
identity of its binned source product. Corrected WIC/SI13 belong here because
they make the energy/flux calculation inspectable without repeating binning.

`PrecipitationImage` is implemented as the initial Product-2 boundary. It
matches only the sensors required by each method, records source indices, and
regrids SI data onto the WIC grid. Because both grids are regular and exactly
nested in xi/eta, the implementation uses explicit four-cell bilinear weights
rather than a general scattered-point triangulation. All four source values
must be finite. Targets in the outer half of a physical boundary cell use that
nearest SI cell, while internal gaps and targets outside the physical grid
remain NaN. Variance is propagated with squared bilinear weights internally;
boundary targets inherit the source-cell uncertainty. Each array is
interpolated independently;
the interpolation stage does not impose matching finite support between signal
and uncertainty. Sample counts remain only in Product 1.

The class accepts an injected physics function for controlled comparisons and
otherwise selects the corresponding shared icPhysics precipitation function.
No equations are copied into `PrecipitationImage`, and it contains no
Hall/Pedersen calculation.

### 3. Conductance forward model

Read a precipitation product and apply a named forward model to
`(E0, Fe)` and their uncertainty information. Robinson is the initial model,
but it should be replaceable without rerunning sensor binning or precipitation
inference.

The corresponding calculation must also exist as a small public importable
function, independent of orbit files and orchestration. A simple scientific
API should accept arrays and return ordinary arrays or a dictionary containing
Hall, Pedersen, and propagated uncertainty. Avoid factories or a generalized
plugin framework until a second model demonstrates what interface is actually
shared.

Conductance products should identify their precipitation source and forward
model. This yields the dependency chain:

```text
binned FUV
    -> precipitation/image_ratio
    -> precipitation/zhang_paxton
        -> conductance/robinson
        -> conductance/another_model
```

Changing a downstream model then regenerates only its own stage.

`ConductanceImage` remains suitable for this stage after it is changed to
consume a precipitation product rather than three `BinnedImage` objects. A
fourth `SplineImage` product remains explicitly deferred.

### NetCDF boundary

Each product class should have its own `to_nc` implementation and later its
own corresponding reader because each schema has different scientific
variables and provenance. Reusing the method name is useful; sharing one large
generic serializer is not. A small icBuilder-only helper may write repeated
time, grid, weight, Kp, source-product, units, and product-type metadata.
Avoid a base product class or inheritance hierarchy.

## Shared use by icAnalyzer

**Confirmed dependency constraint:** icAnalyzer must not depend on the full
icBuilder package. icBuilder requires pipeline packages such as `fuvpy` that
are not needed in the VAE environment and should not propagate into
icAnalyzer.

The numerical proton correction, camera-response conversion, energy/flux
inference, and conductance forward models should therefore live in the new
`icPhysics` package below both projects, rather than in scripts,
`ConductanceImage` methods, or the icBuilder distribution. `icBuilder` and
`icAnalyzer` can then call exactly the same functions for orbit processing and
each VAE or Monte-Carlo realization.

The intended dependency direction is:

```text
ZhangPaxton2008
        |
        v
icPhysics
        |                         |
        v                         v
icBuilder                     icAnalyzer
(fuvpy, ApexPy, I/O)          (VAE/analysis stack)
```

Keep `fuvpy`, ApexPy coordinate processing, file reading, NetCDF writing,
orbit iteration, Kp acquisition/matching, and multiprocessing in icBuilder.
The shared functions should have lightweight NumPy/SciPy-level dependencies
and make units, masks, and uncertainty inputs explicit. They may accept Kp as
a numerical input, but should not download or independently assign Kp to
IMAGE frames; icBuilder should serialize the matched value for icAnalyzer.

The shared package should remain small and scientific: ordinary public
functions, no file-product classes, framework, plugin system, or dependence on
either consuming project. ZhangPaxton2008 remains its own focused package;
the shared layer may depend on it rather than copying its equations.

## Probabilistic extension

The longer-term proposal is to train a joint WIC/SI12 VAE and pass sampled FUV
realizations plus sampled Zhang--Paxton E0 through the same precipitation and
conductance functions. The resulting Hall/Pedersen ensemble is the primary
probabilistic product; means, joint covariance blocks, low-rank forms, and
Gaussian approximations are derived representations.

The proposed common log-E0 amplitude perturbation is a transparent rank-one
model, but Zhang--Paxton latitude spread is not measured event-to-event model
error. Until calibrated, describe it as a profile-spread-derived
model-discrepancy or sensitivity prior, not a formal Zhang--Paxton predictive
uncertainty. Preserve the current arithmetic-mean nominal collapse, use
spherical-area weighting, avoid over-weighting wide MLT slices, retain Kp
dependence, and explicitly choose whether the nominal curve is the ensemble
mean or median.

The VAE likelihood must likewise treat current WIC/SI12 uncertainty estimates
as estimates requiring validation, distinguish latent auroral signal from
observation noise, and respect the channels' different effective resolutions.

## Scientific gates retained

- Confirm whether Zhang--Paxton mean energy is compatible with the energy
  definition used by the camera response and conductance model.
- Decide the Zhang--Paxton oval threshold and latitude domain before using its
  profile spread probabilistically.
- Define the assumed proton energy and its uncertainty explicitly.
- Decide how zeros, background, and missing support are represented in
  log-conductance statistics.
- Preserve the Hall--Pedersen cross-covariance rather than calculating only
  separate Hall and Pedersen covariances.
- Do not treat a standard VAE encoder as a calibrated partial-image posterior;
  sparse conditioning needs an explicit masked-observation inference method.

## Next design decision

The `PrecipitationImage` API, NetCDF schema, and method-specific orbit builder
are now implemented. Continue the smallest prototype as:

1. load a precipitation product through `icReader.load()`;
2. generate a Robinson conductance file from that precipitation file without
   repeating binning or precipitation inference;
3. call the same precipitation and Robinson functions from a small
   icAnalyzer-side test.

No full-corpus regeneration or VAE redesign should begin until this boundary
and the variables owned by each stage are agreed.
