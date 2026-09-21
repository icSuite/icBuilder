# Detector CS Post-processing Implementation Plan

Last reviewed: 2026-09-20
Status: Implemented; corrected full-orbit set 0085, 0086, 0261, and 0968 verified

## Objective

Create the fixed-grid post-processing stage for detector Product 2 and Product
3. The stage consumes completed `precipitation_detector` and
`conductance_detector` orbit files and writes analysis-ready
`precipitation_cs` and `conductance_cs` products on one explicitly frozen
46-by-46 Cubed-Sphere grid.

The two outputs share the same detector-footprint mapping and reduction code.
They remain independent scientific products: `conductance_cs` is reduced
directly from `conductance_detector`, never recalculated from
`precipitation_cs`.

## Settled decisions

- This is post-processing. Detector Products 2 and 3 remain the canonical
  results of the physical retrieval and conductance calculation.
- The initial detector-CS grid is 46 by 46, matching the current binned FUV
  and conductance products.
- Freeze the grid edges, projection, radius, reference height, coordinates,
  and a durable `grid_id`; do not let the installed secsy version silently
  choose the grid dimensions.
- Build one WIC-detector-footprint mapping per frame and reuse it for all
  Product-2 and Product-3 fields needed for that frame.
- Write separate `precipitation_cs` and `conductance_cs` files with their own
  product identity, source identity, schema version, and restart validation.
- Use detector Product 3 as the source of binned Hall and Pedersen conductance.
  Do not apply Robinson to binned E0 and Fe.
- The legacy 36-by-36 Zhang--Paxton lookup and old bin-first conductance path
  are outside this implementation. They require a later redesign but do not
  block the image-ratio detector-CS products.

## Proposed file structure

```text
icbuilder/detectorcs.py
icbuilder/precipitationcs.py
icbuilder/conductancecs.py
scripts/pipeline/make_detector_cs_orbit_files.py
tests/test_detectorcs.py
tests/test_precipitationcs.py
tests/test_conductancecs.py
```

`detectorcs.py` owns only the common geometry and numerical reducers.
`precipitationcs.py` and `conductancecs.py` keep their field contracts and
NetCDF serialization explicit. The combined orbit runner coordinates both
products so their common mapping is calculated once.

## Frozen grid contract

Add a detector-CS grid constructor without changing the unresolved legacy
Zhang--Paxton grid path. The constructor uses explicit edges and validates only
that the resulting shape is `(46, 46)`. The NetCDF files store the grid
coordinates; no coordinate hash or exact byte-level coordinate comparison is
part of the contract.

## Shared reduction engine

For each frame:

1. read the WIC detector MLAT/MLT geometry;
2. call the existing footprint-overlap geometry code once;
3. retain the sparse overlap-area matrix and target-cell area;
4. apply field-specific reducers grouped by common validity masks;
5. discard the mapping before moving to the next frame; and
6. write the small CS result arrays after the complete orbit is reduced.

The common module should expose direct functions rather than a generic plugin
framework. Suggested operations are:

```text
make_detector_cs_mapping(mlat, mlt, grid)
reduce_area_mean(values, valid, mapping)
reduce_measurement_variance(variance, valid, mapping)
reduce_covariance(covariance, valid, mapping)
reduce_flag_fraction(flag, valid, mapping)
reduce_support(valid, mapping, cell_area)
```

Where several continuous fields share a validity mask, reduce them as one
matrix block so the overlap denominator is calculated once.

## Reduction semantics

### Continuous fields

For overlap areas `a_i` and finite valid source values `x_i`, store the
renormalized overlap mean:

```text
x_cs = sum(a_i x_i) / sum(a_i)
```

This applies to detector-derived intensive quantities including E0, Fe, Fp,
P, H, and continuous diagnostic fields. Do not clip the resulting mean inside
the generic reducer. Physical clipping, if required, belongs to the source
product contract.

### Measurement uncertainty

For independent detector-pixel variances `v_i`, propagate the uncertainty of
the overlap mean as:

```text
v_cs = sum(a_i^2 v_i) / sum(a_i)^2
d_cs = sqrt(v_cs)
```

Input `d*` fields must therefore be squared before reduction. This is not the
within-cell spatial spread. Store any spatial-spread diagnostic under a
separate, explicit name and never substitute it for propagated measurement
uncertainty.

The independence assumption is incomplete because SI coregistration induces
correlation between neighbouring WIC pixels. The required covariance is not
available in the current detector products. Record this limitation in both CS
products; do not silently claim complete uncertainty.

### Covariance

Propagate `varE0Fe` using the same squared normalized overlap weights:

```text
cov_cs = sum(a_i^2 cov_i) / sum(a_i)^2
```

Do not introduce unmodelled Hall--Pedersen covariance.

### Support, masks, quality, and flags

- Store valid covered area divided by target-cell area as `coverage`.
- Store the number of valid detector footprints intersecting each cell as
  `source_count`.
- Retain values for any positive valid overlap during the first implementation
  and expose coverage so downstream selections are explicit. Revisit a minimum
  coverage threshold only from edge tests and diagnostics.
- Store quality as an area-weighted diagnostic. Do not multiply it into the
  physical spatial mean unless a separately justified estimator is adopted.
- Store clipping and similar Boolean state as both an area fraction and an
  any-contributor flag where interpretation benefits from both.
- Keep central validity separate from uncertainty validity.

Constant-field reproduction and area-integral accounting are separate tests.
Renormalized means must reproduce a constant field; coverage must reveal
partial observation rather than making it look complete.

## Product-2 CS contract

The first `precipitation_cs` schema should retain:

- Product-2 time and source-frame identity;
- fixed grid coordinates and identity;
- E0, Fe, Fp, Ep/Ep_model, their available uncertainty, and varE0Fe;
- selected corrected-channel or ratio diagnostics needed to audit the
  retrieval, after a field-by-field relevance review;
- method-valid and uncertainty-valid coverage;
- method quality, source count, coverage, and proton-energy clipping fraction;
- precipitation, proton-energy, and count-uncertainty provenance; and
- exact detector Product-2 source identity and software/schema provenance.

E0 and Fe are reduced from detector E0 and Fe. They must not be recalculated
from spatially averaged counts or ratios.

## Product-3 CS contract

The first `conductance_cs` schema should retain:

- Product-3 time and the same fixed grid identity;
- P, H, dP, and dH reduced directly from detector Product 3;
- conductance central-valid and uncertainty-valid coverage;
- method quality, source count, and coverage;
- the binned E0/Fe state needed to interpret P/H, either copied from the common
  in-memory Product-2 reduction or independently reduced with the same mapping;
- conductance and upstream precipitation/proton provenance; and
- exact detector Product-3 source identity plus its Product-2 linkage.

The output must state explicitly that averaging detector conductance is not
equivalent to evaluating Robinson from averaged precipitation.

## Orbit runner

`make_detector_cs_orbit_files.py` should:

- default to producing both CS products for `IR_hardy/robinson`;
- accept `--base-input` and `--base-output`, with output defaulting to input;
- support repeated `--orbit` selection and orbit-level `--workers`;
- allow one requested output when only one downstream representation is
  wanted, while keeping both as the default;
- validate that Product 2 and Product 3 have identical time axes and detector
  geometry and that Product 3 identifies the selected Product 2;
- build the common frame mapping only once when both outputs are pending;
- use same-directory partial files and atomic replacement;
- classify outputs as missing, invalid, mismatched, or complete;
- treat a changed source schema, source identity, grid identity, or reduction
  contract as invalid or mismatched rather than silently skipping it; and
- remain safe under Python multiprocessing spawn/forkserver semantics.

Suggested layout:

```text
precipitation_cs/image_apex_130km_46x46_v1/IR_hardy/or_XXXX.nc
conductance_cs/image_apex_130km_46x46_v1/IR_hardy/robinson/or_XXXX.nc
```

## Verification

### Unit tests

1. The fixed grid is exactly 46 by 46.
2. A constant detector field remains constant wherever coverage is positive.
3. A hand-calculated two-footprint example verifies mean, coverage, count,
   variance, covariance, and flag fractions.
4. Missing values use field-specific denominators and do not contaminate other
   fields.
5. Measurement uncertainty and within-cell spatial spread remain distinct.
6. Partial edge coverage is reported rather than hidden by renormalization.
7. Product-2 and Product-3 reduction in one frame invokes the mapping builder
   once.

### Product tests

1. Synthetic detector Products 2 and 3 round-trip into self-describing CS
   products with the expected dimensions, units, masks, and provenance.
2. Old schema, wrong grid, changed source, truncated variables, and mismatched
   Product-2/Product-3 pairs are rejected.
3. Atomic restart skips complete products and regenerates invalid products.
4. A deliberately nonuniform synthetic field demonstrates that direct binned
   detector conductance differs from Robinson applied to binned E0/Fe; the
   stored Product 3 must match the former.
5. Serial and two-worker orbit execution produce identical files apart from
   recorded execution metadata.

### Full-orbit gate

Run orbit 0085 first and require:

- both output shapes equal `(122, 46, 46)`;
- exact agreement with a direct frame-by-frame application of the common
  reducer;
- no central values outside source-field physical ranges apart from expected
  floating-point tolerance;
- uncertainty support consistent with the corrected detector products;
- explicit accounting of uncovered/partially covered cells;
- a successful immediate restart skip; and
- recorded runtime and peak memory.

Then run orbits 0086, 0261, and 0968 with two workers. Audit coverage,
contributor counts, finite support, uncertainty support, and output size before
starting the server corpus.

## Implementation sequence

1. Freeze the 46-by-46 grid and add coordinate-identity tests.
2. Implement and hand-test the common mapping and reducers.
3. Implement explicit Product-2 CS loading, field reduction, serialization,
   and status validation.
4. Implement Product-3 CS using the same frame mapping and reducers.
5. Add the combined restart-safe, parallel orbit runner.
6. Complete synthetic integration and nonlinear-ordering tests.
7. Pass the orbit-0085 execution/performance gate.
8. Pass the four-orbit parallel gate and update the live handoff.
9. Only then prepare the server command and corpus run.

## Deferred work

- Redesign of the Zhang--Paxton method for arbitrary detector geometry or the
  new 46-by-46 CS representation.
- Replacement or removal of the legacy 36-by-36 lookup.
- Migration of the old bin-first Product-2/Product-3 infrastructure.
- `fuv_cs`, unless the joint FUV VAE input design requires it.
- Full detector-to-CS covariance induced by SI coregistration.
- Hall--Pedersen cross-covariance and probabilistic/ensemble products.
- Changes to icReader, the VAE, spline, covariance, or sparse-reconstruction
  consumers until the two new CS schemas are verified.

The Zhang--Paxton redesign will eventually affect production choices for
precipitation and conductance, but it is not a prerequisite for producing the
image-ratio `precipitation_cs` and Robinson `conductance_cs` products defined
here.

## Implementation result

The approved structure is implemented in `icbuilder/detectorcs.py`,
`icbuilder/precipitationcs.py`, `icbuilder/conductancecs.py`, and
`scripts/pipeline/make_detector_cs_orbit_files.py`. The grid constructor in
`icbuilder/grids.py` supplies the 47 xi/eta edges that define the 46-by-46 grid
and checks only the resulting shape.

The paired runner builds one footprint mapping per detector frame and reuses it
for both outputs. It preserves explicit field-specific masks, propagates
measurement variance and E0--Fe covariance with squared overlap weights,
stores valid and uncertainty coverage separately, and records clipping as
fraction plus any-contributor flag. It supports separate input/output bases,
atomic restart, orbit selection, and real two-worker multiprocessing.

Ten focused detector-CS tests pass, including hand-calculated reducers,
grid shape, source-pair rejection, schema/restart behavior, actual
two-worker execution, and a nonlinear fixture proving that Product-3 CS is
the mean of detector conductance rather than Robinson evaluated on binned
precipitation. The full repository suite has 126 passes and the same four
deferred legacy 36-by-36/Zhang--Paxton failures.

The first orbit-0085 implementation was scientifically correct but read every
compressed field one frame at a time. Profiling orbit 0086 isolated about 2.2
seconds per frame in repeated NetCDF decompression while footprint construction
took only about 0.02 seconds. The final implementation builds all frame
mappings, then reads each compressed variable once for the orbit. Optimized
and original outputs are identical for every Product-2 CS and Product-3 CS
variable on orbits 0086 and 0261.

The complete corrected four-orbit gate (0085, 0086, 0261, and 0968) finished
in 39.60 seconds with two workers and 1,614,780-KB reported peak RSS. The 893
frames produced 1,477,851 central-valid CS cells, all with propagated
uncertainty. Product-2 CS totals 90.84 MB and Product-3 CS totals 52.85 MB.
Median valid coverage by orbit is 0.988, 0.994, 0.956, and 1.000; 8,276 cells
have coverage below 0.01 and 28,228 below 0.1. Positive overlap remains stored
without a minimum threshold so those edge slivers are explicit. All four files
skip together on restart in 0.33 seconds.

Direct orbit-0085 frame checks at indices 0, 61, and 121 match the shared
reducer; output E0/Fe/P/H remain inside their detector-source ranges. The
local implementation and performance gates are complete. Server corpus
execution remains the next operational step after the code is committed and
the corrected detector Products 2 and 3 are available there.

### Audit of the user's completed four-orbit pipeline tree

The generated files under
`/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test` were inspected
directly on 2026-09-20 rather than inferred from the scratch performance run.
All four Product-3, Product-2 CS, and Product-3 CS orbit pairs have complete
restart status, matching time/index axes, matching source identities, and the
frozen grid identity. Detector Product 3 copies E0, Fe, their uncertainty, and
covariance exactly from Product 2. Its P, H, dP, and dH agree with a fresh
Robinson evaluation to float32 precision. The two CS products have identical
central and uncertainty masks, precipitation state, coverage, and contributor
counts. Their central and uncertainty fields are finite exactly on their
declared masks and all conductance values are nonnegative.

Seven CS frames are empty across the four orbits. Each corresponding detector
Product-2 frame also has zero valid pixels because SI12 or SI13 is absent, so
this is upstream missing support rather than loss during CS mapping. Visual
checks at each orbit's maximum integrated Hall frame preserve the detector
auroral structure after reduction.

The audit also exposed two science-acceptance issues that are not file-format
or CS-reducer failures. Hardy mean energy falls below the 0.47-keV camera
response boundary for most clipped pixels: 52.1--66.2 percent of detector-valid
pixels are clipped in these orbits, while only 571--2,872 pixels per orbit are
clipped above 46.7 keV. The clipped value is used in the central retrieval and
Hardy energy uncertainty remains unmodelled. In addition, `varE0Fe` is zero
everywhere and 78--97 percent of central-valid CS cells have dP>P or dH>H.
This is consistent with the very weak/zero-flux pixels admitted by the current
method mask and the unstable image ratio, but means the products are largely
uncertainty dominated. Whether clipped-energy and low-signal retrievals should
remain central-valid, be flagged more strongly, or be excluded is the next
scientific decision; it must not be conflated with successful production.
