# Detector-First Product Architecture

Last reviewed: 2026-09-03

Status: The detector-first direction is accepted. Choices labelled proposed
or open below still require a decision or a focused test. This document
describes the intended architecture. Experimental `current_fuvpy_v1`
`fuv_detector` and image-ratio `precipitation_detector` slices now implement
the first two detector-space boundaries. Product 3 and the candidate
Zhang--Paxton Product 2 still use the current bin-first `modular_pipeline`
path.

## Purpose

The current modular pipeline bins each sensor before precipitation and
conductance are calculated. That makes the numerical grid part of the physics
path and forces SI12 and SI13 through their own Cubed-Sphere products before
they are combined with WIC.

The intended architecture instead keeps the physical calculation on the WIC
detector geometry:

1. correct each camera on its native detector grid;
2. coregister SI12 and SI13 onto the WIC detector pixels;
3. calculate precipitation on the WIC detector geometry;
4. calculate Hall and Pedersen conductance on the same geometry; and
5. bin completed products onto a fixed Cubed-Sphere grid only when a regular
   representation is needed.

For a selected and recorded camera-preprocessing configuration, the detector
products preserve the observation geometry and are the canonical scientific
products. The publication preprocessing configuration is not yet selected.
The Cubed-Sphere products are analysis-ready spatial representations for the
VAE, covariance calculations, sparse reconstruction, and splines.

Storage is not a design constraint. Prefer scientifically clear,
independently inspectable products over avoiding repeated coordinates or
fields.

## Decision boundary

| Settled direction | Still open |
| --- | --- |
| Coregister SI12/SI13 to WIC detector geometry | Publication background and camera-preprocessing branch |
| Calculate Products 1--3 before CS binning | Exact sensor time-match and frame-validity rules |
| Derive Product 2 CS and Product 3 CS independently | Detector radiometry, uncertainty, and induced correlation |
| Use one fixed CS grid for analysis products | CS extent/resolution and field-specific reducers |
| Keep physical equations in icPhysics | Grid-independent Zhang--Paxton lookup and coordinate conventions |

The left column defines the architecture. The right column defines the tests
and decisions needed before implementation or publication.

## Terminology

Two independent labels describe a product:

- the **product level** states what physical quantity it contains;
- the **representation** states which spatial geometry it uses.

The proposed explicit names are:

| Level | Detector representation | Cubed-Sphere representation |
| --- | --- | --- |
| Product 1: FUV observations | `fuv_detector` | `fuv_cs` if needed |
| Product 2: precipitation | `precipitation_detector` | `precipitation_cs` |
| Product 3: conductance | `conductance_detector` | `conductance_cs` |

Names such as Product 2b and Product 3b are useful in conversation, but the
stored product should say what it contains and which geometry it uses.

`detector` means the time-dependent WIC detector geometry. Detector row and
column are array indices, not a fixed physical grid. Every frame must retain
its own geographic and magnetic coordinates and viewing geometry.

`cs` means one named, fixed Cubed-Sphere grid shared by the complete dataset.
Every cell then has the same meaning in every frame, which is required for the
VAE, spatial covariance, and spline model.

## Dependency graph

```text
native WIC -----------+
native SI12 ----------+--> Product 1: fuv_detector
native SI13 ----------+          |
                                  | proton and electron models
                                  v
                         Product 2: precipitation_detector
                                  |
                                  | conductance forward model
                                  v
                         Product 3: conductance_detector

Optional fixed-grid representations:

fuv_detector --------> fuv_cs
precipitation_detector --> precipitation_cs
conductance_detector ---> conductance_cs ---> VAE / covariance / splines
```

The important point is that the two lower Cubed-Sphere products are separate
branches:

```text
precipitation_cs = bin(precipitation_detector)

conductance_detector = forward_model(precipitation_detector)

conductance_cs = bin(conductance_detector)
```

In general,

```text
bin(forward_model(precipitation_detector))
    != forward_model(bin(precipitation_detector))
```

because the conductance conversion is nonlinear. `precipitation_cs` must
therefore not become the source of the canonical `conductance_cs` product.

## Upstream sensor files

The existing per-sensor orbit files remain the input boundary. They contain
the camera images and time-dependent detector geometry supplied by the IMAGE
and fuvpy processing.

Background removal, flat-field treatment, and other camera-specific
corrections should normally happen on each native detector grid before
coregistration. A future joint background method would be a new, explicitly
labelled Product-1 preparation path rather than an undocumented change inside
the precipitation calculation.

These source files are not discarded after Product 1 is written. They retain
the native SI geometry needed to audit or repeat the coregistration.

## Product 1: coregistered FUV observations

### Scientific role

Once its camera-preprocessing configuration has been frozen, Product 1 is the
stable observational boundary between camera preprocessing and replaceable
precipitation models. It contains WIC, SI12, and SI13 on the same WIC detector
pixels without applying a proton-energy or electron-energy model.

Current fuvpy, fixed FUVVIEW, active `p`, and active `p2` background branches
produce materially different ratios. None of the historical alternatives has
yet justified replacing current fuvpy. Each retained preprocessing experiment
must therefore have a separate Product-1 label; no branch is silently the
universal canonical observation.

### Proposed frame support

WIC defines the Product-1 time axis. SI12 and SI13 are matched independently
to each WIC frame within a documented tolerance.

- A missing SI frame does not delete the WIC frame.
- An unmatched channel is represented by NaN values, zero coverage, and a
  source index of `-1`.
- Product 2 decides which channels are required by its selected method.

This keeps Zhang--Paxton frames from being lost merely because SI13 is absent,
while retaining SI13 wherever it exists for the image-ratio method and for
diagnostics.

This WIC-led rule is a proposal, not the behavior of the live pipeline. Before
implementation, define the tolerance, tie handling, whether one SI frame may
be reused by adjacent WIC frames, the timestamp used for Kp, and whether WIC
frames rejected by the present fullness tests remain absent or survive with a
validity flag. Store all three source times so this choice remains auditable.

### Spatial coregistration

WIC is not resampled during detector coregistration. SI12 and SI13 are mapped
independently onto the WIC detector pixels using their projected detector
footprints and overlap areas.

The current detector-coregistration experiment is the candidate starting
point:

1. use each sensor's per-frame geographic footprint geometry;
2. express SI footprint corners in continuous WIC row and column coordinates;
3. intersect each SI footprint with the WIC detector pixels;
4. calculate an overlap-area-weighted SI value in each WIC pixel; and
5. retain mapped coverage and mapping diagnostics.

This is still an approximation because the camera point-spread functions are
unknown. It is nevertheless more defensible than treating a large SI pixel as
a point or interpolating a previously binned SI image onto a finer WIC grid.

The current helper also assumes a uniform top-hat footprint and uses specific
round-trip-error and minimum-coverage thresholds. The footprint model,
thresholds, treatment of partial edge footprints, and reuse of a sensor frame
must be validated and then recorded as Product-1 configuration.

There is a deeper radiometric gate. Area-averaging SI counts assumes that the
stored calibrated count rate behaves as an intensive image quantity over a
detector footprint. Confirm that assumption against the IMAGE response-table
definition, integration time, pixel solid angle, and viewing geometry.

The mapping must not imply that several WIC pixels covered by one SI footprint
are independent SI observations. Retain channel-specific coverage and source
support so later uncertainty and covariance work can account for that shared
information.

### Minimum content

Product 1 should contain, per WIC frame and detector pixel:

- selected WIC counts and a detector-noise uncertainty when one is defined;
- coregistered SI12 and SI13 counts and propagated detector-noise uncertainty
  when one is defined;
- separate WIC, SI12, and SI13 quality weights or validity fields;
- SI12 and SI13 coregistration coverage and source-support diagnostics;
- WIC detector row and column;
- geographic latitude and longitude;
- magnetic latitude and MLT at the agreed reference height;
- SZA, DZA, and any other retained viewing geometry;
- WIC time and matched source indices/times for all three sensors; and
- source-file and camera-processing identifiers.

Kp and proton-corrected counts do not belong in Product 1. They first enter
the precipitation calculation.

The current upstream `PreImage` files do not provide a calibrated
detector-count uncertainty. The existing `BinnedImage.sigma` is within-cell
spatial spread with optional small-sample inflation; it must not be moved into
Product 1 and relabelled as detector noise. Defining or sourcing the
detector-level noise model is an explicit schema and VAE-likelihood gate.

## Product 2: precipitation on the WIC detector geometry

### Scientific role

Product 2 applies a named proton correction and electron-precipitation method
to Product 1. Different methods produce separate products that can coexist.

The two current electron paths are:

- `image_ratio`: requires finite WIC, SI12, and SI13 support;
- `zhang_paxton`: requires finite WIC and SI12 support, while SI13 remains an
  optional diagnostic.

The time axis and detector geometry remain aligned with Product 1. Cells or
frames lacking a required input are invalid for that method rather than being
silently removed from the orbit.

The first experimental implementation now provides the `image_ratio` path.
It preserves every Product-1 WIC frame, applies Hardy or constant proton
energy on the time-dependent detector coordinates, uses the provisional
map-SI12-counts-then-infer-flux order, and stores method validity separately
from the three input-channel validity fields. This is a parity product, not a
candidate publication retrieval. The detector-compatible Zhang--Paxton path
remains the next separate implementation step.

### Provisional calculation order

The visible calculation should remain:

1. match Kp to the WIC frame time;
2. evaluate the named proton-energy model;
3. infer proton flux from SI12;
4. proton-correct WIC and SI13;
5. calculate electron mean energy and energy flux with the named method; and
6. store method-specific diagnostics and validity.

Proton correction is a distinct calculation within Product 2. It should not
be hidden inside either electron-energy method.

Steps 2--4 contain a noncommutation question that must be tested. The current
idea maps SI12 counts to WIC and then evaluates the energy-dependent proton
response at each WIC pixel. An alternative is to infer proton flux on the
native SI12 geometry, map that flux to WIC, and only then calculate the WIC
and SI13 proton contributions. These are equivalent only under restrictive
response and within-footprint-energy assumptions. Product 1 should retain the
native-source link and coregistration information until the order is settled.

### Zhang--Paxton consequence

The current Zhang--Paxton lookup is tied to one 36-by-36 Cubed-Sphere grid.
That is incompatible with a detector-first Product 2 because each WIC frame
has different two-dimensional MLT coordinates.

The replacement should be spatial-grid independent: a collapsed lookup or
callable in `(Kp, MLT)` that can be sampled at arbitrary WIC detector pixels.
The full Zhang--Paxton model remains in `ZhangPaxton2008`; the shared collapsed
calculation and lookup interface belong in `icPhysics`. The final Kp spacing,
periodic MLT spacing, interpolation rule, and table ownership require a
separate design decision before implementation.

The coordinate convention is also a scientific gate. The live pipeline uses
Modified Apex coordinates at 130 km, while the Hardy model was published in
corrected geomagnetic coordinates and Zhang--Paxton has its own stated MLAT
and MLT convention. Record more than a reference height: test and document the
coordinate-system approximation used by both models.

### Minimum content

Product 2 should contain:

- the Product-1 geometry and time information needed to interpret it alone;
- Kp and its matched interval;
- the named proton-energy and proton-flux methods;
- raw and response-limited proton energy, proton flux, and clipping flags;
- proton-corrected WIC and SI13 with their uncertainties;
- electron mean energy and energy flux with their uncertainty fields;
- image ratio and related diagnostics when calculated;
- separate input-channel support plus method-valid support;
- the exact precipitation-method configuration; and
- a reference to the Product-1 source.

The current analytic uncertainty fields may be retained for continuity, but
their known limitations must remain explicit. The planned Monte-Carlo path in
icAnalyzer is not replaced by these fields.

## Product 3: conductance on the WIC detector geometry

### Scientific role

Product 3 applies a named conductance forward model to
`precipitation_detector`. Robinson is the first model, not an architectural
assumption.

Each forward model produces a separately labelled Product 3. Changing the
forward model must not require rerunning camera preprocessing, detector
coregistration, proton correction, or electron precipitation.

### Minimum content

Product 3 should contain:

- the detector geometry and time needed to interpret the file alone;
- Hall and Pedersen conductance;
- retained conductance-uncertainty fields and validity;
- the E0 and Fe state used by the forward model, or an unambiguous reference
  to it;
- precipitation, proton, and conductance method identifiers; and
- a reference to the Product-2 source.

The canonical conductance values are calculated before any Cubed-Sphere
binning.

## Fixed Cubed-Sphere representations

### Why bin late

Late binning separates two decisions that are currently entangled:

- how the instrument observations are converted into physical quantities;
- how those quantities are represented for a particular analysis.

It also avoids pretending that SI13 gained WIC resolution merely because an
SI Cubed-Sphere field was interpolated to a WIC Cubed-Sphere field.

### Mapping rule

Project WIC detector footprints onto the fixed Cubed-Sphere cells. Build the
geometry mapping once per WIC frame and reuse its overlap geometry, while
applying a reducer appropriate to each field. Continuous values, variances,
validity flags, clipping flags, coverage, and source counts must not all be
treated as ordinary spatial means.

Each CS product should retain at least:

- the area-weighted field value;
- valid covered area or fractional coverage;
- the number of contributing WIC detector footprints;
- a support or effective-sample diagnostic where shared SI footprints matter;
- the fixed grid identifier and coordinates; and
- the detector product and binning-method identifiers.

The exact treatment of uncertainty is not yet settled. A simple weighted
scatter is not automatically a measurement uncertainty, especially when one
SI footprint contributes to several WIC pixels. Coverage and source counts do
not reconstruct the induced correlation. Decide whether Product 1 stores a
sparse overlap operator, source-footprint identifiers and fractions, or a
smaller approximation. icBuilder owns that mapping geometry; icAnalyzer may
consume it when correlated observation error matters.

Constant-field reproduction and area-integral conservation are different
tests. A renormalized overlap mean can pass the former while hiding missing
covered area in the latter. Validate and report both.

### Which CS products to store

`precipitation_cs` and `conductance_cs` are expected outputs. Store both
because they answer different scientific questions and storage is available.

`fuv_cs` is optional. It is not needed to create detector-level precipitation
or conductance, but it becomes useful if the planned joint WIC/SI12 VAE is
trained directly on FUV observations. Decide this from the first VAE input
design rather than creating it only for symmetry.

That decision must precede the new schemas and corpus generation if the joint
FUV VAE is part of the first publication. The existing annotations are tied to
an older 36-by-36 grid; their dividing lines should move to a grid-independent
coordinate form before being rasterized onto the new fixed grid.

`conductance_cs` is the normal input to the present Hall/Pedersen VAE,
covariance work, sparse reconstruction, and spline workflow.

## Product identity and file layout

The exact directories are not yet fixed. A readable layout would make both
scientific method and representation visible, for example:

```text
fuv_detector/<preprocessing_label>/or_XXXX.nc

precipitation_detector/<retrieval_label>/or_XXXX.nc
precipitation_cs/<grid_id>/<retrieval_label>/or_XXXX.nc

conductance_detector/<retrieval_label>/<forward_model>/or_XXXX.nc
conductance_cs/<grid_id>/<retrieval_label>/<forward_model>/or_XXXX.nc
```

Labels such as `IR_hardy` and `ZP_hardy` remain useful for human-selected
sensitivity datasets. The file itself must still record the full method and
parameter values; the folder name is not provenance.

Every NetCDF product should have a small common identity block:

- `product_type`;
- `representation = detector` or `cs`;
- `schema_version`;
- source product path or identifier;
- relevant method names and parameters;
- coordinate system and reference height; and
- software version or commit when the publication corpus is generated.

This does not require a generic schema framework or inheritance hierarchy.
Each scientific product should keep its own readable writer and reader.

These are new product identities, not in-place schema updates. They supersede
the current `binned_fuv` schema 1 and `precipitation`/`conductance` schema 2
for the future corpus. Old and new files must not share a `product_type` that
makes them appear interchangeable.

## Code and repository ownership

### icBuilder

Owns the processing workflow:

- sensor-file reading and preprocessing selection;
- WIC-led time matching;
- detector geometry and SI-to-WIC coregistration;
- serialization of whatever overlap/source identity is selected for
  coregistration uncertainty;
- Kp assignment to IMAGE frames;
- detector-to-Cubed-Sphere mapping;
- per-orbit orchestration and restart; and
- NetCDF serialization.

### icPhysics

Owns array-in/array-out physical calculations shared by icBuilder and
icAnalyzer:

- camera response and proton correction;
- Hardy proton energy;
- image-ratio precipitation;
- the grid-independent collapsed Zhang--Paxton interface; and
- Robinson and later conductance forward models.

It should not read orbit files, download Kp, define the IMAGE product layout,
or depend on fuvpy.

### icReader

Owns loading and dispatch for each `product_type` and `representation`. It
should expose the stored arrays without silently regridding or recalculating
physics.

### icAnalyzer

Consumes regular CS products for VAE training, covariance estimation, sparse
reconstruction, and assimilation experiments. It may also call icPhysics
directly for Monte-Carlo samples. A generated VAE realization is an analysis
object rather than a reprocessed detector observation, so applying icPhysics
directly on its fixed grid does not redefine the canonical observational
pipeline.

## Proposed implementation sequence

Implementation should proceed through one representative orbit before a full
corpus rerun.

1. **Freeze the source boundary.** Confirm which native WIC, SI12, and SI13
   fields enter Product 1 and record their background/calibration provenance.
2. **Implement `fuv_detector`.** Reuse the established detector-space
   footprint coregistration and preserve unmatched channels rather than
   dropping WIC frames.
3. **Move Product 2 to detector geometry.** Keep proton correction separate.
   The ratio path may be implemented first as a parity and regression test
   because it needs no new Zhang--Paxton lookup; current evidence does not
   support it as the default publication retrieval.
4. **Replace the grid-bound Zhang--Paxton lookup.** Verify an arbitrary-MLT
   lookup against direct collapse values and resolve its energy and coordinate
   definitions before enabling the candidate publication path.
5. **Move Product 3 to detector geometry.** Apply the existing icPhysics
   conductance function without spatial binning.
6. **Implement one reusable WIC-detector-to-CS mapping.** Apply it separately
   to Product 2 and Product 3 and verify area and constant-field behavior.
7. **Add icReader support.** Load all detector and CS representations without
   compatibility fallbacks for the experimental branch.
8. **Validate one orbit end to end.** Compare current and new products stage
   by stage before selecting the publication configuration or rerunning the
   corpus.

## Acceptance checks

Before the architecture is used for a publication dataset, verify:

- sensor time matches and unmatched-channel behavior on real orbits;
- SI-to-WIC constant-field and covered-area conservation;
- the radiometric meaning of area-averaging calibrated detector counts;
- map-counts-then-correct versus correct/infer-on-SI-then-map proton results;
- detector coregistration against the existing Coumans/Frey diagnostic path;
- identity when detector and target geometry are deliberately made equal;
- fixed-grid cell areas, coordinates, and grid identity across all orbits;
- detector-to-CS conservation for constant and synthetic structured fields;
- Product-3 CS values equal binned detector conductance, not conductance
  recalculated from binned precipitation;
- source links, methods, units, masks, and clipping flags survive NetCDF and
  icReader round trips; and
- one-orbit restart does not mix products made with different configurations.

## Open decisions

The architecture deliberately does not settle:

- the fixed Cubed-Sphere extent, resolution, and permanent `grid_id`;
- whether `fuv_cs` is required for the first VAE publication;
- the grid-independent annotation format and conversion of old annotations;
- the final directory names and schema versions;
- the WIC-led matching tolerance, tie/reuse behavior, and frame-validity rule;
- the selected Product-1 background, reflattening, DZA, and hemisphere rules;
- the detector-count noise model needed by the VAE likelihood;
- the radiometric interpretation and operation order of SI12 coregistration
  and proton inference;
- the detector-footprint, edge, round-trip-error, and coverage rules;
- the periodic MLT resolution and interpolation rule of the new
  Zhang--Paxton lookup;
- the Hardy and Zhang--Paxton magnetic-coordinate convention;
- the exact uncertainty propagation through detector coregistration and CS
  binning;
- whether the SI-to-WIC overlap operator or a smaller correlation descriptor
  is stored;
- the most useful measure of independent/effective SI support after
  coregistration;
- whether a publication product should duplicate all geometry or carry a
  smaller self-contained geometry block plus a source reference; and
- which background, electron-energy, and conductance configurations become
  the published defaults.

These choices should be resolved by focused numerical tests. They are not
reasons to return to early binning.
