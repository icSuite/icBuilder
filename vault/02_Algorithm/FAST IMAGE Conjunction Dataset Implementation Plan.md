# FAST--IMAGE Conjunction Dataset Implementation Plan

Last reviewed: 2026-11-05

Status: Implemented in `scripts/fast/06_generate_conjunction_dataset.py` with
focused synthetic tests and a verified one-orbit server smoke test. The first
year-2000 dataset build remains pending.

## Objective

Build one analysis-ready NetCDF dataset containing all objectively matched
FAST--IMAGE conjunction samples. The first analysis must use the complete,
unannotated conjunction population. Manual annotation is a possible later
selection layer, not part of conjunction discovery or this script.

The calculation must run on the server because the canonical IMAGE detector
Products 2 and 3 are not stored locally. The processed FAST products will be
transferred to the same server before the run.

## Fixed workflow

Follow the established DMSP--IMAGE crossing architecture:

1. use the compact IMAGE precipitation CS product as a temporal and spatial
   support gate;
2. open detector Products 2 and 3 only for CS-supported candidates;
3. match the 130-km FAST footprint to an exact WIC detector footprint;
4. extract the FAST and detector-level IMAGE values for every successful
   match;
5. resolve multiple IMAGE-frame matches so each FAST record appears at most
   once; and
6. write one consolidated, provenance-rich dataset.

Do not replace the CS gate with a direct scan of all detector pixels. Do not
write geometrical coincidences that lack matched IMAGE data.

## Inputs

The script requires four product collections:

- processed FAST NetCDF files from `04_process_fast_archive.py`;
- IMAGE precipitation CS products on the fixed 130-km grid;
- the detector precipitation products from which those CS products were
  derived; and
- the corresponding detector conductance products.

All IMAGE products must be read through the public `icReader` interface. Raw
NetCDF access remains appropriate for the icBuilder-owned FAST products and
for writing the conjunction dataset.

The initial command-line boundary is:

```text
--fast-path
--precipitation-cs-path
--precipitation-detector-path
--conductance-detector-path
--output
--year
--orbit                 optional restricted diagnostic run
--fast-file             optional restricted diagnostic run
--frame-support-seconds default 60
--overwrite
```

Use explicit paths on the server. Do not embed the local `/media/bing/LaCie`
layout in the conjunction calculation.

## Input collection checks

Before matching, index all four collections and fail clearly when their
contracts are incompatible.

For each IMAGE orbit used:

- the CS product must be `precipitation_cs` on
  `image_apex_130km_46x46_v1`;
- the precipitation detector product must be `precipitation_detector`;
- the conductance detector product must be `conductance_detector`;
- the CS product must identify the selected precipitation detector product as
  its source;
- Product 3 must identify the selected Product 2 as its source;
- Products 2 and 3 must have identical frame times and detector geometry;
- the reference height must be 130 km and agree with the FAST footprint
  height; and
- the collection must not mix preprocessing, precipitation, count-source,
  smoothing, proton-energy, or conductance configurations.

Use the hashes already recorded by the IMAGE products to verify source
relationships. A filename match alone is insufficient.

For every FAST file used:

- processing must be complete;
- `footprint_height_km` must be 130;
- timestamps must be strictly increasing and unique; and
- the spectral channel dimension and required field contract must be
  supported.

Files that fail structural or configuration checks are errors, not silently
skipped crossings.

## Matching unit

One output `sample` is one FAST measurement uniquely matched to one IMAGE
frame and one WIC detector pixel.

The primary dataset is pointwise. Contiguous samples can be grouped later for
event-level analysis. Do not introduce a subjective event or annotation
definition into the first builder.

Retain stable source indices:

- FAST source filename and record index;
- IMAGE orbit and frame index;
- WIC source index; and
- detector row and column.

These fields define a reproducible sample identity without depending on row
order in the consolidated file.

## Candidate matching

### 1. Temporal candidates

For each IMAGE frame, select FAST records within the closed interval

```text
[image_time - frame_support, image_time + frame_support]
```

The default `frame_support` is 60 seconds, matching the DMSP crossing
collector. Preserve the signed FAST-minus-IMAGE offset in seconds.

A FAST record can initially be a candidate for adjacent IMAGE frames. Do not
discard those alternatives before the CS and detector support tests; the
nearest frame may lack usable IMAGE support while another frame inside the
allowed interval has a valid match.

### 2. CS support gate

Apply the same support rule as `fetch_all_dmsp_crossings.py` to each temporal
candidate:

1. require finite FAST footprint MLAT and MLT;
2. require the footprint to lie inside the CS grid;
3. find the corresponding CS row and column;
4. require finite, positive WIC coverage;
5. require finite, positive SI13 coverage; and
6. require finite CS `wic_corrected` and `si13_corrected` values.

Only candidates satisfying all six conditions proceed to the detector
products. Save the CS row, column, WIC coverage, and SI13 coverage for
successful samples.

This is intentionally the image-ratio Product-2 support contract used by the
current DMSP workflow. If a future IMAGE retrieval has a different required
sensor set, it requires an explicit new gate tied to that retrieval; do not
quietly weaken or generalize this rule.

### 3. Exact detector-footprint match

For each surviving frame:

1. read detector MLAT and MLT;
2. infer detector quadrilateral footprints with the existing
   `icbuilder.footprints` utilities;
3. use the CS cell and its neighbours to restrict candidate detector
   footprints;
4. test point-in-quadrilateral containment in the CS projection;
5. reject the candidate if no detector footprint contains the FAST point;
6. if several footprints contain it, choose the nearest detector-footprint
   centre; and
7. retain the number of containing footprints and the angular separation from
   the selected pixel centre.

Reuse the established DMSP helpers or extract shared helpers if doing so
removes real duplication without obscuring the calculation. The numerical
footprint semantics must not diverge between DMSP and FAST.

### 4. Exact detector data requirement

After containment, require finite matched detector WIC and SI13 corrected
values at the selected Product-2 pixel. A CS-supported point whose exact
detector pixel lacks that pair is counted as a detector-data miss and is not
written.

Do not additionally require `method_valid`, `conductance_valid`, or valid
uncertainty to define a conjunction. Those flags describe results at a pixel
with matched IMAGE data and must remain available for the first analysis.

### 5. Unique IMAGE match per FAST record

After all support tests, retain at most one match for each FAST source record:

1. smallest absolute time offset;
2. earlier IMAGE time for an exact tie;
3. lower IMAGE orbit;
4. lower frame index; and
5. lower detector row, then column.

Record the number of successful candidate frames before this selection. This
distinguishes an unambiguous conjunction from a sample observed by overlapping
IMAGE frame-support intervals.

## Output contract

Write one atomic NetCDF file with dimensions:

```text
sample
energy_index
```

Sort samples by FAST time, then FAST source filename and record index. The
output must not depend on filesystem enumeration order or worker completion
order.

The output descriptor is provisionally:

```text
product_type = fast_image_conjunctions
schema_version = 1
representation = detector_point_samples
```

The final file is published through a `.partial` file, reopened, validated,
and atomically renamed. An existing current output is skipped only when its
input collection signatures and matching configuration agree. `--overwrite`
permits deliberate replacement.

### Match identity and geometry

Store at least:

- FAST source filename and source-record index;
- FAST time and IMAGE time;
- signed time offset and successful-frame-candidate count;
- IMAGE orbit, IMAGE frame, and WIC source index;
- CS row and column plus WIC/SI13 CS coverage;
- detector row and column;
- detector-footprint match count and centre separation;
- FAST spacecraft geographic position and altitude;
- FAST 130-km footprint geographic latitude/longitude, MLAT, and MLT; and
- selected IMAGE detector geographic latitude/longitude, MLAT, MLT, SZA, and
  DZA.

### FAST scientific fields

Copy the processed FAST quantities rather than recalculating them:

- spacecraft potential, FAST validity, and FAST data-quality flags;
- loss-cone angle, footprint-to-spacecraft magnetic-field ratio, and
  footprint magnetic field;
- conservative, centre, and inclusive Earthward number flux, energy flux,
  and mean energy;
- full-spectrum Earthward moments;
- anti-Earthward energy flux, net diagnostics, directional flux ratio, and
  contributing angle-bin count;
- spectrum-derived Hall and Pedersen conductance;
- Robinson-effective mean energy, flux, and validity;
- full time-dependent energy centres and channel widths; and
- centre-rule Earthward and anti-Earthward pitch-angle-integrated spectra.

Do not copy redundant conservative/inclusive spectra in schema 1. The
centre-rule spectra reproduce the baseline conductance, while the saved
conservative/inclusive scalar moments preserve the established structural
sensitivity.

### IMAGE Product-2 fields

At the selected detector pixel store:

- WIC/SI12/SI13 validity and coverage fields;
- `method_valid` and method quality weight;
- `Ep_model`, `Ep`, `dEp`, `Fp`, and `dFp`;
- corrected WIC, SI12 when available, and SI13 values and uncertainties;
- `R`, `dR`, `E0`, `dE0`, `Fe`, `dFe`, and `varE0Fe`;
- response-range and proton-clipping flags when available; and
- frame Kp and Kp interval start.

Prefix IMAGE variables in the conjunction product, for example `image_E0`
and `image_P`, so FAST and IMAGE quantities cannot be confused.

### IMAGE Product-3 fields

At the same detector pixel store:

- Hall and Pedersen conductance and their uncertainties;
- `conductance_valid`;
- `conductance_uncertainty_valid`; and
- optional Fang domain/clipping diagnostics when present.

Do not duplicate Product-2 precipitation variables merely because Product 3
also carries them. Validate their alignment/source relationship and retain one
canonical copy from Product 2.

### Units and flags

Copy documented units and flag meanings from the source variables where
possible. Add explicit units to all matching coordinates and separations.
Preserve integer flags as integer variables; do not turn unavailable flags
into false values.

## Counts and exclusions

The root metadata must report at least:

- FAST input file and record counts;
- temporal candidate count;
- nonfinite-footprint count;
- outside-CS-grid count;
- CS-coverage/data rejection count;
- detector-footprint miss count;
- exact detector-data miss count;
- ambiguous detector-footprint count;
- successful candidate count before unique selection;
- duplicate-frame candidates removed; and
- final sample count.

These audit counts describe why input records did not become conjunction
rows. Rejected candidates themselves are not stored.

## Provenance

Record:

- all four input roots;
- sorted file-collection signatures based on paths, sizes, and modification
  times;
- the signature basis and selected year/orbit/file restrictions;
- the Product-2, CS, and Product-3 source hashes and common configuration;
- FAST processing configuration and space-weather-index identity;
- temporal support and unique-match tie-breaking rule;
- CS gate definition;
- detector-footprint algorithm;
- reference height and coordinate system; and
- icBuilder, icReader, icPhysics, NumPy, xarray, netCDF4, and secsy versions.

Full hashes of every large server product are not required when the products'
recorded source hashes and collection signatures are retained.

No annotation file or annotation status belongs in schema 1. A later
annotation table must join through the stable source/frame identity without
rewriting the all-conjunction dataset.

## Memory and server execution

The detector products and FAST spectra make a full in-memory yearly join
unsafe. Process IMAGE orbits serially in the first implementation:

1. load only temporally overlapping FAST records;
2. apply the compact CS gate;
3. materialize each required detector field once for the selected orbit;
4. gather matched pixels;
5. append compressed chunks to the `.partial` output; and
6. release the orbit products before continuing.

Use chunking along `sample`; spectral fields use `(sample, energy_index)`
chunks. Start with one process. Add orbit-level parallel extraction only after
profiling shows that server storage and memory can support it and without
allowing multiple workers to write the same NetCDF file.

Because FAST preprocessing is still growing, the safe initial update policy
is to rebuild the consolidated file when the FAST source collection changes.
Do not implement in-place append or checkpoint/resume machinery in the first
version.

The implementation will reasonably exceed 200 lines because it combines four
validated product contracts, two-stage spatial matching, unique-match
resolution, chunked serialization, and provenance. Keep it to the named
script plus one focused test file unless a genuinely shared DMSP/FAST
footprint helper is extracted deliberately.

## Verification

Focused tests must cover:

1. temporal-window edges and deterministic ties;
2. rejection outside the CS grid;
3. rejection for zero/NaN WIC or SI13 CS coverage/data;
4. exact detector containment, miss, and ambiguous-footprint resolution;
5. rejection of an exact detector pixel without finite WIC/SI13 data;
6. retention of matched pixels with false method/conductance/uncertainty
   flags;
7. unique selection when one FAST record matches adjacent IMAGE frames;
8. Product-2/CS/Product-3 source mismatch failures;
9. mixed retrieval/configuration failure;
10. FAST and IMAGE field extraction including the spectral dimension;
11. deterministic ordering and collection signatures;
12. atomic output and restart validation; and
13. empty-result behavior with a valid zero-sample dataset or an explicit,
    documented failure chosen before implementation.

Then run one small real server case, initially one IMAGE orbit, and verify:

- every output sample satisfies the CS gate;
- every output point lies in its saved detector footprint;
- direct source reads reproduce selected FAST, Product-2, and Product-3
  values;
- each FAST source record occurs at most once;
- counts reconcile from input through final output; and
- peak memory remains compatible with serial server execution.

Only after this smoke test should the script process all currently available
year-2000 FAST products.

### Verified server smoke test

Orbit 0085 was run against 22 processed FAST files from 2000-05-15/16, the
local precipitation CS product, and the mounted server detector Products 2
and 3. The scratch output is
`/tmp/icbuilder_fast_image_smoke/or_0085.nc` and contains 3,998 unique samples
in a 3.6-MB file.

The audit chain reconciles exactly: 4,434 temporal candidates comprise 116
outside the CS grid, 181 rejected by CS coverage/data, 69 detector-footprint
misses, 57 exact detector-data misses and 4,011 supported candidates. Thirteen
adjacent-frame duplicates were removed, leaving 3,998 unique FAST source
records. Seven supported candidates had ambiguous detector-footprint
containment and used the recorded nearest-centre rule. No candidate had a
nonfinite FAST footprint.

An independent exhaustive reread verified every saved FAST scalar and
spectrum, every Product-2 and Product-3 field, both saved CS coverage values,
all source indices and all 3,998 point-in-detector-footprint relations. The
configuration-aware restart check then skipped the unchanged output. The
smoke test satisfies the implementation gate; the next action is the
available year-2000 build.

## Explicit non-goals

Schema 1 does not:

- annotate or accept/reject conjunctions;
- fit, calibrate, or score either instrument;
- average FAST records in time;
- aggregate point samples into event summaries;
- recalculate FAST or IMAGE precipitation/conductance;
- include unsupported geometrical coincidences;
- depend on local availability of detector products; or
- process the full 2000--2003 corpus before the year-2000 workflow is checked.
