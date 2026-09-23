# Product 2 Unsubtracted Counts Diagnostic Plan

Last reviewed: 2026-09-23
Status: Implemented and verified on orbit 0085

## Objective

Test whether fuvpy background subtraction prevents the detector IMAGE
WIC/SI13 ratio from reproducing the DMSP electron-energy relation. The test
must change only the selected Product-1 count fields. Time matching, detector
geometry, SI coregistration, Hardy proton energy, SI12 proton correction,
measurement uncertainty, ratio retrieval, and DMSP matching remain identical.

The final comparison uses paired DMSP samples on common detector support and
requires WIC solar zenith angle greater than 105 degrees.

## Product contracts

Product 1 schema 2 retains both count stages:

- `wic_counts`, `si12_counts`, and `si13_counts` contain background-subtracted
  `dgimg`;
- `*_unsubtracted_counts` contain calibrated unsubtracted `img`;
- `*_valid` and `*_unsubtracted_valid` describe their separate detector
  support; and
- variance, geometry, source indices, and coregistration mapping are shared.

Product 2 schema 3 gains one root attribute:

```text
count_source = background_subtracted | unsubtracted
```

Its variable layout remains unchanged. Canonical Product-2 fields such as
`wic_corrected`, `si13_corrected`, `R`, `E0`, and their uncertainties describe
the selected count source. Unsubtracted counts are not duplicated in Product
2 because Product 1 remains their observational source.

## Product-2 calculation

Add `count_source="background_subtracted"` to `PrecipitationDetector`.

For the existing branch:

- read `wic_counts`, `si12_counts`, and `si13_counts`;
- preserve current validity and `dgweight`-based method-quality behavior; and
- verify that output arrays remain unchanged.

For the unsubtracted branch:

- read the three `*_unsubtracted_counts` fields;
- require both ordinary and unsubtracted validity for each channel, giving
  common support with the established Product-1 calculation;
- reuse `wic_variance`, `si12_variance`, and `si13_variance`;
- apply the same Hardy/SI12 proton correction and ratio retrieval; and
- do not use background-fit `dgweight` in `method_quality_weight`. Use one on
  successful method support and NaN elsewhere, with this choice recorded in
  metadata.

The serialized `wic_valid`, `si12_valid`, and `si13_valid` fields describe the
selected Product-2 inputs. Background-fit channel weights may remain as source
diagnostics but are explicitly not method weights for unsubtracted counts.

## Orbit runner

Add:

```text
--count-source background-subtracted
--count-source unsubtracted
```

The first remains the default. Default retrieval directories are:

```text
precipitation_detector/IR_hardy/
precipitation_detector/IR_hardy_unsubtracted/
```

Restart validation compares `count_source` in addition to the existing source
file and proton configuration. Atomic `.partial` publication and orbit-level
parallelism remain unchanged.

## Reader support

icReader must expose the six additive Product-1 fields lazily when present,
without rejecting older schema-2 files that lack them. Product-2 readers expose
`count_source`; older schema-3 files without the attribute are interpreted as
`background_subtracted` because that was the only implemented source.

## DMSP extraction and analysis

`fetch_all_dmsp_crossings.py` adds Product-2 `sza` as `wic_sza`. It otherwise
continues to read canonical corrected-count, ratio, energy, uncertainty, and
validity fields, so the same extractor handles both Product-2 branches by
changing `--image-path` and `--output-path`.

The controlled analysis joins both crossing datasets by:

```text
orbit, image_time, dmsp_sat, dmsp_time
```

It requires validity in both branches and `wic_sza > 105`. Histograms must
report paired sample, crossing, and orbit counts. No correlation- or
RMSE-conditioned selection is used for this background sensitivity test.

## Verification gates

1. Focused unit tests cover both count sources, common support, method weights,
   metadata, output labels, and restart mismatch detection.
2. The default branch reproduces the pre-change Product-2 scientific arrays.
3. One real orbit produces both Product-2 variants with identical time,
   geometry, source indices, Kp, and proton energy.
4. The matching extractor saves the same DMSP keys for both branches and
   preserves WIC SZA.
5. Only after those gates pass should annotated or full-corpus products run.

All five implementation gates passed on 2026-09-23. The orbit-0085
background-subtracted output reproduced every prior scientific array exactly;
the unsubtracted output preserved time, geometry, source indices, Kp, and
proton energy while changing the count-dependent retrieval fields. Both DMSP
extracts contained the same 3,847 exact match keys and identical WIC SZA.
Orbit 0085 is sunlit (SZA 62--74 degrees), so a darker orbit or the full corpus
is required to exercise the final SZA-greater-than-105 comparison.

## Production sequence

After both repositories are installed from their updated revisions, generate
the diagnostic Product 2 without overwriting the canonical branch:

```bash
python scripts/pipeline/make_precipitation_detector_orbit_files.py \
  --base-input ~/IMAGE_FUV/ --base-output ~/IMAGE_FUV/ \
  --count-source unsubtracted --workers 10
```

Run the DMSP extractor once for each Product-2 directory, using different
output directories. Existing background crossing files must be overwritten
once because older files do not contain `wic_sza`:

```bash
python scripts/ratio_relation_validation/fetch_all_dmsp_crossings.py \
  --image-path ~/IMAGE_FUV/precipitation_detector/IR_hardy \
  --output-path ~/IMAGE_FUV/ratio_relation_validation/detector_crossings \
  --workers 4 --overwrite

python scripts/ratio_relation_validation/fetch_all_dmsp_crossings.py \
  --image-path ~/IMAGE_FUV/precipitation_detector/IR_hardy_unsubtracted \
  --output-path ~/IMAGE_FUV/ratio_relation_validation/detector_crossings_unsubtracted \
  --workers 4
```

The final comparison command is:

```bash
python scripts/ratio_relation_validation/ratio_validation_background_comparison.py \
  --background-matches ~/IMAGE_FUV/ratio_relation_validation/detector_crossings \
  --unsubtracted-matches ~/IMAGE_FUV/ratio_relation_validation/detector_crossings_unsubtracted
```

## Implementation files

icBuilder:

- `icbuilder/precipitationdetector.py`
- `scripts/pipeline/make_precipitation_detector_orbit_files.py`
- `scripts/ratio_relation_validation/fetch_all_dmsp_crossings.py`
- `scripts/ratio_relation_validation/ratio_validation.py` or a paired
  background-sensitivity analysis script
- focused tests under `tests/`

icReader:

- `icreader/product.py`
- `icreader/detectorproduct.py`
- `tests/test_detector_first_products.py`
