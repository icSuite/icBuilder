# icReader Read Migration Implementation Plan

Last reviewed: 2026-09-20
Status: Implemented and verified locally; exact icReader commit pin pending

## Objective

Make icReader the single high-level read boundary for NetCDF products written
by icBuilder. Remove duplicate product decoding and schema validation from the
active icBuilder pipeline where the committed icReader interface already
covers the product, without moving raw-source ingestion, serialization, or
scientific calculations out of icBuilder.

The migration must preserve numerical output, provenance, atomic restart
semantics, and the optimized read-once behavior of detector-to-CS reduction.

## Implementation result

The active detector-first path now reads generated products through public,
context-managed icReader objects. This includes Product-1 to Product-2,
Product-2 to Product-3, paired detector-to-CS reduction, detector/CS restart
validation, and the maintained diagnostics listed below. Raw NetCDF remains at
the intended boundaries: fuvpy/native or external inputs, writers, lookup
data, and tests that deliberately construct or corrupt serialized files.

icReader now exposes an immutable `variable_attrs` mapping for every declared
variable, with the same mapping available as `ProductField.attrs`. This lets
restart validation retain units and flag-meaning checks without reaching into
the reader's private NetCDF handle. Detector-to-CS construction opens each
source once and has a regression test proving that every compressed cube is
read no more than once.

The migrated Product-2 and Product-3 orbit-0085 files were exactly equal to
the existing products for every variable. Both CS products were also exactly
equal apart from expected runtime provenance paths/version text. The complete
0085, 0086, 0261, and 0968 CS gate took 37.77 seconds and 1,614,272-KB peak
RSS with two workers, compared with 39.60 seconds and 1,614,780 KB before the
migration. All eight generated files were variable-for-variable and
grid-for-grid identical. Restart recognized all four orbits in 0.42 seconds.

icReader passes all 34 tests. icBuilder passes 127 tests with only the four
known, unrelated 36-by-36 grid/Zhang--Paxton lookup failures. Repository-wide
AST parsing and targeted whitespace checks pass.

One packaging closeout remains intentionally unperformed: the new icReader
metadata API is still an uncommitted worktree change, so icBuilder cannot yet
name its exact commit. `pyproject.toml` temporarily follows icReader's
`modular_pipeline` branch. Commit and push icReader first, replace that branch
reference with the resulting commit hash, then commit icBuilder. No scientific
or product-schema change is waiting on this bookkeeping step.

## Verified starting point

- icBuilder is on `modular_pipeline` at `23ea850` with a heavily dirty
  research worktree. Existing changes and generated products must be
  preserved.
- icReader is on `modular_pipeline` at pushed commit `2324d6c`. Its dispatcher
  supports `binned_fuv`, `precipitation`, `conductance`, `fuv_detector`,
  `precipitation_detector`, `conductance_detector`, `precipitation_cs`, and
  `conductance_cs`.
- The five detector-first readers validate exact product type,
  representation, and schema. Their three-dimensional fields are lazy and
  sliceable, and their product objects are context-managed.
- The two CS readers reconstruct and verify the stored `secsy.CSgrid`.
- The old binned, precipitation, and conductance readers are eager. They
  already serve the active legacy calculation boundary in
  `icbuilder.precipitationimage` and `icbuilder.conductanceimage`, but eager
  loading makes them unsuitable for cheap restart checks over a large corpus.

## Ownership boundary

### icReader should own

- dispatch from the serialized product descriptor;
- product schema, dimensions, required variables, type-aware fill handling,
  and CF-time decoding;
- generated-product metadata access;
- lazy reads of detector and CS data fields;
- reconstruction and validation of stored CS grids.

### icBuilder should retain

- reading fuvpy sensor-orbit files, which are pipeline inputs rather than
  icBuilder products;
- reading HDF orbit indexes, DMSP products, calibration files, and the
  Zhang--Paxton lookup;
- all NetCDF writing and atomic publication;
- source-file hashing, size, and modification-time provenance;
- workflow-specific compatibility checks, such as preprocessing label,
  selected precipitation method, proton model, and expected source identity;
- all scientific calculations and spatial reductions;
- raw NetCDF fixture construction, corruption, and serialization-contract
  assertions in tests.

No icBuilder code should reach through an icReader object to its private
NetCDF handle. If a required operation is not public in icReader, either keep
that bounded raw read or extend icReader deliberately first.

## Migration inventory

| icBuilder location | Current read | Target | Decision |
| --- | --- | --- | --- |
| `icbuilder/precipitationdetector.py` | raw Product-1 `Dataset` plus local mask/time helpers | `icreader.open_product()` returning `FUVDetector` | Replace now |
| `icbuilder/conductancedetector.py` | raw Product-2 `Dataset` plus duplicate schema/time helpers | `PrecipitationDetector` reader | Replace now |
| `icbuilder/precipitationcs.py` | raw detector Product 2 | lazy `PrecipitationDetector` reader | Replace now |
| `icbuilder/conductancecs.py` | raw detector Product 3 | lazy `ConductanceDetector` reader | Replace now |
| `icbuilder/detectorcs.py` | paired raw Product-2/Product-3 reads | paired context-managed icReader products | Replace now while preserving read-once reduction |
| detector and CS pipeline restart validators | raw reads of generated detector/CS outputs | icReader contract validation plus builder-specific source/configuration checks | Replace now |
| `icbuilder/precipitationimage.py` and `icbuilder/conductanceimage.py` | icReader already used for generated inputs | existing icReader calls | Keep |
| legacy binned/precipitation/conductance restart validators | raw generated-product checks | eager legacy readers | Defer because loading full orbit cubes just to decide restart is wasteful |
| `icbuilder/fuvdetector.py`, `PreImage`, and binned input stage | raw fuvpy sensor files | no corresponding icReader product | Keep raw |
| product `to_nc()` methods and lookup writer/loader | raw NetCDF serialization or non-product data | outside reader responsibility | Keep raw |
| focused product diagnostics | mixed raw NetCDF/xarray reads | icReader for generated products only | Migrate after production paths |
| tests | raw fixture writes, mutations, and direct schema checks | mixed | Keep raw where the test is intentionally below the reader boundary |

## Public API convention

Use one import spelling in active consumers:

```python
from icreader import open_product
```

New detector and CS products must be opened with a context manager:

```python
with open_product(filename) as source:
    E0 = source.read("E0")
    valid = source.conductance_valid[:]
```

Use named product attributes for common metadata and `source.attrs` for
serialized attributes that are intentionally not promoted by icReader. Use
`source.shape`, decoded time attributes, `source.time_encoding`, and
`source.read()` rather than dimensions, variables, `num2date`, masked-array
handling, or the private `source._nc` object.

Calculations that require a complete orbit may materialize the necessary
field once. Frame-local diagnostics should slice the lazy field. Avoid adding
an icBuilder compatibility wrapper that recreates the old NetCDF variable
interface; call the small icReader API directly.

## Implementation phases

### Phase 0: freeze the dependency and baseline

1. Record icReader commit `2324d6c` as the minimum tested interface.
2. Add one small icReader interface needed to preserve existing restart
   validation: a read-only mapping of serialized variable attributes for both
   eager and lazy variables. The Product-1 runner currently verifies units and
   `frame_quality.flag_meanings`; dropping those checks or accessing
   icReader's private NetCDF handle is not acceptable. Add this metadata API
   and its focused tests in icReader before pinning the final dependency
   commit.
3. Add an installable, reproducible icReader requirement to icBuilder. Until
   icReader has a release tag or package-index release, use an exact VCS commit
   rather than a mutable branch. Document the editable sibling-repository
   setup used for development.
4. Correct the README statement that icReader is optional: it is already used
   by modular Product-2/Product-3 construction and becomes required by the
   detector pipeline after this migration.
5. Capture the current focused/full test baseline and preserve the four known
   deferred legacy grid/lookup failures as a separate known issue.

### Phase 1: replace detector Product-1 and Product-2 loaders

1. Rewrite `load_fuv_detector()` around `open_product()` and the
   `FUVDetector` contract.
2. Keep only icBuilder-specific checks for the selected preprocessing label
   and authoritative source-time decoding marker. Let icReader own product,
   representation, schema, shape, required-variable, masked-value, and CF-time
   validation.
3. Materialize the fields required by the precipitation calculation before
   leaving the context. Cast them to the same dtypes used by the current
   loader--notably float64 for fields currently read through `_as_array()`--so
   the reader migration does not change numerical evaluation. Preserve
   source-file identity calculation in icBuilder.
4. Rewrite `load_precipitation_detector()` in the same way using the
   `PrecipitationDetector` reader. Retain checks for image-ratio method,
   measurement uncertainty mode, and the source-time marker.
5. Delete the duplicated `_as_array()` and `_read_time()` helpers after parity
   tests cover their former behavior.

These builders already need the complete input orbit for their calculations,
so materializing the required fields does not introduce a new memory cost.

### Phase 2: migrate the optimized detector-to-CS reducer

1. Open Product 2 and Product 3 once, together, through icReader for the
   complete reduction lifetime.
2. Change `validate_detector_pair()` to use public reader attributes and
   decoded time arrays. Keep source-SHA matching and exact frame-identity,
   MLAT, and MLT comparisons.
3. Change `PrecipitationCS` and `ConductanceCS` constructors and reducers to
   accept icReader products. Replace `source.variables[name]` and the local
   `read_variable()` helper with `source.read(name)` or lazy slices.
4. Preserve the performance optimization: build mappings once, and read each
   compressed three-dimensional source variable exactly once per orbit. Do
   not regress to frame-by-frame decompression merely because lazy slicing is
   available.
5. Preserve the current reduction dtypes explicitly. `read_variable()` now
   converts continuous fields to float64 even though their serialized dtype is
   often float32; direct icReader reads preserve the stored dtype. Cast at the
   calculation boundary so summation order and output do not change.
6. Remove `_require_source()` checks that only duplicate icReader's exact
   manifest. Keep workflow constraints such as method, conductance model, and
   uncertainty mode.

The current CS writer copies encoded time values and their CF metadata. After
the migration, store decoded time arrays plus `source.time_encoding`, then
encode them explicitly with `date2num()` during writing. Preserve missing
source times as the variable fill value. A before/after product comparison
must prove that time values, units, calendars, and all other stored variables
remain identical.

### Phase 3: replace detector-first restart reads

Use icReader in:

- `scripts/pipeline/make_fuv_detector_orbit_files.py`;
- `scripts/pipeline/make_precipitation_detector_orbit_files.py`;
- `scripts/pipeline/make_conductance_detector_orbit_files.py`; and
- `scripts/pipeline/make_detector_cs_orbit_files.py`.

Opening the product should determine structural validity. The runner should
then apply only its configuration and provenance checks to distinguish
`mismatch` from `complete`. Preserve the current meanings of `missing`,
`invalid`, `mismatch`, and `complete`, along with `.partial` validation and
atomic `os.replace()` publication.

Retain serialization checks that are not part of the current reader manifest,
including Product-1 field units and frame-quality flag meanings, through the
public variable-metadata mapping added in Phase 0. Do not silently weaken the
restart gate merely because shape validation moved into icReader.

For CS restart checks, rely on icReader's exact grid reconstruction and
coordinate-hash validation rather than maintaining a second manual grid-group
validator. For source matching, continue comparing the stored path, size,
mtime, schema, and/or SHA fields currently required by each runner.

Do not migrate the three legacy restart validators in this phase. Their
icReader implementations eagerly load all image cubes, which would turn a
cheap corpus scan into unnecessary decompression. Revisit them only after the
legacy readers gain a lightweight validation mode or lazy fields.

### Phase 4: migrate maintained diagnostics selectively

Convert generated-product reads in the maintained diagnostics, including:

- `audit_large_ratio_corpus.py`;
- `compare_hardy_fixed_proton_corpus.py`;
- `plot_ratio_vs_dza.py` for its precipitation input;
- `dza_threshold_sensitivity.py`;
- `compare_dmsp_image_ratio_pass.py` for its IMAGE product; and
- `test_ratio_smoothing.py`.

Keep xarray or raw NetCDF for the DMSP side, fuvpy sensor orbits, external
calibration products, and other non-icBuilder inputs. Historical one-off
scripts need not be rewritten merely to eliminate an import; migrate one only
when it is maintained or rerun.

### Phase 5: clean up and document the boundary

1. Remove NetCDF imports only from modules that no longer use them for either
   reads or writes. Writer modules will continue to import `Dataset`.
2. Update the processing-pipeline and detector-architecture notes so generated
   products enter calculations through icReader.
3. Document the context-managed lifetime and the distinction between complete
   orbit reads and frame slices.
4. Record the tested icBuilder/icReader commit pair and server installation
   command.

## Verification

### Contract and parity tests

- For Product 1 and Product 2 loaders, compare every returned field, decoded
  time, dtype, fill convention, shape, and provenance value against the
  pre-migration loader on synthetic fixtures and real orbit 0085.
- Verify missing SI source times remain `None`, validity fields remain Boolean,
  integer counts/indices preserve integer dtype and fill values, and float32
  detector fields are not silently promoted unless the calculation already
  requested float64.
- Verify unsupported schema, method, preprocessing label, uncertainty mode,
  source-time marker, and source identity still fail at the same boundary with
  useful errors.
- Exercise every restart status with complete, structurally invalid,
  configuration-mismatched, source-mismatched, and missing files.

### Numerical output gate

Generate Products 2 and 3 and both CS products before and after the migration
in separate scratch directories. Compare all root attributes, dimensions,
time encodings, grid variables, data variables, masks, counts, and fill
locations. Numerical arrays should be exactly equal; any unavoidable metadata
difference must be explained before acceptance.

Run this first for orbit 0085, then for the existing four-orbit set 0085,
0086, 0261, and 0968. Never run against tracked `example_data` for this gate.

### Performance and resource gate

- Count or instrument reads in the CS path to prove every compressed source
  field is materialized at most once per orbit.
- Compare the four-orbit CS wall time and peak RSS with the verified baseline
  of 39.60 seconds and 1,614,780 KB with two workers. Treat a material
  regression as a failed migration.
- Confirm restart scans remain fast and do not materialize detector cubes.
- Run a two-worker spawn/forkserver scratch test to ensure context-managed
  handles are created inside workers and never passed between processes.

### Repository gate

- Run the focused detector Product-1/Product-2/Product-3 and CS tests.
- Run the full test suite and distinguish the four pre-existing legacy
  grid/lookup failures from new failures.
- Run the package-import subprocess test and repository-wide AST parse.
- Run `git diff --check` on the touched source, tests, README, and vault files.
- Repeat one small UiB-server pipeline run in a clean environment with the
  declared icReader dependency.

## Acceptance criteria

The migration is complete when:

1. active detector Product 2, Product 3, and detector-to-CS construction read
   generated source products only through public icReader APIs;
2. detector-first restart validation uses icReader without decompressing full
   detector cubes;
3. duplicate mask/fill/time/schema helpers are removed from the migrated
   modules;
4. raw reads remain only where the input is not an icBuilder product, the code
   is writing or corrupting NetCDF, or a documented performance limitation of
   an eager legacy reader makes raw validation preferable;
5. pre/post migration outputs agree exactly on the four-orbit gate;
6. CS read-once performance and multiprocessing behavior do not regress;
7. the required icReader revision and installation path are reproducible; and
8. documentation explains the final reader/writer ownership boundary.

## Explicitly deferred

- changing product schemas or scientific calculations;
- changing Hardy clipping, low-signal validity, or E0--Fe covariance;
- rewriting icReader's eager legacy readers as lazy readers;
- replacing raw fuvpy, HDF, DMSP, lookup, or calibration reads;
- replacing NetCDF serialization in icBuilder;
- a blanket conversion of abandoned historical scripts.

This is an interface-consolidation change. It must not be used to alter the
scientific products or hide a numerical difference behind a reader rewrite.
