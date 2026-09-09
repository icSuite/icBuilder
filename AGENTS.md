# Repository guidance

## Purpose and project tracks

`icBuilder` processes IMAGE WIC, SI12, and SI13 observations into
height-integrated Hall and Pedersen conductance estimates with propagated
uncertainties.

Keep these parts of the repository distinct:

- `icbuilder/` contains the reusable image, binning, conductance, and spline
  classes.
- `scripts/make_orbit_h5_files.py`, `scripts/make_orbit_nc_files.py`, and
  `scripts/make_conductance_orbit_files.py` form the primary orbit-processing
  workflow.
- plotting, resolution, spline, and neural-network scripts are secondary or
  exploratory workflows and are not all described accurately by the README.

Do not infer current project status from old dated notes. Check Git and the
live code, then read the current-state and handoff notes.

## Project memory

Project memory is stored in the repository-root `vault/`.
Before substantive work, read:

1. `vault/01_Project/Current State.md`
2. `vault/05_Handoff/Handoff - Latest.md`
3. `vault/01_Project/Project Brief.md`

Read `vault/02_Algorithm/Processing Pipeline.md` when inspecting the live
workflow. Read `vault/02_Algorithm/Detector-First Product Architecture.md`
when changing or discussing the intended publication pipeline. Treat the
older dated notes and images at the vault root as historical evidence, not
current instructions.

Use this source-of-truth order:

1. live code, Git state, tests, and regenerated results;
2. `Current State.md` and `Handoff - Latest.md`;
3. decision and algorithm notes;
4. dated historical notes.

If memory conflicts with live evidence, follow the repository and update the
live notes when the task changes project understanding.

## Working rules

- Preserve unrelated or pre-existing worktree changes.
- Use repository-relative paths in documentation.
- Do not run processing scripts against tracked `example_data/` unless
  overwriting versioned artifacts is explicitly intended.
- Use a copied or scratch `--base` directory for trial runs, and verify that
  all expected output directories exist before starting.
- Begin with a small representative orbit and serial execution before
  increasing worker counts. Some workflows are memory-intensive and have high
  multiprocessing defaults.
- Do not run the full external IMAGE corpus unless the input location,
  environment, runtime, and output destination are explicitly in scope.
- Do not commit the external multi-terabyte corpus, credentials, caches, or
  regenerated products solely because they are near the project vault.
- Treat tracked example data, figures, and historical parameter claims as
  evidence requiring provenance, not automatically current validation.
- Use `icReader` rather than this package when the task is only to read
  generated conductance products.

## Scientific coding style

- Write for a small research group. The expected reader is a student or
  scientist who should be able to follow the calculation from top to bottom.
- Write for a scientist reading the code interactively in an editor. Make
  workflows visually scannable with `#%%` sections where appropriate, blank
  lines between conceptual stages, and short comments that identify each block.
- Keep the main calculation linear and visible from top to bottom. Prefer
  familiar intermediate variables and explicit operations over nested
  expressions, generic plumbing, or compressed control flow.
- Extract a helper when it represents a distinct calculation or removes
  meaningful repetition. Do not hide a few obvious sequential steps behind an
  abstraction merely to shorten the main function.
- Comments may serve as navigational headings even when the underlying Python
  is straightforward. Treat formatter conventions and line-length targets as
  secondary to human readability, while avoiding ambiguous or unwieldy code.
- Use nearby code as the stylistic baseline, not as a ceiling. Do not copy weak
  patterns blindly: identify scientific, numerical, or code choices that could
  be improved, explain the tradeoff, and propose a clearer or safer alternative.
- Adopt improvements when they materially improve correctness,
  reproducibility, clarity, or demonstrated performance. Do not add complexity
  merely because it is conventional in large production systems.
- Let complexity follow the science, numerical method, or actual reuse
  requirements. Prefer direct functions, NumPy arrays, ordinary loops and
  dictionaries, and keep the main calculation visible in execution order.
- Unless current requirements justify them, avoid dataclasses, manager or
  factory classes, generic schemas, version and compatibility frameworks,
  checkpoint/resume machinery, and speculative extension points.
- Retain scientific rigor: make units, coordinates, assumptions, provenance,
  and uncertainty explicit, and add focused tests or reference comparisons for
  consequential calculations.
- Scale packaging, validation, documentation, and abstractions to the code's
  real reuse. A reusable package may justify more structure, but that structure
  should solve a current, explained need.
- If a nominally small feature grows beyond roughly 200 lines or more than two
  new source files, pause and explain why before continuing.

## Verification

No automated test suite or CI workflow was present at the 2026-07-26 review.
A safe syntax-only check that does not write bytecode is:

```bash
python -c "import ast, pathlib; files=list(pathlib.Path('icbuilder').rglob('*.py'))+list(pathlib.Path('scripts').rglob('*.py')); [ast.parse(p.read_text(), filename=str(p)) for p in files]; print(f'parsed {len(files)} files')"
```

For documentation-only changes, also run:

```bash
git diff --check
```

Functional verification requires external dependencies and can overwrite HDF5,
NetCDF, NumPy, spline-factor, or figure outputs. Before running a pipeline,
confirm the environment, use an isolated data copy, and inspect the selected
script's path and worker defaults.

## Automatic memory checkpoints

Project-memory maintenance is a default responsibility. Do not wait for the
user to request a vault update or announce that a session is ending.

Checkpoint after a verified fix or result, a durable implementation or
scientific decision, a changed blocker or next action, and any milestone that
would otherwise leave important understanding only in the conversation.

At a meaningful checkpoint:

1. create `vault/04_Sessions/YYYY-MM-DD.md` only when historical
   detail is worth preserving;
2. rewrite `Current State.md` when verified project state changed;
3. append only durable choices to the decision log;
4. replace obsolete content in `Handoff - Latest.md`;
5. update algorithm notes only when scientific interpretation changed;
6. refresh the handoff's `Portfolio impact` section, using `Central update
   needed: No` when no portfolio-level information changed.

Do not write raw logs, transient speculation, or unchanged state into the
vault. An explicit read-only or no-file-changes request disables automatic
memory writes for that task. Never append new work to an older dated note.
Do not edit the central second brain directly; communicate portfolio changes
through the latest handoff.
