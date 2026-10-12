# Electron-energy model validation

This folder revives the DMSP validation of the Zhang--Paxton electron-energy
model using the current detector-era, per-orbit crossing files.

The workflow deliberately has two stages.

1. `build_energy_comparison.py` streams the crossing files and writes one
   compact record per IMAGE frame, DMSP satellite, and temporal support. DMSP
   mean energy is calculated as `sum(Q) / sum(Q / E)`. The output contains
   definitive Kp, the original Zhang--Paxton model evaluated along every
   one-second DMSP position, and the Zhang--Paxton area-mean, area-median,
   flux-weighted, and latitude-spread values. The pointwise prediction is
   reduced over each window with the same `sum(Q) / sum(Q / E)` operation as
   DMSP. Annotation status is attached but does not filter the data.
2. `evaluate_energy_models.py` compares constant and Zhang--Paxton predictors
   using held-out IMAGE orbits. It reports the all-crossing and manually
   accepted cohorts separately. The score is absolute log10 energy error,
   averaged first within orbit and then over orbits.

New comparison files also contain the uncalibrated full spatial Hardy electron
model, evaluated at each one-second DMSP MLAT/MLT and the matched IMAGE-frame
Kp before reduction. `hardy_pointwise_energy` is the ratio of summed model
energy flux to summed model number flux on the valid DMSP window. Its native
flux has units keV cm^-2 s^-1 sr^-1; it does not replace measured DMSP flux.
The evaluator and plots include Hardy automatically when this field is present;
older comparison files and results retain the previous model set. Rebuild the
comparison file to add Hardy, and use separate output directories for a new
evaluation. All compared models use the same finite positive common sample.
Hardy is not calibrated to DMSP; report its valid sample count alongside scores.
When sample counts are present, the common evaluation requires every valid
DMSP second in a window to have Hardy support. This avoids comparing a
partially covered model window against the complete DMSP window.

The model comparison includes the historical hard-coded 2-keV electron
energy, two constants fitted only on training orbits, the original pointwise
Zhang--Paxton model, three latitude-collapse definitions, and two diagnostic
DMSP calibrations of the area-mean collapse.

Build the reusable comparison data:

```bash
python scripts/electron_energy_model/build_energy_comparison.py \
    --matches /path/to/detector_crossings \
    --annotations data/dmsp_frame_annotations.csv \
    --output data/electron_energy_model/comparison_data.nc
```

Run the default comparison with no absolute DMSP flux floor:

```bash
python scripts/electron_energy_model/evaluate_energy_models.py \
    --input data/electron_energy_model/comparison_data.nc \
    --output data/electron_energy_model/results
```

The flux threshold is an evaluation choice, not a build choice. For the
strong-signal sensitivity discussed during ratio validation, rerun with:

```bash
python scripts/electron_energy_model/evaluate_energy_models.py \
    --minimum-flux 2e11 \
    --output data/electron_energy_model/results_flux_2e11
```

Plot either evaluation independently:

```bash
python scripts/electron_energy_model/plot_energy_model_diagnostics.py \
    --results data/electron_energy_model/results

python scripts/electron_energy_model/plot_energy_model_diagnostics.py \
    --results data/electron_energy_model/results_flux_2e11
```

Each run writes sample histograms, raw Zhang--Paxton/DMSP density plots,
held-out residual histograms, and equal-orbit model scores as PNG and PDF.
The default figures use 120-s support; pass `--support-seconds 66` for the
shorter window.

Annotation values retain their original meaning: `-1` is absent from the
annotation table, `0` is rejected, `1` is accepted, and `2` is deferred
because the crossing samples only the dayside. The `accepted` cohort uses only
status `1`; the `all` cohort does not inspect annotation status.
