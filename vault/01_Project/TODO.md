# Project TODO

Last reviewed: 2026-10-09

## Dynamic IMAGE response and ratio pipeline

Confirmed by the user2026-10-09. Planned work; not implemented.
Dependency: the user hopes Benoit Hubert will provide sensor-dependent,
DZA-dependent response curves. Receipt and delivery commitment are unconfirmed.
Scientific background: [[Audit - Published IMAGE Ratio Limitations]].

- [ ] Obtain and inspect the sensor response curves from Benoit Hubert.
  Identify sensors, units, energy/wavelength axes, DZA convention, exposure and
  calibration normalization, and reference atmosphere. Determine whether these
  are spectral sensitivities or precipitation-energy-to-count responses at a
  fixed atmosphere; the latter alone may not support varying the atmosphere.
- [ ] Set up public GLOW using its supplied excitation physics and benchmark
  the available reference responses. Historical cross-section matching is not
  a prerequisite for the requested modern reconstruction.
- [ ] Reuse locally cached space-weather drivers: daily F10.7 and F10.7a,
  with the selected daily/3-hourly Ap convention documented. Preserve driver
  changes when assigning observations to calculation bins.
- [ ] Build orbit-specific response tables across mean electron energy and
  occupied geographic latitude, longitude, UT and DZA bins, conditioned on
  F10.7, F10.7a and Ap. Compute separate WIC/SI13 count responses, then their
  ratio; retain count responses for subsequent energy-flux estimation.
- [ ] Reuse one MSIS profile for each occupied latitude/longitude/UT/driver
  group across DZA bins. Handle DZA in line-of-sight emission/absorption and
  sensor-response calculations. Avoid per-pixel full model evaluations and
  avoid calculating the entire unoccupied multidimensional grid.
- [ ] Expose latitude, longitude, UT and DZA bin sizes. Benchmark runtime and
  response/retrieval differences against finer bins and representative direct
  calculations before choosing defaults. No bin sizes have been approved.
- [ ] Apply the appropriate response within each orbit at detector level,
  preserving individual measured counts and geometry. Binning here shares
  model evaluations; it does not authorize averaging the observations before
  nonlinear retrieval. Record model, driver, calibration and bin provenance.
- [ ] Validate against the fixed Frey baseline and independent DMSP/FAST
  comparisons before adopting the dynamic relation in production. Identify
  any separate SI12/proton-response work required by the supplied curves.

Next actionable step: receive/locate the response curves and establish what
they contain; then choose the smallest runnable forward-model benchmark.
