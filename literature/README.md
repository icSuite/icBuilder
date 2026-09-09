# IMAGE-FUV literature

Last reviewed: 2026-08-30

This folder contains the papers currently used to check the IMAGE-FUV
processing and the WIC/SI13 electron-energy retrieval. File names are shortened
descriptions; use the DOI below when citing a paper.

## DZA conclusion

- `Ohma_et_al_2024_fuvpy_background_removal.pdf` gives the clearest numerical
  guidance. Its optically thin **dayglow model** uses
  `cos(SZA) / cos(DZA)`. The authors say an upper DZA limit is required and
  give 70--80 degrees as an example because the secant term diverges and
  geolocation becomes unreliable. They do not prescribe 70 degrees. This is not a validation of the
  WIC/SI13 energy ratio below that limit. DOI: `10.26464/epp2023051`.
- `Frey_et_al_2003_IMAGE_FUV_quantitative_interpretation.pdf` states that the
  tabulated electron response was integrated for a nadir observation. It also
  shows that disturbed atmospheric composition can preferentially absorb SI13,
  raise WIC/SI13, and imply unrealistically high electron energies. No DZA
  cutoff is given. DOI: `10.1023/B:SPAC.0000007521.39348.a5`.
- `Meurant_et_al_2003_IMAGE_electron_precipitation.pdf` says the quantitative
  reduction includes line-of-sight and atmospheric absorption in the forward
  modeling. Its published WIC/SI13 response curve is for vertical viewing,
  while apparent emission can be calculated for a defined viewing geometry.
  The paper does not document a separate per-pixel DZA correction of the
  electron ratio, and no DZA cutoff is given. DOI: `10.1029/2002JA009685`.
- `Grocott_et_al_2006_oblique_IMAGE_viewing.pdf` rejects an IMAGE interval once
  its view is too oblique for meaningful mapping, but gives no angle.
  DOI: `10.5194/angeo-24-3365-2006`.
- `Ostgaard_et_al_2018_asymmetric_geospace.pdf` cautions that oblique auroral
  intensities depend on auroral structure and are not equivalent to a simple
  nadir view. It gives no retrieval cutoff.
  DOI: `10.5194/angeo-36-1577-2018`.

The literature therefore supports a near-limb mask, but not a common
`cos(DZA)` correction as a solution to the WIC/SI13 problem. A common factor
cancels in the ratio. The remaining angle dependence is channel-specific:
slant atmospheric absorption, different emission-altitude weighting, auroral
structure along the line of sight, and instrument response.

The present `icBuilder` input step uses `DZA < 75 degrees`. This lies inside
Ohma et al.'s example 70--80 degree range for background-model stability,
but no paper reviewed here establishes 75 degrees as a safe quantitative
WIC/SI13 energy-retrieval boundary.

## Local papers

- `Coumans_et_al_2002_IMAGE_NOAA_proton_correction.pdf` -- simultaneous IMAGE
  and NOAA validation, including proton contamination and line-of-sight
  modeling. DOI: `10.1029/2001JA009233`.
- `Coumans_et_al_2004.pdf` -- conductance reconstruction and comparison with
  particle measurements. DOI: `10.5194/angeo-22-1595-2004`.
- `Frey_et_al_2003_IMAGE_FUV_quantitative_interpretation.pdf` -- central
  quantitative interpretation and response tables.
- `Frey_et_al_2017_UV_calibration.pdf` -- wide-field UV calibration and
  flat-field context. DOI: `10.1002/2016JA022700`.
- `Gasparini_et_al_2024.pdf` -- recent fuvpy-based IMAGE conductance workflow;
  it gives count-based retrieval fallbacks but no DZA rule.
  DOI: `10.1029/2024JA032599`.
- `Gerard_et_al_2001_SI12_viewing_geometry.pdf` -- SI12 response depends on
  viewing direction and proton pitch-angle/energy distribution; useful for the
  proton-correction side of the pipeline. DOI: `10.1029/2001JA900119`.
- `Grocott_et_al_2006_oblique_IMAGE_viewing.pdf` -- practical example of
  rejecting an excessively oblique IMAGE view.
- `Mende_et_al_2000_IMAGE_FUV_system_design.pdf` -- FUV system and physical
  basis of the spectral measurements. DOI: `10.1023/A:1005271728567`.
- `Mende_et_al_2000_IMAGE_SI.pdf` -- SI12/SI13 design and response.
  DOI: `10.1023/A:1005292301251`.
- `Meurant_et_al_2003_IMAGE_electron_precipitation.pdf` -- application and
  validation of the WIC/SI13 energy method.
- `Ohma_et_al_2024_fuvpy_background_removal.pdf` -- method implemented by
  fuvpy for dayglow and residual-background removal.
- `Ostgaard_et_al_2018_asymmetric_geospace.pdf` -- processing provenance for
  the IMAGE dataset used in this project.

## Hardy ion-precipitation model

- `Hardy_ion_model_sources.md` -- bibliography, confirmed model structure,
  reconstruction notes, and implementation status.
- `hardy_et_al_1989_A_statistical_model_of_auroral_ion_precipitation.pdf` --
  original statistical ion model and validation maps.
  DOI: `10.1029/JA094iA01p00370`.
- `hardy_et_al_1991_a_statistical_model_of_auroral_ion_precipitation_2_functional.pdf`
  -- complete functional representation and coefficient table.
  DOI: `10.1029/90JA02451`.
- `Hardy_et_al_1991_auroral_particle_review.pdf` -- public overview with
  Kp=3 ion energy-flux, number-flux, and average-energy maps.
  DOI: `10.5636/jgg.43.Supplement1_337`.
- `NASA_1991_Auroral_Particle_Models_Catalog.pdf` -- NASA-TM-105052,
  contemporary documentation of the AFGL ion-precipitation model and its
  historical FORTRAN distribution.

The primary papers are sufficient for implementation. The 1991 paper supplies
all Fourier coefficients and explicitly defines average energy as integral
energy flux divided by integral number flux. The historical FORTRAN remains
useful only as an independent comparison.

## Relevant papers not stored locally

- Mende et al. (2000), *Far ultraviolet imaging from the IMAGE spacecraft. 2.
  Wideband FUV imaging*. DOI: `10.1023/A:1005227915363`.
- Hubert et al. (2002), *Total electron and proton energy input during auroral
  substorms: Remote sensing with IMAGE-FUV*. Its forward efficiencies include
  viewing geometry. DOI: `10.1029/2001JA009229`.
- Galand and Lummerzheim (2004), *Contribution of proton precipitation to
  space-based auroral FUV observations*. DOI: `10.1029/2003JA010321`.
