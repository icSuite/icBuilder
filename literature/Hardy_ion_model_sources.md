# Hardy auroral-ion model sources

Last reviewed: 2026-08-30

## Primary papers

The Hardy proton model is defined by two papers:

1. Hardy, Gussenhoven, and Brautigam (1989), *A statistical model of
   auroral ion precipitation*, JGR 94(A1), 370--392.
   DOI: https://doi.org/10.1029/JA094iA01p00370
2. Hardy, McNeil, Gussenhoven, and Brautigam (1991), *A statistical model
   of auroral ion precipitation: 2. Functional representation of the
   average patterns*, JGR 96(A4), 5539--5547.
   DOI: https://doi.org/10.1029/90JA02451

Both papers are now stored locally. The 1989 paper provides the statistical
maps and scientific basis; the 1991 paper is the implementation paper.

## Documents stored locally

- `hardy_et_al_1989_A_statistical_model_of_auroral_ion_precipitation.pdf`
  contains the original statistical maps and the processing details used to
  construct them.
- `hardy_et_al_1991_a_statistical_model_of_auroral_ion_precipitation_2_functional.pdf`
  contains the complete Fourier coefficient table, Epstein equations, limiting
  background levels, reconstruction instructions, and validation plots.
- `Hardy_et_al_1991_auroral_particle_review.pdf` provides a readable overview
  and Kp=3 maps of ion energy flux, number flux, and average energy.
- `NASA_1991_Auroral_Particle_Models_Catalog.pdf` is NASA-TM-105052,
  *Solar-terrestrial models and application software*. Section 2-22 documents
  the AFGL ion-precipitation model and its historical FORTRAN availability.

## Confirmed model structure

- Inputs: Kp, magnetic local time, and corrected geomagnetic latitude.
- Domain: 50--90 degrees CGL in 30 latitude bins and 48 half-hour MLT bins.
- Activity: seven Kp categories, represented as 0 through 6.
- Source data: about 26.5 million one-second DMSP F6/F7 SSJ/4 ion spectra from
  30 eV to 30 keV.
- Outputs: integral ion energy flux, integral ion number flux, and average
  energy.
- Functional representation: the latitude profile is fit with a generalized
  Epstein function. Its six MLT-dependent parameters are represented by
  sixth-order Fourier series with 13 coefficients each.
- Mean energy is calculated as integral energy flux divided by integral number
  flux, with the units in the paper giving the result directly in keV.
- For intermediate Kp, interpolate the evaluated log-flux values between Kp
  levels. The paper explicitly warns against interpolating the coefficients.
- The published fit generally agrees with the original statistical maps to
  within 25% for energy and number flux.

## Implementation status

The two primary papers contain enough information to implement and validate
the Kp-driven model. The historical FORTRAN would be useful as an independent
cross-check, but it is no longer required. Coefficients should be transcribed
from the PDF table and checked visually because automated PDF text extraction
introduces minus-sign and character-recognition errors.

Do not substitute Lompe's `hardy()` routine: it implements the separate 1987
electron/conductance model, not the 1989/1991 auroral-ion model.
