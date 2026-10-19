# IMAGE IDL Pointing Metadata

Last reviewed: 2026-10-19
Repository evidence: `modular_pipeline` at `917a50a`

## Server audit: orbits197/365 expose SPIN and unmapped-frame limits

Confirmed from the user's server diagnostic output, not locally reread IDL:

- SI12 orbit197 files `s1220002032117.idl` and `s1220002032119.idl` have
  zero valid stored, reconstructed and common DZA pixels (21:17:32.127 and
  21:19:33.132 UTC on2000-07-21). These frames have no usable viewing geometry;
  this does not establish whole-orbit geometry failure.
- WIC197 repeats22503 in `wic20002030825.idl`/`wic20002030827.idl`,22536
  in `wic20002030828.idl`/`wic20002030829.idl`, and22806 across
  `wic20002032119.idl`/`wic20002032120.idl`/`wic20002032121.idl`.
- WIC365 repeats54609 in `wic20003031225.idl`/`wic20003031227.idl`;
  SI13 repeats54609 in `s1320003031225.idl`/`s1320003031227.idl`.
  The repeated-counter records have different increasing timestamps; SPIN
  alone is not evidence of duplicate images. SI12/365 reports no violations.

Cross-camera counter inconsistency is also visible: at21:19:33 UTC WIC
records22806 while SI12 records22913, only48ms apart. The SI12 record is one
of the unmapped frames, so its other metadata are not independently validated.
Nevertheless, a globally unique shared-SPIN key is not justified for all
reported records. At the audit the FUVPY corrector used SPIN in validation,
segmentation, phase construction and reference lookup; deleting the guard
alone would not resolve those assumptions.

The user approved excluding those two SI12 files, retaining repeated-counter
frames and associating sensors by precise time. The implementation now adds
FUVPY read_idl option skip_unmapped=False, explicitly enabled by the2D native
writer. It skips only files with zero valid stored and reconstructed DZA,
warns their filenames and records skipped_unmapped_files on the dataset.
Partial/inconsistent geometry remains an error; strict loading is the default.

FUVPY now accepts finite integer SPIN without strict increase and associates
each camera frame with the nearest precise WIC observation time, independently
of counter equality. Corrections are evaluated at each sensor's own time
inside a supported WIC segment. Original counters and repeated frames remain.
The existing segmentation on WIC counter steps other than+1 and time gaps
>1.5median cadence is retained. It separates pointing smoothing only; the
temporal background model is not split. Fewer than eleven frames means no
pointing correction for that segment, but ordinary background fitting remains
available. Equidistant matches and unsupported intervals are not extrapolated.

138 FUVPY tests and10 native-writer tests pass. New coverage reproduces absent
geometry, repeated/offset/reset counters and independently fitted segments.
Actual197/365 raw files remain remote: the server rerun is still required.
Full85/968 metadata regressions match committed FUVPY bitwise for every data
variable in all cameras, including shared outlier repairs. The full85 native
writer also preserves every saved data field and both background fits.

## Cross-camera pointing and shared orbit968 position anomaly

Physical expectation: approximately fixed relative camera rotations at a
common spacecraft attitude, with instrument-specific alignment and detector
geometry. SI12 and SI13 are two channels of the same SI instrument. The
imagers use TDI composites and have different time responses (about10 seconds
per Earth point for WIC and5 seconds for SI); IMAGE rotates about3 degrees
per second. Small timestamp/phase differences therefore matter. Primary
sources: Mende et al.2000, system design, pp.252,255,265-267,
https://doi.org/10.1023/A:1005271728567, and the NASA SI abstract,
https://image.gsfc.nasa.gov/publication/abstract/2000_mende_etal_03.txt.

Correction scope: a shared spacecraft-attitude repair should be propagated
through each camera's alignment and timing. A spacecraft-position repair
changes Earth intersections and derived coordinates/SZA/DZA for all cameras,
not the inertial ray orientation itself. A WIC-specific alignment correction
belongs to WIC. A limb-derived rotation alone does not identify which cause
it compensates. The WIC-only trials have not automatically updated SI12/SI13;
each sensor's fields must be explicitly recalculated from the corrected state.

All indexed original IDLs in the four local test orbits were inspected;
three-way common SPIN counts are225/213/271/326 for85/86/261/968. The saved
alignment triplets (azimuth, co-elevation, roll), degrees, are:

| Orbit/segment | WIC | SI12 | SI13 |
| --- | --- | --- | --- |
| 85 and86 | (42.93,0.98,-0.05) | (43.87,-0.72,0) | (42.82,-0.54,0) |
| 261 | (43.23,0.91,-0.05) | (42.77,-0.78,0) | (43.64,-0.60,0) |
| 968 before SPIN25524 | (43.25,1.30,-2.45) | (42.63,-0.96,-2.45) | (43.28,-0.78,-2.45) |
| 968 from SPIN25524 | (43.05,1.15,-2.05) | (42.62,-0.94,-2.45) | (43.15,-0.81,-2.45) |

85/86/261 have constant triplets within each orbit. The968 switch occurs at
files *20012950002.idl, around2001-10-22 00:02:04 UTC: WIC full frame294,
SI full frame293. Recovered get_inst_angles_p.pro:88-100 selects exact or
nearest year/day, and getudf_var.pro:369-400 stores its result or overrides.
Current SI tables reproduce both triplets. Current wic_angles.txt day295 is
(43.18,1.19,-2.45), different from the saved record: the exact historical
WIC table revision/override remains unresolved. The step supports a
calibration-selection change, not a demonstrated physical movement.

Use full camera-to-GCI matrices from check_pointing_reconstruction.rotation
and account for exact time/phase before comparing R_WIC.T@R_SI. Matching SPIN
worked in the four local test orbits but is not a general corpus assumption;
the server audit above demonstrates repeated/inconsistent counters.
Setting all sensors to one common attitude makes relative orientation
constant by construction wherever alignment triplets are constant; that
checks metadata encoding, not absolute or image-derived pointing accuracy.
EPOCH's effective exposure-centroid convention remains uncertain even though
the producer evaluates SPIN_PHASE at an integration midpoint. The WIC-only
rays()/geometry() grid must not be applied directly to saved128x128 SI data.

The known WIC full frame292 anomaly is shared by all three channels at
SPIN25522: wic20012942357.idl, s1220012942357.idl, s1320012942357.idl store
exactly O_GCI=[11730.051,2614.008,32008.720] km. Their times differ by at
most16 ms. Independent spot inspection confirms exact equality. A cubic
through the two neighboring spins each side, excluding this spin, gives
position residuals436.185/436.188/436.138 km for WIC/SI12/SI13. This is a
metadata consistency result, not an independent ephemeris. Cross-camera
agreement cannot detect or repair this shared error. It is distinct from the
alignment-table switch two spins later. No SI correction or background fit
was performed in this check.

## Confirmed contents

Direct `scipy.io.readsav()` inspection of local WIC, SI12 and SI13
`imageinfo` records confirms the following fields. Nine frames (three per
sensor, spanning May 2000 and October 2001) contain every field listed here.

| Fields | Meaning established from recovered FUVVIEW3 source |
| --- | --- |
| `O_GCI` | Spacecraft position vector in geocentric inertial coordinates, km |
| `A_GCI3` | Unit direction of the true spacecraft spin axis in GCI at snapshot time |
| `SPIN_PHASE` | Spin phase angle in degrees |
| `SC_SV` | Spin-axis direction in spacecraft coordinates |
| `INST_AZIMUTH`, `INST_CO_ELEV`, `INST_ROLL` | Camera-to-spacecraft alignment angles, degrees |
| `L0` | GCI look unit vector for a central detector pixel |
| `RA`, `DEC` | Right ascension and declination of that central look, degrees |
| `HFOV`, `VFOV` | Horizontal and vertical field of view, degrees |
| `ANGL_RES_R`, `ANGL_RES_C` | Angular row/column increments used by the pointing calculation, degrees |
| `GLAT`, `GLON` | Per-pixel geographic emission-shell locations |
| `MLAT`, `MLON`, `MLT` | Per-pixel magnetic coordinates/local time |
| `DZA`, `SZA` | Per-pixel spacecraft and solar zenith angles at the emission shell, degrees |
| `EMIS_HGT` | Assumed emission-shell height, 130 km in all nine inspected frames |
| `TIME`, `EPOCH`, `SPIN` | Image time and spin identifier |

Recovered source is outside the repository at
`/home/bing/Downloads/fuview3/fuview3`:

- `fuv_init.pro:263–310` defines the record and comments its fields.
- `fuv_calc_pointing.pro:31–67` saves the pointing outputs.
- `fuv_ptg.pro:255–398` constructs detector rays, rotates them into GCI,
  intersects the emission ellipsoid and returns geographic coordinates.
- `fuv_rotation_matrix.pro` owns the camera/spacecraft/GCI rotations.

The nine `L0` vectors reproduce their stored `RA`/`DEC` to about
1.4e-6 degrees using ordinary Cartesian-to-spherical conversion. This checks
internal consistency, not independent pointing accuracy. The precise inertial
epoch/precession convention has not been audited; do not silently substitute
a modern J2000 transform.

## Concrete example

At approximately 2000-05-16 18:01:16 UTC, files
`example_data/wic_data/wic20001371801.idl`,
`example_data/s12_data/s1220001371801.idl` and
`example_data/s13_data/s1320001371801.idl` contain:

| Sensor | Azimuth | Co-elevation | Roll | Central RA | Central DEC |
| --- | --- | --- | --- | --- | --- |
| WIC | 42.93 | 0.98 | -0.05 | 194.3056 | -43.5805 |
| SI12 | 43.87 | -0.72 | 0 | 191.9331 | -44.5933 |
| SI13 | 42.82 | -0.54 | 0 | 192.2080 | -43.5493 |

All angles above are degrees. The common stored `O_GCI` is approximately
`[32319.01, 7295.098, 39742.012]` km. Camera central looks need not be
identical; these differences alone do not establish a registration error.

## Reconstruction and retention limits

The metadata and recovered source provide inputs for a reconstruction of
per-pixel rays and a different assumed emission height. A later WIC baseline
prototype now reproduces saved SZA/DZA in nine frames to <0.00003 degrees
over DZA<75 pixels. [[IMAGE Pointing - Background Stabilization]] owns this
feasibility result and its limits; no improved pointing or new emission height
has been validated.
No quaternion, attitude covariance, explicit pointing-error estimate or
pointing-quality flag is present in the inspected schemas.

Important SI format detail: actual saved image/geolocation arrays are 128 by
128, while `ROWS`, `COLS` remain 256 and angular increments describe that
256-format pointing grid. `downshift.pro` copies these scalar fields unchanged
while resampling arrays with `congrid(...,128,128)`. A reconstruction must
respect this producer resampling and detector-axis orientation; applying the
stored angular increments directly to the saved 128-grid would be wrong.

The inspected live dependency `../fuvpy/fuvpy/src/load.py:126–183` retains
pixel geographic/magnetic coordinates and SZA/DZA. With `viewing_geometry=True`
it uses `O_GCI`, `LAT_CTR`, `LON_CTR`, `EMIS_HGT` to derive solar-view factors,
but does not retain the original position/attitude/boresight fields as output
variables. For that audit, read original IDLs. In the recovered producer,
`LAT_CTR`/`LON_CTR` are assigned the spacecraft's Earth-fixed direction;
their initialization comments describing FOV centre are misleading.

## Verified example pointing check

User requested checking pointing changes in the example orbits. New standalone
`scripts/debugging/check_example_pointing.py` reads all three sensor manifests
and original IDLs. It computes normalized-vector angular changes with
`atan2(norm(cross), dot)`, intervals/rates, spin-axis changes, wrapped spin-phase
increments, camera-alignment ranges and central off-nadir angle relative to
`-O_GCI`. The first frame has no predecessor. Gaps above 1.5 times the median
cadence are separately flagged and excluded from consecutive-pair summaries.

Run in the existing `tConductance` environment:

```bash
python scripts/debugging/check_example_pointing.py --output figures/debugging/example_pointing_20261016
```

All 129 local IDLs are present: 20 frames per sensor for each of orbits 0085
and 0086, and three per sensor for 0968. These are excerpts, not complete
orbits. All 120 consecutive pairs have approximately two-minute cadence;
none meets the gap criterion. Their order matches all nine example sensor
NetCDFs; the NetCDF times omit the IDLs' fractional seconds.

| Orbit | Sensor | Median look step (degrees) | Maximum look step (degrees) | Maximum transition (zero-based) |
| --- | --- | --- | --- | --- |
| 0085 | WIC | 0.459 | 1.281 | 6 to 7 |
| 0085 | SI12 | 0.361 | 1.345 | 0 to 1 |
| 0085 | SI13 | 0.507 | 1.246 | 0 to 1 |
| 0086 | WIC | 0.345 | 0.965 | 15 to 16 |
| 0086 | SI12 | 0.340 | 0.918 | 9 to 10 |
| 0086 | SI13 | 0.340 | 1.218 | 9 to 10 |
| 0968 | WIC | 0.868 | 0.874 | 0 to 1 |
| 0968 | SI12 | 0.866 | 0.911 | 1 to 2 |
| 0968 | SI13 | 0.868 | 0.909 | 1 to 2 |

0085/0086 exhibit irregular changes, also visible relative to Earth. Maximum
WIC off-nadir steps are 1.241 and 0.918 degrees, respectively. Largest WIC
look steps occur at 2000-05-16 03:07:31.584 and 18:29:12.948 UTC. Across all
sensors/frames, camera azimuth/co-elevation/roll are constant within each
excerpt; spin-axis steps never exceed 0.046 degrees. Look-step magnitudes
track absolute spin-phase increments within 0.015 degrees throughout. Thus
stored spin phase dominates this pointing diagnostic, rather than changes in
mounting angles or a large spin-axis reorientation. Orbit 0968's two steps
are similar; three frames cannot establish a full-orbit stability result.

These are changes in the stored pointing solution. Physical pointing motion
versus upstream timing/phase estimation has not been resolved. No quantitative
background-quality comparison or subtraction rerun was done, so a causal
connection to background failures is unverified. Next: compare the identified
transitions with the user's visibly problematic frames and their raw/corrected
images; an attitude-only diagnostic also misses pure rotation about `L0`.

Outputs: `pointing_examples.png`/`.pdf`, `pointing_frames.csv` and
`pointing_summary.csv` under `figures/debugging/example_pointing_20261016/`.
The PNG was visually reviewed. Angular helper checks include normalized and
scaled vectors, zero/90/180-degree separations, a tiny angle and rejection of
zero vectors. All nine peak pairs independently agree with clipped-arccos
calculations to 1e-9 degrees. Original data/products/backgrounds were untouched.
No processing defaults, model selection or publication status changed.

## Direct image comparison at the pointing transitions

The user requested seeing whether the images themselves move. New
`scripts/debugging/plot_example_pointing_frames.py` reads original IDL `IMAGE`
arrays and the existing example sensor NetCDF `dgimg` and `shimg`, matching
every frame by time at the product's whole-second precision. No image
registration, spatial smoothing or new background fitting is applied.

It plots four neighbouring frames around each sensor's largest transition:
WIC 0085 frames 5–8, WIC 0086 frames 14–17, SI12/SI13 0085 frames 0–3 and
SI12/SI13 0086 frames 8–11. Orange borders identify the peak pair. Detector
axes and colour limits are fixed across each complete sensor/orbit excerpt;
BS/SH share signed limits. Grey denotes unavailable pixels. Original IDL
`IMAGE` is before fuvpy, not an uncalibrated detector-count product.

Outputs under `figures/debugging/example_pointing_20261016/frames/` include
six PNG strips, a six-page `largest_pointing_transitions.pdf`, six full
20-frame GIFs and six two-frame peak-pair blink GIFs. `display_manifest.csv`
records selections and scales. Reproduce with:

```bash
python scripts/debugging/plot_example_pointing_frames.py --pointing figures/debugging/example_pointing_20261016
```

Visual review of both WIC strips confirms leftward movement of the original
Earth image across the highlighted pairs, also visible in the stored
subtractions. This supports a real change of the detector scene; it does not
independently validate the pointing solution or establish why subtraction
quality changes. No general registration estimate or background-quality
correlation was performed. All six GIFs have 20 frames, all blinks have two,
all PNGs reopen, and all six source/product frame sequences match in order.
Original example data/products remain unchanged by this plotting task.

## SZA/DZA motion versus the detector scene

The user clarified that registration of the background coordinates matters
more than pointing changes themselves. Checked the same two WIC peak pairs
directly in original IDLs. `scripts/debugging/check_pointing_coordinate_motion.py`
fits independent translations to the image and the SZA/DZA fields, and writes
contour overlays plus `wic_coordinate_motion.csv` under
`figures/debugging/example_pointing_20261016/coordinates/`.

| Pair | Image row/column shift (pixels) | Joint SZA/DZA row/column shift (pixels) |
| --- | --- | --- |
| 0085, 6 to 7 | +1.781 / -14.801 | +0.325 / -16.128 |
| 0086, 15 to 16 | -0.841 / -6.352 | +0.101 / -11.521 |

Signs indicate where an earlier feature appears in the later detector frame:
positive row is upward in the plots, positive column rightward. The angle
fields move in the same broad direction as the image, but the local rigid
translation approximation gives a roughly 5.2-pixel horizontal discrepancy
in 0086 and a smaller roughly 1.3-pixel discrepancy in 0085. The two cases
also have approximately 1-pixel row discrepancies. These are relative-motion
diagnostics, not a calibrated absolute pointing-error estimate.

Image method: sigma-2 Gaussian smoothing for motion estimation only, rows and
columns 30:190 sampled every two pixels, robust least squares with translation,
brightness gain and offset. Three independent initial horizontal shifts agree;
linear/square-root/log intensity comparisons give similar translations. The
plots themselves show unfiltered original images. Geometry method: match
SZA/DZA together without angle offsets on common interior support, at least
26 pixels from invalid coordinates, with a bounded translation fit. Joint
angle residual RMS is 0.061/0.053 degrees. Original SZA/DZA arrays match every
finite stored-preprocessing pixel within 1e-5 degrees, so no later reader
recalculation explains the discrepancy.

Independent read-only verification used strong-limb-gradient and interior-disk
fits. It found image column shifts -14.73 and -6.23 to -6.41 pixels. DZA
50/60/70/80-degree contour intersections shift approximately -16.3 and -11.65
pixels, respectively. Separate limb row bands agree with the image estimate;
the 0086 result is not driven solely by auroral variability or fill values.
Same-record IMAGE/SZA/DZA, filenames, frame indices, consecutive spin numbers
and approximately 119-second time intervals were independently checked.

Known synthetic translations, including image gain/offset, were recovered
within 0.1 pixel; the angle translation check agrees within 0.005 pixel.
Both PNG overlays were reviewed. Right-panel dashed contours are the earlier
angles translated with the measured image motion; solid contours are the
actual later angles. Their separation visualizes the relative mismatch.

Inference: a frame-dependent image/geometry registration issue is a concrete
candidate contributor to subtraction variability, especially for 0086. This
does not establish a faulty upstream quantity or a causal background result:
natural geometry evolution, approximate translation and morphology remain
limitations. No correction, background refit or production/default change was
made. Next: broaden this comparison before adopting a registration correction
and audit the producer's phase/timing-to-image mapping. Recovered
`getudf_var.pro:87–103` interpolates spin phase at integration midpoint and
permits timing/phase overrides; this is a provenance lead, not a proven bug.

## Pointing provenance and the early-2000 warning

Source trace 2026-10-16: pointing is calculated upstream in the original
FUVVIEW3 ground-processing/export stage. The local fuvpy reader imports saved
pixel coordinates and SZA/DZA; it does not determine those angles by fitting
the image. Recovered source outside this repository is
`/home/bing/Downloads/fuview3/fuview3`.

Confirmed calculation chain:

1. `getudf_var.pro:60–103` reads FUV Level-1/UDF records using
   `fuv_read_udf`. Spin phase is interpolated between `pm1_spinphase` and
   `pm2_spinphase`, using their valid METs, at the midpoint of the recorded
   integration period. The code also calculates a spin-rate-extrapolated phase
   but saves the interpolated phase. Timing/phase override variables exist;
   the main procedure does not declare them as keyword arguments, so their
   effective external interface is not established.
2. `getudf_var.pro:188–208` similarly interpolates `pm1_sv_*`/`pm2_sv_*` for
   the spin-axis vector. For 2000/2001 and early 2002 the recovered code uses
   this older path. Its later branch reads satellite star-tracker quaternions
   from `IMS4BYTE`/`IMS2BYTE` and recalculates the spin axis. Separate seasonal
   fallback branches use fitted `spin_axis_YYYY.sav` files.
3. `get_inst_angles_p.pro` loads camera azimuth/co-elevation/roll from the
   `attitude/wic_angles.txt`, `s12_angles.txt` and `s13_angles.txt` databases.
   It uses an exact date or the nearest listed day, rather than interpolating
   alignment angles in time. For day 137 of 2000, the nearest entries are day 145
   (2000-05-24). All three corresponding angle triples match the saved local
   IDLs within 1e-5 degrees. This supports consistency with these tables, not
   identification of the exact exporter revision.
4. `getudf_var.pro:369–437` supplies those angles, spacecraft position and
   camera geometry to `fuv_calc_pointing.pro` and `fuv_ptg.pro`. Camera rays
   are rotated to GCI, intersected with the emission shell, and converted to
   pixel geolocation and SZA/DZA before the record is saved.

**Confirmed source warning:** `main_scr3.pro:1671–1672` explicitly checks
`year == 2000` and `doy < 157`, displaying
“Warning! Pointing unreliable for this day”. Both 0085 and 0086 are
2000-05-16/day 137, inside this interval (before 2000-06-05). Thus this recovered
software would flag both example orbits. The reason for that particular early
interval is not stated in the warning, and it must not be conflated with the
separately described later star-tracker summer outages. No pointing-warning
flag is present in the inspected exported IDL schemas.

Independent published context: Frey et al. (2004), *Substorm onset observations
by IMAGE-FUV*, §7 describes regular pointing corrections using bright UV stars
and residual errors up to 4 pixels in the spin plane and 2 perpendicular to it.
§10 dates regular IMAGE-FUV operations from 2000-05-19, after these examples.
These are historical limitations, not a direct error model for our exports.
Source: https://doi.org/10.1029/2004JA010607 (paragraph 7 also returned in the
indexed full-paper copy at CiteSeerX).

NASA's prelaunch mission requirements identify onboard spin-axis/phase
determination using an autonomous star tracker and three-axis magnetometer:
https://image.gsfc.nasa.gov/publication/document/dmr/image_dmr_3.html (§1320).
Use this for the mission architecture, not a verified operational accuracy
claim. Colin Forsyth's IMAGE-FUV data tutorial independently lists the later
seasonal pointing limitations, but not the early 2000 warning:
https://www.issibern.ch/teams/dataassimilation/PDFs/forsyth_image_fuv_issi.pdf
(slide 8).

The exact source/database versions, UDF inputs and any manual corrections
used for our saved IDLs remain unverified. The control below is now checked;
broader comparisons and producer phase/timing/calibration provenance remain
open. These source findings do not establish a corpus-wide pointing failure
or a specific hardware/software fault. The source audit changed memory only.

## Earlier short control outside the warning intervals: orbit 0968

The later full local control is owned by [[IMAGE Pointing - Full Local Controls]].
It confirms mismatches elsewhere in both 0261 and 0968, including a roughly
nine-pixel excursion coincident with a stored-position discontinuity. The
two short-control transitions below remain locally consistent; they do not
represent the longer orbit.

The original image/coordinate-motion comparison had only covered 0085/0086;
the earlier 0968 analysis checked look-vector steps only. After the user asked
about an outside-interval control, both adjacent WIC pairs in the three repository-example
0968 frames were checked. They span 2001-10-21 23:37:32.413–23:41:37.486 UTC,
outside the early-2000 and 2001 seasonal warning intervals in the recovered code.

| WIC pair | Image row/column shift (pixels) | Joint SZA/DZA row/column shift (pixels) | Absolute horizontal difference |
| --- | --- | --- | --- |
| 0968, 0 to 1 | -0.053 / -6.782 | -0.218 / -6.562 | 0.220 |
| 0968, 1 to 2 | +0.049 / -6.189 | -0.199 / -6.236 | 0.047 |

The changing aurora occupies the original central image-fit region, so the
control uses the visible dayglow/limb region: rows 100:215, columns 165:248.
The selected region excludes the bright source near row 80. Gaussian sigma 2,
robust gain/offset fitting and the joint-angle calculation are otherwise
unchanged. The initial central-region image fit had larger residuals and
horizontal differences 0.15/0.75 pixels; do not use that evolving-aurora fit
as a precise control. The adopted limb image RMS is 0.030/0.016 in scaled
units; angle RMS is 0.234/0.230 degrees, reflecting nontranslation geometry.

Independent read-only checks confirm the result: strong right-limb gradient
row bands give about -6.7 to -6.8 and -6.1 to -6.3 pixels; dayglow and limb-only
photometric fits remain stable over Gaussian sigma 1/2/3. DZA 50/70/80-degree contours
move approximately -6.58/-6.40/-6.35 and -6.30/-6.12/-6.07 pixels. Variation
with contour level demonstrates natural geometry evolution beyond translation.
No approximately 5-pixel discrepancy is reproduced in these two transitions.
This is a four-minute control, not absolute pointing validation or evidence
that all data outside warning intervals are accurate.

The diagnostic now accepts selected pairs, output directory and image-fit
region, and records the latter in its CSV. Image samples remain inside the
detector during fitting; earlier default-region results are unchanged to
1e-8 pixels. Reproduce in the existing tConductance environment:

```bash
python scripts/debugging/check_pointing_coordinate_motion.py --pointing figures/debugging/example_pointing_20261016 --pairs 968:1 968:2 --image-region 100 215 165 248 --output figures/debugging/example_pointing_20261016/coordinates_0968
```

Both pair overlays, a two-page PDF and `wic_coordinate_motion.csv` are under
that output directory. PNGs were visually reviewed, finite stored/source
angles still match, and the original 0085/0086 output files are preserved.
No source data, calibration, background, production default or model choice
was changed. The longer local controls have now been checked in the linked note.

## Larger local IDL inventory found

Confirmed 2026-10-16 after the user asked whether only 0085/0086/0968 were
available. Earlier pointing work inspected the small repository `example_data`
sample only: 43 IDLs per sensor, comprising 20/20/3 frames for those orbits.
That is not the full locally accessible raw-data inventory.

The existing working directory
`/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test` has four larger
original-IDL sets and matching HDF orbit indices:

| Orbit | Observation dates | WIC IDLs | SI12 IDLs | SI13 IDLs |
| --- | --- | --- | --- | --- |
| 0085 | 2000-05-15 to 2000-05-16 | 236 | 226 | 235 |
| 0086 | 2000-05-16 | 213 | 213 | 213 |
| 0261 | 2000-08-28 | 272 | 271 | 272 |
| 0968 | 2001-10-21 to 2001-10-22 | 328 | 326 | 326 |
| Total | | 1049 | 1036 | 1046 |

All indexed filenames are present in `wic_data`, `s12_data` and `s13_data`;
there are no unindexed local IDLs in those folders. These are source-file
counts, not accepted/preprocessed science-frame counts. The working
`orbitdates.csv` independently lists the four orbit intervals. The standard
`/home/bing/Dropbox/work/data/IMAGE_FUV` currently has downstream products but
no direct raw-IDL sensor folders. No exhaustive whole-filesystem/server search
was performed or needed to establish the additional local availability.

0261 and the longer 0968 sequence provide locally available controls outside
the warning intervals identified in recovered FUVVIEW3. The subsequent full
WIC check is owned by [[IMAGE Pointing - Full Local Controls]]. The inventory
check itself performed no motion test, data transfer or background fit.
