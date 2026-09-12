# KMT-2026-BLG-0521 / OGLE-2026-BLG-0058

## Astrometry of the microlensing target — summary

Prepared 2026-09-12. Solution: `~/KMTdata/Results/v8/IFfinal_<field>.mat`
(variable `IFsys`), tie in `Tie_<field>.mat`.

### 1. The data and the algorithm

**Data.** Two KMTNet cut-outs of the same field, **BLG41** and **BLG01**,
300 x 300 pixels at 0.4 arcsec/pix, in I band, spanning 2016 to 2026 in ten
observing seasons (2020 is a short stub and is dropped). **19 133** and
**19 981** exposures respectively, taken alternately on the same nights, a
median of 12.2 minutes apart. The target sits at the field centre at
I_OGLE = 18.13, RA 17:52:38.09, Dec −31:47:36.1.

**Calibration set.** Stars are selected on magnitude and isolation:
**14 < I < 19** with **no companion brighter than I = 18 within 2 arcsec**, the
companion list taken from the OGLE catalogue, which resolves about six times
more sources than KMT detects. The set is then **iterated once**: a first
solution is computed, and stars whose residual RMS lies more than **2 sigma**
above the running median for their magnitude are dropped and the fit repeated.
The cut is made in log space, so it is multiplicative and does not simply shave
off the faint end. This leaves **260** stars on BLG41 and **235** on BLG01. The
target passes the selection on its own merits. Everything else is dropped at the
input, so all measurement, detrending and quality cuts act on this set alone.

**The fit.** Each exposure gets a full 2-D affine transformation, 6 parameters,
and each star a position and proper motion, 4 parameters; these are solved
together by alternating least squares. Alongside them the fit removes
differential chromatic refraction (fitted in 6 equal-population colour bins), an
annual term, a pixel-phase term, and two SysRem components computed once over
the whole decade. Every star is freely fitted — nothing is held at a catalogue
value — and each is weighted by its own residual scatter, with moving-median
outlier rejection.

**Absolute frame.** A free fit determines only *relative* astrometry, so the
result is tied to Gaia DR3 afterwards. The tie is fitted on stars with
**RUWE < 1.4 and 15 < I < 17** and removes the full six-parameter gauge — a
constant *and* a linear gradient across the field — then is applied to every
source.

---

### 2. Results

**Achieved accuracy**, residual RMS over the decade after position and proper
motion are removed:

| | BLG41 ΔX / ΔY | BLG01 ΔX / ΔY |
|---|---|---|
| calibration stars brighter than I = 17 | **4.80 / 4.83 mas** | **4.75 / 4.95 mas** |
| the target, I = 18.13 | **23.77 / 19.67 mas** | **21.51 / 20.45 mas** |

The frame is therefore good to under **5 mas** on well-measured stars, and the
target — three magnitudes fainter and blended — is measured to about **20 mas**
per epoch, consistent with the running median of the residual RMS at its
magnitude.

**Proper motion of the target**, absolute, after the Gaia tie:

| | μ_α cos δ | μ_δ |
|---|---|---|
| BLG41 | **−1.966 ± 0.347** | **−7.019 ± 0.358** |
| BLG01 | **−1.514 ± 0.227** | **−7.073 ± 0.376** |
| difference between the fields | −0.452 ± 0.415 | +0.054 ± 0.519 |

mas/yr. The two independent fields agree within their errors, at 1.09 and 0.10
sigma. The quoted error is internal, from the scatter of the ten season means
about the fitted line; the external check against Gaia gives about
**1.0 / 0.62 mas/yr** per star, and that is the figure to use in practice.

**No astrometric anomaly is detected at the target.** Its χ² on monthly binned
residuals is 1500 and 1313 for 180 degrees of freedom, against field medians of
χ²/DoF of 6.5 and 6.8 — the target sits inside the body of the distribution,
not in the tail. Whatever astrometric excursion the microlensing event produces
is smaller than the scatter of an ordinary star of this brightness.

@@FIG report_v8_RMS.png | Residual RMS against the OGLE catalogue magnitude, per axis and per field, after position and proper motion have been removed. Blue: all 260 / 235 calibration stars, every one freely fitted. Green: the stars the Gaia tie is fitted on, RUWE below 1.4 and 15 < I < 17. Red: the target. Black: running median. The dashed lines mark the 14 and 19 mag selection limits. Compared with a set that has not been iterated, the 2 sigma cut removes the scattered points that previously sat well above the median relation.

@@FIG report_v8_target_beforeafter.png | The target before and after proper-motion detrending, in both fields. Columns are the fields; the top row keeps the proper motion and the bottom row removes it. Blue circles are the offset in right ascension, red squares in declination, binned in one sidereal month with the standard error within the bin. The proper motion is plainly visible as the declination trend in the top row and is absent from the bottom row; the residual RMS quoted in the lower panels is per bin.

---

## More details

### 3. Object selection, step by step

| step | BLG41 | BLG01 |
|---|---|---|
| raw sources in the MSc file | 1177 | 1187 |
| survive the standard quality cuts | 594 | 621 |
| have an OGLE I magnitude | 540 | 478 |
| in 14 < I < 19 | 529 | 470 |
| ... and isolated: no I < 18 companion within 2.0 arcsec | 287 | 253 |
| ... and surviving the 2 sigma RMS cut | **260** | **235** |

The isolation cut is the dominant one, removing 46% of the magnitude-selected
sources on both fields. The outlier iteration then removes a further **27 and
18** stars — 9% and 7% — spread across all magnitudes rather than concentrated
at the faint end: 6%, 14%, 10% and 3% by magnitude band on BLG41. Their median
residual RMS is 22.4 and 48.5 mas against 10.8 and 10.6 mas for the stars kept.

For the tie only, Gaia is required: 247 of 260 and 207 of 235 sources have a
counterpart, of which 198 and 163 pass RUWE < 1.4, and 104 and 85 also fall in
15 < I < 17.

### 4. The Gaia tie

**The map.** Gaia standard coordinates are fitted directly onto cut-out pixels
from the tie stars — not carried through OGLE, which would inherit its ~110 mas
astrometric error. Median residual **22.1 mas (BLG41)** and **21.8 (BLG01)**,
with scales of 2.5195 and 2.5289 pix/arcsec against the 2.500 expected at
0.4″/pix.

**The gauge.** In the proper motions the gauge freedom is a constant plus a
linear gradient across the field — six parameters. Removing only the constant
leaves a shear of about −0.036 mas/yr per pixel, reproduced at nearly the same
value in both independently reduced fields. With the full gauge removed, the
scatter of our proper motions about Gaia is **0.994 / 0.608 mas/yr** (BLG41) and
**1.049 / 0.641** (BLG01); removing only the constant would leave 3 to 4 times
that. The gauge constant removed is −0.126 / +5.105 and +2.151 / +4.324 mas/yr.

### 5. What the fit does, in detail

Each of the two steps runs **two passes** internally. `IFsys.AnnualEffect` and
`IFsys.PixPhase` read `0` in the saved solutions, which does *not* mean those
effects were ignored: they are fitted and subtracted from the data in pass 1 and
switched off in pass 2, which works on data from which they have been removed.

| | pass 1 | between | pass 2 (final) |
|---|---|---|---|
| per-epoch affine | fitted | | fitted |
| per-source position and motion | fitted | | fitted |
| DCR / chromatic refraction | fitted | | fitted again |
| annual term | fitted and subtracted | | forced off |
| pixel-phase term | fitted and subtracted | | forced off |
| SysRem | | applied | |
| iterative reweighting | 15 iterations | | 6 iterations |

- **Per-epoch affine** — 6 parameters per exposure in `ParE`, stored as the
  deviation from the identity: translation (2), rotation (1), uniform scale (1),
  shear (2). 520 equations for 6 unknowns each epoch. Translation dominates at
  about 37 / 28 mas rms; rotation and scale each move a star at the cut-out edge
  by roughly 10 mas.
- **Per source** — position at `JD0` = 2019-06-01 and proper motion, 4
  parameters in `ParS`, in pixels and pixels per year at 400 mas/pix.
- **DCR** — `[1, sin(pa)·secz, cos(pa)·secz]` plus six higher orders in
  `(secz·sin pa)ⁿ`, `(secz·cos pa)ⁿ` for n = 2, 3, 4: 18 parameters per colour
  bin, 6 equal-population bins. Only the *differential* part is recoverable —
  refraction common to the field is exactly a translation and is absorbed by
  `t_x, t_y` every epoch.
- **Annual term** — a quartic in the phase of the calendar year, 10 parameters,
  fitted globally, one curve shared by every star. Not parallax, which is
  per-source and switched off.
- **Pixel phase** — a quintic in the sub-pixel phase, 10 parameters, global,
  recoverable only because the pipeline stored the per-epoch registration shifts.
- **SysRem** — 2 components over the whole decade, reused rather than refitted
  per season, with point-wise rejection of bad cells. Refitting per season was
  tested and is worse: see section 8.

Timings and convergence:

| | BLG41 | BLG01 |
|---|---|---|
| step 1, decade SysRem | 16.9 min | 15.4 min |
| step 2, final solution | 31.1 min | 29.1 min |
| convergence, last 3 iterations | 0.0002 mas | 0.0047 mas |

No per-season fitting enters the solution; seasons are used only for reporting.

### 6. Comparison stars

Six, chosen to be the same physical objects in both fields — three near the
target's magnitude and three brighter. All six survive the outlier cut and sit
at identical positions in the two fields:

| tag | I | dist [pix] | decade RMS ΔX / ΔY, BLG41 | BLG01 |
|---|---|---|---|---|
| `m18_d04` | 18.14 | 4.5 | 28.86 / 25.64 | 25.24 / 27.09 |
| `m18_d17` | 18.09 | 17.3 | 44.90 / 23.22 | 80.36 / 27.53 |
| `m18_d48` | 18.05 | 48.4 | 18.06 / 15.57 | 13.86 / 15.58 |
| `m16_d24` | 16.68 | 24.4 | **4.14 / 4.31** | **4.13 / 4.32** |
| `m16_d30` | 16.63 | 30.4 | 6.69 / 5.53 | 6.55 / 5.62 |
| `m16_d31` | 16.75 | 30.6 | 8.42 / 7.82 | 11.16 / 7.84 |

mas. The bright trio reproduces between the fields almost exactly — `m16_d24`
gives −3.22 / −7.85 on BLG41 and −3.22 / −7.85 on BLG01 — while the faint trio
does not: `m18_d17` reaches 45 and 80 mas decade RMS on only ~2400 usable
epochs. It is a poor star.

**`m18_d04` is the neighbour blended with the target**, 4.5 pixels away under
~7 pixel seeing. It is a legitimate member of the calibration set — the target
at I = 18.13 is fainter than the I < 18 companion threshold, so it does not
trigger rejection — but its centroid tracks the target's brightening through the
event, and it is not an independent object.

**Coordinates.** Pixel positions are per field, since the cut-outs are offset by
about 1.1 pixels; sky positions come from inverting each field's own map and the
two fields agree on them to 2 to 9 mas.

| star | I | BLG41 pixel | BLG01 pixel | RA (J2000) | Dec (J2000) |
|---|---|---|---|---|---|
| target | 18.13 | 150.4, 150.4 | 151.6, 150.2 | 17:52:38.085 | −31:47:36.21 |
| `m18_d04` | 18.14 | 146.5, 148.3 | 147.6, 148.0 | 17:52:38.210 | −31:47:37.06 |
| `m18_d17` | 18.09 | 134.5, 143.6 | 135.6, 143.4 | 17:52:38.590 | −31:47:39.00 |
| `m18_d48` | 18.05 | 198.6, 145.4 | 199.8, 145.0 | 17:52:36.585 | −31:47:38.10 |
| `m16_d24` | 16.68 | 161.3, 172.2 | 162.4, 172.1 | 17:52:37.746 | −31:47:27.53 |
| `m16_d30` | 16.63 | 120.4, 154.6 | 121.4, 154.5 | 17:52:39.023 | −31:47:34.60 |
| `m16_d31` | 16.75 | 138.7, 178.7 | 139.8, 178.6 | 17:52:38.436 | −31:47:24.98 |

### 7. Chi-squared of the monthly binned residuals

For every star the residuals are binned in sidereal months; each bin contributes
`(mean/standard error)²` in both axes, and `DoF = 2·N_bins − 4`. Values are not
divided by DoF.

| | BLG41 | BLG01 |
|---|---|---|
| stars measured | 260 | 235 |
| median DoF | 182 | 182 |
| **target χ²** | **1500** (DoF 180) | **1313** (DoF 180) |
| target χ²/DoF | 8.3 | 7.3 |
| field median χ²/DoF | 6.5 | 6.8 |

The target is unremarkable. The field median of χ²/DoF near 7 is itself
informative: the standard error of the mean within a month **understates the
true bin uncertainty by about √7 ≈ 2.6**, because the noise within a month is
correlated rather than independent. The effect is multiplicative and nearly
independent of magnitude, and the faintest stars come closest to unity, since
for them genuine photon noise — which *is* independent between exposures —
dominates.

### 8. A residual systematic that survives every correction

Fitting a straight line to each star's residuals **within each season
separately** leaves slopes that are real and unexplained.

| | BLG41 | BLG01 |
|---|---|---|
| median per-star season-slope rms | 2.90 / 2.80 mas | 2.93 / 2.90 mas |
| expected from white noise | 0.99 mas | 0.88 mas |
| **excess over white noise** | **2.9x** | **3.3x** |

They are **not noise**. Matching 181 stars between the fields through OGLE and
correlating each star's ten-season slope pattern gives a median correlation of
**+0.795** in one axis and **+0.754** in the other, with **96 to 97% of stars
positive**. Two independent reductions, independent frames, essentially every
star showing the same season-by-season pattern: about 79% of the slope variance
is shared systematic rather than measurement noise. The target is entirely
typical in this respect, at the 51st to 75th percentile among stars within
0.4 mag of it.

The cause has not been identified. Six candidates were tested and eliminated:

| candidate | test | result |
|---|---|---|
| DCR / colour | per-season slope against colour | correlation below 0.11 in size |
| DCR colour resolution | full refit with 12 and with 20 colour bins instead of 6 | slopes and cross-field correlation unchanged, see below |
| blending | against nearest-neighbour distance | correlation below 0.05 in size |
| field distortion | spatial correlation; quadratic in position | ~0 at all separations; 1.6 to 2.0% of variance |
| pixel phase | slope against within-season phase drift | r ≈ 0.00, 0% of variance |
| SysRem cadence | refit per season instead of per decade | slopes nearly **double**, to 5.1/7.6 and 5.6/7.7 mas |

The colour-bin refits are the direct test of the DCR model. The 18 DCR
parameters are fitted independently in each colour bin, so with 12 or 20 bins
(19 or 11 stars per bin instead of 37) any colour dependence too fine for six
bins would have been absorbed. Nothing changed:

| | 6 bins (v8) | 12 bins | 20 bins |
|---|---|---|---|
| slope rms, BLG41 | 2.90 / 2.80 | 2.98 / 2.77 | 2.92 / 2.80 |
| slope rms, BLG01 | 2.93 / 2.90 | 2.97 / 2.89 | 2.92 / 2.83 |
| cross-field correlation | +0.795 / +0.754 | +0.794 / +0.742 | +0.796 / +0.766 |
| bright-star RMS, BLG41 | 4.80 / 4.83 | 4.74 / 4.83 | 4.74 / 4.83 |
| bright-star RMS, BLG01 | 4.75 / 4.95 | 4.78 / 4.95 | 4.73 / 4.96 |
| target RMS, BLG41 | 23.77 / 19.67 | 23.80 / 19.72 | 23.78 / 19.66 |
| target RMS, BLG01 | 21.51 / 20.45 | 21.47 / 20.40 | 21.53 / 20.34 |
| target PM, BLG41 | −1.966 / −7.019 | −1.879 / −7.020 | −1.909 / −7.116 |
| target PM, BLG01 | −1.514 / −7.073 | −1.684 / −7.048 | −1.595 / −7.105 |

All in mas or mas/yr, X / Y. The cross-field correlation, which is the
quantity that measures the shared systematic, is unchanged to the third decimal,
the bright-star RMS moves by at most 0.06 mas and the target's proper motion by
at most 0.17 mas/yr, or 0.75 of its error. Six colour bins already capture all
the colour dependence there is; the six-bin solution is kept.

The per-season SysRem test is the most informative. SysRem represents residuals as a per-source
coefficient times a per-epoch mode, which is exactly the structure the slopes
appear to have; giving it its own two components in each season — ten times the
freedom — should have absorbed them. Instead the slopes grew and every other
metric degraded by 5 to 11%, because each season's correction is constrained by
about 1900 epochs instead of 19 000 and only 88 to 90% of epochs received any
correction at all. The slopes are therefore **not** a low-rank common mode
within seasons, consistent with their being per-star and spatially uncorrelated.

Since the effect reproduces between two cut-outs of the same stars on the same
nights, it is in the images rather than in the reduction. It is the likely origin
of the χ²/DoF of about 7 in section 7, of the factor 2.6 by which the internal
proper-motion errors are optimistic, and of the floor on what this dataset can
deliver.

### 9. Overlap of the two fields

Stars are identified through the OGLE catalogue; matching the two source lists
directly is ambiguous because a close pair may be split in one field and merged
in the other.

| | BLG41 | BLG01 |
|---|---|---|
| sources | 594 | 621 |
| matched to an OGLE entry | 551 (93%) | 497 (80%) |
| **distinct** OGLE stars behind them | **454** | **394** |

The gap between the last two rows is the blending, measured: 551 detections
correspond to only 454 stars on BLG41, so about 97 detections are duplicates.

| OGLE I | BLG41 | BLG01 | in both | distinct | found in both |
|---|---|---|---|---|---|
| < 15 | 28 | 22 | 22 | 28 | 78.6% |
| 15 – 16 | 54 | 42 | 42 | 54 | 77.8% |
| 16 – 17 | 172 | 149 | 146 | 175 | 83.4% |
| 17 – 18 | 135 | 119 | 117 | 137 | 85.4% |
| 18 – 19 | 55 | 44 | 41 | 58 | 70.7% |
| **all brighter than 19** | **444** | **376** | **368** | **452** | **81.4%** |

About four stars in five are detected in both fields, flat from the brightest
down to I = 18. **There is nothing fainter than I = 19 to compare**: our source
lists span 12.87 to 18.98 in both fields, while OGLE itself reaches 21.40 and
holds 2660 entries fainter than 19. That is the KMT detection limit in this
crowding, not a catalogue cut.

### 10. Files

| file | content |
|---|---|
| `report_v8_RMS.png` | residual RMS against OGLE I, per axis and field |
| `report_v8_target_beforeafter.png` | the target before and after proper-motion detrending |
| `report_v8_pm_vs_gaia.png` | our absolute proper motion against Gaia's, after the tie |
| `report_v8_chi2.png` | χ² distribution of the monthly binned residuals, target marked |
| `report_v8_radec_<field>_<tag>.png` | RA and Dec residuals against time and against each other, colour-coded by time; bins with fewer than 10 epochs and epochs at sec z > 1.3 are dropped |
| `report_v8_motion_<field>_<tag>.png` | 2x2 motion figures, proper motion retained and removed; same binning cuts as the RA/Dec figures |
| `source_motion_v8_<field>[_<tag>].csv` | per-epoch positions for the target and each comparison star |

CSV columns are `JD, X_mas, Y_mas, errX_decade, errY_decade, errX_season,
errY_season, season`, with positions relative to that star's own mean and
**proper motion not removed**. Headers carry the star's magnitude, pixel
position, sky position and absolute proper motion.

In the motion figures the columns are the two axes, **X left and Y right**; the
top row keeps the proper motion and the bottom row removes it. The legend quotes
the motion along the pixel axis and the heading quotes it on the sky. **Pixel X
runs opposite to RA**, so those two carry opposite signs in X.
