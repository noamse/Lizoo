# KMT-2026-BLG-0521 / OGLE-2026-BLG-0058

## Astrometry of the microlensing target — summary

Prepared 2026-09-10. Solution: `~/KMTdata/Results/v6/IFfinal_<field>.mat`
(variable `IFsys`), tie in `Tie_<field>.mat`.

### 1. The data and the algorithm

**Data.** Two KMTNet cut-outs of the same field, **BLG41** and **BLG01**,
300 x 300 pixels at 0.4 arcsec/pix, in I band, spanning 2016 to 2026 in ten
observing seasons (2020 is a short stub and is dropped). **19 133** and
**19 981** exposures respectively, taken alternately on the same nights, a
median of 12.2 minutes apart. The target sits at the field centre at
I_OGLE = 18.13, RA 17:52:38.09, Dec −31:47:36.1.

**Calibration set.** The analysis is restricted to stars selected on magnitude
and isolation alone: **14 < I < 19** with **no companion brighter than I = 18
within 2 arcsec**, the companion list taken from the OGLE catalogue, which
resolves about six times more sources than KMT detects. This leaves **287**
stars on BLG41 and **253** on BLG01. The target passes the selection on its own
merits. Everything else is dropped at the input, so all measurement, detrending
and quality cuts act on this set alone.

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
| calibration stars brighter than I = 17 | **5.15 / 5.06 mas** | **4.93 / 5.12 mas** |
| the target, I = 18.13 | **23.75 / 19.63 mas** | **21.36 / 20.37 mas** |

The frame is therefore good to about **5 mas** on well-measured stars, and the
target — three magnitudes fainter and blended — is measured to about **20 mas**
per epoch, consistent with the running median of the residual RMS at its
magnitude.

**Proper motion of the target**, absolute, after the Gaia tie:

| | μ_α cos δ | μ_δ |
|---|---|---|
| BLG41 | **−1.937 ± 0.335** | **−7.140 ± 0.360** |
| BLG01 | **−1.658 ± 0.256** | **−7.089 ± 0.377** |
| difference between the fields | −0.279 ± 0.421 | −0.052 ± 0.521 |

mas/yr. The two independent fields agree within their errors, at 0.66 and 0.10
sigma. The quoted error is internal, from the scatter of the ten season means
about the fitted line; the external check against Gaia gives about
**1.0 / 0.65 mas/yr** per star, and that is the figure to use in practice.

**No astrometric anomaly is detected at the target.** Its χ² on monthly binned
residuals is 1574 and 1313 for 180 degrees of freedom, against field medians of
the same quantity — the 63rd and 53rd percentile of the population, and the 68th
and 62nd among stars within 0.4 mag of it. Whatever astrometric excursion the
microlensing event produces is smaller than the scatter of an ordinary star of
this brightness.

@@FIG report_v6_RMS.png | Residual RMS against the OGLE catalogue magnitude, per axis and per field, after position and proper motion have been removed. Blue: all 287 / 253 calibration stars, every one freely fitted. Green: the stars the Gaia tie is fitted on, RUWE below 1.4 and 15 < I < 17. Red: the target. Black: running median. The dashed lines mark the 14 and 19 mag selection limits.

@@FIG report_v6_target_beforeafter.png | The target before and after proper-motion detrending, in both fields. Columns are the fields; the top row keeps the proper motion and the bottom row removes it. Blue circles are the offset in right ascension, red squares in declination, binned in one sidereal month with the standard error within the bin. The proper motion is plainly visible as the declination trend in the top row and is absent from the bottom row; the residual RMS quoted in the lower panels is per bin.

---

## More details

### 3. Object selection, step by step

| step | BLG41 | BLG01 |
|---|---|---|
| raw sources in the MSc file | 1177 | 1187 |
| survive the standard quality cuts | 594 | 621 |
| have an OGLE I magnitude | 540 | 478 |
| in 14 < I < 19 | 529 | 470 |
| ... and isolated: no I < 18 companion within 2.0 arcsec | **287** | **253** |

The isolation cut is the dominant one, removing 46% of the magnitude-selected
sources on both fields. For the tie only, Gaia is required: 273 of 287 and 225
of 253 sources have a counterpart, of which 220 and 177 pass RUWE < 1.4, and
118 and 90 also fall in 15 < I < 17.

### 4. The Gaia tie

**The map.** Gaia standard coordinates are fitted directly onto cut-out pixels
from the tie stars — not carried through OGLE, which would inherit its ~110 mas
astrometric error. Median residual **23.4 mas (BLG41)** and **22.1 (BLG01)**,
with scales of 2.5195 and 2.5290 pix/arcsec against the 2.500 expected at
0.4″/pix. Restricting the tie set to 15 < I < 17 improves the map by 35 to 40%
over using all RUWE-clean stars, since bright stars have better positions.

**The gauge.** In the proper motions the gauge freedom is a constant plus a
linear gradient across the field — six parameters. Removing only the constant
leaves a shear of **−0.037 mas/yr per pixel**, visible as `corr(Δμ_α, Y) = −0.75`
and reproduced at nearly the same value in both independently reduced fields.
Removing the full gauge:

| scatter of our proper motions about Gaia | constant only | **full gauge** |
|---|---|---|
| BLG41 | 4.34 / 3.57 mas/yr | **1.02 / 0.64** |
| BLG01 | 3.48 / 3.20 mas/yr | **1.06 / 0.67** |

**What the magnitude cut on the tie set is worth.** Evaluated on the same test
set — all RUWE-clean stars, so neither tie is flattered by its own training
sample — it changes the proper-motion agreement by a few per cent, inconsistent
in sign. The narrower set is better measured per star, 0.57 against 0.98 mas/yr,
and determines the gauge slightly more precisely, 0.053 against 0.066 mas/yr,
but halving the star count nearly cancels that and the gauge is **under 1% of
the total variance**. The cut is kept because it improves the positional map,
not because it improves the motions.

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
  shear (2). 574 equations for 6 unknowns each epoch. Translation dominates at
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
  per season, with point-wise rejection of bad cells.

Timings and convergence:

| | BLG41 | BLG01 |
|---|---|---|
| step 1, decade SysRem | 19.8 min | 16.7 min |
| step 2, final solution | 35.6 min | 31.3 min |
| convergence, last 3 iterations | 0.0028 mas | 0.0007 mas |

No per-season fitting enters the solution; seasons are used only for reporting.

### 6. Comparison stars

Six, chosen to be the same physical objects in both fields — three near the
target's magnitude and three brighter:

| tag | I | dist [pix] | decade RMS ΔX / ΔY, BLG41 | BLG01 |
|---|---|---|---|---|
| `m18_d04` | 18.14 | 4.5 | 29.03 / 27.91 | 24.84 / 26.86 |
| `m18_d17` | 18.09 | 17.3 | 51.68 / 27.85 | 85.47 / 28.47 |
| `m18_d48` | 18.05 | 48.4 | 17.92 / 16.76 | 14.58 / 15.57 |
| `m16_d24` | 16.68 | 24.4 | **4.19 / 4.27** | **4.17 / 4.38** |
| `m16_d30` | 16.63 | 30.4 | 6.91 / 5.85 | 6.43 / 5.55 |
| `m16_d31` | 16.75 | 30.6 | 8.48 / 7.79 | 11.20 / 7.71 |

mas. The bright trio reproduces between the fields almost exactly — `m16_d24`
gives −3.31 / −7.81 on BLG41 and −3.31 / −7.89 on BLG01 — while the faint trio
does not: `m18_d17` gives −6.03 against −10.05 mas/yr and reaches 52 and 85 mas
decade RMS on only ~2400 usable epochs. It is a poor star.

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
| stars measured | 287 | 253 |
| median DoF | 182 | 182 |
| **target χ²** | **1574** (DoF 180) | **1313** (DoF 180) |
| target χ²/DoF | 8.7 | 7.3 |
| field median χ²/DoF | 7.2 | 7.1 |
| target percentile, all stars | 63% | 53% |
| target percentile, stars within 0.4 mag | 68% | 62% |

The target is unremarkable. The field median of χ²/DoF ≈ 7 is itself
informative: the standard error of the mean within a month **understates the
true bin uncertainty by about √7 ≈ 2.6**, because the noise within a month is
correlated rather than independent. The effect is multiplicative and nearly
independent of magnitude — the Spearman correlation of χ²/DoF with I is −0.12
and −0.17 — and the faintest stars come closest to unity, since for them genuine
photon noise, which *is* independent between exposures, dominates.

### 8. Overlap of the two fields

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

### 9. Files

| file | content |
|---|---|
| `report_v6_RMS.png` | residual RMS against OGLE I, per axis and field |
| `report_v6_target_beforeafter.png` | the target before and after proper-motion detrending |
| `report_v6_pm_vs_gaia.png` | our absolute proper motion against Gaia's, after the tie |
| `report_v6_chi2.png` | χ² distribution of the monthly binned residuals, target marked |
| `report_v6_radec_<field>_<tag>.png` | RA and Dec residuals against time and against each other, colour-coded by time |
| `report_v6_motion_<field>_<tag>.png` | 2x2 motion figures, proper motion retained and removed |
| `source_motion_v6_<field>[_<tag>].csv` | per-epoch positions for the target and each comparison star |

CSV columns are `JD, X_mas, Y_mas, errX_decade, errY_decade, errX_season,
errY_season, season`, with positions relative to that star's own mean and
**proper motion not removed**. Headers carry the star's magnitude, pixel
position, sky position and absolute proper motion.

In the motion figures the columns are the two axes, **X left and Y right**; the
top row keeps the proper motion and the bottom row removes it. The legend quotes
the motion along the pixel axis and the heading quotes it on the sky. **Pixel X
runs opposite to RA**, so those two carry opposite signs in X.
