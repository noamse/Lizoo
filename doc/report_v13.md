# KMT-2026-BLG-0521 / OGLE-2026-BLG-0058

## Astrometry of the microlensing target — summary

Prepared 2026-09-14. Solution: `~/KMTdata/Results/v13/IFfinal_<field>.mat`
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
annual term, a pixel-phase term, and six SysRem components computed once over
the whole decade (v8, the previous reference, used two; section 8 gives the
series). Every star is freely fitted — nothing is held at a catalogue
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
| calibration stars brighter than I = 17 | **4.40 / 4.68 mas** | **4.48 / 4.87 mas** |
| the target, I = 18.13 | **19.53 / 19.28 mas** | **20.19 / 19.02 mas** |

The frame is therefore good to under **5 mas** on well-measured stars, and the
target — three magnitudes fainter and blended — is measured to about **20 mas**
per epoch, consistent with the running median of the residual RMS at its
magnitude.

**Proper motion of the target**, absolute, after the Gaia tie:

| | μ_α cos δ | μ_δ |
|---|---|---|
| BLG41 | **−2.271 ± 0.221** | **−6.931 ± 0.363** |
| BLG01 | **−1.960 ± 0.235** | **−7.023 ± 0.447** |
| difference between the fields | −0.311 ± 0.323 | +0.092 ± 0.576 |

mas/yr. The two independent fields agree within their errors, at 0.96 and 0.16
sigma. The quoted error is internal, from the scatter of the ten season means
about the fitted line; the external check against Gaia gives about
**1.0 / 0.60 mas/yr** per star, and that is the figure to use in practice.
That external scatter is magnitude dependent: 0.55 / 0.35 mas/yr for I < 17,
1.8 / 1.0 at 17–18 and 2.5–4.5 / 1.5–2.2 at 18–19, where both KMT blending and
Gaia's own errors grow. At the target's magnitude the disagreement with Gaia to
expect is therefore 1.5–2 mas/yr in declination, not 0.6.

**The blend.** The target has a companion of equal brightness (`m18_d04`,
I = 18.14, Gaia G = 19.69) at 4.5 pixels = 1.8″, itself paired with a G = 20.3
star 1.0″ further out for which Gaia solves no proper motion; under the 3.1″
KMT PSF the three are two blobs, and the calibration set's isolation cut, which
excludes companions brighter than I = 18 within 2″, lets the pair through by
0.14 mag. Gaia gives the target **−2.48 ± 0.47 / −9.13 ± 0.31 mas/yr**, i.e. our
declination motion is 2.1–2.2 mas/yr less negative than Gaia's in both fields
(right ascension agrees to 0.2 and 0.5 mas/yr), and
the companion +1.74 / +0.27, 9.4 mas/yr away from the target in declination.
Three tests show that this offset is **not** the companion pulling our centroid:
the measured separation of the two blobs is 91% of Gaia's, with the shortfall
in the direction of the G = 20.3 star rather than of the target, so the mutual
pull is at most a few per cent and could bias the proper motion by at most
~0.5 mas/yr; the position of the target moves with seeing by 5 mas per pixel
of FWHM, the same as isolated stars, and the seeing has no trend over the
decade; and the relative motion of the two blobs measured directly from the raw
positions, without any frame solution, is +2–3 / +0.4–0.8 mas/yr against Gaia's
+4.2 / +9.4. The 2.1 mas/yr is inside the 1.5–2 mas/yr KMT-minus-Gaia scatter
for stars of this magnitude (26–29% of the stars within 0.5 mag of the target
are further from Gaia in declination than it is), and Gaia's own five-parameter solution for the
target carries 1.0 mas of astrometric excess noise in a 2″ group of three
sources. Which of the two declination values is right cannot be decided from
these data; what the blend does do for certain is add to the per-epoch scatter
of the target, and dilute any astrometric excursion of the source by the flux
fraction of the companion inside the fitting aperture.

**No astrometric anomaly is detected at the target.** Its χ² on monthly binned
residuals is 1272 and 1547 for 180 degrees of freedom, against field medians of
χ²/DoF of 5.1 and 5.2 — the 75th and 82nd percentile of the population and the
75th and 78th among stars within 0.4 mag of it: inside the body of the
distribution, not in the tail. Whatever astrometric excursion the microlensing event produces
is smaller than the scatter of an ordinary star of this brightness.

@@FIG report_v13_RMS.png | Residual RMS against the OGLE catalogue magnitude, per axis and per field, after position and proper motion have been removed. Blue: all 260 / 235 calibration stars, every one freely fitted. Green: the stars the Gaia tie is fitted on, RUWE below 1.4 and 15 < I < 17. Red: the target. Black: running median. The dashed lines mark the 14 and 19 mag selection limits. Compared with a set that has not been iterated, the 2 sigma cut removes the scattered points that previously sat well above the median relation.

@@FIG report_v13_target_beforeafter.png | The target before and after proper-motion detrending, in both fields. Columns are the fields; the top row keeps the proper motion and the bottom row removes it. Blue circles are the offset in right ascension, red squares in declination, binned in one sidereal month with the standard error within the bin. The proper motion is plainly visible as the declination trend in the top row and is absent from the bottom row; the residual RMS quoted in the lower panels is per bin.

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
astrometric error. Median residual **22.0 mas (BLG41)** and **21.6 (BLG01)**,
with scales of 2.5195 and 2.5289 pix/arcsec against the 2.500 expected at
0.4″/pix.

**The gauge.** In the proper motions the gauge freedom is a constant plus a
linear gradient across the field — six parameters. Removing only the constant
leaves a shear of about −0.036 mas/yr per pixel, reproduced at nearly the same
value in both independently reduced fields. With the full gauge removed, the
scatter of our proper motions about Gaia is **0.967 / 0.593 mas/yr** (BLG41) and
**1.034 / 0.606** (BLG01); removing only the constant would leave 3 to 4 times
that. The gauge constant removed is −0.097 / +5.129 and +2.076 / +4.328 mas/yr.

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
- **Annual term** — a quartic in the phase of the calendar year, 10 parameters
  per colour bin, one curve shared by every star of the bin (the same six
  equal-population colour bins as the DCR term). Not parallax, which is
  per-source and switched off.
- **Pixel phase** — a quintic in the sub-pixel phase, 10 parameters, global,
  recoverable only because the pipeline stored the per-epoch registration shifts.
- **SysRem** — 6 components over the whole decade, reused rather than refitted
  per season, with point-wise rejection of bad cells. Refitting per season was
  tested and is worse, and the number of components was chosen from a series
  of 0 to 12: see section 8.

Timings and convergence:

| | BLG41 | BLG01 |
|---|---|---|
| step 1, decade SysRem | 31.7 min | 29.9 min |
| step 2, final solution | 54.8 min | 53.3 min |
| convergence, last 3 iterations | 0.0098 mas | 0.0018 mas |

No per-season fitting enters the solution; seasons are used only for reporting.

### 6. Comparison stars

Six, chosen to be the same physical objects in both fields — three near the
target's magnitude and three brighter. All six survive the outlier cut and sit
at identical positions in the two fields:

| tag | I | dist [pix] | decade RMS ΔX / ΔY, BLG41 | BLG01 |
|---|---|---|---|---|
| `m18_d04` | 18.14 | 4.5 | 25.72 / 24.96 | 25.43 / 26.21 |
| `m18_d17` | 18.09 | 17.3 | 58.55 / 23.99 | 93.54 / 33.84 |
| `m18_d48` | 18.05 | 48.4 | 14.22 / 15.23 | 13.22 / 14.94 |
| `m16_d24` | 16.68 | 24.4 | **4.04 / 4.17** | **3.99 / 4.25** |
| `m16_d30` | 16.63 | 30.4 | 5.26 / 5.33 | 5.17 / 5.31 |
| `m16_d31` | 16.75 | 30.6 | 8.26 / 7.58 | 7.81 / 7.73 |

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
| **target χ²** | **1272** (DoF 180) | **1547** (DoF 180) |
| target χ²/DoF | 7.1 | 8.6 |
| field median χ²/DoF | 5.1 | 5.2 |
| target percentile, all stars | 75% | 82% |
| target percentile, stars within 0.4 mag | 75% | 78% |

The target is unremarkable. The field median of χ²/DoF near 5 is itself
informative: the standard error of the mean within a month **understates the
true bin uncertainty by about √5 ≈ 2.3** (it was 2.6 with two SysRem
components), because the noise within a month is
correlated rather than independent. The effect is multiplicative and nearly
independent of magnitude, and the faintest stars come closest to unity, since
for them genuine photon noise — which *is* independent between exposures —
dominates.

**What the error bars in the figures are.** Each plotted point is the mean
position of the frames in one sidereal-month bin. Its *standard error* (s.e.)
is σ/√N, with σ the scatter of the individual frames inside the bin and N their
number — the uncertainty of the mean, not the spread of the frames. For the
target a typical bin has σ ≈ 21–23 mas and N ≈ 75, hence an s.e. of about
2.5 mas, which is what the bar would be if the frames were independent. They
are not, so the bars drawn are the s.e. multiplied by an empirical factor,
√(χ²/DoF) of that star's bin means about zero, computed per star and per axis
and quoted in each panel title (1.9–3.4 for the target, 1.5–2 for the
best-measured comparison stars): the standard error the frame count suggests,
inflated to match the scatter the bins actually show.

### 8. A residual systematic that survives every correction

Fitting a straight line to each star's residuals **within each season
separately** leaves slopes that are real and unexplained.

| | BLG41 | BLG01 |
|---|---|---|
| median per-star season-slope rms | 2.54 / 2.68 mas | 2.45 / 2.66 mas |
| expected from white noise | 0.91 mas | 0.84 mas |
| **excess over white noise** | **2.8x** | **2.9x** |

They are **not noise**. Matching 181 stars between the fields through OGLE and
correlating each star's ten-season slope pattern gives a median correlation of
**+0.729** in one axis and **+0.664** in the other, with **96% of stars
positive** (with two SysRem components these were +0.795 and +0.754). Two
independent reductions, independent frames, essentially every star showing the
same season-by-season pattern: about 70% of the slope variance is shared
systematic rather than measurement noise. The target is entirely typical in
this respect, at the 50th to 75th percentile among stars within 0.4 mag of it.

The cause has not been identified. Eleven candidates were tested; the first
nine were tested on the v8 solution, with two SysRem components, and are quoted
with its numbers; the tenth is the SysRem series that led to v13, the eleventh
was tested on v13.

| candidate | test | result |
|---|---|---|
| DCR / colour | per-season slope against colour | correlation below 0.11 in size |
| DCR colour resolution | full refit with 12 and with 20 colour bins instead of 6 | slopes and cross-field correlation unchanged, see below |
| calibrator contamination | second clipping iteration on the season scale, refit without the 21 + 22 stars with 2σ-high season-offset χ² | per-star RMS unchanged to within 1%, see below |
| season-edge epochs | refit without the first and last 14 nights of every season (5% of the epochs) | RMS on the common epochs unchanged to within 1%; cross-field correlation +0.770 / +0.734, see below |
| high-airmass epochs | refit without sec z > 1.3 (30% of the epochs) | fit 3–4% better for I > 17 on the common epochs; slope excess and cross-field correlation barely lower, see below |
| SysRem rank | 0, 2, 3, 4, 6, 8, 12 decade-wide components | half the slope signal is a rank-2 common mode; the rest erodes by ~4% per component and is not low-rank, see below |
| annual term | v13 refitted with the annual term switched off in both steps and both passes | nothing changes: RMS within ±1%, slopes within ±2%, cross-field +0.704 / +0.670, see below |
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

The season-scale clipping asks whether the frame is being pulled by
calibrators that themselves carry coherent season offsets — the decade RMS used
for the v8 cut averages such offsets away (a 16.75-mag star among the comparison
set has its X residual at +15 mas for two seasons and −15 mas for one, with a
perfectly normal decade RMS). Applying a second 2σ cut, in log space per
magnitude bin, on each star's season-offset χ²/DoF — Σ(season mean / s.e.)²
over seasons and axes — removes 21 stars in BLG41 and 22 in BLG01, two-thirds of
them at 16–17 mag; the target sits at z = +0.9 / +0.6, inside the population.
Refitting on the 239 + 213 survivors changes nothing on the stars common to
both solutions:

| common stars | BLG41 | BLG01 |
|---|---|---|
| RMS, I < 17 (117 / 105 stars) | 4.66 / 4.77 → 4.65 / 4.80 | 4.61 / 4.95 → 4.60 / 4.94 |
| RMS, 17–18 | +1.0% / +0.5% | −1.5% / +1.5% |
| RMS, 18–19 | +0.1% / +2.4% | −0.9% / +1.1% |
| Gaia scatter, RUWE-clean (181 / 146) | 0.967 / 0.558 → 0.934 / 0.572 | 1.067 / 0.629 → 0.982 / 0.661 |
| target PM | −1.966 / −7.019 → −2.003 / −6.997 | −1.514 / −7.073 → −1.770 / −7.024 |

Per-star RMS moves by under 1% with half the stars going each way; the Gaia
scatter improves in X and worsens in Y; the target's proper motion moves by
0.04 and 0.26 mas/yr (0.1σ and 1.1σ), bringing the two fields to within 0.23
mas/yr in X. The apparent 3% gain in the quoted bright-star RMS (4.64 / 4.72 and
4.58 / 4.94) is entirely the removal of the worst stars from the median, not an
improvement of the solution. The per-epoch frame with some 240 stars is
insensitive to these 8–9%, so calibrator contamination is not the origin of
the slopes, and the calibration set is unchanged.

The two epoch cuts ask whether the slopes are made at the season edges or at
high airmass. They are nearly independent selections: the first and last 14
nights of a season hold only 25–50 epochs each at a median sec z of 1.33–1.40,
but frames at sec z > 1.3 are taken at the start and end of every night all
season long, so only 12% of them fall in the edge trim. Both were refitted with
the v8 stars and fit otherwise unchanged (the truncated 2026 season keeps its
final, mid-season nights). Compared on exactly the surviving epochs, so that
dropping the noisiest frames is not mistaken for a better solution:

| on the common epochs | edge trim, BLG41 / BLG01 | sec z ≤ 1.3, BLG41 / BLG01 |
|---|---|---|
| epochs kept | 96% / 96% | 78% / 77% |
| RMS change, I < 17 | +0.6 / +0.1%, −0.1 / −0.4% | −0.6 / −0.2%, 0.0 / +0.1% |
| RMS change, 17–18 | 0.0 / −0.4%, +0.4 / +0.1% | −0.3 / −3.0%, −3.9 / −3.9% |
| RMS change, 18–19 | +0.6 / −0.5%, +1.1 / +0.2% | −0.5 / −1.0%, −1.0 / −1.0% |
| RMS change, target | −0.1 / −0.5%, +2.1 / +1.1% | −2.6 / 0.0%, −3.9 / +2.0% |
| target PM | −1.926 / −6.988, −1.583 / −7.056 | −1.530 / −7.195, −1.583 / −7.187 |
| Gaia scatter | 0.985 / 0.598, 1.039 / 0.648 | 0.995 / 0.596, 1.027 / 0.698 |
| slope rms | 2.68 / 2.56, 2.69 / 2.67 | 2.73 / 2.79, 2.81 / 2.72 |
| excess over white noise | 2.9 / 3.2 | 2.6 / 2.9 |
| cross-field correlation | +0.770 / +0.734 | +0.770 / +0.730 |

The edge trim is a null: on the common epochs nothing changes, and the 7%
lower slope rms is only the shorter season. The airmass cut is a small, real
improvement of the fit for stars fainter than I = 17 and for the target's X
scatter — the differential-refraction and high-airmass PSF model evidently
does worst on faint stars — and it moves the target's proper motion in BLG41
by 0.44 mas/yr in X (1.3σ), bringing the two fields to within 0.05 mas/yr;
that is one number moving within its error, at the cost of 30% of the data.
Either way the slope excess and the cross-field correlation survive with
three-quarters of the epochs at sec z < 1.3, so the shared systematic is made
neither at the season edges nor at high airmass.

**The SysRem series** is what changed the reference from v8 to v13. With the
v8 stars and fit, the number of decade-wide SysRem components was varied:

| components | bright RMS, BLG41 / BLG01 | target RMS | Gaia scatter | target PM, BLG41 | target PM, BLG01 | slope rms | excess | cross-field |
|---|---|---|---|---|---|---|---|---|
| 0 | 6.14/5.99, 6.48/7.15 | 30.1/25.7, 30.9/26.1 | 0.995/0.613, 1.022/0.703 | −2.195 / −6.971 | −1.844 / −7.001 | 5.12/5.72, 5.44/6.11 | 3.9 / 4.3 | +0.882 / +0.900 |
| 2 (v8) | 4.80/4.83, 4.75/4.95 | 23.8/19.7, 21.5/20.5 | 0.994/0.608, 1.049/0.641 | −1.966 / −7.019 | −1.514 / −7.073 | 2.90/2.80, 2.93/2.90 | 2.9 / 3.3 | +0.795 / +0.754 |
| 3 | 4.54/4.76, 4.67/4.91 | 20.0/19.5, 20.4/20.5 | 1.002/0.602, 1.004/0.650 | −2.219 / −6.988 | −1.881 / −7.077 | 2.76/2.79, 3.00/2.89 | 3.0 / 3.4 | +0.762 / +0.731 |
| 4 | 4.47/4.75, 4.50/5.00 | 20.0/19.5, 20.2/19.8 | 0.986/0.598, 1.023/0.678 | −2.232 / −7.043 | −1.868 / −7.149 | 2.65/2.69, 2.64/2.92 | 2.9 / 3.1 | +0.723 / +0.731 |
| **6 (v13)** | 4.40/4.68, 4.48/4.87 | 19.5/19.3, 20.2/19.0 | **0.967/0.593, 1.034/0.606** | −2.271 / −6.931 | −1.960 / −7.023 | 2.54/2.68, 2.45/2.66 | 2.8 / 2.9 | +0.729 / +0.664 |
| 8 | 4.27/4.64, 4.36/4.60 | 19.5/18.2, 20.0/18.9 | 0.966/0.589, 1.037/0.665 | −2.291 / −7.010 | −1.916 / −7.089 | 2.32/2.39, 2.41/2.51 | 2.7 / 2.9 | +0.727 / +0.665 |
| 12 | 4.18/4.32, 4.17/4.44 | 19.3/18.2, 19.9/18.3 | 0.969/0.623, 1.023/0.671 | −2.189 / −7.007 | −1.982 / −7.090 | 2.23/2.15, 2.34/2.20 | 2.6 / 2.9 | +0.656 / +0.599 |

All in mas or mas/yr, X / Y. Three things follow. First, without SysRem the
season slopes double and their cross-field correlation rises to +0.88 / +0.90:
before any correction the season structure is almost entirely a shared common
mode, and the first two components remove half of it — but each further
component removes only ~4% more, so what remains is not low-rank and will not
be SysRem'd away. Second, the external accuracy against Gaia is best at six
components and drifts worse in declination beyond (BLG01 0.606 → 0.665 → 0.671),
the first sign of the components absorbing real per-star motion; the residual
RMS keeps falling with rank, but that is what free parameters do and is not
evidence of a better solution. Third, the target's right-ascension motion is
rank-dependent in a specific way: with 0, 3, 4, 6, 8 and 12 components it sits
at −2.19 to −2.29 (BLG41) and −1.84 to −1.98 (BLG01) mas/yr, and only the
two-component v8 solution gives −1.97 / −1.51 — the second component carried
something the target's residual projected onto. Six components were therefore
adopted: best external accuracy, target motion at its plateau, bright-star RMS
8% / 6% and target X scatter 18% below v8. Declination is unaffected by the
rank, −6.93 to −7.15 throughout.

**The annual term** — a quartic in the phase of the calendar year, fitted in
pass 1 and subtracted before pass 2 — was switched off entirely (v14). Nothing
moves: bright-star RMS 4.38 / 4.69 and 4.48 / 4.83 mas against 4.40 / 4.68 and
4.48 / 4.87, target RMS within 1%, target proper motion −2.211 / −6.932 and
−1.855 / −7.098 against −2.271 / −6.931 and −1.960 / −7.023 (at most 0.4σ),
slope rms 2.49 / 2.73 and 2.47 / 2.69 against 2.54 / 2.68 and 2.45 / 2.66,
cross-field correlation +0.704 / +0.670 against +0.729 / +0.664. The only
consistent direction is the external scatter against Gaia, slightly worse
without the term (1.005 / 0.598 and 1.029 / 0.648 against 0.967 / 0.593 and
1.034 / 0.606), so it is kept. Whatever the annual term describes is already
absorbed by SysRem and the refraction model, and the season structure is not a
mis-modelled annual effect.

The per-season SysRem test is the most informative on the *nature* of the
remainder. SysRem represents residuals as a per-source
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
| `report_v13_RMS.png` | residual RMS against OGLE I, per axis and field |
| `report_v13_target_beforeafter.png` | the target before and after proper-motion detrending |
| `report_v13_pm_vs_gaia.png` | our absolute proper motion against Gaia's, after the tie |
| `report_v13_chi2.png` | χ² distribution of the monthly binned residuals, target marked |
| `report_v13_radec_<field>_<tag>.png` | RA and Dec residuals against time and against each other, colour-coded by time; bins with fewer than 10 epochs and epochs at sec z > 1.3 are dropped |
| `report_v13_motion_<field>_<tag>.png` | 2x2 motion figures, proper motion retained and removed; same binning cuts as the RA/Dec figures |
| `source_motion_v13_<field>[_<tag>].csv` | per-epoch positions for the target and each comparison star |

CSV columns are `JD, X_mas, Y_mas, errX_decade, errY_decade, errX_season,
errY_season, season`, with positions relative to that star's own mean and
**proper motion not removed**. Headers carry the star's magnitude, pixel
position, sky position and absolute proper motion.

In the motion figures the columns are the two axes, **X left and Y right**; the
top row keeps the proper motion and the bottom row removes it. The legend quotes
the motion along the pixel axis and the heading quotes it on the sky. **Pixel X
runs opposite to RA**, so those two carry opposite signs in X.
