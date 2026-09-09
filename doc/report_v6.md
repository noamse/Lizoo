# KMT-2026-BLG-0521 / OGLE-2026-BLG-0058

## The v6 reduction — wide calibration set, freely fitted, Gaia applied afterwards

Prepared 2026-09-09. Solution: `~/KMTdata/Results/v6/IFfinal_<field>.mat`
(variable `IFsys`), with the tie in `Tie_<field>.mat`.

---

### 1. A short description

The analysis is restricted to a calibration set defined by **magnitude and
isolation alone**: 14 < I < 19 with no companion brighter than I = 18 within
2 arcsec. Every one of those stars is freely fitted and every one takes part in
determining the per-epoch frame at every step. **No Gaia quality cut is applied
at selection**; RUWE enters only afterwards, when the finished solution is tied
to Gaia. The target passes the selection on its own merits and needs no special
handling.

This is v5 with the magnitude window widened from 14–17 to 14–19 and the RUWE
cut moved to the end. Both changes help: **v6 measures the target better than
any previous reduction, including the reference one.**

---

### 2. Object selection, step by step

| step | BLG41 | BLG01 |
|---|---|---|
| raw sources in the MSc file | 1177 | 1187 |
| survive the standard quality cuts | 594 | 621 |
| have an OGLE I magnitude | 540 | 478 |
| in 14 < I < 19 | 529 | 470 |
| ... and isolated: no I < 18 companion within 2.0 arcsec | **287** | **253** |
| **= analysis set actually fitted** | **287** | **253** |

The isolation cut is now the dominant one, removing 46% of the magnitude-selected
sources on both fields. The set is more than twice v5's 130 / 100.

**The target is a member.** At I = 18.13 it lies inside the window and passes
the isolation test, so unlike v4 and v5 it is not added by hand and no
"passenger" stars are needed.

For the tie only, Gaia is required: **273 of 287** and **225 of 253** sources
have a Gaia counterpart, and of those **220** and **177** pass RUWE < 1.4.

---

### 3. Gaia, applied only at the end

**The map.** Gaia standard coordinates are fitted directly onto cut-out pixels
from the RUWE-clean stars — not carried through OGLE, which would inherit its
~110 mas astrometric error. Median residual **35.1 mas (BLG41)** and 37.1
(BLG01), with scales of 2.5198 and 2.5291 pix/arcsec against the 2.500 expected
at 0.4″/pix.

**The gauge, and a correction to all earlier work.** A free fit determines
astrometry only up to a gauge: applying a global affine to all positions and
compensating in the per-epoch transformations changes no residual. In the proper
motions this freedom is **a constant plus a linear gradient across the field** —
six parameters, not two.

Every tie made in this project before v6 removed **only the constant**. The
residual shows up plainly in the cross terms: `corr(Δμ_α, Y) = −0.75` on BLG41,
a **shear of −0.037 mas/yr per pixel**, and −0.035 on BLG01 — nearly the same
value in two independently reduced fields, so it is a property of the method,
not noise. Removing the full six-parameter gauge:

| scatter of our proper motions about Gaia | constant only | **full gauge** |
|---|---|---|
| BLG41 | 4.34 / 3.57 mas/yr | **0.98 / 0.66** |
| BLG01 | 3.48 / 3.20 mas/yr | **1.00 / 0.70** |

A factor of three to four. The gauge removed here is a constant of
+0.015 / +5.170 (BLG41) and +2.275 / +4.329 (BLG01) mas/yr, plus the gradient
above. **The absolute proper motions in the earlier reports are affected by
this and should be read as superseded.**

---

### 4. The chain, and every effect fitted in it

Two steps, both on the 287 / 253 source set:

1. **Full-decade SysRem correction** — 2 components on the astrometric
   residuals, computed once and saved.
2. **Final solution**, reusing that correction.

| | BLG41 | BLG01 |
|---|---|---|
| step 1 | 19.8 min | 16.7 min |
| step 2 | 35.6 min | 31.3 min |
| convergence over stars brighter than 17, last 3 iterations | 0.0028 mas | 0.0007 mas |

**Each step runs two passes internally.** `IFsys.AnnualEffect` and
`IFsys.PixPhase` read `0` in the saved solutions, which does **not** mean those
effects were ignored: they are fitted and subtracted from the data in pass 1 and
forced off in pass 2, which works on data from which they have already been
removed.

| | pass 1 | between | pass 2 (final) |
|---|---|---|---|
| per-epoch affine | fitted | | fitted |
| per-source position and motion | fitted | | fitted |
| DCR / chromatic refraction | fitted | | fitted again |
| annual term | fitted and subtracted from the data | | forced off |
| pixel-phase term | fitted and subtracted from the data | | forced off |
| SysRem | | applied | |
| iterative reweighting | 15 iterations | | 6 iterations |

- **Per-epoch affine** — 6 parameters per exposure in `ParE`, stored as the
  deviation from the identity: translation (2), rotation (1), uniform scale (1),
  shear (2). Each epoch's six parameters come from 287 sources, so 574 equations
  for 6 unknowns, against 260 in v5 and 1188 in the reference reduction.
- **Per source** — position at `JD0` = 2019-06-01 and proper motion, 4
  parameters in `ParS`. All free; `ParSFixed` is empty.
- **DCR** — `[1, sin(pa)·secz, cos(pa)·secz]` plus six higher orders in
  `(secz·sin pa)ⁿ`, `(secz·cos pa)ⁿ` for n = 2, 3, 4: 9 per axis, 18 in total,
  fitted in 6 equal-population colour bins. Only the *differential* part is
  recoverable — refraction common to the field is exactly a translation and is
  absorbed by `t_x, t_y` every epoch.
- **Annual term** — a quartic in the phase of the calendar year, 10 parameters,
  fitted globally, one curve shared by every star. Not parallax, which is
  per-source and switched off.
- **Pixel phase** — a quintic in the sub-pixel phase, 10 parameters, global.
- **SysRem** — 2 components over the whole decade, reused rather than refitted
  per season, with point-wise rejection.
- **Weighting** — each source weighted by its own residual scatter, recomputed
  each iteration, with moving-median outlier rejection over 30 epochs.

No per-season fitting enters the solution. Seasons are used only for reporting.

---

### 5. Results

| | BLG41 ΔX / ΔY | BLG01 ΔX / ΔY |
|---|---|---|
| stars brighter than I = 17 | **5.15 / 5.06** | **4.93 / 5.12** |
| the target | **23.75 / 19.63** | **21.36 / 20.37** |

mas. Against the other reductions:

| target residual RMS | BLG41 | BLG01 |
|---|---|---|
| v3, reference, ~600 sources | 24.62 / 21.00 | 25.48 / 21.78 |
| v5, 130 / 100 sources, RUWE-cut | 25.61 / 22.27 | 28.45 / 23.15 |
| **v6, 287 / 253 sources** | **23.75 / 19.63** | **21.36 / 20.37** |

**v6 is the best on every axis in both fields** — 3.5% better than the reference
on BLG41 and 16% better on BLG01. The lesson from v5 is now clear: restricting
the set to bright isolated stars helps, but restricting it *too far* costs more
than it gains. 287 isolated stars beat both 130 isolated stars and ~600
unfiltered ones.

---

### 6. Proper motion of the target, with uncertainties

| | μ_α cos δ | μ_δ |
|---|---|---|
| BLG41 | **−2.017 ± 0.335** | **−7.163 ± 0.360** |
| BLG01 | **−1.781 ± 0.256** | **−7.136 ± 0.377** |
| difference between fields | −0.236 ± 0.421 | −0.028 ± 0.521 |

mas/yr. The two fields agree **within their errors**, at 0.56 and 0.05 sigma.

The quoted error is the **internal** one: the slope uncertainty from the scatter
of the ten season means about the fitted line. It is deliberately not computed
from the per-epoch count, which would give an absurdly small number — the
nightly noise is strongly correlated. Median internal uncertainty over all
sources is 0.127 / 0.114 (BLG41) and 0.123 / 0.112 (BLG01) mas/yr.

**That internal error is still optimistic by a factor of two to three.** The
external check — the scatter of our proper motions about Gaia's, 0.98 / 0.66 and
1.00 / 0.70 mas/yr per star — is larger than the internal estimate would predict
even after Gaia's own 0.15 mas/yr is allowed for. Use **~1.0 / 0.7 mas/yr** as
the realistic per-star uncertainty; the internal figure captures the random part
only and misses systematics that bias a star's slope while leaving its seasons
internally consistent.

Gaia's own value for the target is −2.479 / −9.129 mas/yr. Our μ_α cos δ agrees;
our μ_δ differs by about 2 mas/yr, which is 3 sigma on the external error. The
target is blended with a neighbour 4.5 pix away and its motion here is fitted
straight through a multi-year magnification event, so this is not a clean
comparison.

---

### 7. Comparison stars

Six, chosen to be the same physical objects in both fields — three near the
target's own magnitude and three brighter:

| tag | I | dist [pix] | decade RMS ΔX / ΔY, BLG41 | BLG01 |
|---|---|---|---|---|
| `m18_d04` | 18.14 | 4.5 | 29.03 / 27.91 | 24.84 / 26.86 |
| `m18_d17` | 18.09 | 17.3 | 51.68 / 27.85 | 85.47 / 28.47 |
| `m18_d48` | 18.05 | 48.4 | 17.92 / 16.76 | 14.58 / 15.57 |
| `m16_d24` | 16.68 | 24.4 | **4.19 / 4.27** | **4.17 / 4.38** |
| `m16_d30` | 16.63 | 30.4 | 6.91 / 5.85 | 6.43 / 5.55 |
| `m16_d31` | 16.75 | 30.6 | 8.48 / 7.79 | 11.20 / 7.71 |

mas. The bright trio reproduces between the fields almost exactly — `m16_d24`
gives −3.31 / −7.81 on BLG41 and −3.31 / −7.89 on BLG01 — while the I ≈ 18 trio
does not: `m18_d17` gives −6.03 against −10.05 mas/yr, and its decade RMS reaches
52 and 85 mas with only ~2400 usable epochs. It is a poor star and should not be
used as a control.

`m18_d04` is the neighbour blended with the target, 4.5 pix away under ~7 pix
seeing. It is a legitimate member of the calibration set — the target at
I = 18.13 is fainter than the I < 18 companion threshold, so it does not trigger
rejection — but its centroid tracks the target's brightening through the event,
and it is not an independent object.

---

### 8. Chi-squared of the monthly binned residuals

For every star the residuals are binned in sidereal months; each bin contributes
`(mean/standard error)²` in both axes, and `DoF = 2·N_bins − 4`. Values are not
divided by DoF.

| | BLG41 | BLG01 |
|---|---|---|
| stars measured | 287 | 253 |
| median DoF | 180 | 180 |
| **target χ²** | **1574** (DoF 180) | **1313** (DoF 180) |
| target χ²/DoF | 8.7 | 7.3 |
| **field median χ²/DoF** | **7.2** | **7.1** |

**The target is unremarkable.** Its χ² sits in the body of the distribution, not
in the tail: 8.7 against a field median of 7.2 on BLG41, and 7.3 against 7.1 on
BLG01. Whatever astrometric excursion the event produces is not large enough to
push the target out of the population of ordinary stars.

The field median of χ²/DoF ≈ 7 is itself informative: it says the standard error
of the mean within a month **underestimates the true bin uncertainty by about
√7 ≈ 2.6**, because the noise within a month is correlated rather than
independent. This is the same factor that makes the internal proper-motion
errors of section 6 optimistic.

---

### 9. Overlap of the two fields

Stars are identified through the OGLE catalogue, which resolves about six times
more sources than KMT detects; matching the two source lists directly is
ambiguous because a close pair may be split in one field and merged in the other.

| | BLG41 | BLG01 |
|---|---|---|
| sources | 594 | 621 |
| matched to an OGLE entry | 551 (93%) | 497 (80%) |
| **distinct** OGLE stars behind them | **454** | **394** |

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

---

### 10. Files

| file | content |
|---|---|
| `report_v6_RMS.png` | residual RMS against OGLE I, per axis and field |
| `report_v6_pm_vs_gaia.png` | our absolute proper motion against Gaia's, after the tie |
| `report_v6_chi2.png` | χ² distribution of the monthly binned residuals, target marked |
| `report_v6_radec_<field>_<tag>.png` | RA against time, Dec against time, and RA against Dec colour-coded by time |
| `report_v6_motion_<field>_<tag>.png` | the 2x2 motion figures, proper motion retained and removed |
| `source_motion_v6_<field>.csv` | the target, per epoch |
| `source_motion_v6_<field>_<tag>.csv` | each comparison star, per epoch |

CSV columns are `JD, X_mas, Y_mas, errX_decade, errY_decade, errX_season,
errY_season, season`, with positions relative to that star's own mean and
**proper motion not removed**. Headers carry the star's magnitude, pixel
position, sky position and absolute proper motion.

In the motion figures the columns are the two axes, **X left and Y right**; the
top row keeps the proper motion and the bottom row removes it. The legend quotes
the motion along the pixel axis and the heading quotes it on the sky. **Pixel X
runs opposite to RA**, so those two carry opposite signs in X.
