# KMT-2026-BLG-0521 / OGLE-2026-BLG-0058
## Astrometric analysis — summary

Prepared 2026-09-06. Solution: `~/KMTdata/Results/v3_Final/IFfinal_260058_CTIO_<field>.mat`

---

### 1. Data

| quantity | BLG41 | BLG01 |
|---|---|---|
| exposures | 19 133 | 19 981 |
| after cuts | 17 354 | 17 687 |
| sources fitted | 594 | 621 |
| nights | 1 985 | 2 009 |

KMT CTIO, 2016–2026, I band. Two **overlapping cut-outs** of the same star,
300 × 300 pix at 0.4″/pix (2′ × 2′). They share **no exposures**: the two fields
are observed alternately, a median of 12.2 min apart, on 1 947 common nights, so
they are independent measurements.

V-band frames exist for 2016–2022 (~1 340 per field) and are used only to check
the colour; they are not used in the astrometry.

---

### 2. What **is** in the fit

**Per epoch**
- full affine transformation, 6 parameters: translation, rotation, scale, shear

**Per source**
- position (x₀, y₀) and proper motion (μx, μy) — 4 parameters

**Detrending terms**
- **differential chromatic refraction (DCR)** — `[1, sin(pa)·secz, cos(pa)·secz]`
  plus higher orders in `secz·sin(pa)` and `secz·cos(pa)`, fitted **in 6 colour
  bins** (equal population). 
- **annual term**
- **pixel-phase correction**, from the per-epoch registration shifts
- **SysRem**, 2 components on the astrometric residuals, computed **once over the
  whole decade** and reused by every sub-fit
- iterative reweighting from the residuals, with moving-median outlier rejection

**Convergence.** 15 weighted iterations before SysRem, 6 after. The last three
iterations move the bright-star median by 0.003 mas (BLG41) and 0.001 (BLG01).

---

### 3. What is **not** in the fit

- **parallax** — bulge sources give ≈0.12 mas, far below the noise floor
- **any microlensing or astrometric-deviation model** — the solution is model-free;
  no theoretical curve of any kind has been fitted or subtracted
- **absolute astrometry inside the fit** — the solution itself is *relative*:
  proper motions are in the registered pixel frame, 400 mas/pix is assumed
  rather than measured, and the orientation is not used. A Gaia tie is available
  afterwards (`ml.util.gaiaTie`, see section 7) and gives the orientation, the
  scale and absolute proper motions, but it is applied to the finished solution
  and does not constrain the fit
- acceleration or binary motion
- any correction for the faint-end photometric bias (§6)

---

### 4. Steps

1. **`KMT_pipelineI`** — registration, PSF photometry, source matching →
   `MatchedSources`; cross-match to the OGLE catalogue for I and V, using a
   **per-field** pixel offset (BLG41 229.202/229.283, BLG01 227.542/229.689).
2. **`ml.util.mmsFromMatchedSources`** — magnitudes rebuilt from `FLUX_PSF` (the
   pipeline leaves every `MAG_*` column empty), per-epoch zero point anchored on
   OGLE I, then SysRem photometry; airmass and parallactic angle computed from
   JD and the field centre; colour from OGLE V−I; quality cuts.
3. **Joint fit 2016–2025** → proper motions
4. **Full-decade fit** with those motions held → one SysRem correction for the run.
5. **Ten single-season fits** reusing that correction.
6. **Final global fit, 2016–2026, motions solved.** 

---

### 5. Calibrating stars — selection, and their colour in the plots

Selection (all conditions required):

- **14 < I_OGLE < 19**
- **no companion brighter than I = 18 within 2.0 arcsec**, i.e. **5 KMT pixels**
  at 0.4″/pix. The companion list is the OGLE catalogue, which is at 0.26″/pix
  and resolves about six times more sources than KMT detects, mapped into KMT
  pixels with the per-field offset above.
- a star's own OGLE counterpart is excluded from its companion test.

Surviving: **287 stars in BLG41, 253 in BLG01** — these are the blue points in
Figure 1, and are *not* the smaller set the reference-frame test uses. See the
table below.

**Colour key in the RMS plots**

| colour | meaning |
|---|---|
| grey dots | all field stars |
| **blue circles** | **calibrating stars**, as selected above |
| black line | running median of *all* stars, in 0.5 mag bins |
| red star | the target |

**Two different selections appear in this report, and they must not be
confused.** Both use the isolation rule above; they differ only in the magnitude
window:

| | window | in window | isolated | share of all sources |
|---|---|---|---|---|
| **calibrating stars** — the blue points in Figure 1 | 14 < I < 19 | 529 / 470 | **287 / 253** | 48% / 41% |
| **frame stars** — what `UseRefSources` would fit the frame from | 14 < I < 16 | 87 / 62 | **36 / 35** | 6% / 6% |

BLG41 first, BLG01 second. The blue points are "clean" only in the sense of being
un-blended: the 251 extra stars on BLG41 are the faint half, 16 < I < 19, and
they measure 9.06 / 9.51 mas, essentially the field average of 9.90 / 10.23. The
36 bright ones measure 3.89 / 3.79. Requiring brightness as well as isolation is
what cuts 287 down to 36.

One point to be clear about: in this solution the per-epoch transformation is
fitted from **all** sources, each weighted by its own residual scatter. Neither
of the two sets above defines the frame. `UseRefSources`, which would restrict
the frame fit to the 36 / 35 frame stars — not to the 287 / 253 blue points —
is available but is **off** here.

**Why the frame is fitted from all stars.** This was measured rather than
assumed. `UseRefSources` was run at full scale on both fields, selecting
**14 < I < 16** and no companion brighter than I = 18 within a given radius —
the frame-star row of the table above, not the blue points — everything else at
the defaults:

| run | ref stars | rstd bright X/Y | seasonal wander |
|---|---|---|---|
| **BLG41, frame from all sources** | 594 | **8.85 / 9.23** | **2.21 / 2.57** |
| BLG41, companion radius 5 pix | 36 | 9.06 / 9.73 | 1.97 / 2.80 |
| BLG41, companion radius 3 pix | 48 | 8.67 / 9.91 | 1.96 / 2.82 |
| **BLG01, frame from all sources** | 621 | **6.85 / 7.40** | ~1.5 / 1.3 |
| BLG01, companion radius 5 pix | 10 | **17.27 / 14.55** | 5.74 / 3.55 |

mas. **BLG01 fails outright** — two and a half times worse, with the field's
seasonal wander nearly quadrupled. **BLG41 is roughly neutral**, not clearly
worse: the bright rstd degrades slightly while the event source's seasonal
wander improves by 11% in X.

**That BLG01 run was compromised, and has now been redone.** The test ran on
2026-08-26; the per-field OGLE registration fix (`e89b484`) landed on 2026-09-03.
Until then BLG41's offset was applied to BLG01, holding its OGLE coverage at
23.5% instead of 60.9% — and the reference selection needs an OGLE magnitude, so
it drew from a starved pool: **10 isolated stars, against 35 with the offsets
right**. BLG41, which the bug never touched, had 36 both times.

The whole four-step chain was therefore rerun on both fields with
`UseRefSources` and the correctly registered catalogue, against the same chain
fitting the frame from all sources. Identical inputs, identical convergence
settings:

| field | frame | bright rstd X/Y | target rstd X/Y | bright wander | target wander |
|---|---|---|---|---|---|
| BLG41 | all 594 | **6.414 / 6.879** | 24.62 / 21.00 | 2.62 / 2.77 | 10.75 / 7.00 |
| BLG41 | 36 refs | 6.569 / 6.905 | **24.54 / 20.83** | 2.69 / 2.97 | 11.49 / 7.03 |
| BLG01 | all 621 | 6.475 / **7.007** | 25.48 / 21.78 | 2.83 / 3.06 | 13.63 / 9.13 |
| BLG01 | 35 refs | **6.445** / 7.225 | **25.21** / 21.84 | 2.84 / **3.01** | 13.66 / **9.00** |

mas. **BLG01's failure is gone.** With 35 well-distributed reference stars in
place of 10 it now behaves exactly as BLG41 always did: every difference is
within 3%, and it falls both ways. The 2.5x collapse was an artefact of the
registration bug, not a property of the field or of the method.

**Neither field gains anything.** Across both fields the reference frame moves
the bright residual by at most 3%, improves the target marginally on BLG41 and
BLG01 in X, and degrades it slightly in Y. Nothing here approaches the factor of
three to five the target would need.

The tie to Gaia says the same. Re-tying the reference-frame solutions, the
target's **absolute** proper motion is stable against the frame choice even
though the relative one is not:

| field | relative, all -> refs | absolute, all -> refs |
|---|---|---|
| BLG41 | +2.20 / −2.31 -> +2.50 / −1.82 | −2.92 / −7.25 -> −3.16 / −7.18 |
| BLG01 | −0.66 / −2.90 -> −0.45 / −2.94 | −2.22 / −7.24 -> −2.25 / −7.29 |

mas/yr. The relative motion shifts by up to 0.49 mas/yr, the absolute by at most
0.24 — confirming that the shift is gauge, as restricting the frame to a
different set of stars must produce. But the all-source solutions tie *better*:
their two fields' rotations agree to **0.22 deg** (131.01 and 130.79) where the
reference-frame ones differ by **3.8 deg** (133.24 and 129.42), and their
per-source scatter about the Gaia relation is lower (3.75/3.83 against 4.54/4.29
on BLG41). That is a further, independent argument for the full field.

**Why the subset cannot win here.** The reference stars are individually
excellent — 36 survivors on BLG41 measure 3.89 / 3.79 mas and are detected in
**100%** of epochs, against 9.90 / 10.23 mas and 81% for the field — but they
are few, and two numbers show the cost:

| | ideal frame noise | share of field inverse variance | leverage at target | max leverage |
|---|---|---|---|---|
| BLG41, all 594 sources | 0.295 / 0.310 mas | 100% | 0.0017 | 0.011 |
| BLG41, 36 references | 0.653 / 0.673 mas | 20% | 0.0290 | 0.248 |
| BLG01, all 621 sources | 0.297 / 0.312 mas | 100% | 0.0016 | 0.010 |
| BLG01, 35 references | 0.662 / 0.700 mas | 20% | 0.0337 | 0.416 |
| BLG01, 10 references, the bug | — | — | 0.0761 | **0.939** |

Leverage is `h(r) = a' (Σ a a')⁻¹ a` with `a = [x, y, 1]`: the variance of the
frame shift predicted at a field position, in units of one reference star's
variance. The bugged selection reaches a max leverage near 1 — the frame at the
field corner no better determined than a single star — with its stars spanning
only 55 x 72 pix against 92 x 88 for the full field. That is the mechanism by
which ten clustered stars wrecked the solution, though the precise amplification
to a factor 2.5 was never quantified and, the run having been redone, no longer
needs to be.

Two further reasons, independent of counting:

- **The fit already does this, more gently.** Each source is weighted by its own
  residual scatter — `calculateNee` takes `median(W,1,'omitnan')` per source — so
  a noisy faint star already contributes almost nothing. `RefSrcFlag` replaces
  that soft weighting with a hard cut. The effective number of sources behind
  the weighted frame is 223 (BLG41) and 225 (BLG01), not 594 and 621: the fit
  has already done most of the selecting.
- **The dominant noise is shared.** The two-field comparison put about 74% of
  the nightly noise in the atmosphere, common to every star. Reference stars sit
  under the same atmosphere as the target, so a cleaner frame cannot reach it.
  Only the ~26% reduction noise is attackable, and 26% in variance is 15% in
  amplitude, which caps the possible gain — about the size of the shifts in the
  table above.

The approach is sound and remains implemented. `UseRefSources` stays **off** by
default, now on a clean measurement of both fields rather than one clean field
and one corrupted run.

---

### 6. The event is not confined to one season

With the OGLE parameters t₀ = JD 2461214.2 (2026-06-22), t_E = 231.7 d and
u₀ = 0.499, the magnification at each season's mean epoch is as below.

- **u** is the angular separation between lens and source at that epoch, in
  units of the Einstein radius: `u = sqrt(u₀² + ((t − t₀)/t_E)²)`. It is
  dimensionless, falls to u₀ = 0.499 at closest approach, and is large long
  before and after the event.
- **A** is the **flux** magnification, `A = (u² + 2) / (u·sqrt(u² + 4))`. A = 1
  means unmagnified; A = 2 means twice as much light.
- **Δmag** is the same quantity in magnitudes, `Δmag = −2.5·log₁₀(A)`. The two
  columns are not independent: A = 2.067 is Δmag = −0.788, and −0.788 mag
  inverts to A = 2.066. Where the text below quotes "a magnification of about
  2.1" and "0.8 mag" it is quoting the same number twice.

| season | (t − t₀) [d] | u | A | Δmag |
|---|---|---|---|---|
| 2016 | −3644 | 15.74 | 1.000 | 0.000 |
| 2019 | −2564 | 11.08 | 1.000 | 0.000 |
| 2022 | −1474 | 6.38 | 1.001 | −0.001 |
| 2023 | −1114 | 4.83 | 1.003 | −0.003 |
| 2024 | −754 | 3.29 | 1.012 | −0.013 |
| **2025** | −394 | 1.77 | **1.085** | **−0.089** |
| **2026** | −44 | 0.53 | **2.067** | **−0.788** |

The event is above 1% magnification for **4.4 years**, above 5% for 2.6 years and
above 10% for 2.0 years.

The model predicts 0.80 mag of brightening at the 2026 season median, whereas
the measured magnitudes give 1.085 (18.724 at baseline against 17.639 in 2026).
The 0.29 mag difference is the faint-end photometric bias of section 7, which is
+0.39 mag at I ≈ 18.1 and +0.13 at I ≈ 17.3: a *difference* of two measured
magnitudes inflates by the difference of their biases, 0.80 + 0.26 = 1.06
against 1.085 observed. It is a consistency check on the bias, not a discrepancy
with the model.

Calling 2026 "the event season" and 2016–2025 "pre-event" is therefore wrong:
**2025 is already magnified by 9%** and 2024 by 1.3%.

This matters for the proper-motion baseline. The astrometric deviation of a
microlensing event does not peak with the photometric one: for these parameters
it is largest near u = √2, which falls in **2025** (2.46 mas), with 2024 at
1.91 mas and 2026 back down to 1.63 mas as the source passes closest approach.
So the joint fit that supplies the proper motions, which holds out 2026 alone,
still contains the two seasons carrying the **largest** astrometric signal.

Any proper motion quoted here is therefore fitted through the event, and would
absorb part of a real astrometric excursion. A baseline free of it would have to
stop around 2022, which costs a third of the time span and most of the leverage
on the proper motion. Nothing in this report depends on that choice — no
astrometric model is fitted — but a future attempt to measure the deviation
must not treat 2016–2025 as an event-free baseline.

### 7. Caveats that matter for reading the plots

**The measured magnitudes run faint at the faint end.** `KMT − I_OGLE` is flat to
±0.04 down to I = 17, then rises to **+0.50 by I = 18.75**. The target's own
measured magnitude is 18.70 against its catalogue value of **18.13**. All plots
here therefore use the **OGLE catalogue I** on the magnitude axis. PSF and
aperture photometry bracket OGLE and disagree by ~0.96 mag at I = 18.5, which
points to crowding: OGLE lists 3 783 sources in this cut-out where KMT detects
594. **The astrometry is not affected** — positions come from the PSF centroid,
not the flux, and removing the photometric anchor entirely was measured to cost
0.012 mas.

**The CSV carries two error options; choose deliberately.**

| column | what it is | value |
|---|---|---|
| `errX_decade`, `errY_decade` | the black running-median curve of the RMS plot read at the target's catalogue magnitude I = 18.13, over the whole run. **Constant.** | 21.80 / 23.48 (BLG41), 22.68 / 22.63 (BLG01) |
| `errX_season`, `errY_season` | the target's **own** residual RMS within its observing season. Follows its brightness. | 11.4-33.8, see below |

Neither is a per-epoch measurement uncertainty. The decade value is a population
estimate at a fixed magnitude; the season value is the target's own scatter.

They differ by a factor of two where it matters most. The target is at I ≈ 18.1
for nine seasons and is magnified by 0.80 mag in 2026, so:

| season | target mag | errX / errY, BLG41 | BLG01 |
|---|---|---|---|
| 1 (2016) | 18.73 | 19.28 / 18.22 | 19.32 / 19.51 |
| 4 (2019) | 18.71 | 31.92 / 26.31 | 29.21 / 26.56 |
| 6 (2022) | 18.74 | 33.81 / 27.17 | 38.55 / 31.45 |
| 9 (2025) | 18.61 | 16.85 / 15.13 | 16.76 / 14.75 |
| **10 (2026)** | **17.65** | **12.46 / 11.38** | **11.82 / 11.54** |
| decade constant | — | 21.80 / 23.48 | 22.68 / 22.63 |

Using the constant would overstate the uncertainty at the peak of the event by about a
factor two, and understate it in 2022 by a third.

**The sign of a proper-motion component carries no physical meaning.** The fitted
motion of the target is +2.20 / −2.31 mas/yr in BLG41 and −0.66 / −2.90 in
BLG01: the same star, the same nights, opposite signs in X. This is a gauge
difference between the two frames, not anything about the source. The per-epoch
transformation is free at every epoch, so it absorbs any motion common to the
field, and what survives is defined only relative to the mean motion of each
cut-out's own star population. Those populations differ:

| quantity | BLG41 | BLG01 |
|---|---|---|
| median PM of the field's own sources | +2.128 / +0.183 | −0.261 / −0.617 |

Measured directly on the 383 well-measured stars common to both fields, the
disagreement is the same for every star:

| quantity | X | Y |
|---|---|---|
| PM difference, BLG01 − BLG41, all common stars | −2.41 ± 0.43 | −0.60 ± 0.57 |
| the target | −2.86 | −0.59 |

The target's offset is the population's offset. Fitting the difference as a
constant plus a linear term in position gives, in addition, a relative rotation
of 0.28 deg/yr and a scale drift of 4.7e-3 /yr between the two frames; removing
both leaves 1.14 / 0.35 mas/yr, consistent with per-star measurement error.

So the Y values (−2.31 and −2.90) agree only because the Y gauge offset is small,
and the X values disagree only because the X offset is large. **What is physical
is the differential motion between stars within one field**, unless the frame is
tied to an external one — which it now can be.

**Tied to Gaia, the two fields agree.** `ml.util.gaiaTie` carries a solution onto
the Gaia frame through OGLE, giving the orientation, the plate scale and an
absolute proper motion for every matched source:

| quantity | BLG41 | BLG01 |
|---|---|---|
| KMT sources matched to Gaia | 532 of 594 (90%) | 447 of 621 (72%) |
| median separation | 0.259 arcsec | 0.358 arcsec |
| calibrators (Gaia PM error < 0.4 mas/yr, I < 18) | 426 | 327 |
| fitted rotation of the pixel axes on the sky | **+131.01 deg** | **+130.79 deg** |
| gauge offset removed | −0.36 / +4.41 mas/yr | −2.84 / +3.88 |
| scatter about the Gaia relation | 3.75 / 3.83 mas/yr | 3.38 / 3.51 |

**The target's proper motion:**

| frame | BLG41 | BLG01 | difference |
|---|---|---|---|
| relative (detector x, y), as fitted | +2.196 / −2.310 | −0.660 / −2.900 | 2.856 / 0.590 |
| **absolute (μ_α cos δ, μ_δ), Gaia frame** | **−2.917 / −7.248** | **−2.215 / −7.244** | **0.702 / 0.004** |

mas/yr, taken from `IF.ParS(3:4,Ie)` at 400 mas/pix, target index 383 (BLG41)
and 450 (BLG01). **The two rows are not in the same axes**: the relative values
are in the cut-out's detector x/y, the absolute ones are (μ_α cos δ, μ_δ), since
`gaiaTie` returns `Tie.A \ (PM − Offset)` and so carries the pixel frame into
the Gaia frame through a 131 deg rotation with parity and a 0.92 scale.

The sign disagreement disappears once both are on the same physical
frame, and the two rotations, derived independently from separately reduced
data, agree to 0.22 deg — which is what makes the tie credible rather than a fit
to noise. The relative disagreement itself is gauge, not astrophysics: both
cut-outs have nearly the same orientation (131.01 and 130.79 deg), so the
+2.20 against −0.66 in x is the unconstrained constant-plus-gradient freedom of
each field's own frame, which is exactly what the tie removes.

**Against Gaia's own value for this star.** The target has a Gaia counterpart,
so the absolute motion can be checked directly rather than only field against
field:

| | μ_α cos δ | μ_δ |
|---|---|---|
| Gaia DR3, the target | −2.479 | −9.129 |
| ours, BLG41 | −2.917 (−0.44) | −7.248 (+1.88) |
| ours, BLG01 | −2.215 (+0.26) | −7.244 (+1.89) |

mas/yr, with the difference from Gaia in brackets. μ_α cos δ agrees well. μ_δ sits
1.88 mas/yr away, and by nearly the same amount in both fields — which should
**not** be read as two independent confirmations, because BLG01 reused BLG41's
OGLE-to-Gaia solution and a common systematic in that bridge would reproduce
this pattern exactly. Against the per-source scatter it is about 0.5 sigma and
unremarkable. It is not the event: 1.88 mas/yr over the decade baseline is some
19 mas of accumulated displacement, an order of magnitude beyond the ~2.5 mas
peak deviation the model allows.

**Which error applies.** Three different numbers circulate and they measure
different things:

- **3.4 to 3.8 mas/yr** — scatter about the Gaia relation, *per source*. This is
  the one that applies to a single star's absolute proper motion, the target
  included.
- **0.16 / 0.17 mas/yr** — the error on the tie *transformation*, from 532 and
  447 matched sources. It bounds the frame, not any individual star.
- **~1.0 / 0.7 mas/yr** — the target's own relative proper motion uncertainty,
  from the scatter of its season means about its fitted line.

Reduction choices contribute far less: the pixel-phase on/off variants give
+2.21 / −2.31 and +2.27 / −2.27 (BLG41), −0.65 / −2.89 and −0.71 / −2.94
(BLG01), a spread of about 0.06 mas/yr.

Two cautions. The μ_δ agreement of 0.004 mas/yr between the fields is far better
than either field's own precision and is therefore luck; the honest statement is
that the two agree to within errors of order 1 mas/yr on the mean. And the tie
fixes the **frame, not the precision**: our per-source scatter of 3.4 to 3.8
mas/yr exceeds Gaia's own spread of 2.9 to 3.3, so most of the per-star
difference is ours and nothing in the residual analysis improves.

**Plotted positions are relative to the target's own mean position**, not to the
field centre and not to any absolute reference. Each motion curve has the
target's decade-average position subtracted, so the zero line is that average.
The frame itself is set by the ensemble of 594 (BLG41) and 621 (BLG01) field
stars through the per-epoch transformation, not by the geometric centre of the
cut-out: differences along a curve are meaningful, the absolute level is not,
and the level cannot be compared between the two fields, whose cut-outs have
different origins and independently fitted frames. For scale, the target lies
+0.42, +0.39 pix (+168, +155 mas) from the centre of the BLG41 cut-out and
+1.55, +0.15 pix (+620, +60 mas) from the centre of BLG01; the centring absorbs
both.

**Only sources with an OGLE counterpart appear in the RMS plots** (540 of 594 in
BLG41, 478 of 621 in BLG01), which biases the plotted sample slightly bright.

---

### 8. Where the numbers stand

| quantity | BLG41 | BLG01 |
|---|---|---|
| bright-star residual (I < 17), ΔX/ΔY | 6.414 / 6.879 mas | 6.475 / 7.007 mas |
| calibrating stars, ΔX/ΔY | 8.29 / 8.30 | 7.47 / 8.14 |
| target, whole decade | 24.62 / 21.00 | 25.48 / 21.78 |
| target, 2026 season only | ≈ 11.7 / 11.1 | ≈ 12.2 / 11.7 |

The target's decade figure is dominated by the seasons in which it sits near its
baseline I ≈ 18.1; in 2026 the magnification of about 2.1 brightens it by 0.8 mag
and its residual halves. 2025 is brighter than baseline too, by 0.09 mag. **12 mas is the target's 2026 precision, 24 mas its
decade-average precision** — they are different quantities.

---

### 9. Files in this report

| file | contents |
|---|---|
| `report_RMS_afterPM.png` | residual RMS against OGLE I, per axis, per field — **after position and proper motion removed** |
| `report_motion_<field>_nobin.png` | source motion, unbinned; 2x2 panel, columns X and Y, top row = position with PM **retained**, bottom row = residual |
| `report_motion_<field>_binsid.png` | the same, binned in **one sidereal month** (27.321661 d); error bars are the **RMS within the bin** |
| `source_motion_<field>.csv` | `JD, X_mas, Y_mas, errX_decade, errY_decade, errX_season, errY_season, season` |

Eight motion figures in total: 2 fields × 2 axes × 2 binnings.
