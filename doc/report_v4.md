# KMT-2026-BLG-0521 / OGLE-2026-BLG-0058

## The restricted-set reductions (v4 and v5) — standalone report

Prepared 2026-09-08. Solutions: `~/KMTdata/Results/v4b/IFfinal_<field>.mat` and
`~/KMTdata/Results/v5/IFfinal_<field>.mat` (variable `IFsys` in each). This
document describes **only** these two runs, which share one source set and
differ in exactly one respect: v4 pins the frame to Gaia during the fit, v5
fits everything freely and ties to Gaia afterwards. The decade
reduction that remains the project's reference solution is described in
`KMT260058_astrometry_report.pdf`.

---

### 1. What this run is

A reduction built to a different specification from the reference one, to test
two ideas: restrict the analysis to a clean calibration set, and tie the frame
to Gaia from the start rather than afterwards.

1. **Analysis set.** 14 < I_OGLE < 17, no companion brighter than I = 18 within
   2.0 arcsec, and a Gaia counterpart with **RUWE < 1.4**. The target is added
   explicitly — at I = 18.13 it falls outside the magnitude window. Every other
   source is **dropped at the input**, so all measurements, detrending and
   quality cuts act on this set alone.
2. **The frame is pinned to Gaia.** Each calibration star's position *and*
   proper motion are held at Gaia's values, converted into cut-out pixels. The
   per-epoch affine is therefore fitted against Gaia positions propagated to
   that epoch by Gaia's own motions. The target is left free; pinning it would
   erase the signal. The six-parameter gauge is gone by construction and no
   post-hoc tie is needed.
3. **Passengers.** Three stars at the target's own magnitude (I ~ 18), outside
   the calibration set, are carried along unpinned so that something comparable
   to the target can be plotted in the same frame.

| | BLG41 | BLG01 |
|---|---|---|
| raw sources in the MSc file | 1177 | 1187 |
| surviving the standard quality cuts — **the reference reduction's set** | **594** | **621** |
| **analysis set of these runs** | **130** | **100** |
| of which pinned to Gaia (v4 only) | 126 | 96 |
| free: target + passengers | 4 | 4 |
| epochs | 19 133 | 19 981 |

The row that matters for comparison is the second: the reference reduction fits
**594 and 621** sources, and these runs fit 130 and 100 of them. The raw file
holds about 1180 entries, but half of those are marginal detections that the
standard conversion discards before any reduction sees them — they are found in
a median of **10%** of epochs, against 80% for the sources that survive. Every
source detected in more than half the epochs is kept, without exception. So the
restriction imposed here is 130 out of 594, not 130 out of 1177.

---

### 2. The calibration set

Selection, all conditions required:

- **14 < I_OGLE < 17**
- **no companion brighter than I = 18 within 2.0 arcsec** (5 KMT pixels at
  0.4″/pix), the companion list being the OGLE catalogue, which resolves about
  six times more sources than KMT detects
- **a Gaia DR3 counterpart with RUWE < 1.4**

| | BLG41 | BLG01 |
|---|---|---|
| in the magnitude window, isolated | 144 | 127 |
| with a Gaia counterpart | 144 | 116 |
| **after RUWE < 1.4** | **126** | **96** |

The RUWE cut removes 12% and 17% — consistent with the 87% pass rate for
G < 17 across the whole cone, so crowding is not biasing it heavily.

---

### 3. Gaia

**RUWE is not in the local catsHTM DR3 table.** It carries
`astrometric_chi2_al` and `astrometric_n_good_obs_al` but neither RUWE nor the
`u0(G, BP−RP)` table needed to normalise UWE into it. The true values were
queried from the Gaia archive instead — 9464 sources in a 3 arcmin cone, cached
in `~/KMTdata/GaiaRef/`. This mattered: at G < 17, **87%** pass RUWE < 1.4,
where the unnormalised UWE proxy predicted 74%, so the proxy would have
discarded about 13% of good stars for no reason.

Gaia's proper-motion errors for these stars are **0.070 / 0.044 mas/yr**, some
fifty times better than our internal frame. That is what makes pinning
worthwhile in principle.

**The Gaia-to-pixel map is fitted directly** from the matched pairs, not carried
through OGLE. The OGLE bridge is needed only to identify which star is which;
routing coordinates through it inherits OGLE's ~110 mas astrometric error.

| | through OGLE | **fitted directly** |
|---|---|---|
| BLG41, median offset / scatter | 252 / 159 mas | **23 / 16 mas** |
| BLG01 | 351 / 148 mas | **21 / 20 mas** |

The fitted scale is 2.519 (BLG41) and 2.529 (BLG01) pix/arcsec against the
2.500 expected at 0.4″/pix — an independent check that the map is right. The
remaining ~20 mas is static per star, so it displaces each star's residuals by
a constant and does not affect proper motion or any time-dependent signal.

---

### 4. The chain

With proper motions supplied by Gaia, the joint proper-motion fit and the
per-season fits of the reference reduction are unnecessary. Two steps remain:

1. **Full-decade SysRem correction**, frame pinned, 2 components on the
   astrometric residuals.
2. **Final solution** reusing that correction, with the pixel-phase term on,
   15 weighted iterations before SysRem and 6 after.

Convergence over the frame stars, last three iterations: **0.0009 mas (BLG41)**
and **0.0026 mas (BLG01)**.

---

### 5. Results

| | BLG41 | BLG01 |
|---|---|---|
| frame-star residual RMS, ΔX / ΔY | 6.13 / 5.20 mas | 5.30 / 5.42 mas |
| target residual RMS, ΔX / ΔY | 29.62 / 22.63 mas | 28.73 / 23.25 mas |

For comparison, the reference reduction gives bright-star 6.41 / 6.88 and
6.48 / 7.01, and target **24.62 / 21.00** and **25.48 / 21.78**.

**The calibration stars are measured better and the target worse.** The stars
are bright and isolated, which is what they were selected for; the target is
**14 to 20% worse** than in the reference reduction. The cause is the frame: it
is now built from 126 stars rather than ~600 weighted sources, which roughly
doubles its ideal noise and raises its leverage by a factor of twenty to fifty.
Restricting the analysis set is what costs. The Gaia tie is not at fault.

---

### 6. Proper motions, and whether the two fields agree

Converting the target's motion out of pixel axes with each field's own map:

| | BLG41 | BLG01 | difference |
|---|---|---|---|
| **v4, Gaia-pinned** | −1.799 / −6.780 | −1.158 / −7.147 | **0.640 / 0.367** |
| reference, post-hoc tie | −2.92 / −7.25 | −2.22 / −7.24 | 0.702 / 0.004 |

mas/yr, as (μ_α cos δ, μ_δ). Both agree at about the 1 mas/yr level, which is
the target's own uncertainty, but **pinning the frame did not tighten the
agreement**. The post-hoc tie was already doing this job.

**The passengers agree better than the target does.** Three stars at the
target's magnitude, freely fitted, carried through both fields:

| star | BLG41 | BLG01 | difference |
|---|---|---|---|
| I = 18.14, 4.5 pix from the target | −0.63 / −5.47 | −1.08 / −5.47 | 0.45 / 0.00 |
| I = 17.91, 21 pix | −6.81 / −9.99 | −6.94 / −10.06 | **0.13 / 0.07** |
| I = 18.05, 48 pix | +3.06 / −6.34 | +2.80 / −6.48 | **0.26 / 0.14** |
| **the target**, I = 18.13 | −1.80 / −6.78 | −1.16 / −7.15 | 0.64 / 0.37 |

The two clean passengers reproduce two to five times better than the target.
The 4.5-pix star's apparent agreement is not a good sign — see section 7.

**Sensitivity to the set.** Adding only those three passengers to a set of 127
moved the fitted target motion by **0.42 / 0.63 mas/yr on BLG41** and
0.09 / 0.41 on BLG01. A result that moves by half a mas per year when three
faint stars join the set is not a precise result, and this is the most direct
measurement here of how uncertain the target's motion actually is.

---

### 7. A blended neighbour tracks the event

The passenger 4.5 pixels from the target behaves unlike anything else in the
run. Comparing the 2026 season against the 2016–2022 baseline:

| star | distance | ΔX (BLG01 / BLG41) | ΔY (BLG01 / BLG41) | season scatter |
|---|---|---|---|---|
| target | 0 | +2.7 / +2.6 | +7.1 / +6.4 | 3.6–5.3 |
| **passenger 1** | **4.5 pix** | **−44.1 / −57.8** | **−27.9 / −32.2** | 16–19 |
| passenger 2 | 21 pix | +5.6 / +3.5 | +2.9 / +3.3 | 2.3–3.7 |
| passenger 3 | 48 pix | −3.4 / −5.5 | −2.5 / −2.1 | 2.1–5.3 |

mas. The excursion is **40 to 58 mas at 4.5 pixels and gone by 21 pixels**,
reproduced independently in both fields, and that star's season-to-season
scatter is four times the others' in every season, not only in 2026. With a
seeing FWHM near 7 pixels, a star 4.5 pixels away shares most of its light with
the target: as the target brightens through the event, the blend's centroid
moves. The neighbour moves *away* from the target, consistently in both fields,
which is what PSF fitting a blend tends to do as flux is reattributed to the
brightening component.

**Coordinates.** Pixel positions are per field, since the two cut-outs are
offset by about 1.1 pixels; sky positions are derived by inverting each field's
own Gaia-to-pixel map, and the two fields agree on them to **2 to 9 mas**, which
is an independent check on both maps.

| star | I_OGLE | dist [pix] | BLG41 pixel | BLG01 pixel | RA (J2000) | Dec (J2000) |
|---|---|---|---|---|---|---|
| `target` | 18.13 | 0.0 | 150.4, 150.4 | 151.6, 150.2 | 17:52:38.085 | −31:47:36.21 |
| `cal_d24` | 16.68 | 24.5 | 161.3, 172.2 | 162.4, 172.1 | 17:52:37.746 | −31:47:27.53 |
| `cal_d30` | 16.63 | 30.5 | 120.4, 154.6 | 121.4, 154.5 | 17:52:39.023 | −31:47:34.60 |
| `cal_d43` | 16.94 | 42.9 | 110.9, 166.6 | 111.9, 166.5 | 17:52:39.318 | −31:47:29.89 |
| `out_d04` | 18.14 | 4.5 | 146.4, 148.3 | 147.6, 148.0 | 17:52:38.210 | −31:47:37.06 |
| `out_d21` | 17.91 | 21.3 | 155.3, 171.1 | 156.4, 170.9 | 17:52:37.935 | −31:47:28.01 |
| `out_d48` | 18.05 | 48.5 | 198.6, 145.4 | 199.8, 145.0 | 17:52:36.585 | −31:47:38.10 |

The target's recovered position, 17:52:38.085 −31:47:36.21, sits 0.11 arcsec
from the OGLE EWS coordinates the pipeline was pointed at, which is the accuracy
of the centroid plus the map, not a disagreement.

Three consequences:

- **That star is unusable as a comparison.** It was chosen as the closest at the
  target's magnitude, and being closest is exactly what disqualifies it.
  Passengers 2 and 3 are the meaningful ones.
- **Its good inter-field agreement in section 6 is misleading.** Both fields see
  the same blending artefact from the same event, so agreeing is not evidence
  of a good measurement.
- **It bears directly on the target.** The blend is mutual. If the neighbour's
  centroid moves 40 to 58 mas when the target brightens, the target's centroid
  is being pulled by the same shared flux. The target's own +2.7 / +7.1 mas in
  2026 therefore cannot be read cleanly as astrometric microlensing; some or all
  of it may be blending. That is a systematic at the tens-of-mas level against a
  predicted signal of a few mas.

---

### 8. Overlap of the two fields

Matching the two source lists directly, after fitting the offset between the
cut-outs (−1.11 / +0.24 pix):

| match radius | BLG41 matched | BLG01 matched |
|---|---|---|
| 0.5 pix | 339 / 594 (57%) | 352 / 621 (57%) |
| **1.0 pix** | **546 / 594 (92%)** | **575 / 621 (93%)** |
| 1.5 pix | 576 / 594 (97%) | 609 / 621 (98%) |

These are indeed two observations of nearly the same sky, and the cut-outs span
the same pixel range. The unmatched 7% are faint and marginal: median MAG_PSF
17.18 against 16.87 for matched stars, and a detection fraction of 0.54 against
0.84. They are stars one reduction found and the other did not, or blends split
differently. Restricted to MAG_PSF < 17 the agreement is 94%.

Only 57% agree within 0.5 pix, so two independent reductions place a common
star half a pixel to a pixel apart; 1 pixel is the appropriate radius.

Matching *through* Gaia and OGLE instead gives a misleading 58%, because their
identification coverage over our sources is only 90% / 72% (Gaia) and 91% / 77%
(OGLE) — sources neither catalogue identifies cannot be paired at all and count
as non-overlapping.

---

### 9. v5 — the same set, freely fitted

v4 changes two things at once against the reference reduction: it shrinks the
source set, and it pins the frame. v5 separates them. It uses **exactly the same
source list** as v4b — read from the same file, not re-selected — the same two
steps, the same iteration counts and the same SysRem and pixel-phase settings.
The only difference is that nothing is pinned: every calibration star is freely
fitted and takes part in the frame determination at every step, and Gaia enters
only afterwards, as a post-hoc tie.

So **v3 to v5 measures what restricting the set costs**, and **v5 to v4b
measures what pinning does**.

**Target residual RMS**

| | BLG41 ΔX / ΔY | BLG01 ΔX / ΔY |
|---|---|---|
| **v3** — ~600 sources, free, post-hoc tie | **24.62 / 21.00** | **25.48 / 21.78** |
| **v5** — 130/100 sources, free, post-hoc tie | 25.61 / 22.27 | 28.45 / 23.15 |
| **v4b** — the same sources, pinned to Gaia | 29.62 / 22.63 | 28.73 / 23.25 |

| step | isolates | BLG41 X / Y | BLG01 X / Y |
|---|---|---|---|
| v3 → v5 | restricting the set | +4.0% / +6.0% | +11.7% / +6.3% |
| v5 → v4b | pinning the frame | +15.7% / +1.6% | +1.0% / +0.4% |

Both hurt and neither dominates consistently. Restricting the set costs 4 to 12%
everywhere; pinning costs almost nothing except in BLG41's X, where it costs
16%. The reference reduction remains the best for the target on every axis in
both fields.

**Frame stars.** Over exactly the same 126 / 96 stars:

| | BLG41 | BLG01 |
|---|---|---|
| v5, freely fitted | **4.87 / 4.81** | **4.77 / 4.77** |
| v4b, pinned to Gaia | 6.13 / 5.20 | 5.30 / 5.42 |

4.8 mas is the best frame-star astrometry this project has produced. The
comparison is not entirely fair — the free fit adjusts four parameters per star
and so absorbs part of each star's own systematics, which the pinned version
cannot — but the ~1.3 mas gap is the price of holding stars at catalogue values.

Convergence over the frame stars: 0.0041 mas (BLG41) and 0.0002 (BLG01).

**Inter-field agreement of the target's absolute motion.**

| | difference between the fields |
|---|---|
| v3 | 0.702 / 0.004 |
| v5, `gaiaTie` | 0.072 / 0.123 |
| v5, positional-map tie | 0.363 / 0.294 |
| v4b | 0.640 / 0.367 |

mas/yr. v5 agrees best, but **the headline 0.072 should not be quoted alone**.
`gaiaTie` fits a 2 by 2 matrix from a proper-motion regression, and on these
small sets it returned rotations of **146.34 deg for BLG41 and 108.65 for
BLG01** — 38 deg apart, where the direct positional map puts the two cut-outs
0.126 deg apart. That matrix is fitting noise. Re-tying independently, with the
linear part taken from the positional map and only the gauge constant from the
proper-motion comparison, gives 0.363 / 0.294 instead. The honest statement is
that v5's fields agree to **0.07 to 0.36 in RA and 0.12 to 0.29 in Dec**, and
that the choice of tie method is itself a systematic of 0.2 to 0.6 mas/yr.

The gauge removed from v5 is −0.219 / +4.758 mas/yr (BLG41, 130 calibrators) and
+2.834 / +3.983 (BLG01, 100). With it removed, the seven plotted stars agree
between the fields to 0.25 to 0.78 mas/yr in RA and 0.02 to 0.87 in Dec.

---

### 10. Verdict

Neither v4 nor v5 is adopted as the primary solution. Restricting the analysis
to the calibration set costs 4 to 12% on the target, and pinning the frame to
Gaia costs a further 0 to 16%; the reference reduction is better on every axis
in both fields. What the pair does establish, by differencing, is that **both
choices cost something and the set restriction is the more consistent penalty**
— pinning is nearly free except on BLG41's X axis.

v5 nevertheless produces the best frame stars (4.8 mas) and the best inter-field
agreement of the target's absolute motion of any reduction here, and it does so
while leaving every calibration star free to take part in the frame. If a
restricted-set reduction is wanted, v5 rather than v4 is the one to use.

Its most valuable output is not the astrometry but section 7: the direct
measurement of how far a blended neighbour is dragged by the event, which
applies to the target in every reduction we have made.

---

### 11. Figures and files

| file | content |
|---|---|
| `report_v4_RMS.png` | residual RMS against OGLE I, per axis and field, for this run |
| `report_v4_motion_<field>_target.png` | the target, sidereal-month bins, 2x2 |
| `report_v4_motion_<field>_cal_d{24,30,43}.png` | three calibration stars, proper motion **pinned to Gaia** |
| `report_v4_motion_<field>_out_d{04,21,48}.png` | the three I ~ 18 passengers, freely fitted |
| `report_v5_motion_<field>_*.png` | the same seven stars in **v5**, where nothing is pinned; the quoted sky motions have the gauge removed and are absolute |

The `_d NN` suffix is the star's distance from the target in pixels, and names
the **same physical star in both fields** — the figures are directly comparable
across fields. (An earlier version ranked the comparison stars per field, so
`cal1` denoted different objects in BLG41 and BLG01.)

**Axis convention.** The Gaia-to-pixel map has a negative determinant: it
contains a parity flip, and the pixel X axis runs *opposite* to increasing RA.

| | rotation, flip removed | scale |
|---|---|---|
| BLG41 | −0.317 deg | 2.5194 pix/arcsec |
| BLG01 | −0.444 deg | 2.5290 pix/arcsec |

So the cut-outs are essentially aligned with RA and Dec apart from that X flip,
with a third of a degree of residual rotation, and the two fields differ in
orientation by only 0.126 deg. In each figure the **legend** quotes the motion
along the pixel axis and the **panel heading** quotes it on the sky, so the two
carry opposite signs in X by construction; both are labelled accordingly.

In every motion figure the columns are the two axes, X left and Y right; the top
row keeps the proper motion and the bottom row removes it. **A pinned star's red
line is Gaia's proper motion, not a fitted one**, so its lower panels test the
frame, whereas the target's and the passengers' test the fit.
