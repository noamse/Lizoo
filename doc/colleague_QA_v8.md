# Answers to a colleague's questions on the KMT-2026-BLG-0521 astrometry report

Numbers below are from the v8 solution; the v6e report the questions refer to
differs only in the third digit (χ² 1574 → 1500, etc.). Figure produced for
question 2: `doc/qa_v8_pmdiff_vs_gaia.png`.

---

## 1. The 9 mas monthly scatter, the 20 mas per-point error and the factor 2.6

Short version: the 210 points per month and the 1.4 mas white-noise error are
both too optimistic by about a factor of two, the 9 mas is an unweighted rms
dominated by a quarter of the bins that hold a handful of frames, and once the
two are put on the same footing the inflation is 2.5, i.e. √(χ²/DoF). The noise
is red on every timescale from a night to a season, which is why the
proper-motion error has to be taken from the season-to-season scatter.

### 1a. The target is not measured in every exposure

The target is a blended I = 18.1 star at the detection limit of a 3.1″ PSF.
It is detected and survives the per-epoch quality flags in **7444 of 19 133**
exposures in BLG41 and **7531 of 19 981** in BLG01 — 39% — so a sidereal-month
bin holds a median of **74 / 81** measurements, not 210, and the numbers per
bin range from 3 to about 250. The per-epoch scatter is **23.1 / 20.6 mas** (X,
the two fields). The white-noise error of a monthly mean is therefore
σ_pt/√n ≈ **2.7 / 2.3 mas** (median over bins, 2.35 / 2.07 mas), not 1.4.

### 1b. The 9 mas is the unweighted rms of the bin means

For the target in BLG41 (X axis):

| monthly bins | N | rms of bin means |
|---|---|---|
| all | 93 | 10.1 mas |
| n ≥ 50 frames | 67 | **5.8 mas** |
| n < 50 frames | 26 | 16.7 mas |
| weighted by 1/s.e.² | 93 | 5.5 mas |

BLG01: 9.3 → 3.9 mas for the 68 well-populated bins, 17.0 mas for the 24 sparse
ones, 4.0 mas weighted. The sparse bins are the season edges — February/March
and September/October, when the bulge is only up at sec z > 1.3 and KMTNet
cadence drops to a few frames per night. Their means scatter by 17 mas, and
they dominate an unweighted rms while contributing almost nothing to a χ².
(The report figures now drop bins with fewer than 10 frames and epochs above
sec z = 1.3 for this reason.)

### 1c. Reconciling with χ²/DoF ≈ 7

χ²/DoF is the mean over bins of (mean / s.e.)², each bin normalised by its own
s.e., so it is dominated by the well-populated bins and measures how much
*those* exceed white noise: √(χ²/DoF) = **2.6 / 2.2** for the target. Directly:
5.8 mas / 2.35 mas = **2.5** in BLG41, 3.9 / 2.07 = **1.9** in BLG01. Same
number. The 6.4 in the question is (unweighted rms inflated by the sparse
bins) / (white noise assuming 3× too many frames) — the two factors of ≈ 1.8
multiply.

### 1d. The physical origin: red noise on every timescale

The ratio of the observed scatter of bin means to the white-noise expectation
grows steadily with the bin length. Target, X axis, BLG41 / BLG01, and the
median over 93 / 86 calibration stars with 16 < I < 17:

| bin | frames per bin | target: observed / white | 16–17 mag stars: observed / white |
|---|---|---|---|
| night | 5 / 5 | 1.3 / 1.2 | 1.4 / 1.4 |
| week | 20 / 20 | 2.4 / 2.3 | 2.0 / 2.1 |
| sidereal month | 74 / 81 | 3.8 / 4.1 (unweighted; 2.5 / 1.9 for χ²) | 3.9 / 3.4 (2.8 / 2.7 for χ²) |
| season, about the fitted line | 668 / 738 | 4.5 / 3.2 (RA); 5.4 / 5.8 (Dec) | 7 / 6 |

Frames 12 minutes apart are already correlated at the 30% level in variance
(nightly ratio 1.3–1.4); by a month the correlation halves the effective number
of frames several times over, and by a season the scatter of the mean is 5–7×
white noise. So the correlation is not confined "inside the monthly bins"; the
months within a season are themselves correlated — that is the within-season
slope systematic of section 8, and it is the same for the target and for the
field stars.

### 1e. What the quoted proper-motion error assumes

The ± in the table is the scatter of the **ten season means** about the fitted
straight line, divided by √Σ(t−t̄)² — the standard error of a slope fitted to
ten points. Season means, BLG41, in sky coordinates:

| season | n | RA mean | Dec mean | white s.e. |
|---|---|---|---|---|
| 2016 | 468 | +0.82 | +1.57 | 0.84 / 0.79 |
| 2017 | 541 | +1.28 | −0.78 | 0.95 / 0.84 |
| 2018 | 843 | +3.98 | −5.39 | 0.78 / 0.64 |
| 2019 | 827 | −0.12 | −3.14 | 1.05 / 0.80 |
| 2021 | 937 | −1.78 | −0.87 | 0.75 / 0.66 |
| 2022 | 629 | −9.72 | +0.64 | 1.26 / 1.02 |
| 2023 | 668 | −1.84 | −6.55 | 0.81 / 0.71 |
| 2024 | 778 | −2.02 | −4.35 | 0.72 / 0.61 |
| 2025 | 1049 | +1.79 | −2.79 | 0.53 / 0.41 |
| 2026 | 648 | −3.17 | +5.08 | 0.53 / 0.39 |

Scatter about the line **3.60 / 3.72 mas** (BLG41), **2.35 / 3.91 mas** (BLG01),
against a white-noise s.e. of 0.7–0.8 mas per season — the factor 4.5–5.8 of
the last table row. √Σ(t−t̄)² = 10.37 yr, hence **0.347 / 0.359** and
**0.227 / 0.377 mas/yr**. Your formula σ/T·√(12/N) with σ = 3.6 mas, T = 9.9 yr
and N = 10 gives 0.40 — the same thing with uniform sampling assumed.

The agreement of your version with σ = 9 mas and N = 90 months is partly a
coincidence: it treats the months as independent, and the unweighted 9 mas
happens to be close to 3.6·√9. With the weighted monthly scatter (5.5 mas) the
same formula would give 0.23, i.e. the months inside a season are positively
correlated and independent-month errors would be too small by ~1.5. The
season-scatter error is the honest one because it assumes independence only
between seasons, which is the longest timescale on which we can test it — and
even that may be optimistic, since some of the season-offset pattern
reproduces between fields (section 8) and would then not average down.

This error is what the report calls "internal". The external check against
Gaia is a different thing, and it is magnitude-dependent — see 2.

---

## 2. (μ − μ_Gaia)_α against (μ − μ_Gaia)_δ

Figure: `doc/qa_v8_pmdiff_vs_gaia.png` — every Gaia-matched star of the
calibration set in both fields, same colour code as Figure 1 (blue: RUWE < 1.4,
on which the six-parameter gauge is fitted; green: the tie set, RUWE < 1.4 and
15 < I < 17; grey crosses: RUWE ≥ 1.4; red star: the target, with a circle of
radius internal ⊕ Gaia error). Grey ellipses are 1, 2, 3 robust σ of the blue
population.

Result: the target sits at **Δμ = +0.51 / +2.11** (BLG41) and **+0.96 / +2.06**
(BLG01) mas/yr. Both fields put it at the same place — off in declination by
2.1 mas/yr, which is 3.4σ of the population scatter (0.61 / 0.64) and 4.2–4.5σ
of internal ⊕ Gaia (0.47 / 0.31 for this G = 19.5 star, RUWE 1.09). Our
declination proper motion is **−7.02 / −7.07**, Gaia's is **−9.13 ± 0.31**.

Two things qualify that.

**It is typical for its magnitude.** The scatter against Gaia grows steeply
toward the faint end — far faster than either Gaia's errors or our internal
errors:

| I range | n (BLG41 / BLG01) | rstd Δμ_α | rstd Δμ_δ | median Gaia err_δ | fraction with \|Δμ_δ\| > 2 |
|---|---|---|---|---|---|
| 14–17 | 112 / 90 | 0.55 / 0.55 | 0.35 / 0.39 | 0.10 | 0–2% |
| 17–18 | 57 / 45 | 1.75 / 1.84 | 1.07 / 0.93 | 0.18 | 7–14% |
| 18–19 | 29 / 28 | 4.47 / 2.48 | 2.18 / 1.45 | 0.40 | 39–52% |

Among the 38 clean stars within 0.5 mag of the target the Dec scatter is
1.9 / 1.4 mas/yr, the target is at 1.1σ / 1.5σ of that group, and 26–29% of
them are further from Gaia than it is. In the normalised (Mahalanobis) radius
of the figure the target is at the 82nd / 84th percentile of all clean stars
and the 45th / 55th percentile of its magnitude peers. So: unusual against the
bright stars, ordinary against stars of its own brightness. The 1.0 / 0.62
mas/yr "external error" in the report is the all-star figure and understates
the disagreement with Gaia at I ≈ 18 by a factor of 2–3; at this magnitude
the KMT-minus-Gaia scatter is 1.5–2 mas/yr in Dec and 2.5–4.5 in RA.

**The obvious suspect is the blend, and it does not hold up.** The star
`m18_d04` of the comparison set sits 4.5 pixels = 1.8″ from the target at
nearly equal brightness (I = 18.14 vs 18.13), and is itself paired with a
G = 20.3 star 1.0″ further out for which Gaia solves no proper motion; under
the 3.1″ FWHM KMT PSF the three are two blobs, whereas Gaia resolves all of
them. Gaia DR3 gives the neighbour μ = +1.74 ± 0.59 / +0.27 ± 0.40 mas/yr
(G = 19.69, RUWE 0.93), 9.4 mas/yr away from the target in Dec, so a centroid
that partly followed the neighbour would indeed move less negatively than the
source — a 22% pull would reproduce our −7.0 exactly. But a static pull of
that size is excluded by the data, three ways:

| test | if the centroids mixed at the 22% (target) / 65% (neighbour) level the Dec numbers suggest | measured |
|---|---|---|
| separation of the two blobs against Gaia's 1.97″ | shrinks to ~0.25″ | **1.79–1.80″ = 91%**, and the shortfall points toward the G = 20.3 third star, not toward the target |
| seeing dependence of the target's position | strong | 5 mas per pixel of FWHM, the same as isolated stars (2.5–3); seeing trend over the decade −0.002 px/yr, hence a proper-motion change ≤ 0.15 mas/yr |
| relative motion of the pair from the raw positions, no frame solution involved | Gaia's +4.2 / +9.4 mas/yr times (1 − mixing) | +2.9 / +0.8 (BLG41), +2.0 / +0.35 (BLG01), ± 0.4 / 0.3 |

Our two centroids are pulled toward each other by at most a few per cent,
which can bias the proper motion by at most ~0.5 mas/yr. Using Gaia's
neighbour motion as a prior to "deblend" — μ_target = μ_neighbour(Gaia) −
Δμ_pair(ours) — gives −0.1 / −1.1 and +0.9 / −0.4 mas/yr, which is not a
correction but a demonstration that our neighbour blob (Gaia neighbour plus the
unsolved G = 20.3 star) does not follow Gaia's neighbour at all: we measure it
at −5.7 / −6.4 in Dec against Gaia's +0.27. The relative motion of the pair
differs between us and Gaia by 8–9 mas/yr in Dec.

So the 2.1 mas/yr declination offset of the target from Gaia is real and
unexplained by blending in our measurement. Either Gaia's astrometry of this
2″ group of three sources is off — the target's five-parameter solution has
1.0 mas of astrometric excess noise (significance 1.6) and a goodness-of-fit
of 2.0, the neighbour's is clean, the third source got only a two-parameter
solution from 74 observations — or we carry an unexplained 2 mas/yr systematic
on a faint blended star, which is inside the empirical KMT-minus-Gaia scatter
at I = 18–19 (1.5–2.2 mas/yr in Dec; 26–29% of the target's magnitude peers
are further from Gaia than it is). The data at hand cannot decide between the
two. The two fields agree with each other to 0.05 mas/yr in Dec, but they
share the images and the blend, so that is not an independent check.

What the blend does for certain: it adds to the per-epoch scatter of the
target, and any astrometric excursion of the source appears in our centroid
diluted by the flux fraction of the companion inside the fitting aperture.

This is a new finding from the colleague's suggested plot; a note on it is in
section 2 of the report.

---

## 3. "KMT detects" — where the KMT catalogue comes from

There is no external KMT catalogue. The sources are found by our own pipeline
(`KMT_pipelineI.m`, AstroPack) on the KMTNet CTIO I-band science frames of the
BLG41 and BLG01 fields, cut to 300 × 300 pixels (2′ × 2′ at 0.4″/pix) around
the target. The frames on disk are the standard KMTNet frames as delivered
(`.../KMTdata/Images/260058/KB260058_20*_CTIO_I_BLG41/RAW/`); how they were
obtained from the KMTNet collaboration is a question for the PI, not for the
pipeline.

Per exposure: background and variance maps, an empirical PSF measured from the
frame's own stars, then matched-filter source detection with that PSF in three
passes at S/N thresholds 100, 50 and 20 (bright sources detected, PSF-fitted
and subtracted before the next pass, `imProc.sources.multiIterExtractor`), with
PSF-fit positions and fluxes. Detections closer than 0.5 pixel are merged,
keeping the higher S/N. Every catalogue is shifted onto a reference frame by
cross-correlation (frames with a normalised peak correlation below 0.5 are
dropped), and the catalogues are unified across epochs within a 3-pixel radius,
keeping sources detected in at least 200 exposures. The result — a source list
with per-epoch X, Y, flux and PSF quality for every exposure — is what "the KMT
catalogue" means in the report, and "KMT detects" means detected this way in
≥ 200 frames. Magnitudes and V − I come from matching this list to the OGLE-IV
catalogue of the field.

The comparison with OGLE quantifies the resolution: OGLE, at ~1″ seeing and
0.26″/pix, lists about six times more sources in the same area, and 2660 of
its entries fainter than I = 19 have no KMT counterpart. That is the detection
limit set by a 3.1″ PSF in a bulge field at this density, and it is also why
the calibration set requires no OGLE companion brighter than I = 18 within 2″:
the isolation cut is made against the deeper catalogue precisely because KMT
cannot see the contaminants itself. (As question 2 shows, the target's
neighbour at 1.8″ and I = 18.14 passes that cut by 0.14 mag — the limit is
"brighter than 18" — so the target is in the set on its own merits, blend and
all.)
