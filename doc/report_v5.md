# KMT-2026-BLG-0521 / OGLE-2026-BLG-0058

## The v5 reduction — restricted set, freely fitted, Gaia applied afterwards

Prepared 2026-09-09. Solution: `~/KMTdata/Results/v5/IFfinal_<field>.mat`
(variable `IFsys`; the logical `Frame` in the same file marks the calibration
stars). This document describes **only** v5. The reference decade reduction is
in `KMT260058_astrometry_report.pdf`; the Gaia-pinned variant and the v4/v5
comparison are in `KMT260058_v4_report.pdf`.

---

### 1. What v5 is, in one paragraph

The analysis is restricted to a clean calibration set — bright, isolated, and
astrometrically well behaved in Gaia — plus the target and three comparison
stars at the target's own magnitude. Every one of those stars is **freely
fitted** and every one takes part in determining the per-epoch frame at every
step. Gaia is used only to *select* the set, and then again at the very end, as
a post-hoc tie that fixes the frame's gauge. Nothing is held at a catalogue
value during the fit.

This is the counterpart to v4, which uses the identical source set but pins
each calibration star's position and proper motion to Gaia while fitting.

---

### 2. Object selection, step by step

Both fields, from the raw matched-source file down to what is actually fitted:

| step | BLG41 | BLG01 |
|---|---|---|
| raw sources in the MSc file | 1177 | 1187 |
| survive the standard quality cuts | **594** | **621** |
| have an OGLE I magnitude | 540 | 478 |
| in 14 < I < 17 | 292 | 253 |
| ... and isolated: no I < 18 companion within 2.0 arcsec | 144 | 127 |
| ... and matched to a Gaia DR3 source | 144 | 115 |
| ... and that source has a RUWE | 142 | 111 |
| ... and **RUWE < 1.4** → the calibration set | **126** | **96** |
| + the target + 3 comparison stars | +4 | +4 |
| **= analysis set actually fitted** | **130** | **100** |

Notes on the individual cuts:

- The **quality cuts** are the pipeline's standard ones and are not part of this
  design. They remove sources detected in a median of 10% of epochs, against 80%
  for those retained; every source found in more than half the epochs survives.
  The restriction imposed by v5 is therefore **130 of 594**, not 130 of 1177.
- The **magnitude window** costs about half the OGLE-matched sources. It is the
  single largest cut after the standard ones.
- The **isolation cut** costs another half: 292 → 144 on BLG41. The companion
  list is the OGLE catalogue, which resolves about six times more sources than
  KMT detects, so it finds close pairs that our own source list has merged.
- The **RUWE cut** removes 11% (BLG41) and 14% (BLG01) of what reaches it,
  consistent with the 87% pass rate for G < 17 across the whole cone. Crowding
  is not biasing it heavily.
- BLG01 loses 12 sources at the Gaia-match step where BLG41 loses none. This
  follows its poorer OGLE coverage, 478 of 621 against 540 of 594.
- **The target is not in the calibration set** — at I = 18.13 it falls outside
  the magnitude window — and is added explicitly.

---

### 3. Gaia

Gaia enters twice, and in neither case during the fit.

**Selecting the set.** RUWE is **not** in the local catsHTM DR3 table, which
carries `astrometric_chi2_al` and `astrometric_n_good_obs_al` but neither RUWE
nor the `u0(G, BP−RP)` table needed to normalise UWE into it. The true values
were queried from the Gaia archive — 9464 sources in a 3 arcmin cone, cached in
`~/KMTdata/GaiaRef/`. This mattered: at G < 17, **87%** pass RUWE < 1.4, where
the unnormalised UWE proxy predicted 74%.

**Fixing the gauge, afterwards.** A free fit determines relative astrometry
only: applying any global affine to all source positions, and compensating in
the per-epoch transformations, leaves every residual unchanged. Six parameters
are therefore pure gauge, and the same freedom appears one time-derivative up in
the proper motions as a constant plus a linear gradient across the field. The
tie removes it. The gauge found was:

| | μ_α cos δ | μ_δ | calibrators |
|---|---|---|---|
| BLG41 | −0.219 | +4.758 | 130 |
| BLG01 | +2.834 | +3.983 | 100 |

mas/yr. With it removed, the seven plotted stars agree between the two fields to
0.25 to 0.78 mas/yr in RA and 0.02 to 0.87 in Dec.

**A caution on the tie method.** `ml.util.gaiaTie` fits a 2 by 2 matrix from a
proper-motion regression. On these small sets it is poorly determined: it
returned rotations of **146.34 deg (BLG41) and 108.65 (BLG01)**, 38 deg apart,
where the direct positional Gaia-to-pixel map puts the two cut-outs 0.126 deg
apart. Its linear part is fitting noise. The gauge above was therefore taken the
other way round: the linear part from the positional map, which is solid, and
only the constant from the proper-motion comparison. The two methods differ by
0.2 to 0.6 mas/yr on the target, which is a systematic worth carrying.

---

### 4. The chain

Two steps, on the 130 / 100 source subset:

1. **Full-decade SysRem correction.** 2 components on the astrometric residuals,
   computed once over the whole run and reused.
2. **Final solution**, reusing that correction, with the pixel-phase term on,
   15 weighted iterations before SysRem and 6 after.

The joint proper-motion step and the per-season fits of the reference reduction
are not used: with so few, such well-measured sources the motions are solved
directly in the global fit.

| | BLG41 | BLG01 |
|---|---|---|
| step 1 | 16.3 min | 14.8 min |
| step 2 | 30.0 min | 28.0 min |
| convergence over the calibration stars, last 3 iterations | 0.0041 mas | 0.0002 mas |

---

### 5. Results

| | BLG41 ΔX / ΔY | BLG01 ΔX / ΔY |
|---|---|---|
| **calibration stars** (126 / 96) | **4.87 / 4.81** | **4.77 / 4.77** |
| the target | 25.61 / 22.27 | 28.45 / 23.15 |

mas. **4.8 mas is the best frame-star astrometry this project has produced** —
better than the reference reduction's bright stars (6.41 / 6.88 and 6.48 / 7.01)
and better than the same objects pinned to Gaia in v4 (6.13 / 5.20 and
5.30 / 5.42). Some of that advantage is self-fulfilling: a free fit adjusts four
parameters per star and so absorbs part of each star's own systematics.

The target is **not** improved. The reference reduction measures it at
24.62 / 21.00 and 25.48 / 21.78, better on every axis in both fields. Restricting
the analysis set costs 4 to 12% on the target while improving the calibration
stars, because the per-epoch frame is now built from 126 stars instead of ~600
weighted sources.

**Target proper motion**, absolute, gauge removed:

| | μ_α cos δ | μ_δ |
|---|---|---|
| BLG41 | −1.80 | −6.20 |
| BLG01 | −2.13 | −6.51 |

mas/yr; the two fields differ by 0.33 / 0.31. Using `gaiaTie`'s own matrix
instead gives −1.978 / −6.757 and −1.906 / −6.880, differing by 0.072 / 0.123.
The honest statement is that the fields agree to **0.07 to 0.36 in RA and 0.12
to 0.29 in Dec**, the spread being the choice of tie method.

---

### 6. Comparison stars

Six per field, all freely fitted, all plotted in the same style as the target:

| tag | I | distance from the target | decade residual RMS, BLG41 / BLG01 |
|---|---|---|---|
| `cal_d24` | 16.68 | 24 pix | 4.07 / 4.14 |
| `cal_d30` | 16.63 | 30 pix | 6.73 / 6.83 |
| `cal_d43` | 16.94 | 43 pix | 4.48 / 4.60 |
| `out_d04` | 18.14 | 4.5 pix | **28.61 / 32.37** |
| `out_d21` | 17.91 | 21 pix | 17.36 / 17.54 |
| `out_d48` | 18.05 | 48 pix | 19.12 / 18.59 |

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

mas in X. The three `cal_` stars are members of the calibration set, chosen from
those present in **both** fields so the same physical object carries the same
name everywhere. The three `out_` stars sit at the target's own magnitude and
are deliberately outside the set, to give something comparable to the target.

**`out_d04` is not a usable comparison object.** It lies 4.5 pixels from the
target under a seeing FWHM near 7 pixels, so it shares most of its light with
it. Its residual RMS is 29 to 32 mas against 17 to 19 for the other two stars of
the same brightness, in every season and not only during the event; and its 2026
season is displaced by 40 to 58 mas as the target brightens. Its CSV carries a
warning in the header. `out_d21` and `out_d48` are the honest controls, and their
errors bracket the target's own 25 to 27 mas — the check that the target is not
anomalously noisy for what it is.

---

### 7. Files

| file | content |
|---|---|
| `report_v5_RMS.png` | residual RMS against OGLE I, per axis and field, for this run |
| `report_v5_motion_<field>_target.png` | the target, sidereal-month bins, 2x2 |
| `report_v5_motion_<field>_cal_d{24,30,43}.png` | the three calibration comparison stars |
| `report_v5_motion_<field>_out_d{04,21,48}.png` | the three I ~ 18 comparison stars |
| `source_motion_v5_<field>.csv` | the target, per epoch |
| `source_motion_v5_<field>_<tag>.csv` | each comparison star, per epoch |

CSV columns are `JD, X_mas, Y_mas, errX_decade, errY_decade, errX_season,
errY_season, season`. Positions are relative to that star's own mean, with
**proper motion not removed**. `err*_decade` is that star's own residual RMS over
the whole run and `err*_season` the same within each season; note this differs
from the reference reduction's CSV, where the decade column held a field-star
population estimate — v5 keeps no field stars at I ~ 18, so that curve cannot be
formed here.

In every motion figure the columns are the two axes, **X left and Y right**; the
top row keeps the proper motion and the bottom row removes it. The legend quotes
the motion along the pixel axis and the heading quotes it on the sky with the
gauge removed. **Pixel X runs opposite to RA**, so those two carry opposite signs
in X by construction.
