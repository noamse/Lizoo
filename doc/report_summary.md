# KMT-2026-BLG-0521 / OGLE-2026-BLG-0058
## Astrometric analysis — summary

Prepared 2026-09-06. Solution: `~/KMTdata/Results/v3_Final/IFfinal_260058_CTIO_<field>.mat`
(variable `IFsys`). Code at commit `47c52cf` on branch `dev1`.

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
  bins** (equal population). *This is included.*
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
- **absolute astrometry** — there is no Gaia tie. Proper motions are *relative*,
  in the registered pixel frame; 400 mas/pix is assumed, not measured, and the
  orientation on the sky is unknown
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
3. **Joint fit 2016–2025** → proper motions. The event season is excluded so that
   it cannot help define the motion it is later measured against.
4. **Full-decade fit** with those motions held → one SysRem correction for the run.
5. **Ten single-season fits** reusing that correction.
6. **Final global fit, 2016–2026, motions solved.** *This is the solution behind
   every plot in this report.*

---

### 5. Calibrating stars — selection, and their colour in the plots

Selection (all conditions required):

- **14 < I_OGLE < 19**
- **no companion brighter than I = 18 within 2.0 arcsec**, i.e. **5 KMT pixels**
  at 0.4″/pix. The companion list is the OGLE catalogue, which is at 0.26″/pix
  and resolves about six times more sources than KMT detects, mapped into KMT
  pixels with the per-field offset above.
- a star's own OGLE counterpart is excluded from its companion test.

Surviving: **287 stars in BLG41, 253 in BLG01.**

**Colour key in the RMS plots**

| colour | meaning |
|---|---|
| grey dots | all field stars |
| **blue circles** | **calibrating stars**, as selected above |
| black line | running median of *all* stars, in 0.5 mag bins |
| red star | the target |

One point to be clear about: in this solution the per-epoch transformation is
fitted from **all** sources, each weighted by its own residual scatter. The
calibrating stars are a clean, well-measured subset shown for reference; they do
not exclusively define the frame. (`UseRefSources`, which would restrict the
frame fit to them, is available but is **off** here.)

---

### 6. Caveats that matter for reading the plots

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
for nine seasons and is magnified by about 1.1 mag in 2026, so:

| season | target mag | errX / errY, BLG41 | BLG01 |
|---|---|---|---|
| 1 (2016) | 18.73 | 19.28 / 18.22 | 19.32 / 19.51 |
| 4 (2019) | 18.71 | 31.92 / 26.31 | 29.21 / 26.56 |
| 6 (2022) | 18.74 | 33.81 / 27.17 | 38.55 / 31.45 |
| 9 (2025) | 18.61 | 16.85 / 15.13 | 16.76 / 14.75 |
| **10 (2026)** | **17.65** | **12.46 / 11.38** | **11.82 / 11.54** |
| decade constant | — | 21.80 / 23.48 | 22.68 / 22.63 |

Using the constant would overstate the uncertainty during the event by about a
factor two, and understate it in 2022 by a third.

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

### 7. Where the numbers stand

| quantity | BLG41 | BLG01 |
|---|---|---|
| bright-star residual (I < 17), ΔX/ΔY | 6.414 / 6.879 mas | 6.475 / 7.007 mas |
| calibrating stars, ΔX/ΔY | 8.29 / 8.30 | 7.47 / 8.14 |
| target, whole decade | 24.62 / 21.00 | 25.48 / 21.78 |
| target, 2026 season only | ≈ 11.7 / 11.1 | ≈ 12.2 / 11.7 |

The target's decade figure is dominated by the nine pre-event seasons where it
sits at I ≈ 18.1; in 2026 the event brightens it by about 1.1 mag and its
residual halves. **12 mas is the target's 2026 precision, 24 mas its
decade-average precision** — they are different quantities.

---

### 8. Files in this report

| file | contents |
|---|---|
| `report_RMS_afterPM.png` | residual RMS against OGLE I, per axis, per field — **after position and proper motion removed** |
| `report_motion_<field>_<axis>_nobin.png` | source motion, unbinned; top = position with PM **retained**, bottom = residual |
| `report_motion_<field>_<axis>_binsid.png` | the same, binned in **one sidereal month** (27.321661 d); error bars are the **RMS within the bin** |
| `source_motion_<field>.csv` | `JD, X_mas, Y_mas, errX_decade, errY_decade, errX_season, errY_season, season` |

Eight motion figures in total: 2 fields × 2 axes × 2 binnings.
