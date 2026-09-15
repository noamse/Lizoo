# KMT-2026-BLG-0521 / OGLE-2026-BLG-0058 — the astrometric reduction, step by step

How the reference solution (**v13**) is produced from the raw KMTNet frames,
what every step needs and writes, and which scripts make the report figures and
the diagnostic tests. Written 2026-09-15.

All code is MATLAB with the `~/matlab/Lizoo` package (which extends AstroPack)
unless stated. The scripts named below are collected in **`runKMT/` in this repository (`~/matlab/Lizoo/runKMT/` on the reduction machine)**:

| folder | contents |
|---|---|
| `~/matlab/Lizoo/runKMT/reduction/` | the chain from raw frames to the tied solution, in order |
| `~/matlab/Lizoo/runKMT/report/` | the scripts that make the report figures, CSVs and the PDF |
| `~/matlab/Lizoo/runKMT/diagnostics/` | the tests of section 8 of the report and the answers to the colleague's questions |

They are working scripts, not functions: paths, field names and version tags are
written into them. Each is run with `run('<script>')` from a MATLAB session with
`~/matlab/Lizoo` on the path, or from the shell as
`FIELD=BLG41 matlab -batch "run('~/matlab/Lizoo/runKMT/reduction/v13chain.m')"` where a
script reads its field from the environment (`getenv('FIELD')`). The two fields
are always processed as two independent runs.

Conventions used throughout:

- **Fields**: `BLG41` and `BLG01`, two KMTNet CTIO cut-outs of the same sky,
  300 × 300 pixels at 0.4″/pixel around the target (RA 17:52:38.09,
  Dec −31:47:36.1, I_OGLE = 18.13).
- **Pixel axes**: X runs opposite to right ascension, Y along declination.
  1 pixel = 400 mas.
- **Version tags** (`v3`, `v6`, `v8`, `v13`, …) are directories under
  `~/KMTdata/Results/`; each holds the fit of one reduction variant. The
  reference is `v13`.

---

## 1. Data products, from raw to final

| stage | file(s) | made by |
|---|---|---|
| raw frames | `/bigdata3/projects/KMTdata/Images/260058/KB260058_20*_CTIO_I_<field>/RAW/*I*.fits` | KMTNet |
| OGLE-IV catalogue of the field | `~/matlab/Lizoo/OGLEdata/OB260058/OB160058.mat` | OGLE (X, Y, I, V per star) |
| per-exposure source lists, matched across epochs | `~/KMTdata/Results16_26_v2/KMT_260058_<field>_MSc.mat` | step 1 |
| the same, re-matched to OGLE with per-field offsets (**the input of everything below**) | `~/KMTdata/Results16_26_v3/KMT_260058_<field>_MSc.mat` | step 2 |
| Gaia DR3 RUWE for the 120″ cone | `~/KMTdata/GaiaRef/gaia_dr3_ruwe_cone3arcmin.csv`, `ruwe_cone120.mat` | step 3 |
| first full solution on every source | `~/KMTdata/Results/v3_Final/IFfinal_260058_CTIO_<field>.mat` | step 4 |
| calibration set, first pass (magnitude + isolation) | `~/KMTdata/GaiaRef/prep6_<field>.mat` | step 5 |
| solution on that set | `~/KMTdata/Results/v6/IFfinal_<field>.mat` | step 6 |
| calibration set, second pass (2σ RMS clip) | `~/KMTdata/GaiaRef/prep8_<field>.mat` | step 7 |
| **reference solution** | `~/KMTdata/Results/v13/IFfinal_<field>.mat`, `SysCor_<field>.mat` | step 8 |
| Gaia tie of the reference solution | `~/KMTdata/Results/v13/Tie_<field>.mat` | step 9 |
| comparison stars | `~/KMTdata/Results/v13/cmp.mat` | step 10 |
| report figures, CSVs, PDF | `~/matlab/Lizoo/doc/report_v13_*.png`, `source_motion_v13_*.csv`, `KMT260058_v13_report.pdf` | section 3 |

The MSc files are `MatchedSources` objects (AstroPack): per-epoch × per-source
matrices in `MSc.Data` (`X`, `Y`, `FLUX_PSF`, `MAG_PSF`, `SN`, `PSF_CHI2DOF`, …),
per-source data in `MSc.SrcData` (`I_ogle`, `V_ogle`), the epoch times in
`MSc.JD`, and in `MSc.UserData` the pipeline's per-epoch registration shifts
(`ShiftX`, `ShiftY`), cross-correlation quality (`NormPeakCorr`) and flags
(`FlagGoodEpoch`, `FlagGoodSrc`).

The `IFfinal` files hold `IFsys` (an `IterFit` object: the fitted parameters
`ParE` per epoch, `ParS` per source, `ParHalat` per colour bin, the data and
weights, and methods such as `calculateResiduals`, `calculateRstd`,
`calculateWes`), `Info` (bookkeeping: `SrcInd`, `EpochInd`, `SrcData`, SysRem
provenance) and `Config`.

---

## 2. The reduction chain

### Step 1 — source extraction and cross-epoch matching: `reduction/rerunpipe.m`

Runs `KMT_pipelineI` (in `~/matlab/Lizoo`) on every raw frame of a field:

```matlab
addpath('~/matlab/Lizoo'); addpath('~/matlab/Lizoo/OGLEdata/OB260058');
Temp = '/bigdata3/projects/KMTdata/Images/260058/KB260058_20*_CTIO_I_BLG41/RAW/*I*.fits';
[~, JD, MSc, ~] = KMT_pipelineI([], 'TempName', Temp);
save('~/KMTdata/Results16_26_v2/KMT_260058_BLG41_MSc.mat', 'MSc', 'JD', '-v7.3');
```

What it does per exposure: background and variance maps
(`ThresholdBack = 4000`), an empirical PSF from the frame's own stars
(`PopPSFArgs`: wings threshold 0.15, empirical wings), matched-filter detection
with that PSF in three passes at S/N thresholds **100, 50, 20**
(`imProc.sources.multiIterExtractor`: bright sources are fitted and subtracted
before the next pass), PSF-fit positions and fluxes. Detections closer than
0.5 pixel are merged, keeping the higher S/N. Frames with fewer than 250
sources are dropped. Every catalogue is shifted onto the frame of image
`Iref = 500` by cross-correlation (frames with a normalised peak correlation
below 0.5 are flagged bad), and the catalogues are unified across epochs
within a **3-pixel** radius keeping sources with at least **200** detections
(`pipeline.generic.proc2MatchedSources`, `MinNdet = 200`). `FlagGoodSrc`
marks sources with more than 800 detections and a position scatter below
1 pixel. Finally the sources are matched to the OGLE catalogue (1.5 pixel) for
`I_ogle` and `V_ogle`. Takes several hours per field.

### Step 2 — re-match to OGLE with the right per-field offset: `reduction/rematch.m`

The OGLE-to-cut-out registration is a per-field translation, measured as
`[229.202 229.283]` OGLE pixels for BLG41 and `[227.542 229.689]` for BLG01
(0.26″ OGLE pixels; the conversion is `X = (corrX − off)·0.26/0.4 + 150`).
The pipeline had used one offset for both fields. Only `SrcData.I_ogle` and
`V_ogle` depend on it, so the saved MSc files are re-matched rather than
re-extracted:

```matlab
run('~/matlab/Lizoo/runKMT/reduction/rematch.m')   % Results16_26_v2 -> Results16_26_v3, both fields
```

### Step 3 — Gaia RUWE: `reduction/step0b.m` (first half)

RUWE is not in the catsHTM copy of Gaia DR3. A 3′ cone around the target was
exported from the Gaia archive as `gaia_dr3_ruwe_cone3arcmin.csv` and matched
(0.05″) onto the 120″ `catsHTM.cone_search('GAIADR3', …)` that `ml.util.gaiaTie`
uses, so that `RUWE(GaiaInd)` indexes that cone:

```matlab
run('~/matlab/Lizoo/runKMT/reduction/step0b.m')    % writes ~/KMTdata/GaiaRef/ruwe_cone120.mat
```

(The second half of the script only prints the count of frame stars per field.)

### Step 4 — first full solution, every source: `reduction/v3chain.m`

```bash
FIELD=BLG41 matlab -batch "run('~/matlab/Lizoo/runKMT/reduction/v3chain.m')"
FIELD=BLG01 matlab -batch "run('~/matlab/Lizoo/runKMT/reduction/v3chain.m')"
```

Three calls of `ml.scripts.runIterDetrendMSc` on the full MSc: (1) a solution
with the last season excluded, saved as a proper-motion prior; (2) the
decade-wide SysRem correction with those motions held (`'FixedPM'`); (3) the
final solution with `'PixPhase',true,'NiterWeightsBeforeSys',15,
'NiterWeightsAfterSys',6`. This solution's only later use is to define the
calibration set in step 5 (positions, residual RMS and OGLE identity of every
source). ~1.5 h per field.

### Step 5 — calibration set, first pass: `reduction/v6prep.m`

```matlab
run('~/matlab/Lizoo/runKMT/reduction/v6prep.m')    % writes ~/KMTdata/GaiaRef/prep6_<field>.mat
```

From the v3 solution, `ml.util.selectRefSources` with `'MagRange',[14 19]`,
`'CompanionRadius',5` (pixels = 2.0″), `'CompanionMaxMag',18` and the OGLE
companion catalogue of the field (`ml.util.ogleCompanionCat(…,'Field',f)`):
keep stars with 14 < I_OGLE < 19 that have **no OGLE star brighter than
I = 18 within 2″**. The isolation test is made against OGLE because it
resolves about six times more sources than KMT. The target is added
regardless (it passes anyway). Result: 287 / 253 stars. The file holds
`SrcIdx`, the indices of the kept stars **in the MSc file**.

### Step 6 — solution on the calibration set: `reduction/v6chain.m`

```bash
FIELD=BLG41 matlab -batch "run('~/matlab/Lizoo/runKMT/reduction/v6chain.m')"    # and BLG01
```

Selects `MSc.selectBySrcIndex(prep6.SrcIdx)` and runs the same two-step fit as
the reference (see step 8), with two SysRem components. Output
`~/KMTdata/Results/v6/`.

### Step 7 — calibration set, second pass, 2σ RMS clip: `reduction/v8sel.m`

```matlab
run('~/matlab/Lizoo/runKMT/reduction/v8sel.m')     % writes ~/KMTdata/GaiaRef/prep8_<field>.mat
```

On the v6 solution, for every star the decade residual RMS
`sqrt(rx² + ry²)` from `IFsys.calculateRstd`; in 0.5-mag bins of I_OGLE the
median and robust std of `log10(RMS)` are taken and stars more than **2σ
above the median** are dropped (log space, so the cut is multiplicative and
does not shave the faint end). The target is never dropped. Result:
**260 / 235** stars. Indexing is the one thing to get right here:
`Info.SrcInd` of a solution indexes *that run's input list*, so the MSc index
of v6 source `i` is `prep6.SrcIdx(v6.Info.SrcInd(i))`; the file stores
`SrcIdx` (MSc indices), `keep` and `TargetMsc`, the target's MSc index.

### Step 8 — the reference solution: `reduction/v13chain.m`

```bash
NSYS=6 FIELD=BLG41 matlab -batch "run('~/matlab/Lizoo/runKMT/reduction/v13chain.m')"
NSYS=6 FIELD=BLG01 matlab -batch "run('~/matlab/Lizoo/runKMT/reduction/v13chain.m')"
# or both at once: ~/matlab/Lizoo/runKMT/reduction/run_v13.sh (edit the loop to N=6 only)
```

With `NSYS = 6` the output goes to `~/KMTdata/Results/v13c6/`, which is
copied to `~/KMTdata/Results/v13/`. Two steps:

```matlab
P  = load('~/KMTdata/GaiaRef/prep8_BLG41.mat');
L  = load('~/KMTdata/Results16_26_v3/KMT_260058_BLG41_MSc.mat','MSc');
MS = L.MSc.selectBySrcIndex(P.SrcIdx);

% step A: the decade-wide SysRem correction
[~,~,~,I2] = ml.scripts.runIterDetrendMSc(MS, 'NIterSysRem',6, 'Verbosity',0);
SysCorX=I2.SysCorX; SysCorY=I2.SysCorY; SrcInd=I2.SrcInd(:); JD=MS.JD(I2.EpochInd); JD=JD(:);
save('~/KMTdata/Results/v13/SysCor_BLG41.mat','SysCorX','SysCorY','JD','SrcInd','-v7.3');

% step B: the final solution with that correction applied
Pre = load('~/KMTdata/Results/v13/SysCor_BLG41.mat');
[IFsys,~,IFsysB,Info] = ml.scripts.runIterDetrendMSc(MS, ...
    'SysRemCorrection',Pre, 'PixPhase',true, ...
    'NiterWeightsBeforeSys',15, 'NiterWeightsAfterSys',6, 'Verbosity',0);
Ie = find(P.SrcIdx(Info.SrcInd(:)) == P.TargetMsc);     % the target, by MSc index
save('~/KMTdata/Results/v13/IFfinal_BLG41.mat','IFsys','Info','Config','-v7.3');
```

What `runIterDetrendMSc` does, and the parameters that matter:

- **Input conversion** (`ml.util.mmsFromMatchedSources`): builds the `MMS`
  object, computes sec z and the parallactic angle for every epoch from the
  target's coordinates and CTIO (`Args.RA`, `Args.Dec`, `Args.GeoPos`), drops
  epochs with `MaxSecz = 1.6` or duplicate JD, applies the pipeline's
  `FlagGoodEpoch` and `NormPeakCorr` flags (`ApplyPipelineFlags = true`),
  drops sources with fewer than `MinNdetSrc = 200` detections, a detection
  fraction below `MinDetFrac = 0.3` or a position scatter above
  `MaxSrcStdXY = 1` pixel, epochs with fewer than `MinNsrcEpoch = 50`
  sources, and duplicated sources. Magnitudes are rebuilt from `FLUX_PSF`
  (`ZP0 = 25`) and anchored per epoch on `I_ogle` (`ApplyRefZP`, `perepoch`);
  the colour `V − I` from OGLE is placed in `NColourBins = 6`
  equal-population bins (`ColourMode = 'ownbin'`). Recovers the pixel phase
  from `UserData.ShiftX/Y` when they are aligned with the epochs. Result:
  17 347 / 17 684 epochs and 260 / 235 sources.
- **Pass 1** (`ml.scripts.runIterDetrend`, `NiterNoWeightsBeforeSys = 2`
  unweighted then `NiterWeightsBeforeSys = 15` weighted iterations of
  `IterFit.runIter`): alternating least squares over the per-epoch affine
  transformation (6 parameters in `ParE`: translation, rotation, scale, shear,
  as deviations from the identity), the per-source position at `JD0` and
  proper motion (4 parameters in `ParS`, pixels and pixels/yr), the
  differential chromatic refraction (`HALat = true`, `ChromaicHighOrder =
  true`: `[1, sin(pa)·secz, cos(pa)·secz]` plus `(secz·sin pa)ⁿ`,
  `(secz·cos pa)ⁿ` for n = 2–4, 18 parameters per colour bin, `ParHalat`),
  the **annual term** (`AnnualEffect = true`: a quartic in the phase of the
  year, 10 parameters per colour bin, fitted and subtracted from the data
  every iteration) and the **pixel-phase term** (`PixPhase = true`: a global
  quintic in the sub-pixel phase, subtracted). Weights are each star's own
  residual scatter (`UseWeights = true`) with moving-median outlier
  rejection. Parallax is off (`Plx = false`).
- **SysRem** between the passes: `NIterSysRem` decade-wide components
  (`ml.util.sysRemScriptPart`, each component a per-epoch mode times a
  per-source coefficient, fitted on the weighted residuals of pass 1 and
  subtracted from X and Y); a component larger than `MaxSysRemShift = 1`
  pixel is rejected. In step A it is fitted (`RunSysRem = true`), in step B the
  saved correction is applied (`SysRemCorrection`), matched on JD and MSc
  source index. **Six components** is the reference (v13); v8 used two.
- **Pass 2** (`FinalStep = true`, `NiterWeightsAfterSys = 6`): the same fit on
  the corrected data with the annual and pixel-phase terms forced off (they
  have been removed from the data); DCR is fitted again.

~30 min for step A and ~55 min for step B per field, the two fields in
parallel.

Variants used for the tests (all in `diagnostics/`, same structure): `v12chain`
(no SysRem, `'RunSysRem',false`), `v13chain` with other `NSYS`, `v14chain`
(annual term off, `'AnnualEffect',false` in both steps), `v11chain` (epoch cuts
before the fit, with `UserData`'s per-epoch vectors subset alongside — see
`v11flag.m`), `v9chain` (SysRem per season), `v8b12chain` (12 colour bins,
`'NColourBins',12`), `v10chain` (prep10, the season-scale calibrator clip).

### Step 9 — the Gaia tie: `reduction/v8btie.m`

```bash
RUN=v13 matlab -batch "run('~/matlab/Lizoo/runKMT/reduction/v8btie.m')"   # writes Results/v13/Tie_<field>.mat
```

A free fit determines relative astrometry only. `ml.util.gaiaTie` identifies
the Gaia DR3 counterpart of every source (through the OGLE catalogue, cone of
120″). Then, on the tie stars — **RUWE < 1.4 and 15 < I_OGLE < 17** (104 / 85
stars) — Gaia standard coordinates (propagated to `JD0` with Gaia's own proper
motions, `dt = 2019.4127 − 2016.0`) are fitted linearly onto the cut-out
pixels (`Amap`, a 2 × 2 matrix plus offset, three 3σ-clip iterations), giving
the map from pixel motion to sky motion, `sky = 1000·(Amap \ ParS(3:4,:))`
mas/yr. The gauge freedom of a free fit in the proper motions is a constant
plus a linear gradient across the field; it is fitted on all RUWE-clean stars
(198 / 163) against Gaia — six parameters, `Gauge` — and removed from every
source. Output per field: `Tie`, `Amap`, `PMabs` (absolute proper motions of
every source), `clean` (RUWE < 1.4), `tieSet`, `Gauge`, `J` (Gaia index into
the cone). The scatter of the clean stars about Gaia after the tie —
0.967 / 0.593 and 1.034 / 0.606 mas/yr — is the external check of the
reduction.

### Step 10 — comparison stars: `reduction/v6pick.m`, `reduction/v8pick.m`

Six stars chosen once (`v6pick.m`, on the v6 solution, by Gaia identity so
that they are the same objects in both fields: three at I ≈ 18.1 at 4.5, 17
and 48 pixels from the target, three at I ≈ 16.7 at 24–31 pixels), then
located in a later solution by OGLE identity (`v8pick.m` → `cmp.mat` holding
`Sel.BLG41`, `Sel.BLG01`, `Sel.tags`, `Sel.Ie_<field>`). The v13 solution has
exactly v8's source list (verified index by index), so `Results/v13/cmp.mat`
is a copy of v8's. Any solution with a different list needs the pick re-run
with its version name.

### Step 11 — proper motions with errors: `reduction/v13pm.m`

```matlab
run('~/matlab/Lizoo/runKMT/reduction/v13pm.m')
```

For every star, the position with the proper motion retained is averaged per
observing season (`ml.util.splitEpochGroups`, seasons with ≥ 20 epochs), a
straight line is fitted to the season means, and the error of the slope is
the scatter about that line divided by √Σ(t − t̄)² — the standard error of a
slope fitted to ten points. It is converted to the sky through `Amap`. This is
the ± quoted for the target; it assumes independence between seasons only,
which is the longest timescale on which the red noise (report section 7) can
be tested.

---

## 3. The report: figures, CSVs, PDF — `~/matlab/Lizoo/runKMT/report/`

Run in this order after steps 8–11 (each reads `Results/v13/` and writes into
`~/matlab/Lizoo/doc/`):

| script | writes | what |
|---|---|---|
| `v13f1.m` | `report_v13_RMS.png`, `report_v13_pm_vs_gaia.png` | residual RMS vs I_OGLE per axis and field; our absolute PM against Gaia's after the tie |
| `v13f2.m` | `report_v13_chi2.png`, `Results/v13/chi2.mat` | χ² of the sidereal-month binned residuals of every star, `(mean/s.e.)²` summed over bins and axes, DoF = 2·N_bins − 4 |
| `v13f3.m` | `report_v13_radec_<field>_<tag>.png` (6) | RA and Dec residuals against time and against each other, colour = time, for the target and `m16_d24`, `m16_d30` |
| `v13figs.m` | `report_v13_motion_<field>_<tag>.png` (14) | 2 × 2 motion figures per star: X and Y, PM retained and removed |
| `v13e_ba.m` | `report_v13_target_beforeafter.png` | the target in both fields before and after PM detrending |
| `v13csv.m` | `source_motion_v13_<field>[_<tag>].csv` (14) | per-epoch positions of the target and each comparison star |
| `mkreport_v13.py` | `~/matlab/Lizoo/runKMT/report/report_v13.html` | Markdown (`doc/report_v13.md`) + figures → HTML |

then

```bash
cd ~/matlab/Lizoo/runKMT/report && python3 mkreport_v13.py
soffice --headless --convert-to pdf:writer_web_pdf_Export report_v13.html --outdir ~/matlab/Lizoo/runKMT/report
cp ~/matlab/Lizoo/runKMT/report/report_v13.pdf ~/matlab/Lizoo/doc/KMT260058_v13_report.pdf
```

(LibreOffice hangs on `page-break-before:always` and on ragged tables, and
aborts on a stale `.~lock` file — the builder avoids the first two.)

Binning conventions common to the figures: sidereal-month bins (27.321661 d),
bins with fewer than **10 epochs** and epochs at **sec z > 1.3** are dropped
(they are the sparse season edges, whose means scatter by 17 mas against
5–6 for the rest); error bars are the standard error of the bin mean,
σ/√N, **multiplied by** √(χ²/DoF) of that star's bin means about zero, per
star and per axis, quoted in each panel title. Residual RMS quoted per star is
`IterFit.calculateRstd` (robust std over the decade after position and
proper motion are removed).

Two generic versions for any other solution: `runfigs.m` (`RUN=v11z
FLD=BLG41`) makes the target's motion figure, `runrms.m` (`RUN=v11z`) the RMS
figure, both into a scratch directory rather than `doc/`.

CSV columns: `JD, X_mas, Y_mas, errX_decade, errY_decade, errX_season,
errY_season, season`, positions relative to the star's own mean with the
**proper motion not removed**; headers carry magnitude, pixel and sky
position and the absolute proper motion.

---

## 4. Diagnostics — `~/matlab/Lizoo/runKMT/diagnostics/`

Every script reads a solution by version tag (`RUN=<tag>` in the environment
where present, otherwise hard-wired to v8/v13) and prints its result; the
numbers in section 8 of the report and in `doc/colleague_QA_v8.md` come from
these.

**The within-season slope systematic**

| script | what it measures |
|---|---|
| `slopeorigin.m` (`RUN=`) | per star, a straight line fitted to the residuals within each season; the rms of the ten slopes against the white-noise expectation; correlations with magnitude, colour, nearest-neighbour distance, radius; writes `Results/slopeQ.mat` |
| `slopecross.m` | matches the stars of the two fields through OGLE and correlates each star's season-slope pattern between the fields — the test that the slopes are real (median +0.73 / +0.66, 96% positive in v13) |
| `varimp.m`, `multireg.m`, `seeingtest.m` | per-star regression of the residuals on the observing conditions (fwhm, secz, secz·sin pa, secz·cos pa, PSF shape, χ²) and the leave-one-out importance of each; explains 24–36% of the slope variance but not the cross-field part |
| `slopepix.m` | the slopes against the within-season drift of each star's sub-pixel phase |
| `bintest.m`, `v8b12chain.m` | the colour-bin count (`NColourBins`) propagation and the 12/20-bin refits |
| `v9chain.m` | SysRem refitted per season and stitched |
| `v12chain.m`, `run_v12.sh` | no SysRem |
| `v13chain.m` (in `reduction/`), `run_v13b.sh` | the SysRem-rank series 3, 4, 6, 8, 12 |
| `v14chain.m`, `run_v14.sh` | the annual term off |

**Epoch and calibrator selections**

| script | what |
|---|---|
| `v11diag.m` | what the first/last 14 nights of each season contain, against sec z cuts |
| `v11flag.m`, `v11chain.m`, `run_v11.sh`, `v11cmp.m` | the epoch-cut refits (`MODE=trim` or `secz130`) and their comparison with v8 on the common epochs |
| `seasonsel_diag.m`, `v10sel.m`, `v10chain.m`, `v10cmp.m`, `v10gaia.m` | the season-scale calibrator clip (statistic B: season-offset χ²/DoF), the prep10 refit and its comparison with v8 on the common stars |
| `seczdist.m` | airmass distribution of the target's epochs and the effect of sec z caps on its monthly bins |
| `outdiag.m`, `outdiag2.m` | what the outlying monthly bins are made of (frames, nights, airmass, seeing; whether the field shares them) |

**The error budget and the Gaia comparison (colleague's questions)**

| script | what |
|---|---|
| `qa_noise.m` | rms of bin means against bin length (night, week, month, season) versus σ_pt/√N — the red-noise ladder |
| `qa_season.m` | the season means of the target and the reproduction of the quoted PM error |
| `qa_pmdiff.m`, `qa_pmdiff_std.m`, `qa_pmdiff_plain.m` | (μ − μ_Gaia)_α against (μ − μ_Gaia)_δ for every star, with robust-σ, plain-σ, or no contours |
| `qa_pop.m` | the scatter against Gaia by magnitude band |
| `qa_blend.m`, `blend1.m` … `blend4.m` | the target's 1.8″ companion: Gaia sources within 4″, blob separation against Gaia, relative motion of the pair from raw positions, seeing dependence of both centroids, the single-detection epochs, the control on other close pairs, the seeing-trend PM refit |
| `v13cmpstars.m` | decade RMS of the comparison stars in a solution |

---

## 5. Where the numbers in the report come from

| report item | script / output line |
|---|---|
| bright-star and target RMS (section 2) | `v13chain.m` log line `STEP B`, or `IFsys.calculateRstd` |
| target PM ± (section 2) | `v13pm.m` |
| Gaia scatter, tie counts, map residual, gauge (sections 3–4) | `v8btie.m` with `RUN=v13` |
| χ² table (section 7) | `v13f2.m` and `Results/v13/chi2.mat` |
| comparison-star table (section 6) | `v13cmpstars.m`; positions from `ParS`, sky through `Amap` |
| slope statistics, cross-field correlation (section 8) | `slopeorigin.m`, `slopecross.m` with `RUN=v13` |
| candidate tests (section 8) | the variant chains and comparison scripts listed in section 4 |
| field overlap counts (section 9) | from the MSc source lists; unchanged between versions |
