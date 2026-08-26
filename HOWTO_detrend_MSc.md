# Detrending a KMT_pipelineI MatchedSources

How to run the astrometric detrending on the `MSc` object that `KMT_pipelineI`
produces, instead of going through `ml.scripts.runAstrometryField`, which needs
the per-epoch `AstCat` files and the `msMatch` matching step.

## Normal use

```matlab
addpath('~/matlab/Lizoo');
OutFile = ml.scripts.runAstrometryMSc( ...
    '~/KMTdata/Results16_26_v2/KMT_260058_BLG41_MSc.mat', ...
    'EventNum',260058, 'Site','CTIO', 'Field','BLG41', 'Verbosity',1);
```

About 36 min for ~17,000 epochs and ~600 sources. Everything below is optional.

## The three functions

| function | does | use it when |
|---|---|---|
| `ml.scripts.runAstrometryMSc` | convert, detrend, optional Gaia step, save | normal use |
| `ml.scripts.runIterDetrendMSc` | convert and detrend, returns objects, saves nothing | you want the objects in the workspace |
| `ml.util.mmsFromMatchedSources` | convert `MatchedSources` to `MMS` only | inspecting or debugging the conversion |

Each is a thin layer over the next, so arguments for the inner ones are passed
through: `runAstrometryMSc(..., 'runIterDetrendMScArgs',{...})` and
`runIterDetrendMSc(..., 'mmsFromMatchedSourcesArgs',{...})`.

## Input

A `MatchedSources` object, or the path of a `.mat` file holding one. Nothing
else is needed: the converter rebuilds everything the fit requires that
`KMT_pipelineI` does not store.

- magnitudes, every `MAG_*` column being NaN, from `FLUX_PSF` with a per-epoch
  zero point anchored on `SrcData.I_ogle`, then `SysRemPhotometry`
- `secz`, `pa`, `ha`, `alt` from `MSc.JD`, the cut-out coordinates and the site
- colour from `SrcData.V_ogle - I_ogle`
- pixel phase, but only from files written after commit `d1dfc59`, which stores
  the registration shifts in `MSc.UserData`

## Parameters worth knowing

Defaults are sensible; these are the ones to reach for.

**Selection** — `mmsFromMatchedSourcesArgs`

| parameter | default | note |
|---|---|---|
| `MinDetFrac` | `0.3` | least fraction of epochs a source must be detected in. The fit needs a well filled matrix: below roughly 70% fill the sysrem step and the airmass solve both misbehave |
| `MinNdetSrc` | `200` | absolute floor, acts together with `MinDetFrac` |
| `MaxSrcStdXY` | `1` | pix, drops wandering or blended sources |
| `MaxSecz` | `1.6` | airmass cut, removes about 9% of epochs |
| `MinNsrcEpoch` | `50` | epochs too empty to constrain their own affine |

**Colour** — `ColourMode`, default `ownbin`

`ownbin` puts the sources without an OGLE colour into the lowest colour bins
rather than mixing them in with real colours. `fill` gives them the median
colour, `fillcmd` a colour predicted from their magnitude, `restrict` keeps
only sources that have one. Measured on BLG41, `ownbin` is about 15% better
than `fill`; on BLG01 the four are within 3% of each other. `restrict` is
dangerous when few sources have colours: on BLG01 it left 168 sources and the
result was clearly the worst.

**Reference stars for the frame** — top level on both scripts

By default the per-epoch transformation is fitted from every surviving source,
each weighted by the median of its weights over epochs, so faint and blended
stars help define the frame that everything is then measured against.
`UseRefSources` fits it from a small clean subset instead, while still solving
the source parameters for every star.

| parameter | default | note |
|---|---|---|
| `UseRefSources` | `false` | switch the whole thing on |
| `RefMagRange` | `[14 16]` | OGLE I window the frame stars come from |
| `RefCompanionRadius` | `5` | drop a candidate with a bright companion this close [pix]. The seeing is about 7.7 pix, so anything inside this radius is thoroughly blended |
| `RefCompanionMaxMag` | `18` | companions fainter than this are ignored |
| `RefCompanionCat` | `[]` | `[X Y Mag]` matrix or the path of an OGLE `.mat`. Worth supplying: the matched source list merges close pairs that OGLE resolves, and without it the cut keeps roughly twice as many stars |

```matlab
ml.scripts.runAstrometryMSc(MScFile, 'EventNum',260058, 'Field','BLG41', ...
    'UseRefSources',true, 'RefMagRange',[14 16], 'RefCompanionRadius',5, ...
    'RefCompanionCat','~/matlab/Lizoo/OGLEdata/OB260058/OB160058.mat');
```

On BLG41 this keeps 36 stars of 594, measured to 4.2/4.0 mas and detected in
every epoch, against 14.0/14.8 mas and 81% for the field as a whole. On BLG01
it keeps only 10, OGLE covering 29% of its sources, which is thin for a six
parameter fit: `Info.RefSrc` records how many survived each criterion.

**Fit** — `runIterDetrendMScArgs`

| parameter | default | note |
|---|---|---|
| `PixPhase` | `false` | needs `UserData`; errors clearly if it is not there |
| `ChromaicHighOrder` | `true` | set `false` if the fit will not converge. The high order airmass columns are strongly collinear and on a sparse matrix the `bicg` solve inside `runIter` stops converging |
| `NiterWeightsBeforeSys` | `10` | first pass |
| `NiterWeightsAfterSys` | `4` | final pass |
| `RunSysRem` | `true` | |

**Coordinates** — needed only for a different target, since they set the
airmass and parallactic angle. `RA`, `Dec` default to the KB260058 cut-out
centre and are taken from `MSc.UserData` when it has them. `GeoPos` defaults to
CTIO.

**Output** — `OutputDir` defaults to `~/KMTdata/Results/AstrometryMSc/`. It is
created and tested for writability before the fit starts, so a bad path fails
in a second rather than after half an hour. Note that `/bigdata3/projects` is
mounted read-only and cannot be written to. Set `PerSourcesTargetPath` to also
write the per-source csv tables; it is empty by default.

## Output

`AstrometryMSc_<event>_<site>_<field>.mat`, holding one `File` struct. About
1 GB, since it carries three objects each with a large `Data`.

| field | contents |
|---|---|
| `IFsys` | the `IterFit` solution. `ParS(1:2,:)` positions, `ParS(3:4,:)` proper motions |
| `IFsysBeforeSysRem` | the first pass solution |
| `Obj` | the `MMS` handed to the final pass |
| `Info` | conversion and run record: `SrcInd` and `EpochInd` back into the input, what each cut removed, colour statistics, `SysRemApplied`, `PixPhaseAvailable` |
| `IndForPhotRefernce` | source nearest pixel (150,150) |
| `FieldCenterDeg` | cut-out centre |
| `ParScalibrated`, `GaiaTable`, `DeltaPM_KMT_GAIA` | empty unless the Gaia step ran and succeeded |

Useful afterwards:

```matlab
[RstdX, RstdY] = IFsys.calculateRstd;                 % already in mas
Mu = 400 .* sqrt(IFsys.ParS(3,:).^2 + IFsys.ParS(4,:).^2);   % mas/yr
ml.scripts.IterFitToPerSourceFormat(IFsys, '~/KMTdata/Events/kmt260058/');
```

## Two things to be aware of

**The proper motions are relative.** `GaiaCalib` is off by default and does not
currently solve for these 300x300 pixel cut-outs: `astrometryCore` finds no
pattern match, the field being 2 by 2 arcmin in the bulge with Gaia about twice
as deep as the KMT source list. So `ParS(3:4,:)` stays in the registered pixel
frame, its orientation on the sky is unknown, and the conversion to mas assumes
exactly 400 mas/pix rather than measuring it.

**Do not feed it a sparse matrix.** The chain was tuned on roughly 83% filled
matrices. Keep the fill above about 70%, which the default cuts do, or both the
sysrem step and the airmass solve fail silently rather than with an error.
