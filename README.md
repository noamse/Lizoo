# Ultra-Trim `runEventBatch` Chain for KMTNet

This repository contains a **minimal subset** of a larger photometry and astrometry pipeline, centered on the script **`+ml/+scripts/runEventBatch.m`**. It has been trimmed to include only the functions directly required for running the batch event reduction workflow.

## Purpose

The code was originally developed to process data from the **Korean Microlensing Telescope Network (KMTNet)**.  
Its main goals are:

- **Automated batch reduction** of microlensing event fields.  
- **Astrometric calibration** against reference catalogs (e.g., Gaia).  
- **Iterative detrending and model fitting** using custom algorithms (`@IterFit`).  
- **PSF photometry pipeline** (`@ImRed`) tuned for crowded stellar fields.  

This version is published for reproducibility and demonstration: it allows researchers to see how events are processed in KMTNet-like datasets, and to reuse core routines in their own workflows.

## Workflow Overview

1. **Batch driver:** `ml.scripts.runEventBatch`  
   Iterates over KMTNet fields/events and launches processing for each.

2. **Field-level processing:** `ml.scripts.processKMTEventField`  
   - Loads KMTNet reference catalogs.  
   - Runs the photometry/astrometry pipeline (`ImRed.runPipe`).  
   - Produces matched source catalogs and detrended light curves.

3. **Astrometry step:** `ml.scripts.runAstrometryField` and `gaiaAstrometryKMT`  
   Aligns KMTNet star catalogs to external references (Gaia DR3).

4. **Iterative detrending:** `ml.scripts.runIterDetrend`  
   Applies the `IterFit` class to remove systematics and refine event light curves.

The chain is modular: each step can be run independently, or managed in bulk via `runEventBatch`.

## External Dependencies

To apply this pipeline on real KMTNet images, you will also need:
- **Image / catalog classes:** `AstroImage`, `AstroCatalog` (not included here).  
- **Reference catalogs:** KMTNet reference catalog, Gaia DR3.  
- **Optional visualization:** DS9 utilities (`ds9.*`), if plotting is enabled.

The repository includes stubs for Telegram notifications (`+ut/sendTelegram*.m`) so the code can run without network access or secrets.

## Usage

1. Clone this repo and add it to your MATLAB path:
   ```matlab
   addpath(genpath('/path/to/NoamFun_runEventBatch_ultratrim'));
   ```

2. Prepare your data directory with KMTNet-format images and reference catalogs.

3. Run a batch:
   ```matlab
   ml.scripts.runEventBatch('RootDir','/path/to/kmtnet/data','Fields',{'BLG01'});
   ```

The output will include reduced catalogs, detrended light curves, and astrometric solutions aligned to Gaia.

## License

MIT (see `LICENSE`).
