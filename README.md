# Lizoo: KMTNet Event Batch Reduction Pipeline

**Lizoo** is a MATLAB workflow for batch processing of event fields from the Korean Microlensing Telescope Network (KMTNet). It performs image-level photometry, astrometric calibration, and iterative detrending to produce science-grade catalogs and light curves suitable for downstream analysis.

> **Note.** This repository contains a **slim package** of the code used and described in the accompanying methods paper: https://arxiv.org/abs/2507.11615

## What is Lizoo?

Lizoo provides a practical, end-to-end path from KMTNet images to calibrated, detrended time-series and matched source catalogs. It coordinates:
- **PSF photometry and image processing** via an `ImRed` pipeline (PSF construction, iterative PSF photometry, metadata and product generation).
- **Astrometric calibration** against external reference catalogs (e.g., Gaia DR3).
- **Iterative detrending / model refinement** using the `IterFit` framework for removing systematics in crowded fields.

The code layout follows MATLAB packages and classes, with the top-level driver orchestrating field-by-field processing.

## Workflow Overview

1. **Batch driver** — `ml.scripts.runEventBatch`  
   Iterates over KMTNet fields and invokes the full reduction for each field.

2. **Field processing** — `ml.scripts.processKMTEventField`  
   Loads reference catalogs, calls `ImRed.runPipe` for PSF photometry, and prepares matched catalogs and photometric series.

3. **Astrometry** — `ml.scripts.runAstrometryField`, `ml.scripts.gaiaAstrometryKMT`  
   Aligns source catalogs to Gaia to establish a common astrometric frame.

4. **Detrending** — `ml.scripts.runIterDetrend`  
   Applies `IterFit` to model and mitigate systematics, refining light curves for scientific analysis.

Each stage can be run independently or scheduled end-to-end via the batch driver.

## Dependencies

- MATLAB (R2021a or newer recommended)
- **Classes**: `AstroImage`, `AstroCatalog` (not bundled; required by `ImRed` for IO/representation)
- **Catalogs**: KMTNet reference catalog and Gaia DR3 (paths configured in your environment)
- **Optional**: DS9 utilities (`ds9.*`) if you enable visualization
- **Notifications**: Telegram functions are **stubbed** in `+ut/` to avoid external secrets

## Quick Start

1. Add the repository to your MATLAB path (the folder containing `NoamFun/`):
   ```matlab
   addpath(genpath('/path/to/lizoo'));
   ```

2. Prepare input data in KMTNet format (images and reference catalogs).

3. Run the batch processor:
   ```matlab
   ml.scripts.runEventBatch('RootDir','/path/to/kmtnet/data','Fields',{'BLG01'});
   ```

Outputs include reduced source catalogs, detrended light curves, and astrometric solutions aligned to Gaia.

## Reproducibility and Scope

This distribution focuses on the **core methods** used for the KMTNet pipeline as described in the paper above. Some components not required for the documented workflow have been omitted to keep the package compact and easy to audit.

## Citation

If Lizoo contributes to your work, please cite the associated paper:  
- Segev, N. et al., *Lizoo: KMTNet Event Batch Reduction Pipeline*, arXiv:2507.11615.

## License

MIT (see `LICENSE`).

## Contact

Issues and feature requests are welcome via GitHub. For scientific inquiries related to the KMTNet processing described in the paper, please refer to the corresponding author listed on the arXiv manuscript.
