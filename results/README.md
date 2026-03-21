# KMTNet-Gaia DR3 Cross-Matched Astrometric and Photometric Catalogue

**Version:** 1.0.0  
**Date:** March 2026  
**Author:** [Your Name / Noam Se]  
**Institution:** Weizmann Institute of Science  
**Contact:** [Your Email Address]  

[![DOI](https://zenodo.org/badge/DOI/[Insert-DOI-Here].svg)](https://doi.org/[Insert-DOI-Here])

---

## 1. Overview
This repository contains a comprehensive astrometric and photometric catalogue of stellar sources in the Galactic Bulge, combining multi-epoch observations from the Korea Microlensing Telescope Network (KMTNet) with the Gaia Data Release 3 (DR3) reference frame. 

The internal KMTNet spatial coordinates and proper motions have been geometrically calibrated to the Gaia DR3 reference frame using an M-estimator Sample Consensus (MSAC) algorithm. This robust 2D affine transformation isolates a high-fidelity "Core Belt" of reference stars to correct for instrumental and atmospheric distortions, ensuring precise kinematic alignment.

For full details on the pipeline, calibration methodology, and error estimation, please refer to:
> **[Insert Title of Your Thesis or Paper]** (2026). *[Insert Journal or University Name]*. Link: [Insert Link if available]

---

## 2. File Structure
The dataset is provided in two primary, machine-readable CSV files:

1. `KMT_Gaia_Crossmatch_Master.csv`: The main catalogue containing individual source kinematics, photometry, and cross-matched Gaia parameters.
2. `KMT_Gaia_Field_Metadata.csv`: A summary table containing the field-level metadata, observational baselines, and the derived Affine calibration coefficients.

---

## 3. Data Dictionary: Master Catalogue (`KMT_Gaia_Crossmatch_Master.csv`)

| Column Name | Units | Data Type | Description |
| :--- | :--- | :--- | :--- |
| `FieldName` | -- | String | KMTNet Bulge field designation (e.g., BLG42). |
| `FieldID` | -- | Integer | Internal processing index for the field. |
| `RA` | deg | Float | Calibrated Equatorial Right Ascension (J2000). |
| `Dec` | deg | Float | Calibrated Equatorial Declination (J2000). |
| `Ref_Epoch` | yr | Float | Reference epoch of the KMTNet field in Julian Years. |
| `Nobs` | -- | Integer | Total number of successful epochs utilized in the astrometric fit. |
| `pmra` | mas/yr | Float | Calibrated KMTNet proper motion in RA ($\mu_{\alpha *} \equiv \mu_\alpha \cos \delta$). |
| `pmra_err` | mas/yr | Float | Empirical KMTNet formal proper motion uncertainty in RA. |
| `pmdec` | mas/yr | Float | Calibrated KMTNet proper motion in Dec ($\mu_\delta$). |
| `pmdec_err` | mas/yr | Float | Empirical KMTNet formal proper motion uncertainty in Dec. |
| `I` | mag | Float | Calibrated apparent magnitude in the Cousins $I$-band. |
| `V-I` | mag | Float | KMTNet $V-I$ color index. |
| `pmra_gaia` | mas/yr | Float | Gaia DR3 reference proper motion in RA. |
| `pmra_gaia_err` | mas/yr | Float | Formal uncertainty of the Gaia DR3 proper motion in RA. |
| `pmdec_gaia` | mas/yr | Float | Gaia DR3 reference proper motion in Dec. |
| `pmdec_gaia_err` | mas/yr | Float | Formal uncertainty of the Gaia DR3 proper motion in Dec. |
| `RpGaia` | mag | Float | Gaia DR3 broad-band $G_{RP}$ magnitude. |
| `Bp-Rp` | mag | Float | Gaia DR3 integrated color index ($G_{BP} - G_{RP}$). |
| `Ggaia` | mag | Float | Gaia DR3 broad-band $G$ magnitude. |
| `rmsX` | mas | Float | Positional root-mean-square (RMS) scatter in the focal plane (X-axis). |
| `rmsY` | mas | Float | Positional root-mean-square (RMS) scatter in the focal plane (Y-axis). |
| `GaiaPlx` | mas | Float | Gaia DR3 absolute stellar parallax. |
| `GaiaPlxErr` | mas | Float | Gaia DR3 formal parallax uncertainty. |

*Note: Missing or non-applicable values are represented as `NaN`.*

---

## 4. Usage Notes
* **Proper Motion Formulation:** All `pmra` columns in this dataset represent the true arc projection on the sky ($\mu_{\alpha *} = \mu_\alpha \cos \delta$), strictly matching the Gaia DR3 convention.
* **Epochs:** Because KMTNet fields were initiated at slightly different times, the `Ref_Epoch` varies slightly by field. Time-dependent kinematic propagation should account for this specific Julian Year rather than assuming a universal J2016.0 epoch.

---

## 5. Acknowledgements & Citation
If you use this dataset in your research, please cite the original publication:
```bibtex
@phdthesis{se2026kmtnet,
  author       = {[Se, Noam]}, 
  title        = {[Insert Exact Thesis/Paper Title Here]},
  school       = {Weizmann Institute of Science},
  year         = 2026,
  doi          = {[Insert-DOI-Here]}
}
