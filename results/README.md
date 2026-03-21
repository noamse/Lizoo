# KMTNet-Gaia DR3 Cross-Matched Astrometric and Photometric Catalogue

**Version:** 1.0.0  
**Date:** March 2026  
**Author:** Noam Segev  
**Institution:** Weizmann Institute of Science  
**Contact:** [noam.segev@weizmann.ac.il](mailto:noam.segev@weizmann.ac.il)  

---

## 1. Overview
This repository contains a comprehensive astrometric and photometric catalogue of stellar sources in the Galactic Bulge, combining multi-epoch observations from the Korea Microlensing Telescope Network (KMTNet) with the Gaia Data Release 3 (DR3) reference frame. 

The internal KMTNet spatial coordinates and proper motions have been geometrically calibrated to the Gaia DR3 reference frame using an M-estimator Sample Consensus (MSAC) algorithm. This robust 2D affine transformation isolates a high-fidelity **"Core Belt"** of reference stars ($14.5 < I < 17.5$ mag) to correct for instrumental and atmospheric distortions, ensuring precise kinematic alignment.

For full details on the pipeline, calibration methodology, and error estimation, please refer to:
> **Towards sub-milliarcsecond astrometric precision using high-cadence seeing-limited imaging** (2026). *Monthly Notices of the Royal Astronomical Society (MNRAS)*. [https://doi.org/10.1093/mnras/staf2234]

---

## 2. File Structure
The dataset is provided in two primary, machine-readable CSV files:

1. `data/KMT_Gaia_Crossmatch_Master.csv`: Individual source kinematics, photometry, and Gaia cross-match parameters.
2. `data/KMT_Gaia_Field_Metadata.csv`: Field-level metadata, summary statistics, and derived Affine calibration coefficients.

---

## 3. Data Dictionaries

### 3.1 Master Source Catalogue (`KMT_Gaia_Crossmatch_Master.csv`)

| Column Name | Units | Description |
| :--- | :--- | :--- |
| `FieldName` | -- | KMTNet Bulge field designation (e.g., BLG42). |
| `FieldID` | -- | Internal processing index for the field. |
| `RA`, `Dec` | deg | Calibrated Equatorial coordinates (J2000). |
| `Ref_Epoch` | yr | Reference epoch of the field in Julian Years. |
| `Nobs` | -- | Total number of epochs utilized in the astrometric fit. |
| `pmra` | mas/yr | Calibrated KMTNet proper motion ($\mu_{\alpha *} \equiv \mu_\alpha \cos \delta$). |
| `pmra_err` | mas/yr | Empirical formal proper motion uncertainty in RA. |
| `pmdec` | mas/yr | Calibrated KMTNet proper motion ($\mu_\delta$). |
| `pmdec_err` | mas/yr | Empirical formal proper motion uncertainty in Dec. |
| `I` | mag | Calibrated apparent magnitude in Cousins $I$-band. |
| `V-I` | mag | KMTNet $V-I$ colour index. |
| `pmra_gaia`, `pmdec_gaia` | mas/yr | Gaia DR3 reference proper motions. |
| `RpGaia`, `Ggaia` | mag | Gaia DR3 broad-band magnitudes. |
| `Bp-Rp` | mag | Gaia DR3 integrated colour index. |
| `rmsX`, `rmsY` | mas | Positional RMS scatter in the focal plane. |
| `GaiaPlx`, `GaiaPlxErr` | mas | Gaia DR3 parallax and uncertainty. |

### 3.2 Field Metadata Summary (`KMT_Gaia_Field_Metadata.csv`)

| Column Name | Units | Description |
| :--- | :--- | :--- |
| `CenterRA`, `CenterDec` | deg | Coordinates of the KMTNet field centre. |
| `Nepochs` | -- | Total number of available epochs in the raw data. |
| `Nused` | -- | Number of epochs used in the final astrometric fit. |
| `SpatialCalibCoefs` | -- | 6 affine coefficients $[a, b, c, d, t_x, t_y]$ for spatial calibration. |
| `PMCalibCoefs` | -- | 6 affine coefficients for proper motion calibration. |
| `KMTGAIAPMRSTDRA` | mas/yr | Robust std. dev. of RA residuals against Gaia. |
| `KMTGAIAPMRSTDDec` | mas/yr | Robust std. dev. of Dec residuals against Gaia. |
| `SourceInd` | -- | Range of indices in the Master Catalogue belonging to this field. |

---

## 4. Usage Example (Python)
To link the metadata to the master catalogue and filter by field:

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

# 1. Load the Master Catalogue
# Ensure the path matches your repository structure
file_path = 'data/KMT_Gaia_Crossmatch_Master.csv'
df = pd.read_csv(file_path)

# 2. Apply requested filters
# I < 16.5 and proper motion uncertainty < 0.2 mas/yr
mask = (df['I'] < 16.5) & (df['pmdec_err'] < 0.2)
filtered_df = df[mask].dropna(subset=['pmdec', 'pmdec_gaia'])

# 3. Extract data for plotting
kmt_pmdec = filtered_df['pmdec']
gaia_pmdec = filtered_df['pmdec_gaia']

# 4. Calculate Statistics
corr, _ = pearsonr(kmt_pmdec, gaia_pmdec)
rms_diff = np.sqrt(np.mean((kmt_pmdec - gaia_pmdec)**2))

# 5. Create the Plot
plt.figure(figsize=(8, 8))
plt.scatter(gaia_pmdec, kmt_pmdec, s=5, alpha=0.4, color='midnightblue', label='Stellar Sources')

# Add a 1:1 reference line
lims = [
    np.min([plt.xlim(), plt.ylim()]),  # min of both axes
    np.max([plt.xlim(), plt.ylim()]),  # max of both axes
]
plt.plot(lims, lims, 'r--', alpha=0.75, zorder=0, label='1:1 Reference')

# Formatting
plt.title(f'Proper Motion Comparison: KMTNet vs. Gaia DR3 ($\mu_\delta$)\n'
          f'Filters: $I < 16.5$, $\sigma_{{\mu}} < 0.2$ mas/yr', fontsize=12)
plt.xlabel('Gaia DR3 $\mu_\delta$ [mas yr$^{-1}$]', fontsize=11)
plt.ylabel('Calibrated KMTNet $\mu_\delta$ [mas yr$^{-1}$]', fontsize=11)

# Add text box with statistics
stats_text = f'N = {len(filtered_df):,}\nPearson $r$ = {corr:.4f}\nRMS Diff = {rms_diff:.3f} mas/yr'
plt.gca().text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, 
               fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.legend(loc='lower right')
plt.grid(True, linestyle=':', alpha=0.6)
plt.tight_layout()

# Save and Show
plt.savefig('pmdec_comparison.png', dpi=300)
plt.show()
