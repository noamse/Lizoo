# KMTNet-Gaia DR3 Cross-Matched Astrometric and Photometric Catalogue

**Version:** 1.0.0  
**Date:** March 2026  
**Author:** Noam Segev  
**Institution:** Weizmann Institute of Science  
**Contact:** [noam.segev@weizmann.ac.il](mailto:noam.segev@weizmann.ac.il)  

---

## 1. Overview
This repository contains the primary data products of the **"KMT Astrometry"** project (also known internally as the **Michael Lansing** project). [cite_start]This work represents a core component of my PhD thesis at the Weizmann Institute of Science[cite: 1].

The catalogue provides a comprehensive astrometric and photometric record of stellar sources in the Galactic Bulge. [cite_start]It successfully bridges high-cadence, seeing-limited ground observations from the **Korea Microlensing Telescope Network (KMTNet)** with the precision of the **Gaia Data Release 3 (DR3)** reference frame[cite: 1].

### Methodology & Calibration
[cite_start]To achieve high-fidelity kinematic alignment, we utilized a robust 2D affine transformation powered by an **M-estimator Sample Consensus (MSAC)** algorithm[cite: 1]. Key features of the pipeline include:
* [cite_start]**The "Core Belt":** A high-fidelity reference set of stars ($14.5 < I < 17.5$ mag) used to isolate and correct for instrumental and atmospheric distortions[cite: 1].
* [cite_start]**Astrometric Precision:** The internal KMTNet spatial coordinates and proper motions are geometrically calibrated to the Gaia DR3 frame to ensure sub-milliarcsecond precision[cite: 1].
* **Collaboration:** This pipeline and the resulting catalogue were developed in collaboration with **Yossi Shvartzvald** and **Krzysztof Rybicki (Steve)**.

[cite_start]For full details on the pipeline, calibration methodology, and error estimation, please refer to the thesis [cite: 1] or the associated publication:
> **Towards sub-milliarcsecond astrometric precision using high-cadence seeing-limited imaging** (2026). *Monthly Notices of the Royal Astronomical Society (MNRAS)*. [cite_start][https://doi.org/10.1093/mnras/staf2234] [cite: 1]

---

## 2. File Structure
[cite_start]The dataset is provided in two primary, machine-readable CSV files[cite: 1]:

1. [cite_start]`data/KMT_Gaia_Crossmatch_Master.csv`: Individual source kinematics, photometry, and Gaia cross-match parameters[cite: 1].
2. [cite_start]`data/KMT_Gaia_Field_Metadata.csv`: Field-level metadata, summary statistics, and derived Affine calibration coefficients[cite: 1].

---

## 3. Data Dictionaries

### 3.1 Master Source Catalogue (`KMT_Gaia_Crossmatch_Master.csv`)

| Column Name | Units | Description |
| :--- | :--- | :--- |
| `FieldName` | -- | [cite_start]KMTNet Bulge field designation (e.g., BLG42)[cite: 1]. |
| `FieldID` | -- | [cite_start]Internal processing index for the field[cite: 1]. |
| `RA`, `Dec` | deg | [cite_start]Calibrated Equatorial coordinates (J2000)[cite: 1]. |
| `Ref_Epoch` | yr | [cite_start]Reference epoch of the field in Julian Years[cite: 1]. |
| `Nobs` | -- | [cite_start]Total number of epochs utilized in the astrometric fit[cite: 1]. |
| `pmra` | mas/yr | [cite_start]Calibrated KMTNet proper motion ($\mu_{\alpha *} \equiv \mu_\alpha \cos \delta$)[cite: 1]. |
| `pmra_err` | mas/yr | Empirical formal proper motion uncertainty in RA[cite: 1]. |
| `pmdec` | mas/yr | Calibrated KMTNet proper motion ($\mu_\delta$)[cite: 1]. |
| `pmdec_err` | mas/yr | Empirical formal proper motion uncertainty in Dec[cite: 1]. |
| `I` | mag | Calibrated apparent magnitude in Cousins $I$-band[cite: 1]. |
| `V-I` | mag | KMTNet $V-I$ colour index[cite: 1]. |
| `pmra_gaia`, `pmdec_gaia` | mas/yr | Gaia DR3 reference proper motions[cite: 1]. |
| `RpGaia`, `Ggaia` | mag | Gaia DR3 broad-band magnitudes[cite: 1]. |
| `Bp-Rp` | mag | Gaia DR3 integrated colour index[cite: 1]. |
| `rmsX`, `rmsY` | mas | Positional RMS scatter in the focal plane[cite: 1]. |
| `GaiaPlx`, `GaiaPlxErr` | mas | Gaia DR3 parallax and uncertainty[cite: 1]. |

### 3.2 Field Metadata Summary (`KMT_Gaia_Field_Metadata.csv`)

| Column Name | Units | Description |
| :--- | :--- | :--- |
| `CenterRA`, `CenterDec` | deg | Coordinates of the KMTNet field centre[cite: 1]. |
| `Nepochs` | -- | Total number of available epochs in the raw data[cite: 1]. |
| `Nused` | -- | Number of epochs used in the final astrometric fit[cite: 1]. |
| `SpatialCalibCoefs` | -- | 6 affine coefficients $[a, b, c, d, t_x, t_y]$ for spatial calibration[cite: 1]. |
| `PMCalibCoefs` | -- | 6 affine coefficients for proper motion calibration[cite: 1]. |
| `KMTGAIAPMRSTDRA` | mas/yr | Robust std. dev. of RA residuals against Gaia[cite: 1]. |
| `KMTGAIAPMRSTDDec` | mas/yr | Robust std. dev. of Dec residuals against Gaia[cite: 1]. |
| `SourceInd` | -- | Range of indices in the Master Catalogue belonging to this field[cite: 1]. |

---

## 4. Usage Example (Python)
To link the metadata to the master catalogue and filter by field[cite: 1]:

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

# 1. Load the Master Catalogue
file_path = 'data/KMT_Gaia_Crossmatch_Master.csv'
df = pd.read_csv(file_path)

# 2. Apply filters (e.g., bright stars with low PM uncertainty)
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

lims = [np.min([plt.xlim(), plt.ylim()]), np.max([plt.xlim(), plt.ylim()])]
plt.plot(lims, lims, 'r--', alpha=0.75, zorder=0, label='1:1 Reference')

plt.title(f'Proper Motion Comparison: KMTNet vs. Gaia DR3 ($\mu_\delta$)\n'
          f'Filters: $I < 16.5$, $\sigma_{{\mu}} < 0.2$ mas/yr', fontsize=12)
plt.xlabel('Gaia DR3 $\mu_\delta$ [mas yr$^{-1}$]', fontsize=11)
plt.ylabel('Calibrated KMTNet $\mu_\delta$ [mas yr$^{-1}$]', fontsize=11)

stats_text = f'N = {len(filtered_df):,}\nPearson $r$ = {corr:.4f}\nRMS Diff = {rms_diff:.3f} mas/yr'
plt.gca().text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, 
               fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.legend(loc='lower right')
plt.grid(True, linestyle=':', alpha=0.6)
plt.tight_layout()
plt.show()
