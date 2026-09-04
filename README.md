# ProFPS — a plasma proteomic signature of frailty pace

Code and model coefficients for **"A Proteomic Signature of Frailty Pace Captures
Functional Aging and Predicts Age-Related Diseases and Mortality"** (manuscript
under review).

ProFPS is a 112-protein plasma signature trained to predict the **pace of frailty** —
the rate of longitudinal decline in body weight, gait and handgrip strength — in the
Framingham Heart Study (FHS) Offspring cohort, and evaluated without refitting in
FHS Generation 3 and the UK Biobank.

## Contents

| File | Description |
|---|---|
| `train_ProFPS.R` | Full training pipeline: Olink QC → imputation and inverse-normal transformation → participant-specific slopes (LME) → frailty pace (PCA) → elastic net → final model |
| `ProFPS_112_coefficients.csv` | The published model: 112 non-zero protein coefficients, keyed by UniProt accession and Olink assay name |

## Computing ProFPS in your own data

**No code from this repository is required to apply ProFPS.** The signature is a
weighted sum of protein levels, so the coefficient table is sufficient: multiply
each protein by its coefficient and sum across the 112 proteins.

Protein values must be preprocessed the same way as in training: Olink NPX,
k-nearest-neighbour imputation of missing values (k ≈ √n), then a rank-based
inverse-normal transformation applied per protein. Effect estimates in the paper
correspond to a 1-SD higher ProFPS after within-cohort standardization.

The 112 coefficients are also reported in Supplementary Table 1 of the manuscript.

## Software

Analyses were run in R version 4.5.2 with the following packages:
`glmnet` 4.1-10, `nlme` 3.1-168, `caret` 7.0-1, `haven` 2.5.5, `tidyr` 1.3.1
and `dplyr` 1.1.4.

## Data availability

Individual-level data are not included in this repository and cannot be
redistributed by the authors.

- **Framingham Heart Study (Offspring and Gen 3)**: available upon request through
  the NHLBI Biologic Specimen and Data Repository Information Coordinating Center
  (BioLINCC), <https://biolincc.nhlbi.nih.gov/studies/framoffspring/> and
  <https://biolincc.nhlbi.nih.gov/studies/gen3/>.
- **UK Biobank**: available through the application process described at
  <https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/>. This research was
  conducted under UK Biobank Application Number 42614.

## Citation

The associated manuscript is currently under review. This section will be updated
with the full citation and DOI once the paper is published.

## License

MIT — see [LICENSE](LICENSE).
