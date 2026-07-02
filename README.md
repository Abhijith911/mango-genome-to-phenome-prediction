# 🧬 Genome-to-Phenome Prediction of Seed Thickness in Mango

**An Explainable Genome-to-Phenome Framework Integrating GWAS, Machine Learning, SHAP Explainability, and Candidate Gene Analysis for Seed Thickness Prediction in *Mangifera indica* L.**

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)]()
[![R](https://img.shields.io/badge/R-4.x-blue)]()
[![Streamlit](https://img.shields.io/badge/Streamlit-App-FF4B4B)]()

---

## 📖 Overview

This repository presents an end-to-end **genome-to-phenome prediction framework** for **seed thickness in mango**, built around a single non-negotiable design principle: the train/test split happens **before** any marker discovery, not after.

A large share of published GWAS-assisted machine learning pipelines run association testing on the full sample collection and only split into train/test afterward — which means the "held-out" test set has already influenced which SNPs get selected as features. This repository avoids that leakage entirely, giving performance numbers that reflect genuine generalization to unseen accessions rather than an inflated, leakage-boosted score.

The framework integrates:

- 🌱 Genome-Wide Association Study (GWAS) via FarmCPU
- 🤖 Comparative Machine Learning (5 algorithms × 6 feature sets)
- 📊 SHAP Explainability
- 🧬 Candidate Gene Analysis (±50 kb)
- 🌐 Interactive 6-tab Streamlit Web Application

---

## 🚀 Key Features

- **Leakage-free design** — GWAS restricted to the training partition only
- **FarmCPU GWAS** via the `rMVP` R package, with 3 PCs for population-structure correction
- **Five ML algorithms** benchmarked across six SNP feature-set sizes (30 combinations total)
- **SHAP explainability** — summary, dependence, and waterfall plots via TreeExplainer
- **Candidate gene analysis** within ±50 kb of significant loci, annotated against the CATAS Mindica 2.1 reference
- **Interactive Streamlit app** exposing every stage of the pipeline to non-specialist users
- **Fully reproducible** — fixed random seed (42) for the train/test partition

---

## 🔬 Pipeline Architecture

```
 Mango Accessions (161)     Phenotypic Data (Seed Thickness)     Genotypic Data (135,079 SNPs)
      │                              │                                    │
      └──────────────────────────────┼────────────────────────────────────┘
                                      ▼
                    Leakage-Free Train/Test Split (128 / 33)
                                      │
                                      ▼
                    GWAS with FarmCPU (Training Set Only)
                                      │
                                      ▼
                    Rank SNPs by p-value → Significant SNP ID
                                      │
                                      ▼
              Nested Feature Sets: Top-10 / 100 / 500 / 1000 / 2000 / Significant
                                      │
                                      ▼
                    SNP Encoding (0/1/2) + Median Imputation
                                      │
                                      ▼
        Machine Learning: LightGBM · ElasticNet · Random Forest · XGBoost · SVR
                                      │
                                      ▼
           Model Evaluation (5-fold CV + Independent Test Set)
                                      │
                                      ▼
                 SHAP Explainability (Summary / Dependence / Waterfall)
                                      │
                                      ▼
                    Candidate Gene Analysis (±50 kb window)
                                      │
                                      ▼
                 Interactive Streamlit Web Application (6 tabs)
```

---

## 📊 Dataset

| Attribute | Value |
|---|---|
| Species | *Mangifera indica* L. |
| Accessions | 161 |
| SNPs | 135,079 (high-confidence, MAF > 0.05) |
| Trait | Seed thickness (mm) |
| Phenotype range | 10.15 – 23.52 mm (mean 15.8, SD 2.08) |
| Reference genome | CATAS Mindica 2.1 (`GCF_011075055.1`) |
| Sequencing | Illumina NovaSeq 6000, 150 bp paired-end |
| Source | Eltaher et al., *BMC Genomics* (2025), [DOI: 10.1186/s12864-025-11278-6](https://doi.org/10.1186/s12864-025-11278-6) |

> The original genotype and phenotype data are available through the supplementary material of the source publication and are **not redistributed** in this repository.

---

## 🧬 GWAS Results

FarmCPU, run exclusively on the 128-accession training partition with 3 principal components as covariates, identified **two genome-wide significant SNPs**:

| SNP | Chromosome | Position | MAF | Effect (mm) | P-value |
|---|---|---|---|---|---|
| `NC_058138.1_2620063_T_C` | 2 | 2,620,063 | 0.202 | −1.990 | 1.12 × 10⁻⁷ |
| `NC_058147.1_2507594_C_T` | 11 | 2,507,594 | 0.205 | +1.680 | 2.83 × 10⁻⁷ |

The opposing effect directions at comparable allele frequencies suggest seed thickness in this panel is shaped by at least two loci of moderate, oppositely directed effect rather than a single dominant gene.

**Notable candidate genes within ±50 kb:**

- **Chr 2** — F-box domain genes (SCF ubiquitin ligase complex), a disease-resistance protein, tryptophan–tRNA ligase
- **Chr 11** — `KAN4` (KANADI4) transcription factor, multiple cytochrome P450 monooxygenases (phenylpropanoid/lignin biosynthesis), an allene oxide synthase, a DOF-family zinc-finger transcription factor

These genes are plausible — though not experimentally validated — regulators of seed-coat lignification and integument development.

---

## 🤖 Machine Learning Results

Five regression algorithms were evaluated across six SNP feature-set sizes (30 combinations total), with 5-fold cross-validation on the training set and a single final evaluation on the held-out test set.

| Feature Set | Model | Train R² | CV R² | Test R² | RMSE (mm) | MAE (mm) |
|---|---|---|---|---|---|---|
| Top-1000 | **LightGBM** | 0.982 | 0.329 | **0.433** | **1.783** | 1.385 |
| Top-100 | **ElasticNet** | 0.698 | 0.405 | 0.411 | 1.817 | 1.435 |
| Top-500 | LightGBM | 0.978 | 0.286 | 0.393 | 1.843 | 1.517 |
| Top-1000 | XGBoost | 0.989 | 0.236 | 0.386 | 1.855 | 1.460 |
| Top-100 | Random Forest | 0.917 | 0.450 | 0.376 | 1.869 | 1.539 |
| Top-100 | SVR | 0.699 | 0.420 | 0.303 | 1.976 | 1.482 |

**Key takeaway:** LightGBM (Top-1000) wins on raw test accuracy but shows a large train–CV gap, indicating overfitting risk. ElasticNet (Top-100) achieves nearly identical test accuracy with far better CV–test agreement, making it the more *reliable* model when only training-time cross-validation is available to guide model selection.

Both models are deployed in the Streamlit app for side-by-side comparison.

---

## 📊 Explainable AI (SHAP)

Model interpretation was performed on the best-performing model (LightGBM, Top-1000) using the TreeExplainer algorithm.

- **Summary Plot** — ranks SNPs by mean absolute SHAP value and shows the direction of their effect
- **Dependence Plot** — shows how genotype dosage (0/1/2) at a given SNP relates to its SHAP contribution
- **Waterfall Plot** — decomposes a single accession's prediction into per-SNP contributions

Interestingly, the two most influential SNPs by SHAP value (`NC_058152.1_13686895_C_T` and `NC_058138.1_2664305_G_A`) do **not** coincide with the genome-wide significant FarmCPU markers — illustrating that SHAP surfaces sub-threshold markers with real predictive value inside the multivariate model, complementing rather than replacing the univariate GWAS signal.

---

## 🌐 Streamlit Web Application

The deployed application provides six interactive modules:

| Tab | Description |
|---|---|
| **LightGBM Prediction** | Predicts seed thickness from a Top-1000 SNP genotype |
| **ElasticNet Prediction** | Predicts seed thickness from a Top-100 SNP genotype |
| **GWAS Results** | Manhattan plot, QQ plot, PCA scatter, significant-SNP table |
| **Sequence Explorer** | ±10 bp reference context around each significant SNP |
| **Candidate Gene Explorer** | Full annotated gene list within ±50 kb, filterable by region |
| **SHAP Explainability** | Summary, dependence, and waterfall plots + sortable importance table |

```bash
streamlit run app/app.py
```

---

## ⚙️ Installation

```bash
# Clone the repository
git clone https://github.com/Abhijith911/mango-genome-to-phenome-prediction.git
cd mango-genome-to-phenome-prediction

# Install dependencies
pip install -r requirements.txt

# Launch the app
streamlit run app/app.py
```

---

## 📂 Repository Structure

```
mango-genome-to-phenome-prediction/
├── app/            # Streamlit web application
├── python/         # ML pipeline, SHAP, candidate gene scripts
├── models/         # Trained model artifacts
├── data/           # Processed data (raw genotype/phenotype not redistributed)
├── figures/        # Generated plots (Manhattan, QQ, SHAP, PCA, etc.)
├── results/        # Model performance tables, GWAS summary statistics
├── paper/          # Full IEEE manuscript
└── reference/       # CATAS Mindica 2.1 reference files (GTF, GBFF, FASTA)
```

---

## 🛠️ Technologies Used

`Python` · `R` · `Streamlit` · `LightGBM` · `XGBoost` · `SHAP` · `Scikit-learn` · `rMVP` (FarmCPU) · `Pandas` · `NumPy` · `Matplotlib` · `OpenPyXL` · `pyfaidx`

---

## 📄 Research Paper

The complete manuscript, *"An Explainable Genome-to-Phenome Framework Integrating Genome-Wide Association Study, Machine Learning, and Candidate Gene Analysis for Seed Thickness Prediction in Mango (Mangifera indica L.)"*, is available in the [`paper/`](./paper) directory.

### Citation

Please cite the original dataset publication:

> Eltaher, S., Li, J., Freeman, B., Singh, S., & Ali, G. S. (2025). A genome-wide association study identified SNP markers and candidate genes associated with morphometric fruit quality traits in mangoes. *BMC Genomics*, 26(1), 120.

---

## 👨‍💻 Author

**Abhijith Santhosh**
B.Tech Electronics & Computer Engineering, Saintgits College of Engineering, Kottayam, Kerala, India

[GitHub](https://github.com/Abhijith911) · [LinkedIn](https://linkedin.com/in/abhijithsanthosh2005)

---

## ⭐ Acknowledgements

This work builds on the publicly available mango genomic resources published by **Eltaher et al. (BMC Genomics, 2025)** and was developed as part of undergraduate research in genome-to-phenome prediction. Reference genome data (CATAS Mindica 2.1) accessed via NCBI.
