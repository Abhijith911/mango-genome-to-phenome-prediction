# 🧬 Genome-to-Phenome Prediction of Seed Thickness in Mango

> **An Explainable Genome-to-Phenome Framework Integrating Genome-Wide Association Study (GWAS), Machine Learning, SHAP Explainability, and Candidate Gene Analysis for Seed Thickness Prediction in Mango (*Mangifera indica* L.)**

---

## 📖 Overview

This repository presents an end-to-end **Genome-to-Phenome Prediction** framework developed for predicting **seed thickness in mango** using genomic information.

The framework integrates:

- 🌱 Genome-Wide Association Study (GWAS)
- 🤖 Machine Learning
- 📊 SHAP Explainability
- 🧬 Candidate Gene Analysis
- 🌐 Interactive Streamlit Web Application

Unlike conventional workflows, this project performs the **train-test split before GWAS**, preventing data leakage and ensuring unbiased model evaluation.

---

## 🚀 Key Features

- Leakage-free GWAS + Machine Learning pipeline
- FarmCPU GWAS using **rMVP**
- Comparison of multiple ML algorithms
- SHAP Explainability
- Candidate Gene Analysis (±50 kb)
- Interactive Streamlit Web Application
- End-to-end reproducible workflow

---

# 📊 Dataset

- **Species:** *Mangifera indica* L.
- **Accessions:** 161
- **SNPs:** 135,079
- **Trait:** Seed Thickness
- **Reference Genome:** CATAS Mindica 2.1

**Dataset Source**

Eltaher et al.

BMC Genomics (2025)

DOI:
https://doi.org/10.1186/s12864-025-11278-6

> The original genotype and phenotype datasets are available through the supplementary material of the above publication and are **not redistributed in this repository**.

---

# 🔬 Methodology

```text
Raw Genotype + Phenotype
           │
           ▼
Train/Test Split (80/20)
           │
           ▼
FarmCPU GWAS (Training Set Only)
           │
           ▼
Rank SNPs by p-value
           │
           ▼
Top SNP Feature Sets
           │
           ▼
SNP Encoding (0 / 1 / 2)
           │
           ▼
Median Imputation
           │
           ▼
Machine Learning
           │
           ▼
Model Evaluation
           │
           ▼
SHAP Explainability
           │
           ▼
Candidate Gene Analysis
           │
           ▼
Streamlit Web Application
```

---

# 🤖 Machine Learning Models

The following algorithms were evaluated:

- LightGBM
- ElasticNet
- XGBoost
- Random Forest
- Support Vector Regression (SVR)

The final Streamlit application deploys:

- ✅ LightGBM (Top1000 SNPs)
- ✅ ElasticNet (Top100 SNPs)

---

# 📈 Best Results

| Model | Feature Set | Test R² | RMSE |
|--------|-------------|---------|------|
| LightGBM | Top1000 | **0.433** | **1.783** |
| ElasticNet | Top100 | **0.411** | **1.817** |

---

# 📊 Explainable AI

Model interpretation was performed using **SHAP (SHapley Additive exPlanations)**.

Generated outputs include:

- SHAP Summary Plot
- SHAP Dependence Plot
- SHAP Waterfall Plot

These visualizations explain how individual SNPs contribute to seed thickness prediction.

---

# 🧬 Candidate Gene Analysis

Significant GWAS loci were explored using ±50 kb windows around lead SNPs.

Genes associated with:

- Disease resistance
- Cytochrome P450 proteins
- F-box proteins
- Transcription factors

were identified and discussed in the accompanying paper.

---

# 🌐 Streamlit Web Application

The application provides six interactive modules:

1. Seed Thickness Prediction
2. GWAS Results
3. Sequence Explorer
4. Candidate Gene Analysis
5. Model Comparison
6. SHAP Explainability

---

# 📂 Repository Structure

```
app/
python/
models/
data/
figures/
results/
paper/
reference/
```

---

# ⚙️ Installation

Clone the repository

```bash
git clone https://github.com/Abhijith911/mango-genome-to-phenome-prediction.git
```

Install dependencies

```bash
pip install -r requirements.txt
```

Launch Streamlit

```bash
streamlit run app/app.py
```

---

# 📄 Research Paper

The complete IEEE journal manuscript is available in the **paper/** directory.

---

# 🛠️ Technologies Used

- Python
- R
- Streamlit
- LightGBM
- SHAP
- Scikit-learn
- rMVP
- XGBoost
- Pandas
- NumPy
- Matplotlib
- OpenPyXL

---

# 👨‍💻 Author

**Abhijith Santhosh**

B.Tech Electronics & Computer Engineering

Saintgits College of Engineering

Kottayam, Kerala, India

GitHub:
https://github.com/Abhijith911

LinkedIn:
https://linkedin.com/in/abhijithsanthosh2005

---

# 📚 Citation

If you use this work for academic purposes, please cite the accompanying publication.

---

## ⭐ Acknowledgements

This work was inspired by the publicly available mango genomic resources published by **Eltaher et al. (BMC Genomics, 2025)** and was developed as part of undergraduate research in Genome-to-Phenome prediction.