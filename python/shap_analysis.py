# =============================================================================
# SHAP Analysis for LightGBM Top1000 Model
# Run this AFTER compare_models.py
#
# HOW TO RUN:
#   python shap_analysis.py
#
# WHAT IT GENERATES:
#   1. shap_summary_plot.png       - Feature importance ranked by SHAP
#   2. shap_dependence_plot.png    - Top SNP effect on predictions
#   3. shap_waterfall_plot.png     - Single accession explanation
#   4. SHAP_Values.xlsx            - Full SHAP values table
# =============================================================================

import pandas as pd
import numpy as np
import joblib
import shap
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for saving figures
import json
import warnings
warnings.filterwarnings("ignore")

# -------------------------------------------------------
# 1. Load model and features
# -------------------------------------------------------
print("\n[1/6] Loading LightGBM Top1000 model...")

model    = joblib.load("LightGBM_Top1000.pkl")
imputer  = joblib.load("imputer_Top1000.pkl")

with open("features_Top1000.json", "r") as f:
    features = json.load(f)

print(f"  Model loaded: LightGBM")
print(f"  Features: {len(features)} SNPs")

# -------------------------------------------------------
# 2. Load and encode training data
# -------------------------------------------------------
print("\n[2/6] Loading and encoding training data...")

def encode_snp_column(series, snp_name):
    parts = snp_name.split("_")
    if len(parts) < 5:
        return series
    ref = parts[3]
    alt = parts[4]
    het_codes = {"R", "Y", "S", "W", "K", "M"}

    def convert(x):
        if pd.isna(x):
            return np.nan
        x = str(x)
        if x == ref:
            return 0
        elif x == alt:
            return 2
        elif x in het_codes:
            return 1
        return np.nan

    return series.apply(convert)

def encode_dataset(df):
    df = df.copy()
    snp_cols = [c for c in df.columns
                if c not in ["Accession", "SeedThickness"]]
    for col in snp_cols:
        df[col] = encode_snp_column(df[col], col)
    return df

# Load train data
train_df = pd.read_csv("SeedThickness_Top1000_Train.csv")
train_df = encode_dataset(train_df)

X_train = train_df.drop(columns=["Accession", "SeedThickness"])
y_train = train_df["SeedThickness"]
accessions_train = train_df["Accession"].values

# Reindex to match feature list
X_train = X_train.reindex(columns=features)
X_train = X_train.apply(pd.to_numeric, errors="coerce")

# Impute
X_train_imp = pd.DataFrame(
    imputer.transform(X_train),
    columns=features
)

print(f"  Training samples: {X_train_imp.shape[0]}")
print(f"  Features: {X_train_imp.shape[1]}")

# -------------------------------------------------------
# 3. Compute SHAP values
# -------------------------------------------------------
print("\n[3/6] Computing SHAP values (this may take 2-5 minutes)...")

explainer   = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_train_imp)

print(f"  SHAP values shape: {shap_values.shape}")
print(f"  Base value: {explainer.expected_value:.4f}")

# -------------------------------------------------------
# 4. Summary Plot (Feature Importance)
# -------------------------------------------------------
print("\n[4/6] Generating SHAP Summary Plot...")

plt.figure(figsize=(10, 8))
shap.summary_plot(
    shap_values,
    X_train_imp,
    max_display   = 20,
    show          = False,
    plot_size     = (10, 8)
)
plt.title(
    "SHAP Summary Plot — LightGBM Top1000\nSeed Thickness Prediction",
    fontsize = 13,
    fontweight = "bold",
    pad = 15
)
plt.tight_layout()
plt.savefig("shap_summary_plot.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: shap_summary_plot.png")

# -------------------------------------------------------
# 5. Dependence Plot (Top SNP)
# -------------------------------------------------------
print("\n[5/6] Generating SHAP Dependence Plot...")

# Find top SNP by mean absolute SHAP value
mean_abs_shap = np.abs(shap_values).mean(axis=0)
top_snp_idx   = np.argmax(mean_abs_shap)
top_snp_name  = features[top_snp_idx]

print(f"  Top SNP: {top_snp_name}")
print(f"  Mean |SHAP|: {mean_abs_shap[top_snp_idx]:.4f}")

# Find second most important SNP for color interaction
second_snp_idx  = np.argsort(mean_abs_shap)[-2]
second_snp_name = features[second_snp_idx]

plt.figure(figsize=(9, 6))
shap.dependence_plot(
    top_snp_idx,
    shap_values,
    X_train_imp,
    interaction_index = second_snp_idx,
    show  = False,
    ax    = plt.gca()
)
plt.title(
    f"SHAP Dependence Plot: {top_snp_name}\nColoured by {second_snp_name}",
    fontsize = 11,
    fontweight = "bold"
)
plt.xlabel(f"Genotype Encoding: {top_snp_name}\n(0=REF/REF, 1=REF/ALT, 2=ALT/ALT)")
plt.ylabel("SHAP Value (Impact on Seed Thickness Prediction)")
plt.tight_layout()
plt.savefig("shap_dependence_plot.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: shap_dependence_plot.png")

# -------------------------------------------------------
# 6. Waterfall Plot (Single Accession)
# -------------------------------------------------------
print("\n[6/6] Generating SHAP Waterfall Plot...")

# Pick the accession closest to mean seed thickness
y_pred      = model.predict(X_train_imp)
mean_pred   = y_pred.mean()
sample_idx  = int(np.argmin(np.abs(y_pred - mean_pred)))
sample_name = accessions_train[sample_idx]

print(f"  Selected accession: {sample_name}")
print(f"  Predicted: {y_pred[sample_idx]:.4f} mm")
print(f"  Actual:    {y_train.iloc[sample_idx]:.4f} mm")

# Build Explanation object for waterfall
explanation = shap.Explanation(
    values        = shap_values[sample_idx],
    base_values   = explainer.expected_value,
    data          = X_train_imp.iloc[sample_idx].values,
    feature_names = features
)

plt.figure(figsize=(10, 7))
shap.waterfall_plot(
    explanation,
    max_display = 15,
    show        = False
)
plt.title(
    f"SHAP Waterfall Plot — Accession: {sample_name}\n"
    f"Predicted: {y_pred[sample_idx]:.3f} mm  |  "
    f"Actual: {y_train.iloc[sample_idx]:.3f} mm",
    fontsize = 11,
    fontweight = "bold"
)
plt.tight_layout()
plt.savefig("shap_waterfall_plot.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: shap_waterfall_plot.png")

# -------------------------------------------------------
# 7. Export SHAP values to Excel
# -------------------------------------------------------
print("\n[7/7] Saving SHAP values to Excel...")

shap_df = pd.DataFrame(
    shap_values,
    columns = features
)
shap_df.insert(0, "Accession", accessions_train)
shap_df.insert(1, "Actual_SeedThickness", y_train.values)
shap_df.insert(2, "Predicted_SeedThickness", y_pred)

# Summary table: top 20 SNPs by mean |SHAP|
summary_df = pd.DataFrame({
    "SNP"           : features,
    "Mean_Abs_SHAP" : mean_abs_shap,
    "Mean_SHAP"     : shap_values.mean(axis=0)
}).sort_values("Mean_Abs_SHAP", ascending=False).reset_index(drop=True)

summary_df["Rank"] = summary_df.index + 1

with pd.ExcelWriter("SHAP_Values.xlsx", engine="openpyxl") as writer:
    summary_df.to_excel(
        writer, sheet_name="Feature_Importance", index=False
    )
    shap_df.to_excel(
        writer, sheet_name="All_SHAP_Values", index=False
    )

print("  Saved: SHAP_Values.xlsx")

# -------------------------------------------------------
# Summary
# -------------------------------------------------------
print("\n========================================")
print(" SHAP ANALYSIS COMPLETE")
print("  shap_summary_plot.png")
print("  shap_dependence_plot.png")
print("  shap_waterfall_plot.png")
print("  SHAP_Values.xlsx")
print("\n Top 10 Most Important SNPs:")
print(summary_df[["Rank","SNP","Mean_Abs_SHAP"]].head(10).to_string(index=False))
print("========================================\n")