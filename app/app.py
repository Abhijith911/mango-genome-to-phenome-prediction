import streamlit as st
import pandas as pd
import numpy as np
import joblib
from sklearn.impute import SimpleImputer
import json
import os

# ==========================
# LOAD FEATURE LISTS
# ==========================
with open("features_Top1000.json", "r") as f:
    lgbm_features = json.load(f)

train2 = pd.read_csv("SeedThickness_Top100_Train.csv")
elastic_features = [c for c in train2.columns if c not in ["Accession", "SeedThickness"]]

# ==========================
# PAGE CONFIG
# ==========================
st.set_page_config(
    page_title="Mango Seed Thickness Predictor",
    page_icon="🌱",
    layout="wide"
)

st.title("🌱 Mango Seed Thickness Prediction System")

# ==========================
# LOAD MODELS
# ==========================
@st.cache_resource
def load_models():
    lgbm           = joblib.load("LightGBM_Top1000.pkl")
    elastic        = joblib.load("ElasticNet_Top100.pkl")
    lgbm_imputer   = joblib.load("imputer_Top1000.pkl")
    elastic_imputer= joblib.load("imputer_Top100.pkl")
    return lgbm, elastic, lgbm_imputer, elastic_imputer

lgbm_model, elastic_model, lgbm_imputer, elastic_imputer = load_models()

# ==========================
# LOAD PHENOTYPE + GENOTYPE
# ==========================
@st.cache_data
def load_data():
    snp_df = pd.read_excel("Table S2.xlsx", header=2)
    snp_transposed = snp_df.set_index("rs#").T
    snp_transposed.index = snp_transposed.index.astype(str).str.strip().str.lower()

    pheno_df = pd.read_excel(
        "Supplementary Tables.xlsx",
        sheet_name=" Table S1",
        header=2
    )
    pheno_df = pheno_df.set_index("Accession ID")
    pheno_df.index = pheno_df.index.astype(str).str.strip().str.lower()
    return snp_transposed, pheno_df

snp_transposed, pheno_df = load_data()
varieties = sorted(set(snp_transposed.index) & set(pheno_df.index))

# Auto-detect seed thickness column
SEED_COL = next(
    (c for c in pheno_df.columns if "seed" in c.lower() and "thick" in c.lower()),
    None
)
if SEED_COL is None:
    st.error("❌ Could not find seed thickness column in Table S1.")
    st.stop()

# ==========================
# LOAD GWAS & ANNOTATIONS
# ==========================
@st.cache_data
def load_gwas_data():
    signals = pd.read_csv("SeedThickness.FarmCPU_signals.csv")
    signals.columns = ["SNP", "Chr", "Pos", "MAF", "Effect", "SE", "Pvalue"]
    signals["Chrom"] = signals["SNP"].apply(lambda x: "_".join(x.split("_")[:2]))
    return signals

@st.cache_data
def load_annotations():
    if os.path.exists("Gene_Annotations.xlsx"):
        return pd.read_excel("Gene_Annotations.xlsx", sheet_name="Candidate_Genes")
    return None

signals_df   = load_gwas_data()
annotations_df = load_annotations()

# ==========================
# SNP ENCODING HELPERS
# ==========================
def encode_snp_column(series, snp_name):
    parts = snp_name.split("_")
    if len(parts) < 5:
        return series
    ref = parts[3]
    alt = parts[4]
    het = {"R","Y","S","W","K","M"}

    def convert(x):
        if pd.isna(x): return np.nan
        x = str(x)
        if x == ref:        return 0
        elif x == alt:      return 2
        elif x in het:      return 1
        return np.nan

    return series.apply(convert)

def decode_genotype(numeric_val, snp_name):
    """Convert 0/1/2 back to biological genotype string"""
    parts = snp_name.split("_")
    if len(parts) < 5:
        return "N/A"
    ref = parts[3]
    alt = parts[4]
    mapping = {0: f"{ref}/{ref}", 1: f"{ref}/{alt}", 2: f"{alt}/{alt}"}
    try:
        return mapping.get(int(numeric_val), "N/A")
    except:
        return "N/A"

def get_raw_genotype(variety, snp_name):
    """Get raw IUPAC genotype call for a variety at a SNP"""
    if snp_name in snp_transposed.columns:
        return str(snp_transposed.loc[variety, snp_name])
    return "N/A"

def get_sample(variety, feature_list):
    raw = snp_transposed.loc[[variety]].copy()
    raw = raw.reindex(columns=feature_list)
    encoded = {}
    for col in feature_list:
        if col in raw.columns:
            encoded[col] = encode_snp_column(raw[col], col)
        else:
            encoded[col] = np.nan
    df = pd.DataFrame(encoded, index=raw.index)
    df = df.reindex(columns=feature_list)
    return df

def preprocess(df, feature_list, fitted_imputer):
    df = df.reindex(columns=feature_list)
    df = df.apply(pd.to_numeric, errors="coerce")
    arr = fitted_imputer.transform(df)
    return pd.DataFrame(arr, columns=feature_list)

# ==========================
# GENOME SEQUENCE HELPER
# ==========================
@st.cache_resource
def load_genome():
    try:
        from pyfaidx import Fasta
        fna_file = "GCF_011075055.1_CATAS_Mindica_2.1_genomic.fna"
        if os.path.exists(fna_file):
            return Fasta(fna_file)
        return None
    except Exception:
        return None

genome = load_genome()

def get_sequence_context(chrom, pos, window=10):
    """Extract +-window bp around SNP position from reference genome"""
    if genome is None:
        return None, None
    try:
        start = max(0, pos - window - 1)
        end   = pos + window
        seq   = str(genome[chrom][start:end]).upper()
        mid   = window
        left  = seq[:mid]
        ref_base = seq[mid]
        right = seq[mid+1:]
        return left, ref_base, right
    except Exception as e:
        return None, None, None

# ==========================
# TABS
# ==========================
tab1, tab2, tab3, tab4, tab5 ,tab6 = st.tabs([
    "🌱 LightGBM (Top1000)",
    "📊 ElasticNet (Top100)",
    "🔬 GWAS Results",
    "🧬 Sequence Explorer",
    "🧫 Candidate Genes",
    "📈 SHAP Analysis"
])

# ==========================
# TAB 1 - LIGHTGBM
# ==========================
with tab1:
    st.header("🌱 LightGBM Model (Top1000 SNPs)")

    variety = st.selectbox("Select Variety", varieties, key="lgbm")

    if st.button("Predict with LightGBM"):
        sample = get_sample(variety, lgbm_features)
        X      = preprocess(sample, lgbm_features, lgbm_imputer)
        pred   = lgbm_model.predict(X)[0]
        actual = float(pheno_df.loc[variety, SEED_COL])

        col1, col2, col3 = st.columns(3)
        col1.metric("Predicted (mm)", f"{pred:.4f}")
        col2.metric("Actual (mm)",    f"{actual:.4f}")
        col3.metric("Error (mm)",     f"{abs(actual - pred):.4f}")

# ==========================
# TAB 2 - ELASTICNET
# ==========================
with tab2:
    st.header("📊 ElasticNet Model (Top100 SNPs)")

    variety = st.selectbox("Select Variety", varieties, key="elastic")

    if st.button("Predict with ElasticNet"):
        sample      = get_sample(variety, elastic_features)
        X           = preprocess(sample, elastic_features, elastic_imputer)
        pred        = elastic_model.predict(X)[0]
        actual      = float(pheno_df.loc[variety, SEED_COL])
        actual_value= actual if pd.notna(actual) else 0

        col1, col2, col3 = st.columns(3)
        col1.metric("Predicted (mm)", f"{pred:.4f}")
        col2.metric("Actual (mm)",    f"{actual_value:.4f}")
        col3.metric("Error (mm)",     f"{abs(actual_value - pred):.4f}")

# ==========================
# TAB 3 - GWAS RESULTS
# ==========================
with tab3:
    st.header("🔬 GWAS Results — Seed Thickness (FarmCPU)")

    st.subheader("Genome-wide Significant SNPs")
    st.dataframe(
        signals_df[["SNP","Chr","Pos","MAF","Effect","Pvalue"]]
        .style.format({"Pvalue": "{:.2e}", "Effect": "{:.4f}", "MAF": "{:.3f}"}),
        use_container_width=True
    )

    st.subheader("Manhattan Plot")
    manhattan_candidates = [
        "stoneThickness.FarmCPU.Rectangular-Manhattan.StoneThickness.jpg",
        "SeedThickness.FarmCPU.Rectangular-Manhattan.jpg",
        "SeedThickness_Manhattan.jpg",
        "Manhattan.jpg"
    ]
    manhattan_found = False
    for img in manhattan_candidates:
        if os.path.exists(img):
            st.image(img, use_container_width=True)
            manhattan_found = True
            break
    if not manhattan_found:
        st.info("📁 Place your Manhattan plot image in the project folder to display it here.")

    st.subheader("QQ Plot")
    qq_candidates = [
        "SeedThickness.FarmCPU.QQ.jpg",
        "QQ.jpg",
        "SeedThickness.FarmCPU.QQplot.jpg"
    ]
    qq_found = False
    for img in qq_candidates:
        if os.path.exists(img):
            st.image(img, width=900)
            qq_found = True
            break
    if not qq_found:
        st.info("📁 Place your QQ plot image in the project folder to display it here.")


    st.subheader("PCA Plot")
    qq_candidates = [
        "SeedThickness.PCA_2D.jpg",
    ]
    qq_found = False
    for img in qq_candidates:
        if os.path.exists(img):
            st.image(img, width=900)
            qq_found = True
            break
    if not qq_found:
        st.info("📁 Place your PCA plot image in the project folder to display it here.")

# ==========================
# TAB 4 - SEQUENCE EXPLORER
# ==========================
with tab4:
    st.header("🧬 Sequence Explorer")
    st.write("View DNA sequence context around each significant SNP for any mango variety.")

    variety_seq = st.selectbox("Select Variety", varieties, key="seq")

    if st.button("Explore Sequences"):
        for _, snp_row in signals_df.iterrows():
            snp_id = snp_row["SNP"]
            chrom  = snp_row["Chrom"]
            pos    = int(snp_row["Pos"])
            parts  = snp_id.split("_")
            ref    = parts[3] if len(parts) >= 5 else "?"
            alt    = parts[4] if len(parts) >= 5 else "?"

            st.markdown(f"---")
            st.markdown(f"### SNP: `{snp_id}`")

            col1, col2, col3 = st.columns(3)
            col1.metric("Chromosome", chrom)
            col2.metric("Position",   f"{pos:,}")
            col3.metric("P-value",    f"{snp_row['Pvalue']:.2e}")

            # Raw genotype call for this variety
            raw_call = get_raw_genotype(variety_seq, snp_id)

            # Encoded genotype
            encoded_series = encode_snp_column(
                pd.Series([raw_call], name=snp_id), snp_id
            )
            encoded_val = encoded_series.iloc[0]
            decoded     = decode_genotype(encoded_val, snp_id)

            gc1, gc2, gc3 = st.columns(3)
            gc1.metric("Raw Call",        raw_call)
            gc2.metric("Encoded (0/1/2)", str(int(encoded_val)) if not np.isnan(encoded_val) else "N/A")
            gc3.metric("Genotype",        decoded)

            # Effect direction
            effect = snp_row["Effect"]
            if not np.isnan(encoded_val):
                direction = ""
                if encoded_val == 2:
                    direction = f"ALT homozygous — effect: {'↑ thicker' if effect > 0 else '↓ thinner'} ({effect:+.3f} mm)"
                elif encoded_val == 1:
                    direction = f"Heterozygous — partial effect ({effect/2:+.3f} mm estimated)"
                elif encoded_val == 0:
                    direction = "REF homozygous — baseline"
                if direction:
                    st.info(f"🧬 **Phenotypic effect:** {direction}")

            # DNA Sequence context from reference genome
            st.markdown("**Reference Genome Sequence Context (±10 bp):**")
            if genome is not None:
                left, ref_base, right = get_sequence_context(chrom, pos, window=10)
                if left is not None:
                    st.code(f"Reference:  {left}[{ref}]{right}")
                    st.code(f"Alternate:  {left}[{alt}]{right}")
                    if decoded != "N/A":
                        alleles = decoded.split("/")
                        st.code(f"This variety ({variety_seq}): {left}[{'/'.join(alleles)}]{right}")
                else:
                    st.warning("Could not fetch sequence for this position.")
            else:
                st.info("📁 Reference genome not loaded — sequence context unavailable.")

# ==========================
# TAB 5 - CANDIDATE GENES
# ==========================
with tab5:
    st.header("🧫 Candidate Genes")
    st.write("Genes within ±50kb of each significant SNP, with functional annotations.")

    if annotations_df is None:
        st.error("❌ Gene_Annotations.xlsx not found. Run setup_gene_annotations.py first.")
    else:
        for _, snp_row in signals_df.iterrows():
            snp_id = snp_row["SNP"]
            st.markdown(f"---")
            st.markdown(f"### SNP: `{snp_id}`")

            col1, col2 = st.columns(2)
            col1.metric("P-value", f"{snp_row['Pvalue']:.2e}")
            col2.metric("Effect",  f"{snp_row['Effect']:+.4f} mm")

            snp_genes = annotations_df[annotations_df["SNP"] == snp_id].copy()
            snp_genes = snp_genes.sort_values("Distance_bp")

            if snp_genes.empty:
                st.warning("No candidate genes found near this SNP.")
                continue

            st.markdown(f"**{len(snp_genes)} candidate genes found within ±50kb:**")

            for _, gene in snp_genes.iterrows():
                region_icon = {
                    "Within gene": "🎯",
                    "Upstream"   : "⬆️",
                    "Downstream" : "⬇️"
                }.get(gene["Region"], "📍")

                with st.expander(
                    f"{region_icon} {gene['Gene_ID']}  —  {gene['Region']}  —  {int(gene['Distance_bp']):,} bp away"
                ):
                    info1, info2 = st.columns(2)
                    info1.markdown(f"**Gene ID:** `{gene['Gene_ID']}`")
                    info1.markdown(f"**Region:** {gene['Region']}")
                    info1.markdown(f"**Distance:** {int(gene['Distance_bp']):,} bp")
                    info2.markdown(f"**Gene Start:** {int(gene['Gene_Start']):,}")
                    info2.markdown(f"**Gene End:** {int(gene['Gene_End']):,}")
                    info2.markdown(f"**Chromosome:** {gene['Chrom']}")

                    st.markdown("**Functional Description:**")
                    desc = gene.get("Description", "Not available")
                    if desc and desc != "Not found":
                        st.success(f"📋 {desc}")
                    else:
                        st.warning("No description available in local genome files.")

                    # Biological hypothesis based on description
                    desc_lower = str(desc).lower()
                    hypothesis = None
                    if "disease resistance" in desc_lower:
                        hypothesis = "Disease resistance proteins may influence seed coat development and thickness through stress-response pathways."
                    elif "f-box" in desc_lower:
                        hypothesis = "F-box proteins regulate protein degradation via ubiquitin pathway — may control cell growth timing during seed development."
                    elif "tryptophan" in desc_lower or "ligase" in desc_lower:
                        hypothesis = "Aminoacyl-tRNA ligases are essential for protein synthesis — variants may affect translational efficiency in developing seed tissue."
                    elif "cytochrome" in desc_lower:
                        hypothesis = "Cytochrome P450 enzymes are involved in biosynthesis of cell wall components and lignin — directly relevant to endocarp/seed thickness."
                    elif "uncharacterized" in desc_lower:
                        hypothesis = "Novel uncharacterized gene — proximity to a significant SNP makes it a candidate for functional characterization in mango seed development."

                    if hypothesis:
                        st.markdown("**Biological Hypothesis:**")
                        st.info(f"💡 {hypothesis}")

        st.markdown("---")
        st.markdown("**Full Candidate Gene Table:**")
        display_cols = ["SNP", "Gene_ID", "Distance_bp", "Region", "Description"]
        st.dataframe(
            annotations_df[display_cols].sort_values(["SNP","Distance_bp"]),
            use_container_width=True
        )

        if st.download_button(
            "⬇️ Download Gene Annotations",
            data=open("Gene_Annotations.xlsx", "rb").read(),
            file_name="Gene_Annotations.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        ):
            pass

with tab6:
    st.header("📈 SHAP Analysis")

    st.subheader("SHAP Summary Plot")
    if os.path.exists("shap_summary_plot.png"):
        st.image("shap_summary_plot.png", width=900)
    else:
        st.warning("shap_summary_plot.png not found.")

    st.subheader("SHAP Dependence Plot")
    if os.path.exists("shap_dependence_plot.png"):
        st.image("shap_dependence_plot.png", width=900)
    else:
        st.warning("shap_dependence_plot.png not found.")

    st.subheader("SHAP Waterfall Plot")
    if os.path.exists("shap_waterfall_plot.png"):
        st.image("shap_waterfall_plot.png", width=900)
    else:
        st.warning("shap_waterfall_plot.png not found.")

    if os.path.exists("SHAP_Values.xlsx"):
        xls = pd.ExcelFile("SHAP_Values.xlsx")

        if "Feature_Importance" in xls.sheet_names:
            st.subheader("Top 20 Important SNPs")
            fi = pd.read_excel(xls, sheet_name="Feature_Importance")
            st.dataframe(fi.head(20), use_container_width=True)

        if "All_SHAP_Values" in xls.sheet_names:
            with st.expander("View SHAP Values"):
                shap_df = pd.read_excel(xls, sheet_name="All_SHAP_Values")
                st.dataframe(shap_df, use_container_width=True)

        with open("SHAP_Values.xlsx", "rb") as f:
            st.download_button(
                "⬇️ Download SHAP Values",
                f.read(),
                file_name="SHAP_Values.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
    else:
        st.warning("SHAP_Values.xlsx not found.")
