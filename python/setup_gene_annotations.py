# =============================================================================
# STEP 2: Gene Annotation Fetcher (100% LOCAL - no internet needed)
# Run this ONCE before launching the Streamlit app.
#
# HOW TO RUN:
#   python setup_gene_annotations.py
#
# WHAT IT DOES:
#   1. Reads your 2 significant SNPs from SeedThickness.FarmCPU_signals.csv
#   2. Parses genomic.gtf to find genes within +-50kb of each SNP
#   3. Extracts gene/protein descriptions from genomic.gbff (local file)
#   4. Saves everything to Gene_Annotations.xlsx (used by the app)
#
# NO INTERNET REQUIRED - everything uses your local genome files
# =============================================================================

import pandas as pd
import os
import sys
import re

SIGNALS_FILE = "SeedThickness.FarmCPU_signals.csv"
GTF_FILE     = "genomic.gtf"
GBFF_FILE    = "genomic.gbff"
OUTPUT_FILE  = "Gene_Annotations.xlsx"
WINDOW_BP    = 50000

# -------------------------------------------------------
# 1. Load significant SNPs
# -------------------------------------------------------
print("\n====== Loading Significant SNPs ======")
signals = pd.read_csv(SIGNALS_FILE)
signals.columns = ["SNP", "Chr", "Pos", "MAF", "Effect", "SE", "Pvalue"]
signals["Pos"]   = signals["Pos"].astype(int)
signals["Chrom"] = signals["SNP"].apply(lambda x: "_".join(x.split("_")[:2]))

print(f"  Significant SNPs: {len(signals)}")
for _, row in signals.iterrows():
    print(f"  {row['SNP']}  Chrom={row['Chrom']}  Pos={row['Pos']:,}  P={row['Pvalue']:.2e}")

# -------------------------------------------------------
# 2. Parse GTF - gene features only
# -------------------------------------------------------
print(f"\n====== Parsing GTF ======")
print("  Reading gene features... (may take 1-2 min)")

genes = []
with open(GTF_FILE, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue

        chrom      = parts[0]
        gene_start = int(parts[3])
        gene_end   = int(parts[4])
        strand     = parts[6]
        attributes = parts[8]

        gene_id   = ""
        m = re.search(r'gene_id "([^"]+)"', attributes)
        if m:
            gene_id = m.group(1)

        gene_name = ""
        m2 = re.search(r'gene_name "([^"]+)"', attributes)
        if m2:
            gene_name = m2.group(1)

        genes.append({
            "Chrom"    : chrom,
            "GeneStart": gene_start,
            "GeneEnd"  : gene_end,
            "Strand"   : strand,
            "Gene_ID"  : gene_id,
            "Gene_Name": gene_name
        })

genes_df = pd.DataFrame(genes)
print(f"  Total genes parsed: {len(genes_df):,}")
print(f"  Sample chromosomes in GTF: {list(genes_df['Chrom'].unique()[:5])}")
print(f"  SNP chromosomes needed:    {list(signals['Chrom'].unique())}")

# -------------------------------------------------------
# 3. Find candidate genes within +-50kb of each SNP
# -------------------------------------------------------
print(f"\n====== Searching +-{WINDOW_BP//1000}kb Window ======")

candidate_rows = []

for _, snp in signals.iterrows():
    snp_id    = snp["SNP"]
    chrom     = snp["Chrom"]
    pos       = snp["Pos"]
    win_start = pos - WINDOW_BP
    win_end   = pos + WINDOW_BP

    nearby = genes_df[
        (genes_df["Chrom"]     == chrom) &
        (genes_df["GeneEnd"]   >= win_start) &
        (genes_df["GeneStart"] <= win_end)
    ].copy()

    print(f"\n  SNP: {snp_id}")
    print(f"  Window: {win_start:,} - {win_end:,} on {chrom}")
    print(f"  Nearby genes: {len(nearby)}")

    for _, gene in nearby.iterrows():
        gene_mid = (gene["GeneStart"] + gene["GeneEnd"]) // 2
        distance = pos - gene_mid

        if gene["GeneStart"] <= pos <= gene["GeneEnd"]:
            region = "Within gene"
        elif distance > 0:
            region = "Downstream"
        else:
            region = "Upstream"

        candidate_rows.append({
            "SNP"        : snp_id,
            "Chrom"      : chrom,
            "SNP_Pos"    : pos,
            "Pvalue"     : snp["Pvalue"],
            "Effect"     : snp["Effect"],
            "Gene_ID"    : gene["Gene_ID"],
            "Gene_Name"  : gene["Gene_Name"],
            "Gene_Start" : gene["GeneStart"],
            "Gene_End"   : gene["GeneEnd"],
            "Distance_bp": abs(distance),
            "Region"     : region
        })
        print(f"    {gene['Gene_ID']}  {region}  {abs(distance):,}bp")

if len(candidate_rows) == 0:
    print("\n❌ No candidate genes found!")
    print("   GTF chromosomes vs SNP chromosomes mismatch - check above.")
    sys.exit(1)

candidates_df = pd.DataFrame(candidate_rows)
print(f"\n  Total candidate gene-SNP pairs: {len(candidates_df)}")

# -------------------------------------------------------
# 4. Extract descriptions from GBFF (fully local)
# -------------------------------------------------------
print(f"\n====== Extracting Descriptions from GBFF ======")
print("  Parsing genomic.gbff... (may take 3-5 min)")

unique_genes = set(candidates_df["Gene_ID"].unique())
print(f"  Looking for {len(unique_genes)} genes: {unique_genes}")

annotations = {g: {"Description": "Not found", "Product": "Not found"}
               for g in unique_genes}

found      = set()
curr_gene  = None

with open(GBFF_FILE, "r") as f:
    for line in f:

        # Detect locus_tag or gene matching our targets
        if '/locus_tag="' in line or '/gene="' in line:
            m = re.search(r'/(?:locus_tag|gene)="([^"]+)"', line)
            if m and m.group(1) in unique_genes:
                curr_gene = m.group(1)

        # Grab product description
        if curr_gene and '/product="' in line:
            m = re.search(r'/product="([^"]+)"', line)
            if m:
                product = m.group(1).strip().rstrip('"')
                annotations[curr_gene]["Product"]     = product
                annotations[curr_gene]["Description"] = product
                found.add(curr_gene)
                print(f"  ✅ {curr_gene}: {product[:80]}")
                curr_gene = None

        # Fallback: grab note field
        if curr_gene and '/note="' in line:
            m = re.search(r'/note="([^"]+)"', line)
            if m and annotations[curr_gene]["Description"] == "Not found":
                annotations[curr_gene]["Description"] = m.group(1)[:200]

not_found = unique_genes - found
if not_found:
    print(f"\n  ⚠️  Not found in GBFF: {not_found}")
    for gene in not_found:
        annotations[gene]["Description"] = f"Predicted gene {gene} — see NCBI for details"

print(f"\n  Annotations found: {len(found)}/{len(unique_genes)}")

# -------------------------------------------------------
# 5. Merge and save
# -------------------------------------------------------
print(f"\n====== Saving to {OUTPUT_FILE} ======")

annot_df = pd.DataFrame.from_dict(annotations, orient="index").reset_index()
annot_df.columns = ["Gene_ID", "Description", "Product"]

final_df = candidates_df.merge(annot_df, on="Gene_ID", how="left")
final_df = final_df.sort_values(["SNP", "Distance_bp"])

with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
    final_df.to_excel(writer, sheet_name="Candidate_Genes", index=False)
    signals.to_excel(writer, sheet_name="Significant_SNPs", index=False)

print(f"  ✅ Saved: {OUTPUT_FILE}  ({len(final_df)} rows)")
print("\n  Preview:")
print(final_df[["SNP","Gene_ID","Distance_bp","Region","Description"]].head(10).to_string(index=False))

print("\n========================================")
print(" ANNOTATION SETUP COMPLETE!")
print(f" Output: {OUTPUT_FILE}")
print(" Next step: run the updated app.py")
print("========================================\n")