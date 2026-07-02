# =============================================================================
# STEP 1: One-time Genome Index Setup
# Run this ONCE before launching the Streamlit app.
# It creates a fast lookup index for the reference genome.
#
# HOW TO RUN:
#   Open terminal/command prompt in your mango_gwas folder and run:
#   python setup_genome_index.py
#
# WHAT IT CREATES:
#   genome.fai  ->  fast index for sequence lookups
#
# REQUIRES:
#   pip install biopython pyfaidx
# =============================================================================

import os
import sys

# -------------------------------------------------------
# 1. Check required files exist
# -------------------------------------------------------
GENOME_FILE = "GCF_011075055.1_CATAS_Mindica_2.1_genomic.fna"
GTF_FILE    = "genomic.gtf"
GBFF_FILE   = "genomic.gbff"
PROTEIN_FILE = "protein.faa"

print("\n====== Checking required files ======")
all_ok = True
for f in [GENOME_FILE, GTF_FILE, GBFF_FILE, PROTEIN_FILE]:
    exists = os.path.exists(f)
    status = "✅ Found" if exists else "❌ MISSING"
    print(f"  {status}: {f}")
    if not exists:
        all_ok = False

if not all_ok:
    print("\n❌ Some files are missing! Make sure all 4 files are in this folder.")
    sys.exit(1)

print("\n✅ All files found!")

# -------------------------------------------------------
# 2. Install dependencies if needed
# -------------------------------------------------------
print("\n====== Checking Python packages ======")
try:
    import pyfaidx
    print("  ✅ pyfaidx already installed")
except ImportError:
    print("  Installing pyfaidx...")
    os.system(f"{sys.executable} -m pip install pyfaidx")
    import pyfaidx

try:
    from Bio import SeqIO
    print("  ✅ biopython already installed")
except ImportError:
    print("  Installing biopython...")
    os.system(f"{sys.executable} -m pip install biopython")
    from Bio import SeqIO

# -------------------------------------------------------
# 3. Build FASTA index (.fai file)
# -------------------------------------------------------
print(f"\n====== Indexing {GENOME_FILE} ======")
print("  This may take 1-3 minutes on first run...")

try:
    fasta = pyfaidx.Fasta(GENOME_FILE, build_index=True)
    chromosomes = list(fasta.keys())
    print(f"  ✅ Index created! Found {len(chromosomes)} sequences")
    print(f"  First 5 chromosome IDs:")
    for c in chromosomes[:5]:
        print(f"    {c}")
    fasta.close()
except Exception as e:
    print(f"  ❌ Error indexing genome: {e}")
    sys.exit(1)

# -------------------------------------------------------
# 4. Quick test — fetch a small sequence
# -------------------------------------------------------
print("\n====== Testing sequence fetch ======")
try:
    fasta = pyfaidx.Fasta(GENOME_FILE)
    test_chrom = chromosomes[0]
    test_seq = str(fasta[test_chrom][100:120])
    print(f"  Test fetch from {test_chrom}[100:120]: {test_seq}")
    print(f"  ✅ Sequence fetch working!")
    fasta.close()
except Exception as e:
    print(f"  ❌ Sequence fetch failed: {e}")
    sys.exit(1)

# -------------------------------------------------------
# 5. Quick check on GTF
# -------------------------------------------------------
print("\n====== Checking GTF file ======")
try:
    with open(GTF_FILE, "r") as f:
        lines = [f.readline() for _ in range(10)]
    gene_lines = [l for l in lines if "\tgene\t" in l]
    print(f"  First 10 lines read successfully")
    print(f"  Sample line: {lines[-1][:100]}")
    print(f"  ✅ GTF file readable!")
except Exception as e:
    print(f"  ❌ GTF read failed: {e}")

# -------------------------------------------------------
# Done
# -------------------------------------------------------
print("\n========================================")
print(" SETUP COMPLETE!")
print(" Index file created: genome.fai")
print(" Next step: run setup_gene_annotations.py")
print("========================================\n")