#!/bin/bash
# =============================================
# Build a custom mouse reference genome (mm10) with EGFP transgene
# Usage: bash scripts/bash/00_build_reference.sh
# Note: Ensure sufficient disk space (> 100 GB) and memory (> 120 GB)
# =============================================

set -e
set -u
set -o pipefail

# =============================================
# User-configurable parameters
# =============================================

# Path to Cell Ranger executable (modify to match your server environment)
CELLRANGER_CMD="/home/mzk/miuhuimin/software/cellranger/cellranger-10.0.0/bin/cellranger"

# Reference genome download URL (10x mouse mm10-2020-A)
REFERENCE_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"
REFERENCE_TAR="refdata-gex-mm10-2020-A.tar.gz"
REFERENCE_DIR="refdata-gex-mm10-2020-A"

# Custom reference output directory (original name + _with_EGFP)
CUSTOM_DIR="${REFERENCE_DIR}_with_EGFP"

# Resource settings (adjust based on your server configuration)
THREADS=20
MEMORY_GB=120

# Final reference index name
GENOME_NAME="myrefdata_base_on_mm10-2020-A"

# =============================================
# Main script (no modifications needed below this line)
# =============================================

# Get the project root directory (two levels up from this script)
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

echo "[$(date)] Working directory: $PROJECT_DIR"

# ---------- Step 1: Download reference genome (skip if already exists) ----------
if [ -f "$REFERENCE_TAR" ]; then
    echo "[$(date)] Reference genome archive already exists. Skipping download."
else
    echo "[$(date)] Downloading reference genome from: $REFERENCE_URL ..."
    wget "$REFERENCE_URL"
fi

# Verify file integrity
echo "[$(date)] Verifying MD5 checksum ..."
md5sum "$REFERENCE_TAR"

# ---------- Step 2: Extract reference genome (skip if already extracted) ----------
if [ -d "$REFERENCE_DIR" ]; then
    echo "[$(date)] Reference genome already extracted. Skipping extraction."
else
    echo "[$(date)] Extracting reference genome ..."
    tar -zxvf "$REFERENCE_TAR"
fi

# ---------- Step 3: Create custom reference directory ----------
mkdir -p "$CUSTOM_DIR"
echo "[$(date)] Copying reference genome to custom directory: $CUSTOM_DIR"
cp -r "$REFERENCE_DIR"/* "$CUSTOM_DIR"/

# ---------- Step 4: Unzip GTF annotation file ----------
# Note: 10x reference genes.gtf is usually compressed as .gz
if [ -f "$CUSTOM_DIR/genes/genes.gtf" ]; then
    echo "[$(date)] genes.gtf already unzipped. Skipping."
else
    echo "[$(date)] Unzipping genes.gtf.gz ..."
    gunzip "$CUSTOM_DIR/genes/genes.gtf.gz"
fi

# ---------- Step 5: Create EGFP sequence file ----------
echo "[$(date)] Creating EGFP.fa ..."
cat > EGFP.fa <<'EOF'
>EGFP
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAAAGCGGCCGCGACTTTAGAATTCAATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCTCCTTTTACGCTATGTGGATACGCTGCTTTAATGCCTTTGTATCATGCTATTGCTTCCCGTATGGCTTTCATTTTCTCCTCCTTGTATAAATCCTGGTTGCTGTCTCTTTATGAGGAGTTGTGGCCCGTTGTCAGGCAACGTGGCGTGGTGTGCACTGTGTTTGCTGACGCAACCCCCACTGGTTGGGGCATTGCCACCACCTGTCAGCTCCTTTCCGGGACTTTCGCTTTCCCCCTCCCTATTGCCACGGCGGAACTCATCGCCGCCTGCCTTGCCCGCTGCTGGACAGGGGCTCGGCTGTTGGGCACTGACAATTCCGTGGTGTTGTCGGGGAAATCATCGTCCTTTCCTTGGCTGCTCGCCTGTGTTGCCACCTGGATTCTGCGCGGGACGTCCTTCTGCTACGTCCCTTCGGCCCTCAATCCAGCGGACCTTCCTTCCCGCGGCCTGCTGCCGGCTCTGCGGCCTCTTCCGCGTCTTCGCCTTCGCCCTCAGACGAGTCGGATCTCCCTTTGGGCCGCCTCCCCGCGGGTGGCATCCCTGTGACCCCTCCCCAGTGCCTCTCCTGGCCCTGGAAGTTGCCACTCCAGTGCCCACCAGCCTTGTCCTAATAAAATTAAGTTGCATCATTTTGTCTGACTAGGTGTCCTTCTATAATAT
EOF

# ---------- Step 6: Create EGFP GTF annotation ----------
echo "[$(date)] Creating EGFP.gtf ..."
printf "EGFP\tunknown\texon\t1\t720\t.\t+\t.\tgene_id \"EGFP\"; transcript_id \"EGFP\"; gene_name \"EGFP\";\n" > EGFP.gtf

# ---------- Step 7: Append EGFP to reference genome ----------
echo "[$(date)] Appending EGFP sequence to genome.fa ..."
cat EGFP.fa >> "$CUSTOM_DIR/fasta/genome.fa"

echo "[$(date)] Appending EGFP annotation to genes.gtf ..."
cat EGFP.gtf >> "$CUSTOM_DIR/genes/genes.gtf"

# ---------- Step 8: Verify the append operations ----------
echo "[$(date)] Verifying genome.fa (last 5 lines):"
tail -n 5 "$CUSTOM_DIR/fasta/genome.fa"

echo "[$(date)] Verifying genes.gtf (last 5 lines):"
tail -n 5 "$CUSTOM_DIR/genes/genes.gtf"

# ---------- Step 9: Build reference index with Cell Ranger ----------
echo "[$(date)] Building reference index with Cell Ranger ..."
echo "[$(date)] This may take 1-2 hours. Please be patient ..."

$CELLRANGER_CMD mkref \
    --genome="$GENOME_NAME" \
    --fasta="$CUSTOM_DIR/fasta/genome.fa" \
    --genes="$CUSTOM_DIR/genes/genes.gtf" \
    --nthreads=$THREADS \
    --memgb=$MEMORY_GB

echo "[$(date)] Reference index build complete! Output directory: ./$GENOME_NAME"

# ---------- Step 10: Clean up temporary files (optional) ----------
# Uncomment the following lines to remove intermediate files
# rm -f EGFP.fa EGFP.gtf
# echo "[$(date)] Temporary files cleaned up."

echo "[$(date)] All done!"