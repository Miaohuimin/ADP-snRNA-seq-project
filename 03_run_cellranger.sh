#!/bin/bash
# =============================================
# Run Cell Ranger count for ADP samples
# Usage: bash scripts/bash/02_run_cellranger.sh
# Note: Run in tmux session to prevent interruption
# =============================================

set -e
set -u
set -o pipefail

# =============================================
# User-configurable parameters
# =============================================

# Cell Ranger executable path (same as used in build_reference script)
CELLRANGER_CMD="/home/mzk/miuhuimin/software/cellranger/cellranger-10.0.0/bin/cellranger"

# Path to custom reference genome (built in step 1)
REFERENCE_DIR="/data/miuhuimin/myrefdata_base_on_mm10-2020-A"

# Project root directory
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

# Raw data directory (where FASTQ files are stored)
RAW_DATA_DIR="${PROJECT_DIR}/data/raw"

# Output directory for Cell Ranger results
OUTPUT_DIR="${PROJECT_DIR}/data/processed"

# Resource settings
CORES=8
MEMORY_GB=120

# Chemistry version
CHEMISTRY="SC3Pv4"

# =============================================
# Sample configuration
# =============================================

# Define samples
SAMPLES=("ADP_C" "ADP_N" "ADP_H")

# =============================================
# Main script
# =============================================

echo "[$(date)] Working directory: $PROJECT_DIR"
echo "[$(date)] Using Cell Ranger: $CELLRANGER_CMD"
echo "[$(date)] Reference genome: $REFERENCE_DIR"

# Check if reference genome exists
if [ ! -d "$REFERENCE_DIR" ]; then
    echo "[$(date)] ERROR: Reference genome not found at: $REFERENCE_DIR"
    echo "[$(date)] Please run 00_build_reference.sh first or check the path."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Loop through each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo "[$(date)] Processing sample: $SAMPLE"
    
    # Define output ID
    OUTPUT_ID="${SAMPLE}_EGFP_output"
    OUTPUT_PATH="${OUTPUT_DIR}/${OUTPUT_ID}"
    
    # Check if sample has already been processed (resume capability)
    if [ -d "$OUTPUT_PATH/outs/filtered_feature_bc_matrix" ]; then
        echo "[$(date)] $SAMPLE already processed. Skipping."
        continue
    fi
    
    # Check if FASTQ files exist
    FASTQ_PATH="${RAW_DATA_DIR}/${SAMPLE}"
    if [ ! -d "$FASTQ_PATH" ] || [ -z "$(ls -A "$FASTQ_PATH" 2>/dev/null)" ]; then
        echo "[$(date)] WARNING: No FASTQ files found for $SAMPLE at $FASTQ_PATH"
        echo "[$(date)] Please run 01_download_data.sh first."
        continue
    fi
    
    # Run Cell Ranger count
    echo "[$(date)] Running Cell Ranger count for $SAMPLE ..."
    echo "[$(date)] This may take several hours. Please be patient."
    
    $CELLRANGER_CMD count \
        --id="$OUTPUT_ID" \
        --create-bam=false \
        --chemistry="$CHEMISTRY" \
        --transcriptome="$REFERENCE_DIR" \
        --fastqs="$FASTQ_PATH" \
        --sample="$SAMPLE" \
        --localcores="$CORES" \
        --localmem="$MEMORY_GB" \
        --disable-ui
    
    echo "[$(date)] $SAMPLE processing complete!"
    echo "[$(date)] Output: $OUTPUT_PATH"
done

# Summary of output sizes
echo "[$(date)] =========================================="
echo "[$(date)] Summary of output sizes:"
for SAMPLE in "${SAMPLES[@]}"; do
    OUTPUT_PATH="${OUTPUT_DIR}/${SAMPLE}_EGFP_output"
    if [ -d "$OUTPUT_PATH" ]; then
        SIZE=$(du -sh "$OUTPUT_PATH" 2>/dev/null | cut -f1)
        echo "[$(date)] $SAMPLE: $SIZE"
    else
        echo "[$(date)] $SAMPLE: Not found"
    fi
done

echo "[$(date)] All samples processed!"
echo "[$(date)] To monitor progress, use: watch -n 10 du -sh ${OUTPUT_DIR}/*"