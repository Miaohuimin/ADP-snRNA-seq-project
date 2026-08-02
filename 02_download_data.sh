#!/bin/bash
# =============================================
# Download ADP scRNA-seq FASTQ files
# Usage: bash scripts/bash/01_download_data.sh
# Note: Replace the download URLs with your own valid links
# =============================================

set -e
set -u
set -o pipefail

# =============================================
# User-configurable parameters
# =============================================

# Project root directory (auto-detected)
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

# Raw data directory
RAW_DATA_DIR="${PROJECT_DIR}/data/raw"

# =============================================
# Sample configuration
# =============================================

# Define samples and their download URLs
# Note: These URLs are placeholders. Replace with actual valid links.
# IMPORTANT: Do not upload URLs containing access keys or temporary signatures!
declare -A SAMPLES=(
    ["ADP_C"]="
        https://example.com/path/to/ADP_C_S1_L007_R1_001.fastq.gz
        https://example.com/path/to/ADP_C_S1_L007_R2_001.fastq.gz
    "
    ["ADP_N"]="
        https://example.com/path/to/ADP_N_S1_L007_R1_001.fastq.gz
        https://example.com/path/to/ADP_N_S1_L007_R2_001.fastq.gz
    "
    ["ADP_H"]="
        https://example.com/path/to/ADP_H_S1_L005_R1_001.fastq.gz
        https://example.com/path/to/ADP_H_S1_L005_R2_001.fastq.gz
        https://example.com/path/to/ADP_H_S1_L006_R1_001.fastq.gz
        https://example.com/path/to/ADP_H_S1_L006_R2_001.fastq.gz
        https://example.com/path/to/ADP_H_S1_L001_R1_001.fastq.gz
        https://example.com/path/to/ADP_H_S1_L001_R2_001.fastq.gz
    "
)

# =============================================
# Main script
# =============================================

echo "[$(date)] Working directory: $PROJECT_DIR"

# Create directories for each sample
for SAMPLE in "${!SAMPLES[@]}"; do
    mkdir -p "${RAW_DATA_DIR}/${SAMPLE}"
    echo "[$(date)] Created directory: ${RAW_DATA_DIR}/${SAMPLE}"
done

# Download files for each sample
for SAMPLE in "${!SAMPLES[@]}"; do
    echo "[$(date)] Processing sample: $SAMPLE"
    
    # Read URLs for this sample (split by newline)
    while IFS= read -r URL; do
        # Skip empty lines
        [ -z "$URL" ] && continue
        
        # Extract filename from URL
        FILENAME=$(basename "$URL" | sed 's/?.*//')
        OUTPUT_PATH="${RAW_DATA_DIR}/${SAMPLE}/${FILENAME}"
        
        # Check if file already exists
        if [ -f "$OUTPUT_PATH" ]; then
            echo "[$(date)] $FILENAME already exists. Skipping download."
        else
            echo "[$(date)] Downloading: $FILENAME"
            wget -c "$URL" -O "$OUTPUT_PATH"
        fi
    done <<< "${SAMPLES[$SAMPLE]}"
done

# Generate MD5 checksums for verification
echo "[$(date)] Generating MD5 checksums ..."
cd "$RAW_DATA_DIR"
find . -name "*.fastq.gz" -type f -exec md5sum {} \; > "${PROJECT_DIR}/results/logs/md5_checksums.txt"

echo "[$(date)] MD5 checksums saved to: results/logs/md5_checksums.txt"
echo "[$(date)] All downloads complete!"

# Note: Verify MD5 checksums against the values provided by the sequencing provider