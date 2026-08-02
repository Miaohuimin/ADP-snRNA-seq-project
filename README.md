# snRNA-seq Analysis Pipeline

This repository contains the analysis scripts for the **Circadian thermal control via a segregated SCN ADP RPa circuit** project.

## Data Availability

The raw sequence data have been deposited in the Genome Sequence Archive (Genomics, Proteomics & Bioinformatics 2025) in National Genomics Data Center (Nucleic Acids Res 2026), China National Center for Bioinformation / Beijing Institute of Genomics, Chinese Academy of Sciences (GSA: CRA047474) that are publicly accessible at https://ngdc.cncb.ac.cn/gsa.

## Requirements

- R (>= 4.0)
- Cell Ranger (>= 7.0)
- Bash shell environment with `wget`, `md5sum`, `tar`

## Repository Structure

```text
.
├── scripts/
│   ├── bash/                              # Upstream processing scripts
│   │   ├── 01_build_reference.sh          # Build custom reference genome (mm10 + EGFP)
│   │   ├── 02_download_data.sh            # Download FASTQ files
│   │   └── 03_run_cellranger.sh           # Run Cell Ranger quantification
│   └── R/                                 # Downstream analysis scripts (Seurat)
│       ├── 01_quality_control.R           # QC, filtering, SCTransform
│       ├── 02_remove_doublets.R           # Doublet removal with DoubletFinder
│       ├── 03_integrate_clustering_annotate.R      # Harmony integration & annotation
│       ├── 04_unstimulated_EGFP_neurons_analysis.R # Unstimulated EGFP+ neuron analysis
│       ├── 05_combined_EGFP_neurons_subclustering.R # EGFP+ neuron sub-clustering
│       ├── 06_reference_mapping.R         # Reference mapping
│       ├── 07_differential_expression.R   # DEG: ADP_H/ADP_C vs ADP_N (clusters 0+19)
│       ├── 08_Neuron_DEG_GO.R             # Neuronal DEG + GO enrichment
│       └── 09_AllCellTypes_DEGs.R # Volcano plots for all cell types
├── config/                                # Configuration files
├── data/                                  # Data directory (raw/processed, gitignored)
└── results/                               # Output figures and tables
```

## Usage

### 1. Upstream analysis (Bash scripts)

All Bash scripts should be executed in the project root directory.

#### 1.1 Build a custom reference genome with EGFP transgene

Before running Cell Ranger quantification, you need to build a custom reference genome (mm10-2020-A) containing the EGFP sequence.

```ba
# Make the script executable (once)
chmod +x scripts/bash/01_build_reference.sh

# Run the script
bash scripts/bash/01_build_reference.sh
```

**Notes**:

- Update the `CELLRANGER_CMD` variable in the script to match your server's Cell Ranger installation path.
- Ensure sufficient disk space (> 100 GB) and memory (> 120 GB).
- The output index directory is `./myrefdata_base_on_mm10-2020-A`.
- If you already have `refdata-gex-mm10-2020-A` downloaded, the script will automatically skip download and extraction.

#### 1.2 Download FASTQ files

Download raw sequencing data from the sequencing provider.

```bas
bash scripts/bash/02_download_data.sh
```

**Important Notes**:

- **Replace all placeholder URLs** in `scripts/bash/02_download_data.sh` with your own valid download links before running.
- **Do NOT upload** URLs containing temporary access keys (e.g., `OSSAccessKeyId`, `Signature`, `Expires`) to public repositories due to security risks.
- After download, verify MD5 checksums against the values provided by the sequencing provider. The script generates a checksum file at `results/logs/md5_checksums.txt`.

#### 1.3 Run Cell Ranger quantification

After the reference genome is built and FASTQ files are ready, run Cell Ranger count for all samples.

```ba
bash scripts/bash/03_run_cellranger.sh
```

**Notes**:

- It is **strongly recommended** to run this script inside a `tmux` or `screen` session to prevent interruption from network disconnection.
- The script will automatically skip samples that have already been processed (based on output directory existence).
- Output will be saved to `data/processed/` with the suffix `_EGFP_output`.
- To monitor progress during running, use: `watch -n 10 du -sh data/processed/*`
- Quality metrics for all samples are summarized in `results/tables/quality_metrics_summary.csv`.

------

### 2. Downstream analysis (R scripts)

After Cell Ranger quantification, perform downstream analysis with Seurat.
All R scripts should be executed **in the order shown below**.

#### 2.1 Quality control and filtering

Read H5 files, calculate QC metrics, filter cells, and apply SCTransform normalization.

```ba
Rscript scripts/R/01_quality_control.R
```

**Notes**:

- Ensure all required R packages are installed (Seurat, ggplot2, tidyverse, patchwork, future).
- The H5 files are expected to be located in `data/processed/`.
- Outputs: `ADPC_data_filtered.rds`, `ADPN_data_filtered.rds`, `ADPH_data_filtered.rds`, plus their SCTransform versions.
- A session info file is saved as `session_info_01_QC.txt`.

#### 2.2 Remove doublets with DoubletFinder

Remove doublets using the DoubletFinder algorithm.

**Important**: This script must be run in a dedicated conda environment:

```ba
# Create the environment (once)
conda create -n r-doubletfinder r-base r-seurat -c conda-forge -y

# Activate and run
conda activate r-doubletfinder
Rscript scripts/R/02_remove_doublets.R
```

**Notes**:

- The script expects SCTransform-normalized RDS files from step 2.1 in `data/processed/`.
- pK values for each sample were pre-determined via `paramSweep` and are hard-coded in the script.
- Outputs: `ADP_C_clean.rds`, `ADP_H_clean.rds`, `ADP_N_clean.rds` saved to `data/processed/`.

#### 2.3 Integration, clustering and annotation

Integrate ADP samples with Harmony, annotate cell types, and generate publication-ready figures.

```ba
Rscript scripts/R/03_integrate_clustering_annotate.R
```

**Notes**:

- Requires `*_clean.rds` files from step 2.2 in the working directory.
- All figures are automatically saved to `results/figures/`.
- Key output files: `UMAP_by_celltype.pdf`, `DotPlot_marker_genes.pdf`, and `iLISI_FeaturePlot.pdf`.

#### 2.4 Unstimulated EGFP+ neuron sub-clustering and differential expression

Extract EGFP+ neurons from ADP_N sample, perform differential expression, sub-clustering, and visualization.

```ba
Rscript scripts/R/04_unstimulated_EGFP_neurons_analysis.R
```

**Notes**:

- Requires the integrated object `ADP_harmony_combined_with_celltype.rds` from step 2.3.
- Outputs: volcano plot, dot plot, and UMAP of EGFP+ neurons.
- All figures saved to `results/figures/` and tables to `results/tables/`.

#### 2.5 EGFP+ neuron sub-clustering and differential expression

Perform sub-clustering analysis on EGFP+ neurons, identify cluster-specific marker genes, and generate visualization plots.

```ba
Rscript scripts/R/05_combined_EGFP_neurons_subclustering.R
```

**Notes**:

- Requires the integrated object `ADP_harmony_combined_with_celltype.rds` from step 2.3.
- EGFP+ cells are identified based on EGFP expression (threshold ≥ 1 UMI).
- **Cluster 3 was removed** from the final analysis due to enrichment of mitochondrial genes, indicating low-quality or apoptotic cells.
- Sub-clustering was performed at resolution 1.2, resulting in 23 clusters after removal of cluster 3.

#### 2.6 Reference mapping of EGFP+ neuron clusters

Project specific EGFP+ neuron clusters (Cluster 0 and Cluster 19) onto the control reference UMAP space.

```ba
Rscript scripts/R/06_reference_mapping.R
```

**Notes**:

- Requires `ADP_E_pos.rds` (from step 2.5) and `ADP_N_pos.rds` (from step 2.4) in the working directory.
- Outputs: UMAP projection plot, pie charts, and mapping score tables.

#### 2.7 Differential expression: ADP_H/ADP_C vs ADP_N

Perform differential expression analysis comparing ADP_H and ADP_C (clusters 0+19) against ADP_N (cluster 0), followed by GO enrichment analysis.

```ba
Rscript scripts/R/07_differential_expression.R
```

**Notes**:

- Requires `ADP_E_pos.rds` and `ADP_N_pos.rds` from previous steps.
- Compares ADP_H (cluster 0+19) vs ADP_N (cluster 0) and ADP_C (cluster 0+19) vs ADP_N (cluster 0).
- Outputs: DEG tables, volcano plots, and GO enrichment results.

#### 2.8 Neuronal DEG and GO enrichment (all samples combined)

Perform differential expression analysis comparing ADP_H and ADP_C against ADP_N in EGFP+ neurons (all samples aggregated), followed by GO enrichment of down-regulated genes.

```bash
Rscript scripts/R/08_Neuron_DEG_GO.R
```

**Notes**:

- Requires `ADP_harmony_combined_with_celltype.rds` from step 2.3.
- Compares ADP_H vs ADP_N and ADP_C vs ADP_N in EGFP+ neurons.
- GO enrichment uses `clusterProfiler` with simplified GO terms (cutoff = 0.7).

#### 2.9 Volcano plots for all major cell types

Generate volcano plots for all major cell types (Astrocyte, OPCs, Oligodendrocyte, Microglia, Mural cell, and Neuron) comparing transcriptional responses to heat (ADP_H vs ADP_N) and cold (ADP_C vs ADP_N) stimulation.

```ba
Rscript scripts/R/09_AllCellTypes_VolcanoPlots.R
```

**Notes**:

- This script automatically checks for and creates all required cell type subset objects (`ADP_Astrocyte_pos`, `ADP_OPCs_pos`, `ADP_Oligo_pos`, `ADP_Microglia_pos`, `ADP_Mural_pos`, `ADP_Neuron_pos`) from `ADP_combined_pos`.
- Each subset is independently re-normalized using `SCTransform` with mitochondrial percentage regression.
- Differential expression analysis is performed using `FindMarkers` (Wilcoxon rank-sum test; |log₂FC| > 0.25, p < 0.05).
- **Top 5 upregulated and downregulated genes** are automatically labeled on each volcano plot based on:
  - Largest absolute log₂ fold change
  - Smallest p-value
- Volcano plots are displayed interactively in RStudio for **manual Y-axis adjustment** to achieve optimal visualization before saving.
- **All DEG tables are automatically saved as CSV files** in `results/tables/` for each cell type and stimulation condition.

**Output files**:

| Type  | File pattern                               | Description                                       |
| :---- | :----------------------------------------- | :------------------------------------------------ |
| Table | `results/tables/DEG_{CellType}_H_vs_N.csv` | DEG results for heat stimulation (ADP_H vs ADP_N) |
| Table | `results/tables/DEG_{CellType}_C_vs_N.csv` | DEG results for cold stimulation (ADP_C vs ADP_N) |

**Cell types processed**: Astrocyte, OPCs, Oligodendrocyte, Microglia, Mural cell, Neuron

------

## Contact

For questions or issues, please open an issue on GitHub or contact the corresponding author.