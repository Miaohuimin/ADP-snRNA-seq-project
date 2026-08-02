#!/usr/bin/env Rscript
# =============================================
# Quality control, filtering, and SCTransform normalization
# Usage: Rscript scripts/R/01_quality_control.R
# =============================================

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  library(patchwork)
  library(future)
})

# =============================================
# User-configurable parameters
# =============================================

# Set working directory to project root
# Note: Change this to your actual project path
# Alternatively, run this script in RStudio with the project opened
# setwd("/path/to/your/project")

# Define input and output paths
H5_DIR <- "./data/processed"      # Directory containing filtered_feature_bc_matrix.h5 files
OUTPUT_DIR <- "./data/processed"  # Where to save filtered Seurat objects

# QC thresholds
MIN_FEATURES <- 400
MAX_FEATURES <- 7000
MAX_MT <- 10
MAX_RIBO <- 20
MAX_APOPTOSIS <- 0.2
MAX_RBC <- 0.1

# SCTransform parameters
REGRESS_VARS <- c("percent.mt")
SCT_VERBOSE <- FALSE

# =============================================
# Create output directory if it doesn't exist
# =============================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================
# 1. Read data (10x H5 format)
# =============================================

cat("[", Sys.time(), "] Reading H5 files ...\n")

ADP_C_data <- Read10X_h5(file.path(H5_DIR, "ADP_C_filtered_feature_bc_matrix.h5"))
ADP_N_data <- Read10X_h5(file.path(H5_DIR, "ADP_N_filtered_feature_bc_matrix.h5"))
ADP_H_data <- Read10X_h5(file.path(H5_DIR, "ADP_H_filtered_feature_bc_matrix.h5"))

cat("ADP_C: ", dim(ADP_C_data), "\n")
cat("ADP_N: ", dim(ADP_N_data), "\n")
cat("ADP_H: ", dim(ADP_H_data), "\n")

# =============================================
# 2. Create Seurat objects
# =============================================

cat("[", Sys.time(), "] Creating Seurat objects ...\n")

ADPC_data <- CreateSeuratObject(
  counts = ADP_C_data,
  project = "ADP_C",
  min.cells = 3,
  min.features = MIN_FEATURES
)

ADPN_data <- CreateSeuratObject(
  counts = ADP_N_data,
  project = "ADP_N",
  min.cells = 3,
  min.features = MIN_FEATURES
)

ADPH_data <- CreateSeuratObject(
  counts = ADP_H_data,
  project = "ADP_H",
  min.cells = 3,
  min.features = MIN_FEATURES
)

# =============================================
# 3. Calculate quality metrics
# =============================================

cat("[", Sys.time(), "] Calculating QC metrics ...\n")

# Function to add QC metrics to a Seurat object
add_qc_metrics <- function(obj) {
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj$percent.ribo <- PercentageFeatureSet(obj, pattern = "^Rps|^Rpl")
  return(obj)
}

ADPC_data <- add_qc_metrics(ADPC_data)
ADPN_data <- add_qc_metrics(ADPN_data)
ADPH_data <- add_qc_metrics(ADPH_data)

# =============================================
# 4. Define gene sets for additional filtering
# =============================================

cat("[", Sys.time(), "] Defining gene sets ...\n")

# Apoptosis-related genes (mouse)
apoptosis_genes <- c("Casp3", "Casp7", "Casp8", "Casp9", "Bax", "Bak1")

# Red blood cell marker genes (mouse)
rbc_genes <- c(
  "Hba-a1", "Hba-a2", "Hbb-b1", "Hbb-b2", "Hbb-bs",
  "Hbb-bt", "Hbe1", "Hbg1", "Hbg2", "Alas2", "Ahsp"
)

# Early immediate genes (for reference, not used in filtering)
ieg_genes <- c(
  "Arc", "Atf3", "Cebpd", "Dusp1", "Dusp5", "Dusp6",
  "Egr1", "Egr3", "Fos", "Fosb", "Fosl1", "Gadd45b",
  "Jun", "Junb", "Jund", "Nr4a1", "Nr4a2", "Nr4a3",
  "Srf", "Tnfaip3", "Zfp36", "Myc"
)

# Save gene lists for future use
saveRDS(apoptosis_genes, file = file.path(OUTPUT_DIR, "apoptosis_genes.rds"))
saveRDS(rbc_genes, file = file.path(OUTPUT_DIR, "rbc_genes.rds"))
saveRDS(ieg_genes, file = file.path(OUTPUT_DIR, "ieg_genes.rds"))

# =============================================
# 5. Calculate additional QC metrics
# =============================================

cat("[", Sys.time(), "] Calculating additional QC metrics ...\n")

# Function to add percentage metrics for gene sets
add_gene_set_metrics <- function(obj, genes, col_name) {
  genes_present <- intersect(genes, rownames(obj))
  obj <- PercentageFeatureSet(obj, features = genes_present, col.name = col_name)
  return(obj)
}

# ADP_C
ADPC_data <- add_gene_set_metrics(ADPC_data, apoptosis_genes, "percent.apoptosis")
ADPC_data <- add_gene_set_metrics(ADPC_data, rbc_genes, "percent.rbc")
ADPC_data <- add_gene_set_metrics(ADPC_data, ieg_genes, "percent.ieg")

# ADP_N
ADPN_data <- add_gene_set_metrics(ADPN_data, apoptosis_genes, "percent.apoptosis")
ADPN_data <- add_gene_set_metrics(ADPN_data, rbc_genes, "percent.rbc")
ADPN_data <- add_gene_set_metrics(ADPN_data, ieg_genes, "percent.ieg")

# ADP_H
ADPH_data <- add_gene_set_metrics(ADPH_data, apoptosis_genes, "percent.apoptosis")
ADPH_data <- add_gene_set_metrics(ADPH_data, rbc_genes, "percent.rbc")
ADPH_data <- add_gene_set_metrics(ADPH_data, ieg_genes, "percent.ieg")

# =============================================
# 6. Filter cells
# =============================================

cat("[", Sys.time(), "] Filtering cells ...\n")

# Define filtering function
filter_cells <- function(obj, sample_name) {
  cat("  Before filtering (", sample_name, "): ", ncol(obj), " cells\n", sep = "")
  
  obj_filtered <- subset(obj,
                         subset = nFeature_RNA >= MIN_FEATURES &
                           nFeature_RNA <= MAX_FEATURES &
                           percent.mt <= MAX_MT &
                           percent.ribo <= MAX_RIBO &
                           percent.apoptosis <= MAX_APOPTOSIS &
                           percent.rbc <= MAX_RBC
  )
  
  cat("  After filtering (", sample_name, "): ", ncol(obj_filtered), " cells\n", sep = "")
  cat("  Removed: ", ncol(obj) - ncol(obj_filtered), " cells (",
      round((ncol(obj) - ncol(obj_filtered)) / ncol(obj) * 100, 1), "%)\n", sep = "")
  
  return(obj_filtered)
}

ADPC_data_filtered <- filter_cells(ADPC_data, "ADP_C")
ADPN_data_filtered <- filter_cells(ADPN_data, "ADP_N")
ADPH_data_filtered <- filter_cells(ADPH_data, "ADP_H")

# Save filtered objects
cat("[", Sys.time(), "] Saving filtered Seurat objects ...\n")
saveRDS(ADPC_data_filtered, file = file.path(OUTPUT_DIR, "ADPC_data_filtered.rds"))
saveRDS(ADPN_data_filtered, file = file.path(OUTPUT_DIR, "ADPN_data_filtered.rds"))
saveRDS(ADPH_data_filtered, file = file.path(OUTPUT_DIR, "ADPH_data_filtered.rds"))

# =============================================
# 7. SCTransform normalization
# =============================================

cat("[", Sys.time(), "] Running SCTransform ...\n")

# Set parallelization settings (adjust based on your system)
plan("sequential")
options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB

ADPC_sct <- SCTransform(
  ADPC_data_filtered,
  vars.to.regress = REGRESS_VARS,
  verbose = SCT_VERBOSE
)

ADPN_sct <- SCTransform(
  ADPN_data_filtered,
  vars.to.regress = REGRESS_VARS,
  verbose = SCT_VERBOSE
)

ADPH_sct <- SCTransform(
  ADPH_data_filtered,
  vars.to.regress = REGRESS_VARS,
  verbose = SCT_VERBOSE
)

# =============================================
# 8. Save SCTransform-normalized objects
# =============================================

cat("[", Sys.time(), "] Saving SCTransform objects ...\n")

saveRDS(ADPC_sct, file = file.path(OUTPUT_DIR, "ADPC_filtered_SCT.rds"))
saveRDS(ADPN_sct, file = file.path(OUTPUT_DIR, "ADPN_filtered_SCT.rds"))
saveRDS(ADPH_sct, file = file.path(OUTPUT_DIR, "ADPH_filtered_SCT.rds"))

# =============================================
# 9. Summary
# =============================================

cat("\n[", Sys.time(), "] ==========================================\n")
cat("QC Summary:\n")
cat("  ADP_C: ", ncol(ADPC_sct), " cells retained\n")
cat("  ADP_N: ", ncol(ADPN_sct), " cells retained\n")
cat("  ADP_H: ", ncol(ADPH_sct), " cells retained\n")
cat("==========================================\n")

# Record session information
session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(OUTPUT_DIR, "session_info_01_QC.txt"))

cat("[", Sys.time(), "] All done!\n")