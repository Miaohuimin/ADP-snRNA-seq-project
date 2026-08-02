#!/usr/bin/env Rscript
# =============================================
# Remove doublets using DoubletFinder
# Usage: Rscript scripts/R/02_remove_doublets.R
# Note: Run this script in a dedicated conda environment:
#       conda create -n r-doubletfinder r-base r-seurat -c conda-forge -y
#       conda activate r-doubletfinder
#       Rscript scripts/R/02_remove_doublets.R
# =============================================

# =============================================
# Load required packages
# =============================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
})

# =============================================
# User-configurable parameters
# =============================================

# Set working directory to project root
# Note: Change this to your actual project path if running outside RStudio
# setwd("/path/to/your/project")

# Input and output directories
DATA_DIR <- "./data/processed"

# DoubletFinder parameters
PCs_USED <- 1:20
PN_VALUE <- 0.25          # pN: artificial doublet proportion (usually 0.25)
SCT_ASSAY <- TRUE         # Use SCT assay
DOUBLET_RATE_PER_1000 <- 0.008  # 10x estimated doublet rate (0.8% per 1000 cells)

# Memory limit
options(future.globals.maxSize = 2000 * 1024^2)  # 2 GB

# =============================================
# Create output directory if it doesn't exist
# =============================================

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================
# Define samples and their manually determined pK values
# (Based on paramSweep results from original analysis)
# =============================================

samples <- list(
  ADP_C = list(file = "ADPC_filtered_SCT.rds", pK = 0.09),
  ADP_H = list(file = "ADPH_filtered_SCT.rds", pK = 0.17),
  ADP_N = list(file = "ADPN_filtered_SCT.rds", pK = 0.01)
)

# =============================================
# Function to run DoubletFinder on a single sample
# =============================================

run_doubletfinder <- function(obj, pK_value, sample_name) {
  
  cat("\n[", Sys.time(), "] Processing sample:", sample_name, "\n")
  cat("  Cells before DoubletFinder:", ncol(obj), "\n")
  
  # Step 1: Run paramSweep to find optimal pK (for validation)
  cat("  Running paramSweep ...\n")
  sweep.res <- paramSweep(obj, PCs = PCs_USED, sct = SCT_ASSAY)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  
  # Step 2: Check if pN=0.25 exists in sweep results
  pN25 <- subset(sweep.stats, pN == PN_VALUE)
  if (nrow(pN25) == 0) {
    cat("  WARNING: pN =", PN_VALUE, "not found in sweep results. Using default pK =", pK_value, "\n")
  } else {
    # Sort by BCreal to confirm best pK
    pN25 <- pN25[order(pN25$BCreal, decreasing = TRUE), ]
    cat("  Top pK values from paramSweep:\n")
    print(head(pN25[, c("pK", "BCreal")]))
    cat("  Using manually selected pK =", pK_value, "\n")
  }
  
  # Step 3: Calculate expected doublet rate
  nCells <- ncol(obj)
  doublet_rate <- DOUBLET_RATE_PER_1000 * (nCells / 1000)
  nExp_poi <- round(nCells * doublet_rate)
  cat("  Expected doublets:", nExp_poi, "(", round(doublet_rate * 100, 2), "% )\n")
  
  # Step 4: Run DoubletFinder
  cat("  Running DoubletFinder with pK =", pK_value, "...\n")
  obj <- doubletFinder(
    obj,
    PCs = PCs_USED,
    pN = PN_VALUE,
    pK = pK_value,
    nExp = nExp_poi,
    sct = SCT_ASSAY
  )
  
  # Step 5: Find the classification column (auto-detected)
  doublet_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  if (length(doublet_col) == 0) {
    stop("ERROR: No DF.classifications column found in metadata!")
  }
  cat("  Classification column:", doublet_col, "\n")
  
  # Step 6: Print doublet classification summary
  cat("  Doublet classification summary:\n")
  print(table(obj[[doublet_col]]))
  
  # Step 7: Filter out doublets
  cells_to_keep <- obj[[doublet_col]] == "Singlet"
  obj_clean <- obj[, cells_to_keep]
  cat("  Cells after removing doublets:", ncol(obj_clean), "\n")
  cat("  Removed:", ncol(obj) - ncol(obj_clean), "cells (",
      round((ncol(obj) - ncol(obj_clean)) / ncol(obj) * 100, 1), "%)\n")
  
  return(obj_clean)
}

# =============================================
# Main loop: Process all samples
# =============================================

cat("[", Sys.time(), "] Starting DoubletFinder pipeline ...\n")
cat("==========================================\n")

results <- list()

for (sample_name in names(samples)) {
  
  # Read the SCTransform-normalized RDS file
  input_file <- file.path(DATA_DIR, samples[[sample_name]]$file)
  
  if (!file.exists(input_file)) {
    cat("[", Sys.time(), "] WARNING: File not found:", input_file, "- Skipping\n")
    next
  }
  
  cat("[", Sys.time(), "] Loading:", input_file, "\n")
  obj <- readRDS(input_file)
  
  # Run DoubletFinder
  obj_clean <- run_doubletfinder(
    obj = obj,
    pK_value = samples[[sample_name]]$pK,
    sample_name = sample_name
  )
  
  # Save cleaned object
  output_file <- file.path(DATA_DIR, paste0(sample_name, "_clean.rds"))
  saveRDS(obj_clean, file = output_file)
  cat("  Saved to:", output_file, "\n")
  
  results[[sample_name]] <- obj_clean
}

# =============================================
# Summary
# =============================================

cat("\n[", Sys.time(), "] ==========================================\n")
cat("DoubletFinder pipeline completed!\n")
cat("Summary of cleaned cells per sample:\n")

for (sample_name in names(results)) {
  if (!is.null(results[[sample_name]])) {
    cat("  ", sample_name, ":", ncol(results[[sample_name]]), "cells\n")
  }
}

# Save session information
session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(DATA_DIR, "session_info_02_doubletfinder.txt"))

cat("[", Sys.time(), "] All done!\n")