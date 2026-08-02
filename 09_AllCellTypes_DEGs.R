#!/usr/bin/env Rscript
# =============================================
# Volcano plots for all major cell types
# Generates volcano plots for 6 cell types under heat and cold stimulation
# Usage: Rscript scripts/R/09_AllCellTypes_DEGs.R
# =============================================

# =============================================
# 1. Load packages
# =============================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# =============================================
# 2. Check if ADP_combined_pos exists
# =============================================

if (!exists("ADP_combined_pos")) {
  stop("Object ADP_combined_pos not found. Please run previous analysis first.")
}

cat("[", Sys.time(), "] ADP_combined_pos loaded successfully.\n")
cat("  Total cells:", ncol(ADP_combined_pos), "\n")
cat("  Cell types:\n")
print(table(ADP_combined_pos$cell_type))

# =============================================
# 3. Create all cell type subset objects if missing
# =============================================

cat("[", Sys.time(), "] Checking and creating cell type objects...\n")

# ===== Astrocyte =====
if (!exists("ADP_Astrocyte_pos")) {
  cat("  ADP_Astrocyte_pos not found, creating from ADP_combined_pos...\n")
  ADP_Astrocyte_pos <- subset(ADP_combined_pos, subset = cell_type == "Astrocyte")
  if (ncol(ADP_Astrocyte_pos) == 0) {
    stop("No cells with cell_type == 'Astrocyte' found.")
  }
  DefaultAssay(ADP_Astrocyte_pos) <- "RNA"
  ADP_Astrocyte_pos <- SCTransform(ADP_Astrocyte_pos, 
                                   vars.to.regress = "percent.mt", 
                                   variable.features.n = 3000,
                                   return.only.var.genes = FALSE,
                                   verbose = FALSE)
  cat("    ADP_Astrocyte_pos created:", ncol(ADP_Astrocyte_pos), "cells\n")
} else {
  cat("  ADP_Astrocyte_pos already exists:", ncol(ADP_Astrocyte_pos), "cells\n")
}

# ===== OPCs =====
if (!exists("ADP_OPCs_pos")) {
  cat("  ADP_OPCs_pos not found, creating from ADP_combined_pos...\n")
  ADP_OPCs_pos <- subset(ADP_combined_pos, subset = cell_type == "OPCs")
  if (ncol(ADP_OPCs_pos) == 0) {
    stop("No cells with cell_type == 'OPCs' found.")
  }
  DefaultAssay(ADP_OPCs_pos) <- "RNA"
  ADP_OPCs_pos <- SCTransform(ADP_OPCs_pos, 
                              vars.to.regress = "percent.mt", 
                              variable.features.n = 3000,
                              return.only.var.genes = FALSE,
                              verbose = FALSE)
  cat("    ADP_OPCs_pos created:", ncol(ADP_OPCs_pos), "cells\n")
} else {
  cat("  ADP_OPCs_pos already exists:", ncol(ADP_OPCs_pos), "cells\n")
}

# ===== Oligodendrocyte =====
if (!exists("ADP_Oligo_pos")) {
  cat("  ADP_Oligo_pos not found, creating from ADP_combined_pos...\n")
  ADP_Oligo_pos <- subset(ADP_combined_pos, subset = cell_type == "Oligodendrocyte")
  if (ncol(ADP_Oligo_pos) == 0) {
    stop("No cells with cell_type == 'Oligodendrocyte' found.")
  }
  DefaultAssay(ADP_Oligo_pos) <- "RNA"
  ADP_Oligo_pos <- SCTransform(ADP_Oligo_pos, 
                               vars.to.regress = "percent.mt", 
                               variable.features.n = 3000,
                               return.only.var.genes = FALSE,
                               verbose = FALSE)
  cat("    ADP_Oligo_pos created:", ncol(ADP_Oligo_pos), "cells\n")
} else {
  cat("  ADP_Oligo_pos already exists:", ncol(ADP_Oligo_pos), "cells\n")
}

# ===== Microglia =====
if (!exists("ADP_Microglia_pos")) {
  cat("  ADP_Microglia_pos not found, creating from ADP_combined_pos...\n")
  ADP_Microglia_pos <- subset(ADP_combined_pos, subset = cell_type == "Microglia")
  if (ncol(ADP_Microglia_pos) == 0) {
    stop("No cells with cell_type == 'Microglia' found.")
  }
  DefaultAssay(ADP_Microglia_pos) <- "RNA"
  ADP_Microglia_pos <- SCTransform(ADP_Microglia_pos, 
                                   vars.to.regress = "percent.mt", 
                                   variable.features.n = 3000,
                                   return.only.var.genes = FALSE,
                                   verbose = FALSE)
  cat("    ADP_Microglia_pos created:", ncol(ADP_Microglia_pos), "cells\n")
} else {
  cat("  ADP_Microglia_pos already exists:", ncol(ADP_Microglia_pos), "cells\n")
}

# ===== Mural cell =====
if (!exists("ADP_Mural_pos")) {
  cat("  ADP_Mural_pos not found, creating from ADP_combined_pos...\n")
  ADP_Mural_pos <- subset(ADP_combined_pos, subset = cell_type == "Mural cell")
  if (ncol(ADP_Mural_pos) == 0) {
    stop("No cells with cell_type == 'Mural cell' found.")
  }
  DefaultAssay(ADP_Mural_pos) <- "RNA"
  ADP_Mural_pos <- SCTransform(ADP_Mural_pos, 
                               vars.to.regress = "percent.mt", 
                               variable.features.n = 3000,
                               return.only.var.genes = FALSE,
                               verbose = FALSE)
  cat("    ADP_Mural_pos created:", ncol(ADP_Mural_pos), "cells\n")
} else {
  cat("  ADP_Mural_pos already exists:", ncol(ADP_Mural_pos), "cells\n")
}

# ===== Neuron =====
if (!exists("ADP_Neuron_pos")) {
  cat("  ADP_Neuron_pos not found, creating from ADP_combined_pos...\n")
  ADP_Neuron_pos <- subset(ADP_combined_pos, subset = cell_type == "Neuron")
  if (ncol(ADP_Neuron_pos) == 0) {
    stop("No cells with cell_type == 'Neuron' found.")
  }
  DefaultAssay(ADP_Neuron_pos) <- "RNA"
  ADP_Neuron_pos <- SCTransform(ADP_Neuron_pos, 
                                vars.to.regress = "percent.mt", 
                                variable.features.n = 3000,
                                return.only.var.genes = FALSE,
                                verbose = FALSE)
  cat("    ADP_Neuron_pos created:", ncol(ADP_Neuron_pos), "cells\n")
} else {
  cat("  ADP_Neuron_pos already exists:", ncol(ADP_Neuron_pos), "cells\n")
}

# =============================================
# 4. Verify all objects were created successfully
# =============================================

cat("\n[", Sys.time(), "] Verifying all cell type objects...\n")

required_objs <- c("ADP_Astrocyte_pos", "ADP_OPCs_pos", "ADP_Oligo_pos", 
                   "ADP_Microglia_pos", "ADP_Mural_pos", "ADP_Neuron_pos")

all_exist <- TRUE
for (obj in required_objs) {
  if (exists(obj) && !is.null(get(obj))) {
    cat("  ", obj, "OK -", ncol(get(obj)), "cells\n")
  } else {
    cat("  ", obj, "FAILED - does not exist or is NULL\n")
    all_exist <- FALSE
  }
}

if (!all_exist) {
  stop("Some cell type objects failed to create. Please check cell_type column.")
}

cat("All cell type objects ready!\n\n")

# =============================================
# 5. Define cell types and object list
# =============================================

cell_types <- c("Astrocyte", "OPCs", "Oligo", "Microglia", "Mural", "Neuron")
obj_list <- list(
  Astrocyte = ADP_Astrocyte_pos,
  OPCs      = ADP_OPCs_pos,
  Oligo     = ADP_Oligo_pos,
  Microglia = ADP_Microglia_pos,
  Mural     = ADP_Mural_pos,
  Neuron    = ADP_Neuron_pos
)

# =============================================
# 6. Define parameters
# =============================================

exclude_pattern <- "^(mt-|MT-|Gm|Rps|Rpl|Rp[ls])"
user_genes <- c("Fos", "Per1", "Jun", "Nr4a3", "Homer1", "Gadd45g")

colors_heat <- c("Up" = "#9F044D", "Down" = "#9c929b", "Not significant" = "gray80")
colors_cold <- c("Up" = "#0F99B2", "Down" = "#9c929b", "Not significant" = "gray80")

x_axis_limits <- c(-5, 5)
y_axis_limits <- c(0, 10)

# =============================================
# 7. Create output directories
# =============================================

dir.create("./results/figures/volcano_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/tables", recursive = TRUE, showWarnings = FALSE)

# =============================================
# 8. Loop through each cell type
# =============================================

for (ct in cell_types) {
  cat("\n==================== Processing cell type:", ct, " ====================\n")
  obj <- obj_list[[ct]]
  
  if (!"orig.ident" %in% colnames(obj[[]])) {
    warning(paste(ct, "object missing orig.ident column, skipping"))
    next
  }
  
  DefaultAssay(obj) <- "SCT"
  obj <- PrepSCTFindMarkers(obj)
  
  # ----- Heat stimulation: ADP_H vs ADP_N -----
  if (sum(obj$orig.ident == "ADP_H") >= 3) {
    cat("  Heat DEG analysis (ADP_H vs ADP_N)...\n")
    deg_heat <- FindMarkers(
      object = obj,
      ident.1 = "ADP_H",
      ident.2 = "ADP_N",
      group.by = "orig.ident",
      test.use = "wilcox",
      logfc.threshold = 0,
      min.pct = 0.1,
      min.cells.group = 3,
      only.pos = FALSE
    )
    
    deg_heat <- deg_heat %>%
      mutate(
        log10_p = -log10(p_val),
        sig = case_when(
          avg_log2FC > 0.25 & p_val < 0.05 ~ "Up",
          avg_log2FC < -0.25 & p_val < 0.05 ~ "Down",
          TRUE ~ "Not significant"
        ),
        is_excluded = grepl(exclude_pattern, rownames(deg_heat), ignore.case = TRUE)
      )
    
    # Select top genes for labeling
    deg_up <- deg_heat %>% filter(sig == "Up", !is_excluded)
    top_up_lfc <- deg_up %>% slice_max(abs(avg_log2FC), n = 5)
    top_up_p   <- deg_up %>% slice_min(p_val, n = 5)
    
    deg_down <- deg_heat %>% filter(sig == "Down", !is_excluded)
    top_down_lfc <- deg_down %>% slice_max(abs(avg_log2FC), n = 5)
    top_down_p   <- deg_down %>% slice_min(p_val, n = 5)
    
    label_genes_heat <- unique(c(
      rownames(top_up_lfc), rownames(top_up_p),
      rownames(top_down_lfc), rownames(top_down_p)
    ))
    
    user_genes_sig <- user_genes[user_genes %in% rownames(deg_heat) & 
                                   deg_heat[user_genes, "sig"] %in% c("Up", "Down")]
    label_genes_heat <- unique(c(user_genes_sig, label_genes_heat))
    
    cat("\n  [Heat] Up |log2FC| top5:", 
        if(nrow(top_up_lfc)>0) paste(rownames(top_up_lfc), collapse=", ") else "none", "\n")
    cat("  [Heat] Up p-value top5:", 
        if(nrow(top_up_p)>0) paste(rownames(top_up_p), collapse=", ") else "none", "\n")
    cat("  [Heat] Down |log2FC| top5:", 
        if(nrow(top_down_lfc)>0) paste(rownames(top_down_lfc), collapse=", ") else "none", "\n")
    cat("  [Heat] Down p-value top5:", 
        if(nrow(top_down_p)>0) paste(rownames(top_down_p), collapse=", ") else "none", "\n")
    cat("  [Heat] Final labeled genes:", paste(label_genes_heat, collapse=", "), "\n")
    
    deg_heat <- deg_heat %>%
      mutate(gene_label = if_else(rownames(deg_heat) %in% label_genes_heat & sig %in% c("Up", "Down"),
                                  rownames(deg_heat), NA_character_))
    
    p_heat <- ggplot(deg_heat, aes(x = avg_log2FC, y = log10_p, color = sig)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = colors_heat, name = "Expression") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", alpha = 0.5) +
      geom_text_repel(aes(label = gene_label), size = 3, box.padding = 0.3, max.overlaps = 30, show.legend = FALSE) +
      coord_cartesian(xlim = x_axis_limits, ylim = y_axis_limits) +
      labs(title = paste(ct, ": ADP_H vs ADP_N"),
           x = "Log2 Fold Change", y = "-Log10(P-value)") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5))
    
    print(p_heat)
    
    # Export DEG table (Heat)
    deg_heat_export <- deg_heat %>%
      dplyr::select(-is_excluded) %>%
      tibble::rownames_to_column("gene")
    filename_heat_table <- paste0("./results/tables/DEG_", ct, "_H_vs_N.csv")
    write.csv(deg_heat_export, filename_heat_table, row.names = FALSE)
    cat("  DEG table saved:", filename_heat_table, "\n")
    
    if (interactive()) {
      readline("Press Enter to view next plot...")
    }
    
  } else {
    cat("  Heat group has <3 cells, skipping\n")
  }
  
  # ----- Cold stimulation: ADP_C vs ADP_N -----
  if (sum(obj$orig.ident == "ADP_C") >= 3) {
    cat("  Cold DEG analysis (ADP_C vs ADP_N)...\n")
    deg_cold <- FindMarkers(
      object = obj,
      ident.1 = "ADP_C",
      ident.2 = "ADP_N",
      group.by = "orig.ident",
      test.use = "wilcox",
      logfc.threshold = 0,
      min.pct = 0.1,
      min.cells.group = 3,
      only.pos = FALSE
    )
    
    deg_cold <- deg_cold %>%
      mutate(
        log10_p = -log10(p_val),
        sig = case_when(
          avg_log2FC > 0.25 & p_val < 0.05 ~ "Up",
          avg_log2FC < -0.25 & p_val < 0.05 ~ "Down",
          TRUE ~ "Not significant"
        ),
        is_excluded = grepl(exclude_pattern, rownames(deg_cold), ignore.case = TRUE)
      )
    
    deg_up <- deg_cold %>% filter(sig == "Up", !is_excluded)
    top_up_lfc <- deg_up %>% slice_max(abs(avg_log2FC), n = 5)
    top_up_p   <- deg_up %>% slice_min(p_val, n = 5)
    
    deg_down <- deg_cold %>% filter(sig == "Down", !is_excluded)
    top_down_lfc <- deg_down %>% slice_max(abs(avg_log2FC), n = 5)
    top_down_p   <- deg_down %>% slice_min(p_val, n = 5)
    
    label_genes_cold <- unique(c(
      rownames(top_up_lfc), rownames(top_up_p),
      rownames(top_down_lfc), rownames(top_down_p)
    ))
    
    user_genes_sig <- user_genes[user_genes %in% rownames(deg_cold) & 
                                   deg_cold[user_genes, "sig"] %in% c("Up", "Down")]
    label_genes_cold <- unique(c(user_genes_sig, label_genes_cold))
    
    cat("\n  [Cold] Up |log2FC| top5:", 
        if(nrow(top_up_lfc)>0) paste(rownames(top_up_lfc), collapse=", ") else "none", "\n")
    cat("  [Cold] Up p-value top5:", 
        if(nrow(top_up_p)>0) paste(rownames(top_up_p), collapse=", ") else "none", "\n")
    cat("  [Cold] Down |log2FC| top5:", 
        if(nrow(top_down_lfc)>0) paste(rownames(top_down_lfc), collapse=", ") else "none", "\n")
    cat("  [Cold] Down p-value top5:", 
        if(nrow(top_down_p)>0) paste(rownames(top_down_p), collapse=", ") else "none", "\n")
    cat("  [Cold] Final labeled genes:", paste(label_genes_cold, collapse=", "), "\n")
    
    deg_cold <- deg_cold %>%
      mutate(gene_label = if_else(rownames(deg_cold) %in% label_genes_cold & sig %in% c("Up", "Down"),
                                  rownames(deg_cold), NA_character_))
    
    p_cold <- ggplot(deg_cold, aes(x = avg_log2FC, y = log10_p, color = sig)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = colors_cold, name = "Expression") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", alpha = 0.5) +
      geom_text_repel(aes(label = gene_label), size = 3, box.padding = 0.3, max.overlaps = 30, show.legend = FALSE) +
      coord_cartesian(xlim = x_axis_limits, ylim = y_axis_limits) +
      labs(title = paste(ct, ": ADP_C vs ADP_N"),
           x = "Log2 Fold Change", y = "-Log10(P-value)") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5))
    
    print(p_cold)
    
    # Export DEG table (Cold)
    deg_cold_export <- deg_cold %>%
      dplyr::select(-is_excluded) %>%
      tibble::rownames_to_column("gene")
    filename_cold_table <- paste0("./results/tables/DEG_", ct, "_C_vs_N.csv")
    write.csv(deg_cold_export, filename_cold_table, row.names = FALSE)
    cat("  DEG table saved:", filename_cold_table, "\n")
    
    if (interactive() && ct != tail(cell_types, 1)) {
      readline("Press Enter to view next cell type...")
    }
    
  } else {
    cat("  Cold group has <3 cells, skipping\n")
  }
}

cat("\n==================== All done! Processed", length(cell_types), "cell types ====================\n")