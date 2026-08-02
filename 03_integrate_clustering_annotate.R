# =============================================
# Integration, clustering, annotation and visualization
# Usage: Rscript scripts/R/04_integrate_annotate.R
# =============================================

# Clear environment and load packages
rm(list=ls())
library(Seurat)
library(harmony)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(clustree)
library(tidydr)
library(lisi)

# =============================================
# User-configurable parameters
# =============================================

# Set working directory (modify to your actual path)
# setwd("/data/miuhuimin/R_analysis")

# Create output directories
dir.create("./results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/tables", recursive = TRUE, showWarnings = FALSE)

# =============================================
# 1. Load clean data (after DoubletFinder)
# =============================================

cat("[", Sys.time(), "] Loading clean data ...\n")

ADP_N_clean <- readRDS("ADP_N_clean")
ADP_C_clean <- readRDS("ADP_C_clean")
ADP_H_clean <- readRDS("ADP_H_clean") 

cat("  ADP_N:", ncol(ADP_N_clean), "cells\n")
cat("  ADP_C:", ncol(ADP_C_clean), "cells\n")
cat("  ADP_H:", ncol(ADP_H_clean), "cells\n")

# =============================================
# 2. Integrate data with Harmony
# =============================================

cat("[", Sys.time(), "] Integrating data ...\n")

# Prepare integration features
ADP_list <- list(ADP_N_clean, ADP_C_clean, ADP_H_clean)
for (i in seq_along(ADP_list)) {
  DefaultAssay(ADP_list[[i]]) <- "SCT"
}

integ_features <- SelectIntegrationFeatures(object.list = ADP_list, nfeatures = 3000)

# Merge objects
ADP_combined <- merge(
  ADP_N_clean,
  y = c(ADP_C_clean, ADP_H_clean),
  add.cell.ids = c("ADPN", "ADPC", "ADPH"),
  project = "ADP_combined"
)

DefaultAssay(ADP_combined) <- "SCT"
VariableFeatures(ADP_combined) <- integ_features

# Run PCA
ADP_combined <- RunPCA(
  ADP_combined,
  assay = "SCT",
  layer = "data",
  npcs = 50,
  features = VariableFeatures(ADP_combined),
  verbose = FALSE
)

# Run Harmony
ADP_combined <- RunHarmony(
  ADP_combined,
  group.by.vars = "orig.ident",
  reduction = "pca",
  reduction.save = "harmony"
)

# UMAP and clustering based on Harmony
ADP_combined <- RunUMAP(
  ADP_combined,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

ADP_combined <- FindNeighbors(
  ADP_combined,
  reduction = "harmony",
  dims = 1:30
)

ADP_combined <- FindClusters(
  ADP_combined,
  resolution = 0.1
)

# =============================================
# 2a. Initial UMAP visualizations (saved)
# =============================================

cat("[", Sys.time(), "] Generating initial UMAP plots ...\n")

# UMAP by sample
p1 <- DimPlot(ADP_combined, group.by = "orig.ident", shuffle = TRUE)
print(p1)
ggsave("./results/figures/UMAP_by_sample.pdf", plot = p1, width = 8, height = 6)

# UMAP split by sample
p2 <- DimPlot(ADP_combined, group.by = "orig.ident", split.by = "orig.ident")
print(p2)
ggsave("./results/figures/UMAP_split_by_sample.pdf", plot = p2, width = 15, height = 5)

# UMAP by cluster
p3 <- DimPlot(ADP_combined, group.by = "seurat_clusters", label = TRUE)
print(p3)
ggsave("./results/figures/UMAP_by_cluster.pdf", plot = p3, width = 8, height = 6)

# Save integrated object
saveRDS(ADP_combined, file = "ADP_combined_SCT_harmony.rds")

# =============================================
# 3. Find all markers for each cluster
# =============================================

cat("[", Sys.time(), "] Finding marker genes ...\n")

ADP_combined <- PrepSCTFindMarkers(ADP_combined)

markers <- FindAllMarkers(
  ADP_combined,
  assay = "SCT",
  features = VariableFeatures(ADP_combined),
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  densify = TRUE,
  verbose = TRUE
)

cat("  Marker genes found:", nrow(markers), "\n")

# Save marker table
write.csv(markers, file = "./results/tables/cluster_markers_all.csv", row.names = FALSE)

# =============================================
# 4. Cell type annotation
# =============================================

cat("[", Sys.time(), "] Annotating cell types ...\n")

# Define marker genes for each cell type
usual_markers <- list(
  "Astrocyte" = c("Slc4a4", "Aqp4", "Agt", "Slc1a3", "Ntsr2"),
  "Oligodendrocyte" = c("Plp1", "Olig1", "Mog", "Mobp"),
  "OPC" = c("Gpr17", "Pdgfra", "Cspg4", "Ermn"),
  "Microglia" = c("Runx1", "Cx3cr1", "Siglech", "C1qc"),
  "Endothelial cell" = c("Slco1a4", "Cldn5", "Flt1", "Ly6a", "Adgrl4"),
  "Pericyte" = c("Igfbp7", "Rgs5", "Vtn", "Slc38a11"),
  "Neuron" = c("Slc17a6", "Stmn2", "Syt1", "Rbfox3", "Snap25", "Gad1", "Slc32a1"),
  "EGFP" = c("EGFP")
)

# Dot plot for initial marker validation
p_dot_init <- DotPlot(ADP_combined, features = unique(unlist(usual_markers))) + 
  RotatedAxis()
print(p_dot_init)
ggsave("./results/figures/DotPlot_initial_markers.pdf", plot = p_dot_init, width = 12, height = 6)

# Assign cell type labels (adjust cluster numbers based on your data)
celltype_labels <- c(
  "Neuron",           # 0
  "Neuron",           # 1
  "Neuron",           # 2
  "Oligodendrocyte",  # 3
  "Astrocyte",        # 4
  "Neuron",           # 5
  "OPCs",             # 6
  "Neuron",           # 7
  "Microglia",        # 8
  "Neuron",           # 9
  "Neuron",           # 10
  "Mural cell",       # 11
  "Neuron",           # 12
  "Neuron"            # 13
)

ADP_combined$cell_type <- factor(
  ADP_combined$seurat_clusters,
  levels = 0:13,
  labels = celltype_labels
)

# =============================================
# 5. UMAP visualization
# =============================================

cat("[", Sys.time(), "] Generating UMAP plots ...\n")

# ----- 5a. Seurat DimPlot (by cell type) -----
p_celltype <- DimPlot(
  ADP_combined,
  group.by = "cell_type",
  reduction = "umap",
  shuffle = TRUE,
  label = TRUE
)
print(p_celltype)
ggsave("./results/figures/UMAP_by_celltype.pdf", plot = p_celltype, width = 8, height = 6)

# ----- 5b. ggplot UMAP (by cell type) -----
# Define colors for cell types (6 categories)
mycolor_celltype <- c("#c8b8d4", "#729fbc", "#c9bd5b", "#105666", "#e4c7b5", "#80c5b3")

umap_data <- Embeddings(ADP_combined, "umap")
meta_data <- ADP_combined@meta.data
df_plot <- data.frame(umap_data, meta_data)

p_gg <- ggplot() +
  geom_point(
    data = df_plot,
    aes(x = UMAP_1, y = UMAP_2, color = cell_type),
    size = 0.3,
    alpha = 0.4,
    shape = 16
  ) +
  scale_color_manual(values = mycolor_celltype) +
  theme_minimal()
print(p_gg)
# ===== SAVE: ggplot UMAP =====
ggsave("./results/figures/UMAP_by_celltype_ggplot.pdf", plot = p_gg, width = 8, height = 6)

# ----- 5c. ggplot UMAP faceted by sample -----
# Define colors for samples (3 categories)
mycolor_sample <- c("#9c929b", "#ca8e98", "#c9bd5b")

p_facet <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.3, alpha = 0.4, shape = 16) +
  scale_color_manual(values = mycolor_celltype) +
  facet_wrap(~ orig.ident, ncol = 3) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
print(p_facet)
# ===== SAVE: Facet UMAP =====
ggsave("./results/figures/UMAP_facet_by_sample.pdf", plot = p_facet, width = 15, height = 5)

# Show sample cell counts
table(ADP_combined$orig.ident)

# =============================================
# 6. iLISI score (integration quality metric)
# =============================================

cat("[", Sys.time(), "] Calculating iLISI scores ...\n")

# Extract PCA coordinates and batch info
pca_data <- Embeddings(ADP_combined, reduction = "pca")
batch_data <- ADP_combined$orig.ident

# Calculate iLISI
ilisi_scores <- compute_lisi(
  pca_data,
  data.frame(batch = batch_data),
  label_colnames = "batch"
)
ADP_combined$iLISI <- ilisi_scores$batch

# Summary statistics
cat("  iLISI median:", median(ADP_combined$iLISI), "\n")
cat("  iLISI mean:", mean(ADP_combined$iLISI), "\n")
cat("  iLISI summary:\n")
print(summary(ADP_combined$iLISI))

# Save iLISI scores to table
ilisi_df <- data.frame(
  cell_id = colnames(ADP_combined),
  iLISI = ADP_combined$iLISI,
  sample = ADP_combined$orig.ident
)
write.csv(ilisi_df, file = "./results/tables/iLISI_scores.csv", row.names = FALSE)

# ----- 6a. FeaturePlot iLISI (Seurat) -----
p_ilisi_feature <- FeaturePlot(
  ADP_combined,
  features = "iLISI",
  reduction = "umap",
  label = TRUE,
  order = TRUE
) + scale_color_viridis_c(option = "plasma")
print(p_ilisi_feature)
# ===== SAVE: FeaturePlot iLISI =====
ggsave("./results/figures/iLISI_FeaturePlot.pdf", plot = p_ilisi_feature, width = 8, height = 6)

# ----- 6b. ggplot iLISI -----
umap_coords <- Embeddings(ADP_combined, reduction = "umap")
df_ilisi <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  iLISI = ADP_combined$iLISI
)

p_ilisi_gg <- ggplot(df_ilisi, aes(x = UMAP_1, y = UMAP_2, color = iLISI)) +
  geom_point(size = 0.3, alpha = 0.4, shape = 16) +
  scale_color_gradient(low = "lightblue", high = "darkblue", name = "iLISI") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  labs(title = "iLISI score (integration of three samples)")
print(p_ilisi_gg)
ggsave("./results/figures/iLISI_UMAP_ggplot.pdf", plot = p_ilisi_gg, width = 8, height = 6)

# Save annotated object
saveRDS(ADP_combined, file = "ADP_harmony_combined_with_celltype.rds")

# =============================================
# 7. Final DotPlot with ordered cell types
# =============================================

cat("[", Sys.time(), "] Generating final dot plots ...\n")

# Define gene order for final dot plot
gene_order <- c(
  "Slc1a3", "Slc4a4", "Ntsr2", "Agt",                    # Astrocyte
  "Plp1", "Mobp", "Mog", "Ermn", "Olig1",                # Oligodendrocyte
  "Gpr17", "Pdgfra", "Cspg4",                            # OPCs
  "C1qc", "Runx1", "Cx3cr1", "Siglech",                  # Microglia
  "Flt1", "Slco1a4", "Mecom", "Vtn",                     # Mural cell
  "Slc17a6", "Stmn2", "Syt1", "Rbfox3", "Snap25", "Gad1", "Slc32a1"  # Neuron
)

cell_order <- c("Astrocyte", "Oligodendrocyte", "OPCs", "Microglia",
                "Mural cell", "Neuron")

ADP_combined$cell_type <- factor(ADP_combined$cell_type, levels = cell_order)

# ----- 7a. DotPlot (vertical) -----
p_dot <- DotPlot(
  ADP_combined,
  features = unique(unlist(gene_order)),
  group.by = "cell_type"
) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradient2(low = "white", high = "#367DB0", midpoint = 0)
print(p_dot)
ggsave("./results/figures/DotPlot_marker_genes.pdf", plot = p_dot, width = 14, height = 6)

# ----- 7b. DotPlot (horizontal) -----
p_dot_h <- DotPlot(
  ADP_combined,
  features = unique(unlist(gene_order)),
  group.by = "cell_type"
) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0.5),
    axis.text.y = element_text(size = 8)
  ) +
  scale_color_gradient2(low = "white", high = "#367DB0", midpoint = 0)
print(p_dot_h)
# ===== SAVE: Horizontal dot plot =====
ggsave("./results/figures/DotPlot_marker_genes_horizontal.pdf", plot = p_dot_h, width = 12, height = 8)

# =============================================
# 8. Save final object and session info
# =============================================

cat("[", Sys.time(), "] Saving final object ...\n")
saveRDS(ADP_combined, file = "ADP_combined_annotated_final.rds")

# Save session information
session_info <- capture.output(sessionInfo())
writeLines(session_info, "./results/session_info_04_integration.txt")

# =============================================
# 9. Summary
# =============================================

cat("\n[", Sys.time(), "] ==========================================\n")
cat("Integration pipeline completed!\n")
cat("  Total cells:", ncol(ADP_combined), "\n")
cat("  Number of clusters:", length(unique(ADP_combined$seurat_clusters)), "\n")
cat("  Cell types:", paste(unique(ADP_combined$cell_type), collapse = ", "), "\n")
cat("  iLISI median:", median(ADP_combined$iLISI), "\n")
cat("  Figures saved to: ./results/figures/\n")
cat("  Tables saved to: ./results/tables/\n")
cat("==========================================\n")
cat("[", Sys.time(), "] All done!\n")