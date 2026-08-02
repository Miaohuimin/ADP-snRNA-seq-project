#!/usr/bin/env Rscript
# =============================================
# EGFP+ neuron sub-clustering and differential expression
# Usage: Rscript scripts/R/04_unstimulated_EGFP_neurons_analysis.R
# =============================================

# =============================================
# 1. Load packages and set paths
# =============================================

rm(list=ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  library(patchwork)
  library(clustree)
  library(ggrepel)
  library(future)
})

# Set working directory (modify to your actual path)
# setwd("/data/miuhuimin/R_analysis")

# Create output directories
dir.create("./results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/tables", recursive = TRUE, showWarnings = FALSE)

# Set memory limit
options(future.globals.maxSize = 5 * 1024^3)  # 5 GB
plan("sequential")

# =============================================
# 2. Define color palette
# =============================================

mycolor22 <- c(
  "#4E79A7", "#83b5b5", "#839958", "#9887bc", "#deb956", "#bfc5d5",
  "#c8b8d4", "#7c7a7d", "#c1d09d", "#344B5C", "#C98F96",
  "#e08ea4", "#9dbdd2", "#779ebd", "#bdbb55",
  "#d0c9b0", "#b696b6", "#C8DEF9", "#80c1c4", "#5d3c11", "#7c7a7d", "#FDDF91"
)

# =============================================
# 3. Load data and extract neurons
# =============================================

cat("[", Sys.time(), "] Loading data ...\n")

ADP_combined <- readRDS("ADP_harmony_combined_with_celltype.rds")

# Extract all neurons
ADP_neurons <- subset(ADP_combined, subset = cell_type %in% c("Neuron"))
cat("  Total neurons:", ncol(ADP_neurons), "\n")

# Extract ADP_N sample
ADP_N <- subset(ADP_neurons, subset = orig.ident %in% c("ADP_N"))
cat("  ADP_N neurons:", ncol(ADP_N), "\n")

# =============================================
# 4. Add EGFP counts and classify cells
# =============================================

cat("[", Sys.time(), "] Classifying EGFP+ cells ...\n")

egfp_gene <- "EGFP"
threshold <- 1

# Function to add EGFP metadata
add_egfp_metadata <- function(obj, threshold = 1) {
  egfp_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[egfp_gene, ]
  obj$EGFP_counts <- egfp_counts
  obj$EGFP_positive <- ifelse(obj$EGFP_counts >= threshold, "Positive", "Negative")
  return(obj)
}

ADP_N <- add_egfp_metadata(ADP_N, threshold)
cat("  ADP_N EGFP+ cells:", sum(ADP_N$EGFP_positive == "Positive"), "\n")
cat("  ADP_N EGFP- cells:", sum(ADP_N$EGFP_positive == "Negative"), "\n")

# =============================================
# 5. Differential expression: EGFP+ vs EGFP-
# =============================================

cat("[", Sys.time(), "] Running differential expression ...\n")

DefaultAssay(ADP_N) <- "RNA"
ADP_N <- SCTransform(ADP_N, vars.to.regress = "percent.mt", verbose = TRUE)
DefaultAssay(ADP_N) <- "SCT"
ADP_N <- PrepSCTFindMarkers(ADP_N)

deg_res <- FindMarkers(
  object = ADP_N,
  ident.1 = "Positive",
  ident.2 = "Negative",
  group.by = "EGFP_positive",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  min.cells.group = 0,
  only.pos = FALSE
)

cat("  DEGs found:", nrow(deg_res), "\n")

# Save DEG results
deg_df <- deg_res %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")
write.csv(deg_df, "./results/tables/ADP_N_EGFP_pos_vs_neg_DEGs.csv", row.names = FALSE)

# =============================================
# 6. Volcano plot
# =============================================

cat("[", Sys.time(), "] Generating volcano plot ...\n")

deg_res <- deg_res %>%
  mutate(
    log10_p = -log10(p_val),
    sig = case_when(
      avg_log2FC > 0.25 & p_val < 0.05 ~ "Up",
      avg_log2FC < -0.25 & p_val < 0.05 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

# User-specified genes to label
user_genes <- c("Ndst4", "Sema3c", "Rftn1", "Nos1", "Dgkk")
genes_to_label <- user_genes[user_genes %in% rownames(deg_res)]

cat("  Genes to label:", paste(genes_to_label, collapse = ", "), "\n")
cat("  Genes not found:", paste(setdiff(user_genes, genes_to_label), collapse = ", "), "\n")

deg_res$gene_label <- if_else(rownames(deg_res) %in% genes_to_label,
                              rownames(deg_res),
                              NA_character_)

colors_volcano <- c("Up" = "#839958", "Down" = "gray80", "Not significant" = "gray80")

volcano <- ggplot(deg_res, aes(x = avg_log2FC, y = log10_p, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = colors_volcano, name = "Expression") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(
    aes(label = gene_label),
    size = 4,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = 30,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  coord_cartesian(xlim = c(-2, 2)) +
  labs(
    title = "Volcano Plot: EGFP+ vs EGFP- (ADP_N)",
    x = "Log2 Fold Change",
    y = "-Log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    plot.title = element_text(hjust = 0.5, color = "black")
  )

print(volcano)
ggsave("./results/figures/Volcano_EGFP_vs_nonEGFP.pdf", plot = volcano, width = 8, height = 6)

# =============================================
# 7. Extract EGFP+ neurons for sub-clustering
# =============================================

cat("[", Sys.time(), "] Extracting EGFP+ neurons ...\n")

ADP_N_pos <- subset(ADP_N, subset = EGFP_positive == "Positive")
cat("  EGFP+ neurons:", ncol(ADP_N_pos), "\n")

# =============================================
# 8. SCTransform on EGFP+ neurons (excluding EGFP)
# =============================================

cat("[", Sys.time(), "] Running SCTransform on EGFP+ neurons ...\n")

DefaultAssay(ADP_N_pos) <- "RNA"

ADP_N_pos <- SCTransform(
  ADP_N_pos,
  vars.to.regress = "percent.mt",
  variable.features.n = 3000,
  return.only.var.genes = FALSE,
  verbose = TRUE
)

# Remove EGFP from variable features
var_genes <- VariableFeatures(ADP_N_pos)
var_genes_no_egfp <- setdiff(var_genes, "EGFP")
VariableFeatures(ADP_N_pos) <- var_genes_no_egfp

cat("  Variable features (excluding EGFP):", length(VariableFeatures(ADP_N_pos)), "\n")

# =============================================
# 9. Dimensionality reduction and clustering
# =============================================

cat("[", Sys.time(), "] Running PCA and clustering ...\n")

DefaultAssay(ADP_N_pos) <- "SCT"

ADP_N_pos <- RunPCA(ADP_N_pos, features = VariableFeatures(ADP_N_pos), verbose = FALSE)
ADP_N_pos <- FindNeighbors(ADP_N_pos, dims = 1:30, verbose = FALSE)

# Test multiple resolutions (0.1 to 0.5)
resolutions <- seq(0.1, 0.5, by = 0.1)
for (res in resolutions) {
  ADP_N_pos <- FindClusters(ADP_N_pos, resolution = res, verbose = FALSE)
}

# =============================================
# 10. Clustree (optional, can be commented out)
# =============================================

p_clustree <- clustree(ADP_N_pos, prefix = "SCT_snn_res.")
print(p_clustree)
ggsave("./results/figures/Clustree_ADP_N_pos.pdf", plot = p_clustree, width = 10, height = 8)

# =============================================
# 11. Find marker genes (resolution = 0.4)
# =============================================

cat("[", Sys.time(), "] Finding marker genes ...\n")

Idents(ADP_N_pos) <- "SCT_snn_res.0.4"
cat("  Clusters at res 0.4:", paste(sort(unique(Idents(ADP_N_pos))), collapse = ", "), "\n")
cat("  Cluster sizes:\n")
print(table(Idents(ADP_N_pos)))

ADP_N_pos <- PrepSCTFindMarkers(ADP_N_pos, verbose = TRUE)

ADP_N_pos_all_markers <- FindAllMarkers(
  object = ADP_N_pos,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0,
  test.use = "wilcox",
  verbose = TRUE
)

cat("  Total marker genes found:", nrow(ADP_N_pos_all_markers), "\n")

# Filter significant markers
sig_markers <- subset(ADP_N_pos_all_markers, p_val_adj < 0.05)

# Top 20 markers per cluster
top20_per_cluster <- sig_markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) %>%
  ungroup()

# Save marker results
write.csv(ADP_N_pos_all_markers, "./results/tables/ADP_N_pos_all_markers.csv", row.names = TRUE)
write.csv(top20_per_cluster, "./results/tables/ADP_N_pos_top20_markers.csv", row.names = FALSE)

# =============================================
# 12. DotPlot
# =============================================

cat("[", Sys.time(), "] Generating dot plot ...\n")

genes_to_plot <- c(
  "Lars2", "Pou2f2", "Ntng1", "Ptprk",
  "Sox6", "Stxbp6", "Lhx6", "Gad1",
  "Reln", "Slc17a6", "Nos1", "Trp73",
  "Meis2", "Rarb", "Wfs1", "Chrm1",
  "Zeb2", "Foxp2", "Oprm1", "Penk",
  "Pde11a", "Lhfp", "Kcnh8", "Lhx8",
  "Prdm16", "Gfra1", "Nfia", "Cacng5",
  "Tshz1", "Chst9", "Lypd1", "Eya2",
  "Prlr", "Greb1", "Dlx1", "Esr2",
  "Chat", "Slc5a7", "Prima1", "Ngfr"
)

# Define cluster labels (adjust based on actual cluster count)
n_clusters <- length(unique(Idents(ADP_N_pos)))
cluster_labels <- paste0("ADP_Neu", sprintf("%02d", n_clusters:1))

# Assign cell_type labels
ADP_N_pos$cell_type <- factor(
  ADP_N_pos$SCT_snn_res.0.4,
  levels = 0:(n_clusters - 1),
  labels = cluster_labels
)

# Order cell types
cell_order <- rev(cluster_labels)
ADP_N_pos$cell_type <- factor(ADP_N_pos$cell_type, levels = cell_order)

# Horizontal dot plot
p_dot <- DotPlot(
  ADP_N_pos,
  features = unique(genes_to_plot),
  group.by = "cell_type"
) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0.5),
    axis.text.y = element_text(size = 8)
  ) +
  scale_color_gradient2(low = "white", high = "#367DB0", midpoint = 0)

print(p_dot)
ggsave("./results/figures/DotPlot_ADP_N_pos.pdf", plot = p_dot, width = 14, height = 8)

# =============================================
# 13. UMAP visualization
# =============================================

cat("[", Sys.time(), "] Generating UMAP plot ...\n")

ADP_N_pos <- FindClusters(ADP_N_pos, resolution = 0.4)

ADP_N_pos <- RunUMAP(
  ADP_N_pos,
  reduction = "pca",
  reduction.name = "umap",
  dims = 1:30,
  spread = 0.3,
  min.dist = 0.1,
  n.neighbors = 30,
  verbose = FALSE
)

# Extract UMAP coordinates
umap_data <- Embeddings(ADP_N_pos, "umap")
meta_data <- ADP_N_pos@meta.data
df_plot <- data.frame(umap_data, meta_data)

# Calculate centroids
centroids <- df_plot %>%
  group_by(seurat_clusters) %>%
  summarise(umap_1 = median(UMAP_1), umap_2 = median(UMAP_2))

# Match cluster labels to centroids
centroids$label <- cluster_labels[as.numeric(as.character(centroids$seurat_clusters)) + 1]

# Use mycolor22 for coloring
n_clusters_umap <- length(unique(df_plot$seurat_clusters))
umap_colors <- mycolor22[1:n_clusters_umap]

# Plot
p_umap <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = seurat_clusters), size = 0.6, alpha = 0.8, shape = 16) +
  scale_color_manual(values = umap_colors) +
  geom_text(
    data = centroids,
    aes(label = label),
    size = 3,
    color = "black",
    fontface = "bold"
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right"
  )

print(p_umap)
ggsave("./results/figures/UMAP_ADP_N_pos.pdf", plot = p_umap, width = 10, height = 8)

# =============================================
# 14. Save final object and session info
# =============================================

cat("[", Sys.time(), "] Saving final object ...\n")
saveRDS(ADP_N_pos, file = "ADP_N_pos.rds")

session_info <- capture.output(sessionInfo())
writeLines(session_info, "./results/session_info_05_EGFP_neurons.txt")

# =============================================
# 15. Summary
# =============================================

cat("\n[", Sys.time(), "] ==========================================\n")
cat("EGFP+ neuron analysis completed!\n")
cat("  Total EGFP+ neurons:", ncol(ADP_N_pos), "\n")
cat("  Number of clusters (res 0.4):", length(unique(ADP_N_pos$seurat_clusters)), "\n")
cat("  DEGs found:", nrow(deg_res), "\n")
cat("  Figures saved to: ./results/figures/\n")
cat("  Tables saved to: ./results/tables/\n")
cat("==========================================\n")
cat("[", Sys.time(), "] All done!\n")