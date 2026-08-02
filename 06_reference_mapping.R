#!/usr/bin/env Rscript
# =============================================
# Reference mapping: project EGFP+ neuron clusters onto control UMAP
# Usage: Rscript scripts/R/06_reference_mapping.R
# =============================================

# =============================================
# 1. Load packages
# =============================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggnewscale)

# =============================================
# 2. Set paths and create output directories
# =============================================

dir.create("./results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/tables", recursive = TRUE, showWarnings = FALSE)

# =============================================
# 3. Read input objects
# =============================================

cat("[", Sys.time(), "] Loading data ...\n")

# EGFP+ neuron sub-clustering result (from script 05)
ADP_E_pos <- readRDS("ADP_E_pos.rds")

# Control reference object (ADP_N EGFP+ neurons, from script 04)
ADP_N_pos <- readRDS("ADP_N_pos.rds")

cat("  ADP_E_pos cells:", ncol(ADP_E_pos), "\n")
cat("  ADP_N_pos cells:", ncol(ADP_N_pos), "\n")

# =============================================
# 4. Prepare reference object (ADP_N_pos)
# =============================================

cat("[", Sys.time(), "] Preparing reference object ...\n")

obj_ctrl <- ADP_N_pos
DefaultAssay(obj_ctrl) <- "SCT"

# Ensure PCA and UMAP are run on reference
if (!"pca" %in% names(obj_ctrl@reductions)) {
  cat("  Running PCA on reference ...\n")
  obj_ctrl <- RunPCA(obj_ctrl, verbose = FALSE)
}

if (!"umap" %in% names(obj_ctrl@reductions)) {
  cat("  Running UMAP on reference ...\n")
  obj_ctrl <- RunUMAP(
    obj_ctrl,
    reduction = "pca",
    reduction.name = "umap",
    dims = 1:30,
    spread = 0.3,
    min.dist = 0.1,
    n.neighbors = 30,
    return.model = TRUE,
    verbose = FALSE
  )
}

cat("  Reference object ready.\n")

# =============================================
# 5. Helper function: map query cluster to reference
# =============================================

map_cluster_to_reference <- function(query_obj, cluster_id, ref_obj, ref_cluster_col = "SCT_snn_res.0.4") {
  
  cat("[", Sys.time(), "] Mapping Cluster", cluster_id, "...\n")
  
  # Subset query cells
  cells_query <- colnames(query_obj)[query_obj$SCT_snn_res.1.2 == as.character(cluster_id)]
  obj_query <- subset(query_obj, cells = cells_query)
  DefaultAssay(obj_query) <- "SCT"
  
  cat("  Query cells:", ncol(obj_query), "\n")
  
  # Find transfer anchors
  anchors <- FindTransferAnchors(
    reference = ref_obj,
    query = obj_query,
    normalization.method = "SCT",
    dims = 1:30,
    reference.reduction = "pca"
  )
  
  # Map query to reference
  obj_query <- MapQuery(
    anchorset = anchors,
    query = obj_query,
    reference = ref_obj,
    refdata = ref_obj[[ref_cluster_col]],
    reduction.model = "umap"
  )
  
  return(obj_query)
}

# =============================================
# 6. Map Cluster 0 and Cluster 19 to reference
# =============================================

# Map Cluster 0
obj_query_0 <- map_cluster_to_reference(
  query_obj = ADP_E_pos,
  cluster_id = 0,
  ref_obj = obj_ctrl
)

# Map Cluster 19
obj_query_19 <- map_cluster_to_reference(
  query_obj = ADP_E_pos,
  cluster_id = 19,
  ref_obj = obj_ctrl
)

# =============================================
# 7. Export mapping scores (NEW)
# =============================================

cat("[", Sys.time(), "] Exporting mapping scores ...\n")

# Extract mapping results for Cluster 0
mapping_0 <- data.frame(
  cell_id = rownames(obj_query_0@meta.data),
  predicted.id = obj_query_0$predicted.id,
  predicted.id.score = obj_query_0$predicted.id.score,
  cluster = "0"
)
write.csv(mapping_0, "./results/tables/mapping_scores_cluster_0.csv", row.names = FALSE)

# Extract mapping results for Cluster 19
mapping_19 <- data.frame(
  cell_id = rownames(obj_query_19@meta.data),
  predicted.id = obj_query_19$predicted.id,
  predicted.id.score = obj_query_19$predicted.id.score,
  cluster = "19"
)
write.csv(mapping_19, "./results/tables/mapping_scores_cluster_19.csv", row.names = FALSE)

cat("  Mapping scores exported to ./results/tables/\n")

# =============================================
# 8. Prepare UMAP coordinates for visualization
# =============================================

cat("[", Sys.time(), "] Preparing UMAP visualization data ...\n")

# Reference UMAP coordinates
ref_umap <- as.data.frame(obj_ctrl@reductions$umap@cell.embeddings)
colnames(ref_umap) <- c("UMAP_1", "UMAP_2")
ref_umap$type <- "Reference"

# Cluster 0 projected coordinates
umap_0 <- as.data.frame(obj_query_0@reductions$ref.umap@cell.embeddings)
colnames(umap_0) <- c("UMAP_1", "UMAP_2")
umap_0$mapping_score <- obj_query_0$predicted.id.score
umap_0$cluster <- "Cluster 0"

# Cluster 19 projected coordinates
umap_19 <- as.data.frame(obj_query_19@reductions$ref.umap@cell.embeddings)
colnames(umap_19) <- c("UMAP_1", "UMAP_2")
umap_19$mapping_score <- obj_query_19$predicted.id.score
umap_19$cluster <- "Cluster 19"

# =============================================
# 9. Plot: Projected UMAP (dual color gradient)
# =============================================

cat("[", Sys.time(), "] Generating UMAP projection plot ...\n")

# Background: reference cells (gray)
p <- ggplot() +
  geom_point(
    data = ref_umap,
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey70",
    size = 0.3,
    alpha = 0.7
  )

# Cluster 0: blue gradient
p <- p +
  geom_point(
    data = umap_0,
    aes(x = UMAP_1, y = UMAP_2, color = mapping_score),
    size = 0.3,
    alpha = 0.85
  ) +
  scale_color_gradient(
    low = "#C8DEF9",
    high = "#4E79A7",
    name = "Cluster 0\nMapping Score"
  ) +
  new_scale_color()

# Cluster 19: red gradient
p <- p +
  geom_point(
    data = umap_19,
    aes(x = UMAP_1, y = UMAP_2, color = mapping_score),
    size = 0.3,
    alpha = 0.85
  ) +
  scale_color_gradient(
    low = "#FFC0E0",
    high = "#9F044D",
    name = "Cluster 19\nMapping Score"
  )

# Theme
p <- p +
  labs(
    x = "UMAP_1",
    y = "UMAP_2",
    title = "Reference Mapping: Cluster 0 & 19 Projected onto Control UMAP"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)
ggsave("./results/figures/UMAP_reference_projection.pdf", plot = p, width = 10, height = 8)

# =============================================
# 10. Pie chart: Mapping distribution of Cluster 0
# =============================================

cat("[", Sys.time(), "] Generating pie charts ...\n")

# Function to create pie chart
create_pie_chart <- function(query_obj, cluster_id, color_palette) {
  
  pred <- query_obj$predicted.id
  df <- as.data.frame(table(pred))
  colnames(df) <- c("Cluster", "Count")
  
  all_clusters <- as.character(0:12)
  df$Cluster <- factor(df$Cluster, levels = all_clusters)
  
  df <- df %>%
    mutate(Proportion = round(Count / sum(Count) * 100, 1))
  
  df <- df %>%
    mutate(label = ifelse(Count > 0, paste0(Cluster, " (", Proportion, "%)"), ""))
  
  p_pie <- ggplot(df, aes(x = "", y = Count, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y", start = 0) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3.5
    ) +
    scale_fill_manual(
      values = color_palette,
      drop = FALSE
    ) +
    labs(title = paste("Mapping Distribution of Cluster", cluster_id)) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p_pie)
}

# Define color palette for reference clusters (0-12)
ref_colors <- c(
  "0"  = "#377EB8", "1"  = "#83b5b5", "2"  = "#839958",
  "3"  = "#9887bc", "4"  = "#deb956", "5"  = "#bfc5d5",
  "6"  = "#c8b8d4", "7"  = "#7c7a7d", "8"  = "#c1d09d",
  "9"  = "#344B5C", "10" = "#C98F96", "11" = "#e08ea4",
  "12" = "#9dbdd2"
)

# Cluster 0 pie chart
p_pie_0 <- create_pie_chart(obj_query_0, 0, ref_colors)
print(p_pie_0)
ggsave("./results/figures/pie_cluster_0_mapping.pdf", plot = p_pie_0, width = 8, height = 6)

# Cluster 19 pie chart
p_pie_19 <- create_pie_chart(obj_query_19, 19, ref_colors)
print(p_pie_19)
ggsave("./results/figures/pie_cluster_19_mapping.pdf", plot = p_pie_19, width = 8, height = 6)

# =============================================
# 11. Combined pie chart (optional)
# =============================================

combined_pie <- p_pie_0 | p_pie_19
print(combined_pie)
ggsave("./results/figures/pie_combined_mapping.pdf", plot = combined_pie, width = 14, height = 6)

# =============================================
# 12. Save session info
# =============================================

session_info <- capture.output(sessionInfo())
writeLines(session_info, "./results/session_info_07_reference_mapping.txt")

# =============================================
# 13. Summary
# =============================================

cat("\n[", Sys.time(), "] ==========================================\n")
cat("Reference mapping completed!\n")
cat("  Cluster 0 mapped cells:", ncol(obj_query_0), "\n")
cat("  Cluster 19 mapped cells:", ncol(obj_query_19), "\n")
cat("  Mapping scores saved to: ./results/tables/\n")
cat("  Figures saved to: ./results/figures/\n")
cat("==========================================\n")
cat("[", Sys.time(), "] All done!\n")