#!/usr/bin/env Rscript
# =============================================
# EGFP+ neuron sub-clustering analysis (after removing cluster 3)
# Usage: Rscript scripts/R/05_combined_EGFP_neurons_subclustering.R
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
  library(tidydr)
  library(ggrepel)
  library(future)
  library(gtools)
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

mycolor22 <- c("#4E79A7","#83b5b5","#839958","#9887bc","#FDDF91","#bfc5d5",
               "#c8b8d4","#7c7a7d","#c1d09d","#EFE9D3","#C98F96",
               "#e08ea4","#9dbdd2","#779ebd","#deb956","#bdbb55",
               "#d0c9b0","#b696b6","#C8DEF9","#80c1c4","#5d3c11","#3D9F3C","#344B5C","#FFC0E0")

# =============================================
# 3. Load data and extract neurons
# =============================================

cat("[", Sys.time(), "] Loading data ...\n")

ADP_combined <- readRDS("ADP_harmony_combined_with_celltype.rds")
ADP_neurons <- subset(ADP_combined, subset = cell_type %in% c("Neuron"))

cat("  Total neurons:", ncol(ADP_neurons), "\n")

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

ADP_neurons <- add_egfp_metadata(ADP_neurons, threshold)
cat("  EGFP+ cells:", sum(ADP_neurons$EGFP_positive == "Positive"), "\n")
cat("  EGFP- cells:", sum(ADP_neurons$EGFP_positive == "Negative"), "\n")

# Extract EGFP+ neurons only
ADP_E_pos <- subset(ADP_neurons, subset = EGFP_positive == "Positive")
cat("  EGFP+ neurons by sample:\n")
print(table(ADP_E_pos$orig.ident))

# =============================================
# 5. SCTransform on EGFP+ neurons (excluding EGFP) 
# =============================================

cat("[", Sys.time(), "] Running SCTransform on EGFP+ neurons ...\n")

DefaultAssay(ADP_E_pos) <- "RNA"

ADP_E_pos <- SCTransform(
  ADP_E_pos,
  vars.to.regress = "percent.mt",
  variable.features.n = 3000,
  return.only.var.genes = FALSE,
  verbose = TRUE
)

# Remove EGFP from variable features
var_genes <- VariableFeatures(ADP_E_pos)
var_genes_no_egfp <- setdiff(var_genes, "EGFP")
VariableFeatures(ADP_E_pos) <- var_genes_no_egfp

cat("  Variable features (excluding EGFP):", length(VariableFeatures(ADP_E_pos)), "\n")

# =============================================
# 6. Dimensionality reduction and clustering
# =============================================

cat("[", Sys.time(), "] Running PCA and clustering ...\n")

DefaultAssay(ADP_E_pos) <- "SCT"

ADP_E_pos <- RunPCA(ADP_E_pos, features = VariableFeatures(ADP_E_pos), verbose = FALSE)
ADP_E_pos <- FindNeighbors(ADP_E_pos, dims = 1:30, verbose = FALSE)

# Test multiple resolutions (0 to 1.2)
resolutions <- seq(0, 1.2, by = 0.1)
for (res in resolutions) {
  ADP_E_pos <- FindClusters(ADP_E_pos, resolution = res, verbose = FALSE)
}

# =============================================
# 7. Clustree
# =============================================

cat("[", Sys.time(), "] Generating clustree ...\n")
p_clustree <- clustree(ADP_E_pos, prefix = "SCT_snn_res.")
print(p_clustree)

# =============================================
# 8. Cluster tree (BuildClusterTree)
# =============================================

Idents(ADP_E_pos) <- "SCT_snn_res.1.2"

# Build cluster tree
ADP_E_pos <- BuildClusterTree(ADP_E_pos, dims = 1:30, reorder = FALSE)

# Save cluster tree using base R graphics
pdf("./results/figures/cluster_tree.pdf", width = 8, height = 6)
PlotClusterTree(ADP_E_pos)
dev.off()

# Also display in R console (optional)
PlotClusterTree(ADP_E_pos)

table(ADP_E_pos$SCT_snn_res.1.2)

# =============================================
# 9. Save initial object
# =============================================

saveRDS(ADP_E_pos, file = "ADP_E_pos.rds")

# =============================================
# 10. Load and prepare for visualization
# =============================================

ADP_E_pos <- readRDS("ADP_E_pos.rds")
head(ADP_E_pos)

# =============================================
# 11. Stacked bar plot: sample proportions per cluster
# =============================================

cat("[", Sys.time(), "] Generating stacked bar plot ...\n")

sample_col <- "orig.ident"
res_cols <- grep("^SCT_snn_res\\.", colnames(ADP_E_pos@meta.data), value = TRUE)
fetch_vars <- c(sample_col, res_cols)

meta_df <- FetchData(ADP_E_pos, vars = fetch_vars, cells = colnames(ADP_E_pos))
colnames(meta_df)[colnames(meta_df) == sample_col] <- "sample"

meta_long <- meta_df %>%
  pivot_longer(cols = all_of(res_cols),
               names_to = "resolution_col",
               values_to = "cluster") %>%
  mutate(resolution = as.numeric(gsub("SCT_snn_res\\.", "", resolution_col))) %>%
  filter(!is.na(sample) & !is.na(cluster))

prop_list <- split(meta_long, meta_long$resolution)

prop_final <- lapply(names(prop_list), function(res_str) {
  res_val <- as.numeric(res_str)
  df <- prop_list[[res_str]]
  tbl <- table(as.character(df$cluster), as.character(df$sample))
  prop_tbl <- prop.table(tbl, margin = 1)
  prop_df <- as.data.frame(prop_tbl, stringsAsFactors = FALSE)
  colnames(prop_df) <- c("cluster", "sample", "prop")
  prop_df$resolution <- res_val
  return(prop_df)
}) %>% bind_rows()

prop_final <- prop_final[is.finite(prop_final$prop), ]
all_clusters <- unique(prop_final$cluster)
prop_final$cluster <- factor(prop_final$cluster,
                             levels = gtools::mixedsort(all_clusters))

sample_colors <- c("ADP_N" = "#E1AF00", "ADP_C" = "#3B9AB2", "ADP_H" = "#F21A00")

p_prop <- ggplot(prop_final, aes(x = cluster, y = prop, fill = sample)) +
  geom_col(position = "fill", width = 0.8) +
  facet_wrap(~ resolution, ncol = 4) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = sample_colors) +
  labs(x = "Cluster", y = "Proportion", fill = "Sample") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom"
  )

print(p_prop)

# =============================================
# 12. Fos expression analysis
# =============================================

cat("[", Sys.time(), "] Analyzing Fos expression ...\n")


if ("stim" %in% colnames(ADP_E_pos@meta.data)) {
  sample_col <- "stim"
} else if ("orig.ident" %in% colnames(ADP_E_pos@meta.data)) {
  sample_col <- "orig.ident"
} else {
  stop("请手动指定样本列名")
}

DefaultAssay(ADP_E_pos) <- "SCT"
fos_data <- FetchData(ADP_E_pos, vars = c("Fos", sample_col, "SCT_snn_res.1.2"))
colnames(fos_data) <- c("Fos", "sample", "cluster")

fos_data <- fos_data[!is.na(fos_data$sample), ]

fos_mean_1.2 <- fos_data %>%
  group_by(cluster, sample) %>%
  summarise(
    mean_fos = mean(Fos, na.rm = TRUE),
    median_fos = median(Fos, na.rm = TRUE),
    sd_fos = sd(Fos, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(cluster, sample)

print(fos_mean_1.2, n = 100)

sample_colors <- c("ADP_N" = "#E1AF00", "ADP_C" = "#3B9AB2", "ADP_H" = "#F21A00")

p_fos_bar <- ggplot(fos_mean_1.2, aes(x = factor(cluster), y = mean_fos, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = sample_colors) +
  labs(x = "Cluster (res 1)", y = "Mean Fos Expression (SCT)", fill = "Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_fos_bar)


# =============================================
# 13. UMAP before removing cluster 3
# =============================================

cat("[", Sys.time(), "] Generating UMAP (before removing cluster 3) ...\n")

DimPlot(ADP_E_pos,
        group.by = "SCT_snn_res.1.2",
        reduction = "umap",
        shuffle = TRUE,
        label = TRUE)

# =============================================
# 14. Remove cluster 3 (mitochondrial gene-enriched low-quality cluster)
# =============================================

cat("[", Sys.time(), "] Removing cluster 3 ...\n")
cat("  Cells before removal:", ncol(ADP_E_pos), "\n")

ADP_E_pos <- subset(ADP_E_pos, subset = SCT_snn_res.1.2 != "3")

cat("  Cells after removal:", ncol(ADP_E_pos), "\n")

# =============================================
# 15. Re-run UMAP after removing cluster 3
# =============================================

cat("[", Sys.time(), "] Re-running UMAP after removal ...\n")

ADP_E_pos <- FindClusters(ADP_E_pos, resolution = 1.2)
ADP_E_pos <- RunUMAP(
  ADP_E_pos,
  reduction = "pca",
  reduction.name = "umap",
  dims = 1:30,
  spread = 0.4,
  min.dist = 0.2,
  n.neighbors = 30,
  verbose = FALSE
)

# =============================================
# 16. ggplot UMAP (your original code, keep as is)
# =============================================
umap_data <- Embeddings(ADP_E_pos, "umap")
meta_data <- ADP_E_pos@meta.data
df_plot <- data.frame(umap_data, meta_data)

centroids <- df_plot %>%
  group_by(seurat_clusters) %>%
  summarise(umap_1 = median(UMAP_1), umap_2 = median(UMAP_2))

p_umap <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = seurat_clusters), size = 0.4, alpha = 0.8, shape = 16) +
  scale_color_manual(values = mycolor22) +
  geom_text(data = centroids, aes(label = seurat_clusters),
            size = 3, color = "black") +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

print(p_umap)

# =============================================
# 17. Subset and re-plot UMAP (your original code)
# =============================================

ADP_E_pos_sub <- subset(ADP_E_pos, subset = SCT_snn_res.1.2 != "3")

set.seed(42)
ADP_E_pos_sub <- RunUMAP(
  ADP_E_pos_sub,
  reduction = "pca",
  dims = 1:30,
  spread = 0.4,
  min.dist = 0.2,
  n.neighbors = 30,
  reduction.name = "umap",
  return.model = TRUE
)

umap_data <- Embeddings(ADP_E_pos_sub, "umap")
meta_data <- ADP_E_pos_sub@meta.data
df_plot <- data.frame(umap_data, meta_data)

df_plot$clusters <- as.factor(df_plot$SCT_snn_res.1.2)
df_plot <- df_plot[df_plot$clusters != "3", ]
df_plot$clusters <- droplevels(df_plot$clusters)

centroids <- df_plot %>%
  group_by(clusters) %>%
  summarise(umap_1 = median(UMAP_1), umap_2 = median(UMAP_2))

mycolor20_ori <- c("#4E79A7","#83b5b5","#839958","#9887bc","#FDDF91","#bfc5d5",
                   "#c8b8d4","#7c7a7d","#c1d09d","#EFE9D3","#C98F96",
                   "#e08ea4","#9dbdd2","#779ebd","#deb956","#bdbb55",
                   "#d0c9b0","#b696b6","#C8DEF9","#80c1c4","#5d3c11",
                   "#7c7a7d","#344B5C","#FFC0E0","#C8DEF9")

my_colors <- mycolor20_ori[-4]

existing_clusters <- levels(df_plot$clusters)
names(my_colors) <- existing_clusters

p_umap_sub <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = clusters), size = 0.6, alpha = 0.8, shape = 16) +
  scale_color_manual(values = my_colors) +
  geom_text(data = centroids, aes(label = clusters), size = 3, color = "black") +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

print(p_umap_sub)


# =============================================
# 18. Stacked dumbbell plot (proportion per cluster, after removal)
# =============================================

cat("[", Sys.time(), "] Generating stacked dumbbell plot ...\n")

meta <- ADP_E_pos@meta.data

if(any(duplicated(names(meta)))) {
  warning("重复列名：", paste(names(meta)[duplicated(names(meta))], collapse = ", "))
  meta <- meta[, !duplicated(names(meta))]
}

cluster_col <- "SCT_snn_res.1.2"
sample_col <- "orig.ident"

prop_df <- meta %>%
  group_by(cluster = .data[[cluster_col]], sample = .data[[sample_col]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

sample_order <- c("ADP_C", "ADP_N", "ADP_H")
prop_df$sample <- factor(prop_df$sample, levels = sample_order)

prop_df$cluster_num <- as.numeric(as.character(prop_df$cluster))
cluster_levels <- sort(unique(prop_df$cluster_num), decreasing = FALSE)
prop_df$cluster <- factor(prop_df$cluster_num, levels = cluster_levels)

prop_df <- prop_df %>%
  arrange(cluster, sample) %>%
  group_by(cluster) %>%
  mutate(
    start = cumsum(lag(prop, default = 0)),
    end = start + prop
  ) %>%
  ungroup()

my_colors_dumbbell <- c("ADP_C" = "#4E79A7",
                        "ADP_N" = "#bfc5d5",
                        "ADP_H" = "#9F044D")

p_dumbbell <- ggplot(prop_df) +
  geom_segment(aes(x = start, xend = end,
                   y = cluster, yend = cluster, color = sample),
               linewidth = 1) +
  geom_point(aes(x = start, y = cluster, color = sample), size = 1.2) +
  geom_point(aes(x = end, y = cluster, color = sample), size = 1.2) +
  scale_color_manual(values = my_colors_dumbbell) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     labels = scales::percent) +
  labs(x = "Cumulative proportion of cells",
       y = "Cluster (Resolution 1.2)",
       color = "Sample") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

print(p_dumbbell)

# =============================================
# 19. Fos dot plot
# =============================================

cat("[", Sys.time(), "] Generating Fos dot plot ...\n")

ieg_genes <- c("Fos")
sct_genes <- rownames(GetAssayData(ADP_E_pos, assay = "SCT", slot = "data"))
ieg_genes_valid <- intersect(ieg_genes, sct_genes)

print(paste("Available IEGs in SCT:", paste(ieg_genes_valid, collapse = ", ")))

expr_matrix <- GetAssayData(ADP_E_pos, assay = "SCT", slot = "data")[ieg_genes_valid, , drop = FALSE]

expr_long <- as.data.frame(t(as.matrix(expr_matrix))) %>%
  tibble::rownames_to_column("cell_id") %>%
  pivot_longer(cols = -cell_id, names_to = "gene", values_to = "expression")

meta <- ADP_E_pos@meta.data[, c("SCT_snn_res.1.2", "orig.ident")]
meta$cell_id <- rownames(meta)
expr_long <- left_join(expr_long, meta, by = "cell_id")

expr_long$cluster <- factor(expr_long$SCT_snn_res.1.2,
                            levels = sort(unique(as.numeric(as.character(expr_long$SCT_snn_res.1.2)))))

ieg_stats <- expr_long %>%
  group_by(gene, cluster, orig.ident) %>%
  summarise(
    mean_expr = mean(expression),
    pct_expressed = sum(expression > 0) / n() * 100,
    .groups = "drop"
  )

write.csv(ieg_stats, "./results/tables/IEG_expression_stats_by_cluster_sample.csv", row.names = FALSE)

# Fos dot plot
fos_mean <- expr_long %>%
  filter(gene == "Fos") %>%
  group_by(cluster, orig.ident) %>%
  summarise(mean_expr = mean(expression), .groups = "drop") %>%
  mutate(cluster_num = as.numeric(as.character(cluster)))

fos_mean <- fos_mean %>%
  mutate(cluster_factor = factor(cluster_num,
                                 levels = sort(unique(cluster_num), decreasing = TRUE)))

sample_colors_fos <- c("ADP_C" = "#4E79A7",
                       "ADP_N" = "#bfc5d5",
                       "ADP_H" = "#9F044D")

p_fos_dot <- ggplot(fos_mean, aes(x = orig.ident, y = cluster_factor)) +
  geom_point(aes(size = mean_expr, color = orig.ident, alpha = mean_expr)) +
  scale_size_continuous(name = "Mean Fos Expression", range = c(1, 8)) +
  scale_alpha_continuous(name = "Mean Fos Expression", range = c(0.4, 1)) +
  scale_color_manual(values = sample_colors_fos, name = "Sample") +
  labs(x = "Sample", y = "Cluster",
       title = "Fos Expression Across Clusters") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )

print(p_fos_dot)

# =============================================
# 20. Combined plot (dumbbell + Fos dot)
# =============================================

cat("[", Sys.time(), "] Generating combined plot ...\n")

offset <- 0.005
combined_plot <- p_dumbbell | p_fos_dot
print(combined_plot)
# =============================================
# 21. Find marker genes for all clusters (excluding EGFP)
# =============================================

cat("[", Sys.time(), "] Finding marker genes ...\n")

DefaultAssay(ADP_E_pos) <- "RNA"
ADP_E_pos <- SCTransform(
  ADP_E_pos,
  vars.to.regress = "percent.mt",
  variable.features.n = 3000,
  return.only.var.genes = FALSE,
  verbose = TRUE
)

var_genes <- VariableFeatures(ADP_E_pos)
var_genes_no_egfp <- setdiff(var_genes, "EGFP")
VariableFeatures(ADP_E_pos) <- var_genes_no_egfp

ADP_E_pos_all_markers <- FindAllMarkers(
  object = ADP_E_pos,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0,
  test.use = "wilcox",
  verbose = TRUE
)

cat("  Total marker genes found:", nrow(ADP_E_pos_all_markers), "\n")

write.csv(ADP_E_pos_all_markers,
          file = "./results/tables/ADP_E_pos_all_markers.csv",
          row.names = TRUE)

# =============================================
# 22. QC metrics per cluster
# =============================================

cat("[", Sys.time(), "] QC metrics per cluster ...\n")

VlnPlot(ADP_E_pos,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "SCT_snn_res.1.2",
        pt.size = 0,
        ncol = 3)

VlnPlot(ADP_E_pos,
        features = "percent.mt",
        group.by = "SCT_snn_res.1.2",
        pt.size = 0) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  labs(title = "Mitochondrial percentage per cluster", y = "% mt reads")


meta <- ADP_E_pos@meta.data
cluster_mt_stats <- meta %>%
  group_by(cluster = SCT_snn_res.1.2) %>%
  summarise(
    n_cells = n(),
    median_mt = median(percent.mt),
    pct_high_mt = mean(percent.mt > 10) * 100
  ) %>%
  arrange(desc(pct_high_mt))

print(cluster_mt_stats)

# =============================================
# 23. Dot plots for marker genes
# =============================================

cat("[", Sys.time(), "] Generating dot plots ...\n")

genes_to_plot <- c(
  "Ntng1","Ptprk","Pou2f2","Lars2",
  "Sox6", "Stxbp6", "Lhx6", "Gad1",
  "Reln", "Slc17a6", "Nos1", "Trp73",
  "Rarb","Meis2",  "Wfs1", "Chrm1",
  "Foxp2","Zeb2", "Oprm1",  "Penk",
  "Lhx8", "Pde11a", "Lhfp", "Kcnh8",
  "Prdm16","Gfra1", "Nfia", "Cacng5",
  "Tshz1", "Chst9", "Lypd1","Eya2",
  "Prlr", "Greb1", "Esr2", "Dlx1",
  "Chat", "Slc5a7", "Prima1", "Ngfr"
)

custom_cluster_order <- unique(c(
  "0", "19","7","9","6", "2","1", "11","5","8","14","4",
  "16", "10","13","21", "15", "12", "20", "22", "18", "17", "23"
))

ADP_E_pos_sub <- subset(ADP_E_pos, subset = SCT_snn_res.1.2 %in% custom_cluster_order)
ADP_E_pos_sub$SCT_snn_res.1.2 <- droplevels(ADP_E_pos_sub$SCT_snn_res.1.2)
ADP_E_pos_sub$SCT_snn_res.1.2 <- factor(
  ADP_E_pos_sub$SCT_snn_res.1.2,
  levels = custom_cluster_order
)

dot_plot <- DotPlot(
  ADP_E_pos_sub,
  features = genes_to_plot,
  group.by = "SCT_snn_res.1.2",
  dot.scale = 6,
  scale = TRUE
) +
  scale_color_gradient2(low = "white", high = "#367DB0") +
  labs(x = NULL, y = "Cluster") +
  theme_bw(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    legend.position = "right"
  )

print(dot_plot)

# =============================================
# 24. Focused dot plot with selected marker genes
# =============================================

genes_to_plot_focused <- c(
  "Lars2","Cmss1","Camk1d","Gphn","Cdk8","Lrrc45",
  "Per1","Nop53","Nr4a3","Fos","Jun","Gadd45g",
  "Prdm16","Nfia","Nr4a2",
  "Nr2e1","Crh","Nts",
  "Esr1","Greb1","Pappa","Lamp5",
  "Slc17a8","Nr2f1",
  "Trpm3","Sim1","Slc18a2"
)
genes_to_plot_focused <- unique(genes_to_plot_focused)

dot_plot_focused <- DotPlot(
  ADP_E_pos_sub,
  features = genes_to_plot_focused,
  group.by = "SCT_snn_res.1.2",
  dot.scale = 4,
  scale = TRUE
) +
  scale_color_gradient2(low = "white", high = "#367DB0") +
  labs(x = NULL, y = "Cluster") +
  theme_bw(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    legend.position = "right"
  )

print(dot_plot_focused)

# =============================================
# 25. Save final object and session info
# =============================================

cat("[", Sys.time(), "] Saving final objects ...\n")
saveRDS(ADP_E_pos, file = "ADP_E_pos_clean.rds")
saveRDS(ADP_E_pos_sub, file = "ADP_E_pos_sub.rds")

session_info <- capture.output(sessionInfo())
writeLines(session_info, "./results/session_info_06_EGFP_subclustering.txt")

# =============================================
# 26. Summary
# =============================================

cat("\n[", Sys.time(), "] ==========================================\n")
cat("EGFP+ neuron sub-clustering completed!\n")
cat("  Final EGFP+ neurons (after removing cluster 3):", ncol(ADP_E_pos), "\n")
cat("  Number of clusters (res 1.2):", length(unique(ADP_E_pos$seurat_clusters)), "\n")
cat("  Marker genes found:", nrow(ADP_E_pos_all_markers), "\n")
cat("  Figures: displayed in R console (not auto-saved)\n")
cat("  Tables saved to: ./results/tables/\n")
cat("==========================================\n")
cat("[", Sys.time(), "] All done!\n")