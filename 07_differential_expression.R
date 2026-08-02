#!/usr/bin/env Rscript
# =============================================
# Differential expression analysis: ADP_H/ADP_C (cluster 0+19) vs ADP_N (cluster 0)
# Usage: Rscript scripts/R/07_differential_expression.R
# =============================================

# =============================================
# 1. Load packages
# =============================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
  library(patchwork)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(tidyverse)
})

# =============================================
# 2. Set paths and create output directories
# =============================================

dir.create("./results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/tables", recursive = TRUE, showWarnings = FALSE)

# =============================================
# 3. Load objects
# =============================================

cat("[", Sys.time(), "] Loading data ...\n")

ADP_E_pos <- readRDS("ADP_E_pos.rds")
ADP_N_pos <- readRDS("ADP_N_pos.rds")

cat("  ADP_E_pos cells:", ncol(ADP_E_pos), "\n")
cat("  ADP_N_pos cells:", ncol(ADP_N_pos), "\n")

# =============================================
# 4. Ensure ADP_E_pos has SCT assay
# =============================================

cat("[", Sys.time(), "] Ensuring SCT assay ...\n")

# Check if SCT assay exists, run if missing
if (!"SCT" %in% names(ADP_E_pos@assays)) {
  cat("  Running SCTransform (first time)...\n")
  DefaultAssay(ADP_E_pos) <- "RNA"
  ADP_E_pos <- SCTransform(
    ADP_E_pos,
    vars.to.regress = "percent.mt",
    variable.features.n = 3000,
    return.only.var.genes = FALSE,
    verbose = TRUE
  )
} else {
  cat("  SCT assay already exists. Re-running to ensure unified model...\n")
  DefaultAssay(ADP_E_pos) <- "RNA"
  ADP_E_pos <- SCTransform(
    ADP_E_pos,
    vars.to.regress = "percent.mt",
    variable.features.n = 3000,
    return.only.var.genes = FALSE,
    verbose = TRUE
  )
}

# View available resolution columns
cat("  Available resolution columns:\n")
print(grep("^SCT_snn_res\\.", colnames(ADP_E_pos@meta.data), value = TRUE))

# =============================================
# 5. Extract ADP_N cluster 0 cells (from ADP_N_pos)
# =============================================

cat("[", Sys.time(), "] Extracting ADP_N cluster 0 cells ...\n")

Idents(ADP_N_pos) <- "SCT_snn_res.0.4"

if (!"0" %in% levels(Idents(ADP_N_pos))) {
  stop("Cluster '0' not found in ADP_N_pos. Check clustering results.")
}

ADP_N_cluster0_cells <- WhichCells(ADP_N_pos, idents = "0")
cat("  ADP_N cluster 0 cells:", length(ADP_N_cluster0_cells), "\n")

# =============================================
# 6. Extract experimental cells (Cluster 0+19 from ADP_H and ADP_C)
# =============================================

cat("[", Sys.time(), "] Extracting experimental cells ...\n")

# Check if resolution 1.2 exists
if (!"SCT_snn_res.1.2" %in% colnames(ADP_E_pos@meta.data)) {
  stop("SCT_snn_res.1.2 not found. Please run clustering first.")
}

# Extract ADP_H cells in clusters 0 and 19
ADP_E_pos$temp_cluster <- ADP_E_pos$SCT_snn_res.1.2
ADP_E_pos$temp_sample <- ADP_E_pos$orig.ident

ADP_H_cells <- ADP_E_pos@meta.data %>%
  filter(temp_cluster %in% c("0", "19"), temp_sample == "ADP_H") %>%
  rownames()
cat("  ADP_H (cluster 0+19):", length(ADP_H_cells), "\n")

ADP_C_cells <- ADP_E_pos@meta.data %>%
  filter(temp_cluster %in% c("0", "19"), temp_sample == "ADP_C") %>%
  rownames()
cat("  ADP_C (cluster 0+19):", length(ADP_C_cells), "\n")

# Check if cells exist
if (length(ADP_H_cells) == 0) {
  stop("No ADP_H cells found in clusters 0 or 19.")
}
if (length(ADP_C_cells) == 0) {
  stop("No ADP_C cells found in clusters 0 or 19.")
}

# =============================================
# 7. Create subset object for differential expression
# =============================================

cat("[", Sys.time(), "] Creating subset for DEG analysis ...\n")

cells_to_keep <- c(ADP_H_cells, ADP_C_cells, ADP_N_cluster0_cells)
cat("  Total cells for DEG:", length(cells_to_keep), "\n")

sub_obj <- subset(ADP_E_pos, cells = cells_to_keep)

# Create group labels
sub_obj$group <- case_when(
  colnames(sub_obj) %in% ADP_H_cells ~ "ADP_H",
  colnames(sub_obj) %in% ADP_C_cells ~ "ADP_C",
  colnames(sub_obj) %in% ADP_N_cluster0_cells ~ "ADP_N",
  TRUE ~ "other"
)

cat("  Group cell counts:\n")
print(table(sub_obj$group))

# Remove "other" cells if any
if (sum(sub_obj$group == "other") > 0) {
  cat("  Removing", sum(sub_obj$group == "other"), "'other' cells...\n")
  sub_obj <- subset(sub_obj, subset = group != "other")
  sub_obj$group <- droplevels(sub_obj$group)
}

Idents(sub_obj) <- "group"

# =============================================
# 8. Differential expression: ADP_H vs ADP_N
# =============================================

cat("[", Sys.time(), "] Running DEG: ADP_H vs ADP_N ...\n")

deg_H <- FindMarkers(
  object = sub_obj,
  ident.1 = "ADP_H",
  ident.2 = "ADP_N",
  assay = "SCT",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1,
  min.cells.group = 3,
  only.pos = FALSE
)

deg_H <- deg_H %>%
  rownames_to_column("gene") %>%
  mutate(
    sig = case_when(
      avg_log2FC > 0.25 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "Not significant"
    ),
    log10_p = -log10(p_val_adj)
  )

cat("  ADP_H vs ADP_N:\n")
cat("    Up:", sum(deg_H$sig == "Up"), "\n")
cat("    Down:", sum(deg_H$sig == "Down"), "\n")
cat("    Not sig:", sum(deg_H$sig == "Not significant"), "\n")

write.csv(deg_H, "./results/tables/DEG_ADP_H1_19_vs_ADP_N0.csv", row.names = FALSE)

# =============================================
# 9. Differential expression: ADP_C vs ADP_N
# =============================================

cat("[", Sys.time(), "] Running DEG: ADP_C vs ADP_N ...\n")

deg_C <- FindMarkers(
  object = sub_obj,
  ident.1 = "ADP_C",
  ident.2 = "ADP_N",
  assay = "SCT",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1,
  min.cells.group = 3,
  only.pos = FALSE
)

deg_C <- deg_C %>%
  rownames_to_column("gene") %>%
  mutate(
    sig = case_when(
      avg_log2FC > 0.25 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "Not significant"
    ),
    log10_p = -log10(p_val_adj)
  )

cat("  ADP_C vs ADP_N:\n")
cat("    Up:", sum(deg_C$sig == "Up"), "\n")
cat("    Down:", sum(deg_C$sig == "Down"), "\n")
cat("    Not sig:", sum(deg_C$sig == "Not significant"), "\n")

write.csv(deg_C, "./results/tables/DEG_ADP_C1_19_vs_ADP_N0.csv", row.names = FALSE)

# =============================================
# 10. Helper function for volcano plots
# =============================================

get_top_genes <- function(df, n = 5) {
  if (nrow(df) == 0) return(character(0))
  df_sig <- df %>% filter(sig != "Not significant")
  if (nrow(df_sig) == 0) return(character(0))
  p_top <- df_sig %>% slice_min(p_val_adj, n = min(n, nrow(df_sig))) %>% pull(gene)
  fc_top <- df_sig %>% slice_max(abs(avg_log2FC), n = min(n, nrow(df_sig))) %>% pull(gene)
  union(p_top, fc_top)
}

# =============================================
# 11. Volcano plot: ADP_H vs ADP_N
# =============================================

cat("[", Sys.time(), "] Generating volcano plot: ADP_H vs ADP_N ...\n")

deg_up <- deg_H %>% filter(sig == "Up")
deg_down <- deg_H %>% filter(sig == "Down")
top_genes_H <- union(get_top_genes(deg_up, 5), get_top_genes(deg_down, 5))

deg_H_labeled <- deg_H %>%
  mutate(gene_label = ifelse(gene %in% top_genes_H, gene, ""))

colors_H <- c("Up" = "#9F044D", "Down" = "#db7ea1", "Not significant" = "gray80")

volcano_H <- ggplot(deg_H_labeled, aes(x = avg_log2FC, y = log10_p, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = colors_H, name = "Expression") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(
    aes(label = gene_label),
    size = 3,
    box.padding = 0.3,
    max.overlaps = 30,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(
    title = "Volcano Plot: ADP_H vs ADP_N",
    x = "Log2 Fold Change",
    y = "-Log10(Adj P-value)"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    plot.title = element_text(hjust = 0.5, color = "black")
  )

print(volcano_H)
ggsave("./results/figures/volcano_ADP_H_vs_ADP_N.pdf", plot = volcano_H, width = 8, height = 8)

# =============================================
# 12. Volcano plot: ADP_C vs ADP_N
# =============================================

cat("[", Sys.time(), "] Generating volcano plot: ADP_C vs ADP_N ...\n")

deg_up <- deg_C %>% filter(sig == "Up")
deg_down <- deg_C %>% filter(sig == "Down")
top_genes_C <- union(get_top_genes(deg_up, 5), get_top_genes(deg_down, 5))

deg_C_labeled <- deg_C %>%
  mutate(gene_label = ifelse(gene %in% top_genes_C, gene, ""))

colors_C <- c("Up" = "#4E79A7", "Down" = "#9DC7DD", "Not significant" = "gray80")

volcano_C <- ggplot(deg_C_labeled, aes(x = avg_log2FC, y = log10_p, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = colors_C, name = "Expression") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(
    aes(label = gene_label),
    size = 3,
    box.padding = 0.3,
    max.overlaps = 30,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(
    title = "Volcano Plot: ADP_C vs ADP_N",
    x = "Log2 Fold Change",
    y = "-Log10(Adj P-value)"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    plot.title = element_text(hjust = 0.5, color = "black")
  )

print(volcano_C)
ggsave("./results/figures/volcano_ADP_C_vs_ADP_N.pdf", plot = volcano_C, width = 8, height = 8)

# =============================================
# 13. GO enrichment analysis (down-regulated genes)
# =============================================

cat("[", Sys.time(), "] Running GO enrichment ...\n")

run_go_simplified <- function(genes, group_name) {
  if (length(genes) < 3) {
    message(paste("  [", group_name, "] Less than 3 genes, skipping.", sep = ""))
    return(NULL)
  }
  
  eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  if (nrow(eg) < 3) {
    message(paste("  [", group_name, "] <3 genes with ENTREZID mapping, skipping.", sep = ""))
    return(NULL)
  }
  
  go <- enrichGO(
    gene = eg$ENTREZID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if (is.null(go) || nrow(go) == 0) {
    message(paste("  [", group_name, "] No significant GO terms.", sep = ""))
    return(NULL)
  }
  
  go_simple <- simplify(go, cutoff = 0.7, by = "p.adjust")
  res_df <- as.data.frame(go_simple)
  res_df$Group <- group_name
  return(res_df)
}

# Extract down-regulated genes
down_genes_H <- deg_H %>% filter(sig == "Down") %>% pull(gene)
down_genes_C <- deg_C %>% filter(sig == "Down") %>% pull(gene)

cat("  Down-regulated genes (H):", length(down_genes_H), "\n")
cat("  Down-regulated genes (C):", length(down_genes_C), "\n")

res_H <- run_go_simplified(down_genes_H, "ADP_H_vs_ADP_N")
res_C <- run_go_simplified(down_genes_C, "ADP_C_vs_ADP_N")

all_go <- bind_rows(res_H, res_C)

if (!is.null(all_go) && nrow(all_go) > 0) {
  write.csv(all_go, "./results/tables/down_regulated_GO_BP_simplified_combined.csv", row.names = FALSE)
  cat("  GO results saved to: ./results/tables/\n")
} else {
  cat("  No significant GO terms found.\n")
}

# =============================================
# 14. GO bar plot (if results exist)
# =============================================

if (!is.null(all_go) && nrow(all_go) > 0) {
  
  cat("[", Sys.time(), "] Generating GO bar plot ...\n")
  
  go_data <- read.csv("./results/tables/down_regulated_GO_BP_simplified_combined.csv")
  
  # Define GO terms of interest
  common_GOs <- c("GO:0007416", "GO:0060078", "GO:0099565", "GO:0001764")
  H_special <- c("GO:0099003", "GO:0016358", "GO:0010959")
  C_special <- c("GO:0050803", "GO:0099072", "GO:0031346")
  all_GOs <- c(common_GOs, H_special, C_special)
  
  # Define display order
  order_df <- data.frame(
    ID = c("GO:0007416", "GO:0060078", "GO:0099565", "GO:0001764",
           "GO:0099003", "GO:0016358", "GO:0010959",
           "GO:0050803", "GO:0099072", "GO:0031346"),
    order = 1:10
  )
  
  # Prepare data
  plot_data <- go_data %>%
    filter(ID %in% all_GOs) %>%
    select(ID, Description, Group, GeneRatio, p.adjust) %>%
    mutate(
      Ratio = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      log10_padj = -log10(p.adjust),
      Group = ifelse(Group == "ADP_H_vs_ADP_N", "H", "C")
    ) %>%
    left_join(order_df, by = "ID")
  
  if (nrow(plot_data) > 0) {
    
    # Split by group
    plot_H <- plot_data %>% filter(Group == "H") %>% arrange(order)
    plot_C <- plot_data %>% filter(Group == "C") %>% arrange(order)
    
    # Set factor levels
    if (nrow(plot_H) > 0) {
      plot_H$Description <- factor(plot_H$Description, levels = unique(plot_H$Description))
    }
    if (nrow(plot_C) > 0) {
      plot_C$Description <- factor(plot_C$Description, levels = unique(plot_C$Description))
    }
    
    # Wrap long text
    if (nrow(plot_H) > 0) {
      plot_H$Description <- stringr::str_wrap(plot_H$Description, width = 30)
    }
    if (nrow(plot_C) > 0) {
      plot_C$Description <- stringr::str_wrap(plot_C$Description, width = 30)
    }
    
    # Plot H
    p_H <- NULL
    if (nrow(plot_H) > 0) {
      p_H <- ggplot(plot_H, aes(x = Ratio, y = Description, fill = log10_padj)) +
        geom_col() +
        scale_fill_gradient(low = "#db7ea1", high = "#9F044D",
                            name = "-log10(p.adjust)") +
        labs(title = "ADP_H vs ADP_N", x = "GeneRatio", y = NULL) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 10, hjust = 0),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        geom_text(aes(label = round(Ratio, 3)),
                  hjust = -0.1, size = 3.5, color = "black") +
        scale_y_discrete(limits = rev(levels(plot_H$Description)))
    }
    
    # Plot C
    p_C <- NULL
    if (nrow(plot_C) > 0) {
      p_C <- ggplot(plot_C, aes(x = Ratio, y = Description, fill = log10_padj)) +
        geom_col() +
        scale_fill_gradient(low = "#9DC7DD", high = "#4E79A7",
                            name = "-log10(p.adjust)") +
        labs(title = "ADP_C vs ADP_N", x = "GeneRatio", y = NULL) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 10, hjust = 0),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        geom_text(aes(label = round(Ratio, 3)),
                  hjust = -0.1, size = 3.5, color = "black") +
        scale_y_discrete(limits = rev(levels(plot_C$Description)))
    }
    
    # Combine plots
    if (!is.null(p_H) && !is.null(p_C)) {
      combined <- p_H + p_C + plot_layout(guides = "collect") &
        theme(legend.position = "right")
    } else if (!is.null(p_H)) {
      combined <- p_H
    } else if (!is.null(p_C)) {
      combined <- p_C
    } else {
      combined <- NULL
    }
    
    if (!is.null(combined)) {
      print(combined)
      ggsave("./results/figures/down_regulated_GO_BP.pdf",
             plot = combined, width = 12, height = 8)
    }
    
  } else {
    cat("  No matching GO terms found for bar plot.\n")
  }
}

# =============================================
# 15. Stacked dumbbell plot (optional - from script)
# =============================================

cat("[", Sys.time(), "] Generating stacked dumbbell plot ...\n")

meta <- ADP_E_pos@meta.data

if (any(duplicated(names(meta)))) {
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
cluster_levels <- sort(unique(prop_df$cluster_num))
prop_df$cluster <- factor(prop_df$cluster_num, levels = cluster_levels)

prop_df <- prop_df %>%
  arrange(cluster, sample) %>%
  group_by(cluster) %>%
  mutate(
    start = cumsum(lag(prop, default = 0)),
    end = start + prop
  ) %>%
  ungroup()

my_colors <- c("ADP_C" = "#4E79A7", "ADP_N" = "#bfc5d5", "ADP_H" = "#9F044D")

p_dumbbell <- ggplot(prop_df) +
  geom_segment(aes(x = start, xend = end, y = cluster, yend = cluster, color = sample),
               linewidth = 1) +
  geom_point(aes(x = start, y = cluster, color = sample), size = 1.2) +
  geom_point(aes(x = end, y = cluster, color = sample), size = 1.2) +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), labels = scales::percent) +
  labs(x = "Cumulative proportion of cells", y = "Cluster (Resolution 1.2)", color = "Sample") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent")
  )

print(p_dumbbell)
ggsave("./results/figures/stacked_dumbbell_ADP_E_pos.pdf", plot = p_dumbbell, width = 8, height = 6)

# =============================================
# 16. Save session info
# =============================================

session_info <- capture.output(sessionInfo())
writeLines(session_info, "./results/session_info_08_differential_expression.txt")

# =============================================
# 17. Summary
# =============================================

cat("\n[", Sys.time(), "] ==========================================\n")
cat("Differential expression analysis completed!\n")
cat("  ADP_H vs ADP_N:\n")
cat("    Up:", sum(deg_H$sig == "Up"), "\n")
cat("    Down:", sum(deg_H$sig == "Down"), "\n")
cat("  ADP_C vs ADP_N:\n")
cat("    Up:", sum(deg_C$sig == "Up"), "\n")
cat("    Down:", sum(deg_C$sig == "Down"), "\n")
cat("  Results saved to: ./results/tables/\n")
cat("  Figures saved to: ./results/figures/\n")
cat("==========================================\n")
cat("[", Sys.time(), "] All done!\n")