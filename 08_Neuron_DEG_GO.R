#!/usr/bin/env Rscript
# =============================================
# Neuronal differential expression and GO enrichment analysis
# ADP_H/ADP_C vs ADP_N in EGFP+ neurons (all samples combined)
# Usage: Rscript scripts/R/08_Neuron_DEG_GO.R
# =============================================

# =============================================
# 1. Load packages
# =============================================

rm(list=ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(tidyverse)
  library(patchwork)
  library(future)
})

# =============================================
# 2. Set paths and create output directories
# =============================================

dir.create("./results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/tables", recursive = TRUE, showWarnings = FALSE)

# Set memory limit
options(future.globals.maxSize = 5 * 1024^3)  # 5 GB
plan("sequential")

# =============================================
# 3. Load data
# =============================================

cat("[", Sys.time(), "] Loading data ...\n")

ADP_combined <- readRDS("ADP_harmony_combined_with_celltype.rds")
cat("  Total cells:", ncol(ADP_combined), "\n")

# =============================================
# 4. Add EGFP metadata and extract EGFP+ neurons
# =============================================

cat("[", Sys.time(), "] Extracting EGFP+ neurons ...\n")

egfp_gene <- "EGFP"
egfp_counts <- GetAssayData(ADP_combined, assay = "RNA", layer = "counts")[egfp_gene, ]
ADP_combined$EGFP_counts <- egfp_counts

threshold <- 1
ADP_combined$EGFP_positive <- ifelse(ADP_combined$EGFP_counts >= threshold,
                                     "Positive", "Negative")

cat("  EGFP+ cells:", sum(ADP_combined$EGFP_positive == "Positive"), "\n")

ADP_combined_pos <- subset(ADP_combined, subset = EGFP_positive == "Positive")

# Extract neurons only
ADP_Neuron_pos <- subset(ADP_combined_pos, subset = cell_type %in% c("Neuron"))
cat("  EGFP+ neurons:", ncol(ADP_Neuron_pos), "\n")
cat("  By sample:\n")
print(table(ADP_Neuron_pos$orig.ident))

# =============================================
# 5. SCTransform on EGFP+ neurons
# =============================================

cat("[", Sys.time(), "] Running SCTransform ...\n")

DefaultAssay(ADP_Neuron_pos) <- "RNA"

ADP_Neuron_pos <- SCTransform(
  ADP_Neuron_pos,
  vars.to.regress = "percent.mt",
  variable.features.n = 3000,
  return.only.var.genes = FALSE,
  verbose = TRUE
)

cat("  Variable features:", length(VariableFeatures(ADP_Neuron_pos)), "\n")

# =============================================
# 6. Differential expression: ADP_H vs ADP_N
# =============================================

cat("[", Sys.time(), "] Running DEG: ADP_H vs ADP_N ...\n")

deg_HN <- FindMarkers(
  object = ADP_Neuron_pos,
  ident.1 = "ADP_H",
  ident.2 = "ADP_N",
  group.by = "orig.ident",
  assay = "SCT",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  min.cells.group = 3,
  only.pos = FALSE
)

deg_HN <- deg_HN %>%
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
cat("    Up:", sum(deg_HN$sig == "Up"), "\n")
cat("    Down:", sum(deg_HN$sig == "Down"), "\n")
cat("    Not significant:", sum(deg_HN$sig == "Not significant"), "\n")

# ========== EXPORT TABLE 1 ==========
write.csv(deg_HN, "./results/tables/Neuron_ADP_H_vs_ADP_N_DEGs.csv", row.names = FALSE)

# =============================================
# 7. Differential expression: ADP_C vs ADP_N
# =============================================

cat("[", Sys.time(), "] Running DEG: ADP_C vs ADP_N ...\n")

deg_CN <- FindMarkers(
  object = ADP_Neuron_pos,
  ident.1 = "ADP_C",
  ident.2 = "ADP_N",
  group.by = "orig.ident",
  assay = "SCT",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  min.cells.group = 3,
  only.pos = FALSE
)

deg_CN <- deg_CN %>%
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
cat("    Up:", sum(deg_CN$sig == "Up"), "\n")
cat("    Down:", sum(deg_CN$sig == "Down"), "\n")
cat("    Not significant:", sum(deg_CN$sig == "Not significant"), "\n")

# ========== EXPORT TABLE 2 ==========
write.csv(deg_CN, "./results/tables/Neuron_ADP_C_vs_ADP_N_DEGs.csv", row.names = FALSE)

# =============================================
# 8. Helper function: get top genes for volcano labeling
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
# 9. Volcano plot: ADP_H vs ADP_N
# =============================================

cat("[", Sys.time(), "] Generating volcano plot: ADP_H vs ADP_N ...\n")

deg_up <- deg_HN %>% filter(sig == "Up")
deg_down <- deg_HN %>% filter(sig == "Down")
top_genes_HN <- union(get_top_genes(deg_up, 5), get_top_genes(deg_down, 5))

deg_HN_labeled <- deg_HN %>%
  mutate(gene_label = ifelse(gene %in% top_genes_HN, gene, ""))

colors_H <- c("Up" = "#9F044D", "Down" = "#db7ea1", "Not significant" = "gray80")

volcano_HN <- ggplot(deg_HN_labeled, aes(x = avg_log2FC, y = log10_p, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = colors_H, name = "Expression") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(aes(label = gene_label), size = 3, box.padding = 0.3,
                  max.overlaps = 30, show.legend = FALSE, na.rm = TRUE) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, max(deg_HN$log10_p, na.rm = TRUE) + 5)) +
  labs(title = "Volcano Plot: ADP_H vs ADP_N (Neurons)",
       x = "Log2 Fold Change", y = "-Log10(Adj P-value)") +
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

print(volcano_HN)

# ========== EXPORT FIGURE 1 ==========
ggsave("./results/figures/Neuron_Volcano_HN.pdf",
       plot = volcano_HN, width = 8, height = 7)

# =============================================
# 10. Volcano plot: ADP_C vs ADP_N
# =============================================

cat("[", Sys.time(), "] Generating volcano plot: ADP_C vs ADP_N ...\n")

deg_up <- deg_CN %>% filter(sig == "Up")
deg_down <- deg_CN %>% filter(sig == "Down")
top_genes_CN <- union(get_top_genes(deg_up, 5), get_top_genes(deg_down, 5))

deg_CN_labeled <- deg_CN %>%
  mutate(gene_label = ifelse(gene %in% top_genes_CN, gene, ""))

colors_C <- c("Up" = "#4E79A7", "Down" = "#9DC7DD", "Not significant" = "gray80")

volcano_CN <- ggplot(deg_CN_labeled, aes(x = avg_log2FC, y = log10_p, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = colors_C, name = "Expression") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(aes(label = gene_label), size = 3, box.padding = 0.3,
                  max.overlaps = 30, show.legend = FALSE, na.rm = TRUE) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, max(deg_CN$log10_p, na.rm = TRUE) + 5)) +
  labs(title = "Volcano Plot: ADP_C vs ADP_N (Neurons)",
       x = "Log2 Fold Change", y = "-Log10(Adj P-value)") +
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

print(volcano_CN)

# ========== EXPORT FIGURE 2 ==========
ggsave("./results/figures/Neuron_Volcano_CN.pdf",
       plot = volcano_CN, width = 8, height = 7)

# =============================================
# 11. GO enrichment function (for down-regulated genes)
# =============================================

run_go_enrichment <- function(deg_df, group_name, exclude_pattern = "^(mt-|MT-|Gm|Rps|Rpl|Rp[ls])") {
  
  down_genes <- deg_df %>%
    filter(sig == "Down") %>%
    filter(!grepl(exclude_pattern, gene, ignore.case = TRUE)) %>%
    pull(gene)
  
  if (length(down_genes) < 3) {
    message("  [", group_name, "] <3 down-regulated genes after filtering, skipping.")
    return(NULL)
  }
  
  cat("  [", group_name, "] Down-regulated genes for GO:", length(down_genes), "\n")
  
  eg <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  if (nrow(eg) < 3) {
    message("  [", group_name, "] <3 genes with ENTREZID mapping, skipping.")
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
    message("  [", group_name, "] No significant GO terms.")
    return(NULL)
  }
  
  # Simplify GO terms
  go_simple <- clusterProfiler::simplify(go, cutoff = 0.7, by = "p.adjust")
  
  return(go_simple)
}

# =============================================
# 12. Run GO enrichment for both comparisons
# =============================================

cat("[", Sys.time(), "] Running GO enrichment ...\n")

go_HN_simple <- run_go_enrichment(deg_HN, "ADP_H_vs_ADP_N")
go_CN_simple <- run_go_enrichment(deg_CN, "ADP_C_vs_ADP_N")

# ========== EXPORT TABLES 3 & 4 ==========
if (!is.null(go_HN_simple)) {
  write.csv(as.data.frame(go_HN_simple),
            "./results/tables/Neuron_HN_downregulated_GO_simplified.csv",
            row.names = FALSE)
  cat("  GO results saved: Neuron_HN_downregulated_GO_simplified.csv\n")
}

if (!is.null(go_CN_simple)) {
  write.csv(as.data.frame(go_CN_simple),
            "./results/tables/Neuron_CN_downregulated_GO_simplified.csv",
            row.names = FALSE)
  cat("  GO results saved: Neuron_CN_downregulated_GO_simplified.csv\n")
}

# =============================================
# 13. GO dot plots
# =============================================

cat("[", Sys.time(), "] Generating GO dot plots ...\n")

if (!is.null(go_HN_simple)) {
  p_dot_HN <- dotplot(go_HN_simple, showCategory = 15,
                      title = "Neuron ADP_H vs ADP_N (BP, simplified)")
  print(p_dot_HN)
  # ========== EXPORT FIGURE 3 ==========
  ggsave("./results/figures/Neuron_GO_dotplot_HN.pdf",
         plot = p_dot_HN, width = 10, height = 8)
}

if (!is.null(go_CN_simple)) {
  p_dot_CN <- dotplot(go_CN_simple, showCategory = 15,
                      title = "Neuron ADP_C vs ADP_N (BP, simplified)")
  print(p_dot_CN)
  # ========== EXPORT FIGURE 4 ==========
  ggsave("./results/figures/Neuron_GO_dotplot_CN.pdf",
         plot = p_dot_CN, width = 10, height = 8)
}

# =============================================
# 14. Custom GO bar plot (selected pathways)
# =============================================

cat("[", Sys.time(), "] Generating custom GO bar plot ...\n")

# Combine both GO results
all_go_list <- list()

if (!is.null(go_HN_simple)) {
  df_H <- as.data.frame(go_HN_simple)
  df_H$Group <- "ADP_H_vs_ADP_N"
  all_go_list[["H"]] <- df_H
}

if (!is.null(go_CN_simple)) {
  df_C <- as.data.frame(go_CN_simple)
  df_C$Group <- "ADP_C_vs_ADP_N"
  all_go_list[["C"]] <- df_C
}

if (length(all_go_list) > 0) {
  
  all_go <- bind_rows(all_go_list)
  
  # ========== EXPORT TABLE 5 ==========
  write.csv(all_go, "./results/tables/Neuron_GO_all_BP_simplified.csv", row.names = FALSE)
  
  # Define selected GO terms
  common_GOs <- c("GO:0050807", "GO:0050803", "GO:0016050", "GO:0007416",
                  "GO:0016358", "GO:0031346", "GO:1990778")
  H_special <- c("GO:0099003", "GO:0006836", "GO:0060078")
  C_special <- c()
  all_GOs <- c(common_GOs, H_special, C_special)
  
  # Define display order
  order_df <- data.frame(
    ID = c("GO:0050807", "GO:0050803", "GO:0016050", "GO:0007416",
           "GO:0016358", "GO:0031346", "GO:1990778",
           "GO:0099003", "GO:0006836", "GO:0060078"),
    order = 1:10
  )
  
  # Prepare plot data
  plot_data <- all_go %>%
    filter(ID %in% all_GOs) %>%
    select(ID, Description, Group, GeneRatio, p.adjust) %>%
    mutate(
      Ratio = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      log10_padj = -log10(p.adjust),
      Group_short = ifelse(Group == "ADP_H_vs_ADP_N", "H", "C")
    ) %>%
    left_join(order_df, by = "ID")
  
  if (nrow(plot_data) > 0) {
    
    # Split by group
    plot_H <- plot_data %>% filter(Group_short == "H")
    plot_C <- plot_data %>% filter(Group_short == "C")
    
    # Create H plot
    p_H_bar <- NULL
    if (nrow(plot_H) > 0) {
      plot_H <- plot_H %>%
        arrange(order) %>%
        mutate(Description = stringr::str_wrap(Description, width = 30)) %>%
        mutate(Description = factor(Description, levels = unique(Description)))
      
      p_H_bar <- ggplot(plot_H, aes(x = Ratio, y = Description, fill = log10_padj)) +
        geom_col() +
        scale_fill_gradient(low = "#db7ea1", high = "#9F044D",
                            name = "-log10(p.adjust)") +
        labs(title = "ADP_H vs ADP_N (Neurons)", x = "GeneRatio", y = NULL) +
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
    
    # Create C plot
    p_C_bar <- NULL
    if (nrow(plot_C) > 0) {
      plot_C <- plot_C %>%
        arrange(order) %>%
        mutate(Description = stringr::str_wrap(Description, width = 30)) %>%
        mutate(Description = factor(Description, levels = unique(Description)))
      
      p_C_bar <- ggplot(plot_C, aes(x = Ratio, y = Description, fill = log10_padj)) +
        geom_col() +
        scale_fill_gradient(low = "#9DC7DD", high = "#4E79A7",
                            name = "-log10(p.adjust)") +
        labs(title = "ADP_C vs ADP_N (Neurons)", x = "GeneRatio", y = NULL) +
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
    if (!is.null(p_H_bar) && !is.null(p_C_bar)) {
      combined_bar <- p_H_bar + p_C_bar + plot_layout(guides = "collect") &
        theme(legend.position = "right")
    } else if (!is.null(p_H_bar)) {
      combined_bar <- p_H_bar
    } else if (!is.null(p_C_bar)) {
      combined_bar <- p_C_bar
    } else {
      combined_bar <- NULL
    }
    
    if (!is.null(combined_bar)) {
      print(combined_bar)
      # ========== EXPORT FIGURE 5 ==========
      ggsave("./results/figures/Neuron_GO_barplot_selected.pdf",
             plot = combined_bar, width = 12, height = 8)
    }
  } else {
    cat("  No matching GO terms found for bar plot.\n")
  }
} else {
  cat("  No GO results available for bar plot.\n")
}

# =============================================
# 15. Save session info
# =============================================

session_info <- capture.output(sessionInfo())
writeLines(session_info, "./results/session_info_09_Neuron_DEG_GO.txt")

# =============================================
# 16. Summary
# =============================================

cat("\n[", Sys.time(), "] ==========================================\n")
cat("Neuron DEG + GO enrichment completed!\n")
cat("  ADP_H vs ADP_N: Up =", sum(deg_HN$sig == "Up"),
    ", Down =", sum(deg_HN$sig == "Down"), "\n")
cat("  ADP_C vs ADP_N: Up =", sum(deg_CN$sig == "Up"),
    ", Down =", sum(deg_CN$sig == "Down"), "\n")
cat("  Tables saved to: ./results/tables/\n")
cat("  Figures saved to: ./results/figures/\n")
cat("==========================================\n")
cat("[", Sys.time(), "] All done!\n")