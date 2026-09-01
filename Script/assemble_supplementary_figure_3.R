#!/usr/bin/env Rscript

if (dir.exists("G:/Work_NJU/PanNET/")) {
  setwd("G:/Work_NJU/PanNET/")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggrepel)
  library(scales)
})

base_dir <- "core_supplementary"
src_dir <- file.path(base_dir, "source_data")
out_dir <- file.path(base_dir, "main_supplementary_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

disease_color <- c(Normal = "#4DBBD5", PanNET = "#E64B35", Islet = "#4DBBD5")

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      legend.key.size = unit(3.2, "mm"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.margin = margin(4, 4, 4, 4)
    )
}

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  suppressMessages(readr::read_csv(path, show_col_types = FALSE))
}

blank_panel <- function(label) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = label, size = 3) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void()
}

pca <- read_csv_safe(file.path(src_dir, "SuppFig3A_pseudobulk_PCA_source.csv"))
de <- read_csv_safe(file.path(src_dir, "SuppFig3C_RNA_DESeq2_source.csv"))
hm <- read_csv_safe(file.path(src_dir, "SuppFig3D_RNA_DESeq2_top_gene_heatmap_source.csv"))
go <- read_csv_safe(file.path(src_dir, "SuppFig3D_unfiltered_km_gene_GO_source.csv"))

p3a <- if (!is.null(pca)) {
  ggplot(pca, aes(PC1, PC2, color = condition, label = donor)) +
    geom_point(size = 1.7) +
    ggrepel::geom_text_repel(size = 2, max.overlaps = 50) +
    scale_color_manual(values = disease_color) +
    labs(x = "PC1", y = "PC2") +
    theme_pub()
} else blank_panel("PCA source missing")

p3b <- if (!is.null(de)) {
  de <- de %>%
    mutate(
      direction = case_when(
        padj < 0.05 & log2FoldChange > 0 ~ "PanNET-up",
        padj < 0.05 & log2FoldChange < 0 ~ "Normal-up",
        TRUE ~ "NS"
      ),
      neg_log10_fdr = -log10(pmax(padj, .Machine$double.xmin))
    )
  ggplot(de, aes(log2FoldChange, neg_log10_fdr)) +
    geom_point(aes(color = direction), size = 0.35, alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) +
    scale_color_manual(values = c("PanNET-up" = "#B40426", "Normal-up" = "#3B4CC0", "NS" = "grey78")) +
    labs(x = "log2FC PanNET vs normal", y = "-log10 FDR") +
    theme_pub()
} else blank_panel("DESeq2 source missing")

p3c <- if (!is.null(hm)) {
  colnames(hm)[1] <- "gene"
  hml <- hm %>% pivot_longer(-gene, names_to = "sample", values_to = "z")
  ggplot(hml, aes(sample, factor(gene, levels = rev(unique(gene))), fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426", midpoint = 0, limits = c(-3, 3), oob = squish) +
    labs(x = NULL, y = NULL, fill = "Z") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank())
} else blank_panel("Heatmap source missing")

p3d <- if (!is.null(go)) {
  go_top <- go %>%
    mutate(
      gene_cluster = factor(gene_cluster, levels = paste0("km", 1:10)),
      fdr_value = ifelse(!is.na(p.adjust), p.adjust, qvalue),
      score = -log10(fdr_value),
      gene_ratio_value = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),
      label = Description
    ) %>%
    filter(!is.na(gene_cluster), is.finite(score), is.finite(gene_ratio_value)) %>%
    group_by(gene_cluster) %>%
    slice_min(order_by = fdr_value, n = 3, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(label = factor(label, levels = rev(unique(label[order(gene_cluster, score)]))))
  ggplot(go_top, aes(gene_ratio_value, label)) +
    geom_point(aes(size = Count, color = score), alpha = 0.9) +
    facet_wrap(~gene_cluster, scales = "free_y", ncol = 5) +
    scale_color_gradient(low = "#4DBBD5", high = "#E64B35") +
    scale_size(range = c(1.0, 4.2)) +
    labs(x = "Gene ratio", y = NULL, color = "-log10 FDR", size = "Genes") +
    theme_pub() +
    theme(axis.text.y = element_text(size = 5), strip.text = element_text(face = "bold", size = 7), legend.position = "right")
} else blank_panel("GO source missing")

sf3 <- plot_grid(p3a, p3b, p3c, p3d, nrow = 2, labels = c("a", "b", "c", "d"), label_size = 9, rel_widths = c(1, 1.1))
ggsave(file.path(out_dir, "Supplementary_Figure_3_RNA_pseudobulk_robustness.pdf"), sf3, width = 12.5, height = 9, device = cairo_pdf, bg = "white", units = "in")
