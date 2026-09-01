#!/usr/bin/env Rscript
if (dir.exists("G:/Work_NJU/PanNET/")) {
  setwd("G:/Work_NJU/PanNET/")
} else if (dir.exists("/mnt/g/Work_NJU/PanNET/")) {
  setwd("/mnt/g/Work_NJU/PanNET/")
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggrepel)
  library(viridis)
  library(scales)
})

base_dir <- "core_supplementary"
src_dir <- file.path(base_dir, "source_data")
tab_dir <- file.path(base_dir, "tables")
out_dir <- file.path(base_dir, "main_supplementary_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_type_colors <- c(
  "Acinar" = "#1f77b4",
  "B" = "#ff7f0e",
  "Beta" = "#279e68",
  "Cycling cell" = "#d62728",
  "Endothelial" = "#aa40fc",
  "Lymphatic Endothelial" = "#8c564b",
  "Macrophage" = "#e377c2",
  "Mast" = "#b5bd61",
  "Monocyte" = "#17becf",
  "NK" = "#aec7e8",
  "Pericyte" = "#ffbb78",
  "Stellate" = "#98df8a",
  "T" = "#ff9896",
  "Myeloid" = "#e377c2",
  "Lymphoid" = "#ff9896"
)
scATAC_ct_color <- c(
  "Myeloid" = "#e377c2",
  "Lymphoid" = "#ff9896",
  "Endothelial" = "#aa40fc",
  "Beta" = "#279e68",
  "Pericyte" = "#ffbb78"
)
patient_colors <- c(P1 = "#4DBBD5", P2 = "#00A087", P3 = "#3C5488", P4 = "#F39B7F", P5 = "#E64B35")
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

save_supp <- function(plot, filename, width, height) {
  ggsave(file.path(out_dir, filename), plot, width = width, height = height, device = cairo_pdf, bg = "white", units = "in")
}

downsample_df <- function(df, n = 35000) {
  if (is.null(df) || nrow(df) <= n) return(df)
  df[sample.int(nrow(df), n), , drop = FALSE]
}

plot_qc_box <- function(df, group, value, ylab, log_y = FALSE) {
  p <- ggplot(df, aes(.data[[group]], .data[[value]], fill = .data[[group]])) +
    geom_boxplot(width = 0.65, outlier.shape = NA, linewidth = 0.25) +
    scale_fill_manual(values = patient_colors, na.value = "grey70") +
    labs(x = NULL, y = ylab) +
    theme_pub() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  if (log_y) p <- p + scale_y_log10(labels = label_number())
  p
}

plot_umap <- function(df, color_col, colors = NULL, title = NULL, n = 40000, size = 0.12) {
  df <- downsample_df(df, n)
  p <- ggplot(df, aes(UMAP1, UMAP2, color = .data[[color_col]])) +
    geom_point(size = size, alpha = 0.75, stroke = 0) +
    labs(x = "UMAP1", y = "UMAP2", title = title) +
    theme_pub() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "right"
    )
  if (!is.null(colors)) p <- p + scale_color_manual(values = colors, na.value = "grey70")
  p
}

plot_dot <- function(df, x, y, size_col, color_col, title = NULL, palette = "viridis") {
  ggplot(df, aes(.data[[x]], .data[[y]])) +
    geom_point(aes(size = .data[[size_col]], color = .data[[color_col]]), alpha = 0.9) +
    scale_size(range = c(0.3, 4.2)) +
    viridis::scale_color_viridis(option = "C") +
    labs(x = NULL, y = NULL, title = title, size = "% cells", color = "Mean") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
}

# Supplementary Figure 1 ------------------------------------------------------
sc_qc <- read_csv_safe(file.path(src_dir, "SuppFig1B_scRNA_QC_source.csv"))
atac_qc <- read_csv_safe(file.path(src_dir, "SuppFig1C_scATAC_QC_source.csv"))
sc_contrib <- read_csv_safe(file.path(src_dir, "SuppFig1E_scRNA_patient_celltype_contribution.csv"))
atac_contrib <- read_csv_safe(file.path(src_dir, "SuppFig1E_scATAC_patient_celltype_contribution.csv"))

p1a <- if (!is.null(sc_qc)) plot_qc_box(sc_qc, "donor_label", "n_UMIs", "UMIs per cell", TRUE) else blank_panel("scRNA QC missing")
p1b <- if (!is.null(sc_qc)) plot_qc_box(sc_qc, "donor_label", "n_genes_by_counts", "Genes per cell", TRUE) else blank_panel("scRNA genes missing")
p1c <- if (!is.null(sc_qc)) plot_qc_box(sc_qc, "donor_label", "pct_counts_mt", "Mitochondrial reads (%)") else blank_panel("scRNA mito missing")
p1d <- if (!is.null(atac_qc)) plot_qc_box(atac_qc, "donor_label", "n_fragment", "Fragments per nucleus", TRUE) else blank_panel("scATAC fragments missing")
p1e <- if (!is.null(atac_qc)) plot_qc_box(atac_qc, "donor_label", "tsse", "TSS enrichment") else blank_panel("TSS missing")
#p1f <- if (!is.null(atac_qc)) plot_qc_box(atac_qc, "donor_label", "frac_dup", "Duplicate fraction") else blank_panel("Duplicate fraction missing")
p1f <- if (!is.null(atac_qc)) {
  ggplot(atac_qc, aes(n_fragment, tsse, color = donor_label)) +
    geom_point(size = 0.25, alpha = 0.35) +
    geom_hline(yintercept = 5, linetype = "dashed", linewidth = 0.25) +
    geom_vline(xintercept = c(1000, 60000), linetype = "dashed", linewidth = 0.25) +
    scale_x_log10(labels = label_number()) +
    scale_color_manual(values = patient_colors) +
    labs(x = "Fragments per nucleus", y = "TSS enrichment") +
    theme_pub() +
    theme(legend.position = "none")
} else blank_panel("2D ATAC QC missing")
# p1h <- if (!is.null(sc_contrib)) {
#   ggplot(sc_contrib, aes(patient_fraction, cell_type)) +
#     geom_col(aes(fill = donor_label), position = "stack", width = 0.75) +
#     scale_fill_manual(values = patient_colors) +
#     labs(x = "Fraction within patient", y = NULL) +
#     theme_pub()
# } else blank_panel("scRNA contribution missing")
# p1i <- if (!is.null(atac_contrib)) {
#   ggplot(atac_contrib, aes(patient_fraction, cell_type)) +
#     geom_col(aes(fill = donor_label), position = "stack", width = 0.75) +
#     scale_fill_manual(values = patient_colors[c("P4", "P5")]) +
#     labs(x = "Fraction within patient", y = NULL) +
#     theme_pub()
# } else blank_panel("scATAC contribution missing")
sf1 <- plot_grid(
  plot_grid(p1a, p1b, p1c, nrow = 1, labels = c("a", "b", "c"), label_size = 9),
  plot_grid(p1d, p1e, p1f, nrow = 1, labels = c("d", "e", "f" ), label_size = 9),
  ncol = 1, rel_heights = c(1, 1, 1)
)
save_supp(sf1, "Supplementary_Figure_1_QC_and_sample_contribution.pdf", width = 5, height = 3)

# Supplementary Figure 2 ------------------------------------------------------
umap_donor <- read_csv_safe(file.path(src_dir, "SuppFig2A_scRNA_donor_label_UMAP_source.csv"))
umap_ct <- read_csv_safe(file.path(src_dir, "SuppFig2A_scRNA_cell_type_UMAP_source.csv"))
markers <- read_csv_safe(file.path(src_dir, "SuppFig2B_marker_dotplot_source.csv"))
cnv <- read_csv_safe(file.path(src_dir, "SuppFig2F_CNV_score_source.csv"))
#p2a <- if (!is.null(umap_donor)) plot_umap(umap_donor, "donor_label", patient_colors, "Patient") else blank_panel("UMAP patient missing")
#p2b <- if (!is.null(umap_ct)) plot_umap(umap_ct, "cell_type", cell_type_colors, "Cell type") else blank_panel("UMAP cell type missing")
#p2c <- if (!is.null(markers)) plot_dot(markers, "gene", "cell_type", "pct_expressed", "mean_expression", "Marker expression") else blank_panel("Marker source missing")
p2a <- if (!is.null(cnv)) {
  cnv <- cnv %>%
    mutate(
      donor_label = recode(
        sample_name,
        "251030T" = "P1",
        "251211T" = "P2",
        "251224T" = "P3",
        "260119T" = "P4",
        "260129T" = "P5",
        .default = sample_name
      ),
      donor_label = factor(donor_label, levels = c("P1", "P2", "P3", "P4", "P5")),
      cell_type = factor(cell_type, levels = names(cell_type_colors)[names(cell_type_colors) %in% unique(cell_type)])
    ) %>%
    filter(!is.na(donor_label), !is.na(cell_type), !is.na(cnv_score))
  ggplot(cnv, aes(donor_label, cnv_score, fill = donor_label)) +
    geom_violin(width = 0.85, scale = "width", trim = TRUE, color = NA, alpha = 0.85) +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.22, fill = "white", alpha = 0.85) +
    facet_wrap(~cell_type, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = patient_colors) +
    labs(x = NULL, y = "CNV score") +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 6.5)
    )
} else blank_panel("CNV source missing")
sf2 <- plot_grid(
  plot_grid(p2a, p2b, nrow = 1, labels = c("a", "b"), label_size = 9, rel_widths = c(1, 1.2)),
  plot_grid(p2c, p2d, nrow = 1, labels = c("c", "d"), label_size = 9, rel_widths = c(1.45, 0.75)),
  ncol = 1, rel_heights = c(1, 1.15)
)
save_supp(sf2, "Supplementary_Figure_2_annotation_lineage_CNV.pdf", 13.5, 9)

# Supplementary Figure 3 ------------------------------------------------------
pca <- read_csv_safe(file.path(src_dir, "SuppFig3A_pseudobulk_PCA_source.csv"))
de <- read_csv_safe(file.path(src_dir, "SuppFig3C_RNA_DESeq2_source.csv"))
hm <- read_csv_safe(file.path(src_dir, "SuppFig3D_RNA_DESeq2_top_gene_heatmap_source.csv"))
go <- read_csv_safe(file.path(src_dir, "SuppFig3D_unfiltered_km_gene_GO_source.csv"))
#p3a <- if (!is.null(pca)) {
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
} else blank_panel("DE source missing")
#p3c <- if (!is.null(hm)) {
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
    slice_min(order_by = fdr_value, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(label = factor(label, levels = rev(unique(label[order(gene_cluster, score)]))))
  ggplot(
    go_top,
    aes(
      x = gene_ratio_value,
      y = label
    )
  ) +
    geom_segment(
      aes(
        x = 0,
        xend = gene_ratio_value,
        y = label,
        yend = label
      ),
      linewidth = 0.35,
      color = "grey70"
    ) +
    geom_point(
      aes(
        size = Count,
        color = score
      ),
      alpha = 0.9
    ) +
    facet_wrap(
      ~gene_cluster,
      scales = "free_y",
      ncol = 2
    ) +
    scale_color_gradient(
      low = "#4DBBD5",
      high = "#E64B35"
    ) +
    scale_size(
      range = c(1.0, 4.2)
    ) +
    scale_x_continuous(
      expand = expansion(
        mult = c(0, 0.08)
      )
    ) +
    labs(
      x = "Gene ratio",
      y = NULL,
      color = "-log10 FDR",
      size = "Genes"
    ) +
    theme_pub() +
    theme(
      axis.text.y = element_text(size = 5),
      strip.text = element_text(
        face = "bold",
        size = 7
      ),
      legend.position = "right"
    )
} else blank_panel("GO source missing")
#sf3 <- plot_grid(p3a, p3b, p3c, p3d, nrow = 2, labels = c("a", "b", "c", "d"), label_size = 9, rel_widths = c(1, 1.1))
save_supp(p3b, "Supplementary_Figure_3_RNA_pseudobulk_robustness.pdf", 3.5, 2.5)
save_supp(p3d, "Supplementary_Figure_3_RNA_km_GO.pdf", 7, 6)

# Supplementary Figure 4 ------------------------------------------------------
score <- read_csv_safe(file.path(src_dir, "SuppFig4_module_score_source.csv"))
cors <- read_csv_safe(file.path(src_dir, "SuppFig4C_per_donor_secretory_synaptic_correlations.csv"))
p4a <- if (!is.null(score)) {
  score %>%
    select(donor, secretory_granule_aucell, Synapse_aucell) %>%
    pivot_longer(-donor, names_to = "module", values_to = "score") %>%
    mutate(module = recode(module, secretory_granule_aucell = "Secretory granule", Synapse_aucell = "Synapse")) %>%
    ggplot(aes(donor, score, fill = donor)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
    facet_wrap(~module, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = patient_colors, na.value = "grey75") +
    labs(x = NULL, y = "AUCell score") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
} else blank_panel("Module score source missing")
# p4b <- if (!is.null(score)) {
#   score %>% filter(!is.na(secretory_granule_aucell), !is.na(Synapse_aucell)) %>% downsample_df(40000) %>%
#     ggplot(aes(secretory_granule_aucell, Synapse_aucell, color = Tissue)) +
#     geom_point(size = 0.15, alpha = 0.35) +
#     geom_smooth(method = "lm", linewidth = 0.35, color = "black", se = FALSE) +
#     scale_color_manual(values = disease_color, na.value = "grey70") +
#     labs(x = "Secretory granule score", y = "Synapse score") +
#     theme_pub()
# } else blank_panel("Scatter source missing")
# p4c <- if (!is.null(cors)) {
#   ggplot(cors, aes(donor, spearman_r, fill = donor)) +
#     geom_col(width = 0.7) +
#     geom_hline(yintercept = 0, linewidth = 0.25) +
#     scale_fill_manual(values = patient_colors, na.value = "grey75") +
#     labs(x = NULL, y = "Spearman r") +
#     theme_pub() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
# } else blank_panel("Correlation source missing")
p4b <- if (!is.null(score)) {
  binned <- score %>%
    select(secretory_granule_aucell, Synapse_aucell) %>%
    filter(!is.na(secretory_granule_aucell), !is.na(Synapse_aucell)) %>%
    mutate(bin = ntile(secretory_granule_aucell, 20)) %>%
    group_by(bin) %>% summarise(secretory = mean(secretory_granule_aucell), synapse = mean(Synapse_aucell), .groups = "drop")
  ggplot(binned, aes(secretory, synapse)) +
    geom_point(size = 1.5, color = "#3C5488") +
    geom_smooth(method = "lm", linewidth = 0.35, color = "black", se = TRUE) +
    labs(x = "Binned secretory score", y = "Binned synapse score") +
    theme_pub()
} else blank_panel("Bin source missing")
p4c <- if (!is.null(score)) {
  binned <- score %>%
    select(translation_aucell, Synapse_aucell) %>%
    filter(!is.na(translation_aucell), !is.na(Synapse_aucell)) %>%
    mutate(bin = ntile(translation_aucell, 20)) %>%
    group_by(bin) %>% summarise(translation = mean(translation_aucell), synapse = mean(Synapse_aucell), .groups = "drop")
  ggplot(binned, aes(translation, synapse)) +
    geom_point(size = 1.5, color = "#3C5488") +
    geom_smooth(method = "lm", linewidth = 0.35, color = "black", se = TRUE) +
    labs(x = "Binned translation score", y = "Binned synapse score") +
    theme_pub()
} else blank_panel("Bin source missing")
p4_top <- plot_grid(
  p4a,
  nrow = 1,
  labels = "a",
  label_size = 9
)

p4_bottom <- plot_grid(
  p4b,
  p4c,
  nrow = 1,
  labels = c("b", "c"),
  label_size = 9,
  rel_widths = c(1, 1)
)

sf4 <- plot_grid(
  p4_top,
  p4_bottom,
  ncol = 1,
  rel_heights = c(1, 1)
)
save_supp(sf4, "Supplementary_Figure_4_secretory_synaptic_axis.pdf", 6, 4)

# # Supplementary Figure 5 ------------------------------------------------------
# atac_group <- read_csv_safe(file.path(src_dir, "SuppFig5A_ATAC_beta_sample_group_UMAP_source.csv"))
# atac_sample <- read_csv_safe(file.path(src_dir, "SuppFig5A_ATAC_beta_sample_name_UMAP_source.csv"))
# peak_status <- read_csv_safe(file.path(src_dir, "SuppFig5E_concordant_peak_status_source.csv"))
# p5a <- if (!is.null(atac_group)) plot_umap(atac_group, "sample_group", c(PanNET = "#E64B35", Normal = "#4DBBD5"), "Beta scATAC group", n = 50000, size = 0.12) else blank_panel("ATAC group UMAP missing")
# p5b <- if (!is.null(atac_sample)) plot_umap(atac_sample, "sample_name", NULL, "Beta scATAC sample", n = 50000, size = 0.12) else blank_panel("ATAC sample UMAP missing")
# p5c <- if (!is.null(peak_status)) {
#   ggplot(peak_status, aes(peak_status, count, fill = peak_status)) +
#     geom_col(width = 0.7) +
#     scale_y_continuous(labels = label_number()) +
#     labs(x = NULL, y = "Peaks") +
#     theme_pub() +
#     theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
# } else blank_panel("Peak status missing")
# sf5 <- plot_grid(p5a, p5b, p5c, nrow = 1, labels = c("a", "b", "c"), label_size = 9, rel_widths = c(1, 1, 0.75))
# save_supp(sf5, "Supplementary_Figure_5_ATAC_accessibility_concordance.pdf", 13, 4.2)

# # Supplementary Figure 6 ------------------------------------------------------
# motif <- read_csv_safe(file.path(tab_dir, "Supplementary_Table_7_chromVAR_motif_activity_partial.csv"))
# p6a <- if (!is.null(motif) && all(c("P4_vs_normal", "P5_vs_normal") %in% names(motif))) {
#   names(motif)[1] <- "motif"
#   lab <- motif %>% mutate(abs_eff = abs(mean_panNET_effect)) %>% arrange(desc(abs_eff)) %>% slice_head(n = 12)
#   ggplot(motif, aes(P4_vs_normal, P5_vs_normal)) +
#     geom_point(aes(color = `TF family`), size = 0.8, alpha = 0.75) +
#     geom_hline(yintercept = 0, linewidth = 0.25, color = "grey60") +
#     geom_vline(xintercept = 0, linewidth = 0.25, color = "grey60") +
#     geom_text(data = lab, aes(label = motif), size = 1.6, check_overlap = TRUE, vjust = -0.5) +
#     labs(x = "P4 vs normal chromVAR effect", y = "P5 vs normal chromVAR effect") +
#     guides(color = 'none') +
#     theme_pub()
# } else blank_panel("Motif concordance missing")
# p6b <- if (!is.null(motif) && "mean_panNET_effect" %in% names(motif)) {
#   names(motif)[1] <- "motif"
#   motif %>%
#     arrange(desc(abs(mean_panNET_effect))) %>% slice_head(n = 35) %>%
#     select(motif, P4_vs_normal, P5_vs_normal) %>%
#     pivot_longer(-motif, names_to = "contrast", values_to = "effect") %>%
#     ggplot(aes(contrast, reorder(motif, effect), fill = effect)) +
#     geom_tile() +
#     scale_fill_gradient2(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426", midpoint = 0) +
#     labs(x = NULL, y = NULL, fill = "Effect") +
#     theme_pub()
# } else blank_panel("Motif heatmap missing")
# sf6 <- plot_grid( p6b, nrow = 1, labels = c("a", "b"), label_size = 9, rel_widths = c(1, 0.9))
# save_supp(sf6, "Supplementary_Figure_6_motif_chromVAR_TF_bridge.pdf", 11, 5.5)

# Supplementary Figure 7 ------------------------------------------------------
filter_counts <- read_csv_safe(file.path(src_dir, "SuppFig7B_GRN_filter_counts_available.csv"))
target_counts <- read_csv_safe(file.path(src_dir, "SuppFig7C_TF_target_counts_source.csv"))
metrics <- read_csv_safe(file.path(src_dir, "SuppFig7C_TF_network_metrics_source.csv"))
tf_go <- read_csv_safe(file.path(src_dir, "SuppFig7H_TF_GO_source.csv"))
if (!is.null(target_counts) && !"TF" %in% names(target_counts)) names(target_counts)[1] <- "TF"
p7a <- if (!is.null(filter_counts)) {
  ggplot(filter_counts, aes(count, reorder(filter_step, count))) +
    geom_col(fill = "#3C5488", width = 0.7) +
    scale_x_continuous(labels = label_number()) +
    labs(x = "Edges", y = NULL) +
    theme_pub()
} else blank_panel("GRN filter counts missing")
p7b <- if (!is.null(target_counts)) {
  ggplot(target_counts, aes(target_count, reorder(TF, target_count))) +
    geom_col(fill = "#8491B4", width = 0.7) +
    labs(x = "Targets", y = NULL) +
    theme_pub()
} else blank_panel("TF targets missing")
p7c <- if (!is.null(metrics)) {
  metrics %>% arrange(desc(weighted_out_degree)) %>% slice_head(n = 25) %>%
    ggplot(aes(weighted_out_degree, reorder(TF, weighted_out_degree))) +
    geom_col(fill = "#E64B35", width = 0.7) +
    labs(x = "Weighted out-degree", y = NULL) +
    theme_pub()
} else blank_panel("TF metrics missing")
p7d <- if (!is.null(tf_go)) {
  tf_go %>% arrange(p.adjust) %>% slice_head(n = 22) %>%
    mutate(label = paste(TF, Description, sep = ": "), score = -log10(p.adjust)) %>%
    ggplot(aes(score, reorder(label, score), fill = TF)) +
    geom_col(width = 0.7) +
    labs(x = "-log10 adjusted P", y = NULL) +
    theme_pub() +
    theme(legend.position = "none")
} else blank_panel("TF GO missing")

p7bc <- plot_grid(
  p7b,
  p7c,
  nrow = 1,
  labels = c("a", "b"),
  label_size = 9
)

p7d_label <- plot_grid(
  p7d,
  labels = "c",
  label_size = 9
)

p7d_center <- plot_grid(
  NULL,
  p7d_label,
  NULL,
  nrow = 1,
  rel_widths = c(0.2, 0.6, 0.2)
)

sf7 <- plot_grid(
  p7bc,
  p7d_center,
  ncol = 1,
  rel_heights = c(1, 1)
)
save_supp(sf7,
          "Supplementary_Figure_7_GRN_filtering_robustness.pdf",
          8,
          6)

# Supplementary Figure 9 ------------------------------------------------------
tme <- read_csv_safe(file.path(src_dir, "SuppFig9C_TME_composition_source.csv"))
if (!is.null(tme) && !"donor_label" %in% names(tme)) {
  if ("sample_name" %in% names(tme)) {
    tme <- tme %>% mutate(donor_label = recode(sample_name, "251030T" = "P1", "251211T" = "P2", "251224T" = "P3", "260119T" = "P4", "260129T" = "P5", .default = sample_name))
  }
}
p9a <- if (!is.null(umap_ct)) {
  non_tumor <- umap_ct %>% filter(!grepl("Beta|Tumou?r|PanNET", cell_type, ignore.case = TRUE))
  plot_umap(non_tumor, "cell_type", cell_type_colors, "Non-tumour compartments", n = 30000, size = 0.18)
} else blank_panel("TME UMAP missing")
p9b <- if (!is.null(tme)) {
  ggplot(tme, aes(donor_label, cell_count, fill = cell_type)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = cell_type_colors, na.value = "grey75") +
    scale_y_continuous(labels = label_number()) +
    labs(x = NULL, y = "Cells") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else blank_panel("TME composition missing")
p9c <- if (!is.null(tme)) {
  ggplot(tme, aes(donor_label, fraction_of_tme, fill = cell_type)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = cell_type_colors, na.value = "grey75") +
    scale_y_continuous(labels = percent_format()) +
    labs(x = NULL, y = "Fraction of TME") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else blank_panel("TME fraction missing")
sf9 <- plot_grid(p9a, p9b, p9c, nrow = 1, labels = c("a", "b", "c"), label_size = 9, rel_widths = c(1.15, 0.9, 0.9))
save_supp(sf9, "Supplementary_Figure_9_TME_composition.pdf", 12, 4.5)

manifest <- data.frame(
  figure = c("Supplementary Figure 1", "Supplementary Figure 2", "Supplementary Figure 3", "Supplementary Figure 4", "Supplementary Figure 5", "Supplementary Figure 6", "Supplementary Figure 7", "Supplementary Figure 9"),
  file = c(
    "Supplementary_Figure_1_QC_and_sample_contribution.pdf",
    "Supplementary_Figure_2_annotation_lineage_CNV.pdf",
    "Supplementary_Figure_3_RNA_pseudobulk_robustness.pdf",
    "Supplementary_Figure_4_secretory_synaptic_axis.pdf",
    "Supplementary_Figure_5_ATAC_accessibility_concordance.pdf",
    "Supplementary_Figure_6_motif_chromVAR_TF_bridge.pdf",
    "Supplementary_Figure_7_GRN_filtering_robustness.pdf",
    "Supplementary_Figure_9_TME_composition.pdf"
  )
)
write.csv(manifest, file.path(out_dir, "assembled_supplementary_figure_manifest.csv"), row.names = FALSE)
