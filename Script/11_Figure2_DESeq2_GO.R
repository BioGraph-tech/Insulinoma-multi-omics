library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
source('../utils/utils.R')
gtf <- rtracklayer::import.gff("../utils/hg38.gtf",format = "gtf") %>% 
  data.frame() 
protein_coding_gene <- gtf %>% dplyr::filter(type == 'gene',gene_type  == 'protein_coding')
saveRDS(protein_coding_gene,'../utils/hg38_protein_coding.Rds')
protein_coding_gene <- readRDS('../utils/hg38_protein_coding.Rds') %>% dplyr::filter(!str_detect(gene_name,'^MT-'))
sce <- readRDS('../Rds/islet_beta_PanNET_intergtated.Rds')

candidate <- table(sce$batch) %>% as.data.frame() %>% dplyr::filter(Freq > 200)
# candidate <- candidate %>% 
#   mutate(sample_group = if_else(str_detect(Var1,'HPAP'),'Normal','PanNET')) %>%
#   dplyr::group_by(sample_group) %>%
#   slice_max(n = 5,order_by = Freq)

sce <- subset(sce,batch %in% candidate$Var1)
sce$batch <- factor(sce$batch,levels = unique(sce$batch))

candidate_gene <- intersect(protein_coding_gene$gene_name,rownames(sce))
sce <- sce[candidate_gene,]
sce <- NormalizeData(sce)
saveRDS(sce,file = './Rds/islet_beta_PanNET_intergtated_filter.Rds')
DotPlot(sce,features = c('INS','IAPP','PCSK1','PCSK2'),group.by = 'Tissue')
pb <- AggregateExpression(
  sce,
  group.by = "batch",
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)


pb_count <- pb$RNA %>% as.matrix()
saveRDS(pb_count,'Rds/scRNA_exp_normal_pannet_counts.Rds')

pd_mtx <- edgeR::cpm(pb$RNA)
saveRDS(pd_mtx,'Rds/scRNA_exp_normal_pannet.Rds')
pd_mtx_stat <- tibble(
  gene = rownames(pd_mtx),
  mean = rowMeans(pd_mtx)
)
pd_mtx_stat_filter <- pd_mtx_stat %>% dplyr::filter(mean > 0.1)

pd_mtx_scale <- scale_mtx(pd_mtx[pd_mtx_stat_filter$gene,])
rk <-  kmeans(pd_mtx_scale,centers = 10)

targetgene_km <- tibble(
  gene =  rownames(pd_mtx_scale),
  km_cluster = rk$cluster
)

meta <- tibble(
  sample = colnames(pd_mtx_scale),
  disease = if_else(str_detect(sample,'HPAP'),'Normal','PanNET')
)

p <- Heatmap(
  pd_mtx_scale[targetgene_km$gene,meta$sample],
  cluster_rows = T,
  cluster_columns = T,
  show_row_names = F,
  show_column_names = T,
  row_split = targetgene_km$km_cluster,
  top_annotation = HeatmapAnnotation(
    disease = meta$disease
  ),
  col = circlize::colorRamp2(
    c(-3, 0, 3),
    c("#3B4CC0", "#F7F7F7", "#B40426")
  ),
  name = 'Z-score'
)
savepng(p,filenames = 'RNA pattern',width = 6,height = 6)
saveRDS(targetgene_km,'Rds/targetgene_km.Rds')
targetgene_km <- readRDS('Rds/targetgene_km.Rds')


library(clusterProfiler)
library(org.Hs.eg.db)
km_gene_GO <- purrr::map(
  sort(unique(targetgene_km$km_cluster)),
  function(i) {
    km_gene <- targetgene_km %>% dplyr::filter(km_cluster == i)
    
    genes <- unique(km_gene$gene)
    genes <- genes[!is.na(genes)]
    genes <- genes[genes != ""]
    genes <- as.character(genes)
    
    clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  },
  .progress = TRUE
)
names(km_gene_GO) <- paste0('km',sort(unique( targetgene_km$km_cluster)))
saveRDS(km_gene_GO,'Rds/km_gene_GO.Rds')
km_gene_GO_plot_list <- imap(km_gene_GO, function(obj,nm){
  dotplot(obj,title =nm)
})
km_gene_GO_plot <- cowplot::plot_grid(plotlist = km_gene_GO_plot_list,nrow = 2)
saveplot(km_gene_GO_plot,filenames = 'Km_GO_plot',width = 20,height = 6)
library(DESeq2)
count_mat <- pb$RNA %>% as.matrix()
meta <- tibble(
  sample = colnames(count_mat),
  disease = if_else(str_detect(sample, "HPAP"), "Normal", "PanNET")
) %>%
  mutate(disease = factor(disease, levels = c("Normal", "PanNET")))
all(meta$sample == colnames(count_mat))
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(count_mat),
  colData = meta,
  design = ~ disease
)
ds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("disease", "PanNET", "Normal"))
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    change = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up_in_PanNET",
      padj < 0.05 & log2FoldChange < -1 ~ "Up_in_Normal",
      TRUE ~ "NS"
    )
  )



res_plot <- res_df %>%
  mutate(
    log10padj = -log10(padj),
    group = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(res_plot, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(aes(color = group), size = 1) +
  scale_color_manual(values = c(
    "Up" = "#B40426",
    "Down" = "#3B4CC0",
    "NS" = "grey80"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    title = "PanNET vs Normal (DESeq2)"
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank()
  )

##NMF-----
top_gene <- read.csv('Table/cNMF_top_genes.csv')
library(clusterProfiler)
library(org.Hs.eg.db)
NMF_gene_GO <- purrr::map(
  colnames(top_gene)[-1],
  function(i) {
    NMF_gene <- top_gene[[i]]

    clusterProfiler::enrichGO(
      gene = NMF_gene,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  },
  .progress = TRUE
)
names(NMF_gene_GO) <- paste0('nmf',str_extract(colnames(top_gene)[-1],'\\d+'))
saveRDS(NMF_gene_GO,'Rds/NMF_gene_GO.Rds')


nmf_gene_GO_plot_list <- purrr::imap(NMF_gene_GO, function(obj,nm){
  dotplot(obj,title =nm,label_format = 50)
})
nmf_gene_GO_plot <- cowplot::plot_grid(plotlist = nmf_gene_GO_plot_list,nrow = 2)
saveplot(nmf_gene_GO_plot,filenames = 'NMF_GO_plot',width = 12,height = 4)


top_gene_df <- top_gene %>% pivot_longer(cols = -X,names_to = 'cNMF_cluster',values_to = 'gene') %>%
  mutate(cNMF = str_extract(cNMF_cluster,'\\d+')) %>%
  dplyr::filter(gene %in% rownames(pd_mtx_scale))



p <- Heatmap(
  pd_mtx_scale[top_gene_df$gene,meta$sample],
  cluster_rows = T,
  cluster_columns = T,
  show_row_names = F,
  show_column_names = T,
  row_split = top_gene_df$cNMF,
  top_annotation = HeatmapAnnotation(
    disease = meta$disease
  ),
  col = circlize::colorRamp2(
    c(-1, 0, 1),
    c("#3B4CC0", "#F7F7F7", "#B40426")
  ),
  name = 'Z-score'
)
#endoplasmic reticulum lumen ↓
#cytoplasmic translation ↓
#regulation of insulin secretion ↑
#immunity ↑
# DNA repair ↑
#  neuron function/synapse assemble ↑
norm_res_df <- res_df %>% dplyr::filter(change  == "Up_in_Normal")

norm_res_df <- targetgene_km %>% dplyr::filter(km_cluster == 3)  %>% 
  left_join(res_df,by = 'gene')

##-----

# Translation score ↓ km1 & km5
# Secretion km8 / Neuroendocrine score ↑ km4
# DNA repair score ↑ km6
# Immune / inflammatory score ↑km3


p <- Heatmap(
  pd_mtx_scale[targetgene_km$gene,meta$sample],
  cluster_rows = F,
  cluster_columns = T,
  show_row_names = F,
  show_column_names = T,
  row_split = targetgene_km$km_cluster,
  top_annotation = HeatmapAnnotation(
    disease = meta$disease
  ),
  col = circlize::colorRamp2(
    c(-3, 0, 3),
    c("#3B4CC0", "#F7F7F7", "#B40426")
  ),
  name = 'Z-score'
)
savepng(p,filenames = 'RNA pattern',width = 6,height = 6)


km_gene_GO <- readRDS('Rds/km_gene_GO.Rds')
targetgene_km <- readRDS('Rds/targetgene_km.Rds')
targetgene_km$km_cluster <- paste0('km',targetgene_km$km_cluster)
targetgene_km$new_clutser <- case_when(targetgene_km$km_cluster %in% c('km1','km5') ~ 'C1',
          targetgene_km$km_cluster %in% c('km8') ~ 'C2',
          targetgene_km$km_cluster %in% c('km4') ~ 'C3',
          targetgene_km$km_cluster %in% c('km6') ~ 'C4',
          targetgene_km$km_cluster %in% c('km3') ~ 'C5',
          targetgene_km$km_cluster %in% c('km7') ~ 'C6',
          targetgene_km$km_cluster %in% c('km2') ~ 'C7',
          targetgene_km$km_cluster %in% c('km9','km10') ~ 'C8'
          )

DEGs_filter <- res_plot %>% dplyr::filter(abs(log2FoldChange) > 1, padj<= 0.05)
targetgene_km_filter <- targetgene_km %>% dplyr::filter(gene %in% DEGs_filter$gene)

p <- Heatmap(
  pd_mtx_scale[targetgene_km_filter$gene,meta$sample],
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = F,
  show_column_names = T,
  row_split = targetgene_km_filter$new_clutser,
  top_annotation = HeatmapAnnotation(
    disease = meta$disease
  ),
  col = circlize::colorRamp2(
    c(-3, 0, 3),
    c("#3B4CC0", "#F7F7F7", "#B40426")
  ),
  name = 'Z-score'
)
savepng(p,filenames = 'RNA pattern 2',width = 6,height = 6)


new_clutser_gene_GO <- purrr::map(
  sort(unique(targetgene_km_filter$new_clutser)),
  function(i) {
    km_gene <- targetgene_km_filter %>% dplyr::filter(new_clutser == i)
    
    genes <- unique(km_gene$gene)
    genes <- genes[!is.na(genes)]
    genes <- genes[genes != ""]
    genes <- as.character(genes)
    
    clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  },
  .progress = TRUE
)
names(new_clutser_gene_GO) <- sort(unique(targetgene_km_filter$new_clutser))
saveRDS(new_clutser_gene_GO,'Rds/new_clutser_gene_GO_gene_GO.Rds')

km_gene_GO_plot_list <- imap(new_clutser_gene_GO, function(obj,nm){
  dotplot(obj,title =nm)
})
km_gene_GO_plot <- cowplot::plot_grid(plotlist = km_gene_GO_plot_list,nrow = 2)
saveplot(km_gene_GO_plot,filenames = 'Km_GO_plot_NEW',width = 20,height = 6)


targetgene_km_filter$new_clutser_anno <- case_when(targetgene_km_filter$new_clutser %in% c('C1','C6') ~ c('C1#GO:0002181;GO:0042254'),
                                            targetgene_km_filter$new_clutser %in% c('C2') ~ 'C2#GO:0034774;GO:0031983',
                                            targetgene_km_filter$new_clutser %in% c('C3') ~ 'C3#GO:0007416;GO:0106027',
                                            targetgene_km_filter$new_clutser %in% c('C4') ~ 'C4#GO:0006302;GO:0006282',
                                            targetgene_km_filter$new_clutser %in% c('C5') ~ 'C5#GO:0002443;GO:0002697',
                                            targetgene_km_filter$new_clutser %in% c('C7') ~ 'C6#GO:0007005;GO:0006754',
                                            targetgene_km_filter$new_clutser %in% c('C8') ~ 'C7#GO:0042391;GO:0006874',
)
targetgene_km_filter$new_clutser_new <- str_extract(targetgene_km_filter$new_clutser_anno,'C\\d+')
targetgene_km_filter$new_clutser_new <- factor(targetgene_km_filter$new_clutser_new,levels = sort(unique(targetgene_km_filter$new_clutser_new)))
saveRDS(targetgene_km_filter,'Rds/targetgene_group_final.Rds')

new_clutser_gene_GO <- purrr::map(
  sort(unique(targetgene_km_filter$new_clutser_new)),
  function(i) {
    km_gene <- targetgene_km_filter %>% dplyr::filter(new_clutser_new == i)
    
    genes <- unique(km_gene$gene)
    genes <- genes[!is.na(genes)]
    genes <- genes[genes != ""]
    genes <- as.character(genes)
    
    clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  },
  .progress = TRUE
)
names(new_clutser_gene_GO) <- sort(unique(targetgene_km_filter$new_clutser_new))
saveRDS(new_clutser_gene_GO,'Rds/new_clutser_gene_GO_gene_GO.Rds')

  
GO_res <- map_dfr(names(new_clutser_gene_GO),function(cluster){
  go_term <- targetgene_km_filter %>% dplyr::filter(new_clutser_new == cluster) %>% pull(new_clutser_anno) %>% unique() %>% str_split_1('#')
  go_term <- go_term[2] %>% str_split_1(';')
  new_clutser_gene_GO[[cluster]]@result %>% as.data.frame() %>%
    dplyr::filter(ID %in% go_term) %>%
    dplyr::select(ID,Description,geneID,qvalue) %>%
    mutate(new_clutser_new = cluster,
           FDR = -log10(qvalue))
  
})
saveRDS(GO_res,'Rds/GO_res.Rds')
saveRDS(pd_mtx_scale,'Rds/pd_mtx_scale_res.Rds')


