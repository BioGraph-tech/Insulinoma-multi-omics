library(tidyverse)
source('../utils.R')
##UMAP
sce_scRNA <- readRDS('../Rds/islet_final_annotated.Rds')
p1 <- show_umap_plot(sce_scRNA,group_by = 'cell_type',label = F,label.size = 2,color = cell_type_colors)
saveplot(p1,filenames = 'UMAP_plot_scRNA',width = 2,height = 2)
sce_scRNA$donor <- as.character(sample_map[as.character(sce_scRNA$sample_name)])
p1.1 <- show_umap_plot(sce_scRNA,group_by = 'donor',label = T,label.size = 2,color = patient_colors)
saveplot(p1.1,filenames = 'UMAP_plot_scRNA_donor',width = 2.5,height = 2)



sce_scATAC <- readRDS('../Rds/scATAC_peak_final_annoated.Rds')
p2 <- show_umap_plot(sce_scATAC,group_by = 'cell_type',label = F,label.size = 2,color = scATAC_ct_color)
saveplot(p2,filenames = 'UMAP_plot_scATAC',width = 2,height = 2)
sce_scATAC$donor <- as.character(sample_map[as.character(sce_scATAC$sample)])
p2.1 <- show_umap_plot(sce_scATAC,group_by = 'donor',label = T,label.size = 2,color = patient_colors)
saveplot(p2.1,filenames = 'UMAP_plot_scATAC_donor',width = 2.5,height = 2)

##Cellular composition------
scRNA_meta <- read.csv('../table/scRNA_PanNET_obs.csv') %>%
  mutate(donor =sample_map[sample_name] )

scRNA_ratio <- scRNA_meta %>%
  dplyr::select(donor,cell_type_level1) %>%
  group_by(donor) %>%
  dplyr::mutate(total_cell = n())%>%
  ungroup() %>%
  group_by(donor,cell_type_level1) %>%
  dplyr::mutate(number_cell = n()) %>% 
  distinct_all() %>%
  dplyr::mutate(ratio = number_cell/total_cell * 100) %>%
  dplyr::select(donor,cell_type_level1,ratio) %>%
  distinct_all() 

p_cell_type_level1_ratio <- ggplot(scRNA_ratio,
                                   aes(x = donor,
                                       y = ratio,
                                       stratum = cell_type_level1,
                                       alluvium = cell_type_level1,
                                       fill = cell_type_level1)) +
  ggalluvial::geom_flow(alpha = 0.4, knot.pos = 0, width = 0.5) +
  ggalluvial::geom_stratum(width = 0.5, color = NA) +
  scale_fill_manual(values = cell_type_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL,
       y = "Cellular composition",
       fill = NULL) +
  theme_test() +
  theme(
    axis.text.x = element_text(size = 8,angle = 0, vjust = 0.5, hjust = 0.5),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8)
  )
saveplot(p_cell_type_level1_ratio,filenames = 'scRNA_ratio',width = 3,height = 1.5)


scATAC_meta <- read.csv('../table/scATAC_PanNET_obs.csv') %>%
  mutate(donor =sample_map[sample] )

scATAC_ratio <- scATAC_meta %>%
  dplyr::select(donor,cell_type) %>%
  group_by(donor) %>%
  dplyr::mutate(total_cell = n())%>%
  ungroup() %>%
  group_by(donor,cell_type) %>%
  dplyr::mutate(number_cell = n()) %>% 
  distinct_all() %>%
  dplyr::mutate(ratio = number_cell/total_cell * 100) %>%
  dplyr::select(donor,cell_type,ratio) %>%
  distinct_all() 

p_cell_type_ratio_scATAC <- ggplot(scATAC_ratio,
                                   aes(x = donor,
                                       y = ratio,
                                       stratum = cell_type,
                                       alluvium = cell_type,
                                       fill = cell_type)) +
  ggalluvial::geom_flow(alpha = 0.4, knot.pos = 0, width = 0.5) +
  ggalluvial::geom_stratum(width = 0.5, color = NA) +
  scale_fill_manual(values = scATAC_ct_color) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL,
       y = "Cellular composition",
       fill = NULL) +
  theme_test() +
  theme(
    axis.text.x = element_text(size = 8,angle = 0, vjust = 0.5, hjust = 0.5),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8)
  )
saveplot(p_cell_type_ratio_scATAC,filenames = 'scATAC_ratio',width = 2,height = 1.5)
##track-----
#devtools::install_github("junjunlab/BioSeqUtils")
library(BioSeqUtils)
library(tidyverse)
import_BigWig <- function (bw_file = NULL, file_name = NULL, chrom = NULL, format = c("bw", 
                                                                     "bedGraph")) 
{
  format <- match.arg(format, choices = c("bw", "bedGraph"))
  bWData <- plyr::ldply(1:length(bw_file), function(x) {
    tmp <- rtracklayer::import(bw_file[x], format = format)
    if (!is.null(chrom)) {
      tmp <- data.frame(tmp) %>% fastplyr::f_filter(seqnames %in% 
                                                      chrom) %>% fastplyr::f_select(-width, -strand)
    }
    else {
      tmp <- data.frame(tmp) %>% dplyr::select(-width, 
                                                   -strand)
    }
    if (is.null(file_name)) {
      if (format == "bw") {
        fixchar <- "/|.bw|.bigwig"
      }
      else {
        fixchar <- "/|.bg|.bedgraph"
      }
      spt <- strsplit(bw_file[x], split = fixchar) %>% 
        unlist()
      sname <- spt[length(spt)]
    }
    else {
      sname <- file_name[x]
    }
    tmp$fileName <- sname
    return(tmp)
  })
  return(bWData)
}

gtf <- rtracklayer::import.gff("../utiils/hg38.gtf",format = "gtf") %>% 
  data.frame() 
coverage_bw_file <- list.files('../track/',full.names = T)
names(coverage_bw_file) <- tools::file_path_sans_ext(basename(coverage_bw_file)) 
coverage_bw <- import_BigWig(coverage_bw_file)


p_track <- trackVisProMax(Input_gtf = gtf,
                     Input_bw = coverage_bw,
                     Input_gene = c("INS","PECAM1",'LYZ','CD3D','PDGFRB'),
                     Intron_line_type = "line",
                     sample_fill_col = scATAC_ct_color,
                     arrow_rel_len_params_list = list(rel_len = 0.2),
                     peak_fill_col = c("#1f77b4"),
                     peak_width = 0.2,
                     trans_topN = 1,
                     add_gene_region_label = TRUE,
                     base_size = 8,
                     trans_exon_arrow_params = list(color = "black",
                                                    arrow = grid::arrow(
                                                      length = unit(2, "mm"),   
                                                      type = "open"           
                                                    )),
                     trans_exon_col_params = list(fill = "black",color = "black")) +
  theme(legend.position = 'none')
saveplot(p_track,'ATAC_track',width = 4,height = 2)
