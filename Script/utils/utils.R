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
  "T" = "#ff9896"
)
scATAC_ct_color <- c(
  "Myeloid" = "#e377c2",
  "Lymphoid" = "#ff9896",
  "Endothelial" = "#aa40fc",
  "Beta" = "#279e68",
  "Pericyte" = "#ffbb78"
)
sample_map <- c(
  "251030T" = "P1",
  "251211T" = "P2",
  "251224T" = "P3",
  "260119T" = "P4",
  "260129T" = "P5"
)
disease_color <- c(
  'Normal'= "#4DBBD5",
  'PanNET' = "#E64B35"
)


show_umap_plot <- function(seurat_obj,group_by = NULL,label = T,label.size = 4,color = cols08){
  if (label) {
    p1 <- DimPlot(seurat_obj,group.by = group_by,label.size = label.size,label = F) +labs(title = NULL) +
      scale_color_manual(values = color)+
      theme_bw() +
      theme(panel.grid = element_blank())
  }else{
    p1 <- DimPlot(seurat_obj,group.by = group_by,label.size = label.size,label = T) +labs(title = NULL) +
      scale_color_manual(values = color)+
      theme_bw() +
      theme(panel.grid = element_blank()) +
      NoLegend()
  }
  return(p1)
}

#source('c://Users/Administrator/Documents/BaiduSyncdisk/Core code/scRNA-seq/function.R')
source('c://Users/armstrong/Desktop/BaiduSyncdisk/Core code/scRNA-seq/function.R')
