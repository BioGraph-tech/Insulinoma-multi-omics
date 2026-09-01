library(GO.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
source('../utils/utils.R')
get_go_gene_table <- function(go_sets, orgdb = org.Hs.eg.db) {
  
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(AnnotationDbi)
  library(GO.db)
  
  # 1️⃣ 展开 GO 列表
  go_long <- enframe(go_sets, name = "main_function", value = "GO") %>%
    unnest_longer(GO)
  
  # 2️⃣ 获取 GO 描述
  go_desc <- tibble(
    GO = unique(go_long$GO),
    Description = Term(unique(go_long$GO))
  )
  
  # 3️⃣ 获取 GO → gene 映射
  go_gene <- AnnotationDbi::select(
    orgdb,
    keys = unique(go_long$GO),
    keytype = "GO",
    columns = c("GO", "ENTREZID", "SYMBOL")
  ) %>%
    rename(Gene = SYMBOL) %>%
    filter(!is.na(Gene))
  
  # 4️⃣ 合并
  res <- go_long %>%
    left_join(go_desc, by = "GO") %>%
    left_join(go_gene, by = "GO") %>%
    distinct(main_function, GO, Description, ENTREZID, Gene) %>%
    arrange(main_function, GO, Gene)
  
  return(res)
}
go_sets <- list(
  translation = "GO:0002181",
  secretory_granule = c("GO:0034774", "GO:0030667", "GO:0045055"),
  DNA_repair = c("GO:0006282"),
  Synapse = c("GO:0007416", "GO:0050808", "GO:0007268")
)
go_gene_long <- get_go_gene_table(go_sets)
writexl::write_xlsx(go_gene_long,'Table/Main_GO_gene.xlsx')
##plot
library(tidyverse)
library(ggpubr)
gene_score <- read.csv('Table/figure2_gene_score.csv')

gene_score$Tissue2 <- case_when(
  gene_score$Tissue == 'PanNET' ~ 'PanNET',
  gene_score$Tissue == 'Islet' ~ 'Normal'
)

p1 <- ggplot(gene_score, aes(secretory_granule_aucell, Synapse_aucell)) +
  ggrastr::geom_point_rast(
    aes(color = Tissue2),
    size = 0.8,
    alpha = 0.2,
    raster.dpi = 600
  ) +
  geom_smooth(
    aes(color = Tissue2),
    method = "lm",
    se = TRUE,
    linewidth = 0.9
  ) +
  stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 3
  ) +
  scale_color_manual(values = disease_color) + 
  facet_wrap(~ Tissue2, scales = "free", ncol = 4) +
  theme_test() +
  labs(
    x = "Secretory granule score",
    y = "Synapse score",
    color = NULL
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )
saveplot(p1,'Secretory_vs_Synapse',width = 4,height = 2)
p2 <- ggplot(gene_score, aes(translation_aucell, Synapse_aucell)) +
  ggrastr::geom_point_rast(
    aes(color = Tissue2),
    size = 0.8,
    alpha = 0.2,
    raster.dpi = 600
  ) +
  geom_smooth(
    aes(color = Tissue2),
    method = "lm",
    se = TRUE,
    linewidth = 0.9
  ) +
  stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 3
  ) +
  scale_color_manual(values = disease_color) + 
  facet_wrap(~ Tissue2, scales = "free", ncol = 4) +
  theme_test() +
  labs(
    x = "Translation score",
    y = "Synapse score",
    color = NULL
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )
saveplot(p2,'Translation_vs_Synapse',width = 4,height = 2)

##phase
ratio_meta <- gene_score %>% 
  dplyr::select(Tissue2,phase) %>%
  group_by(Tissue2) %>%
  dplyr::mutate(total_cell = n())%>%
  ungroup() %>%
  group_by(Tissue2,,phase) %>%
  dplyr::mutate(number_cell = n()) %>% 
  distinct_all() %>%
  dplyr::mutate(ratio = number_cell/total_cell * 100) %>%
  dplyr::select(Tissue2,phase,ratio) %>%
  distinct_all() 


p3 <- ggplot(ratio_meta,
             aes(x = Tissue2,
                 y = ratio,
                 stratum = phase,
                 alluvium = phase,
                 fill = phase)) +
  ggalluvial::geom_flow(alpha = 0.4, knot.pos = 0, width = 0.5) +
  ggalluvial::geom_stratum(width = 0.5, color = NA) +
  #scale_fill_manual(values = ct_color) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL,
       y = "Cell (%)",
       fill = "Phase") +
  theme_test() +
  theme(
    axis.text.x = element_text(size = 6,angle = 0, vjust = 0.5, hjust = 0.5),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8)
  )
#no significant