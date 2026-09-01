library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(grid)



## 1. 对齐 gene 顺序
targetgene_km_filter <- targetgene_km_filter %>%
  filter(gene %in% rownames(pd_mtx_scale))

pd_mtx_scale_sub <- pd_mtx_scale[targetgene_km_filter$gene, ]

row_cluster <- targetgene_km_filter$new_clutser_new
names(row_cluster) <- targetgene_km_filter$gene



col_group <- meta$disease



## 3. 配色
cluster_levels <- sort(unique(row_cluster))
disease_levels <- unique(col_group)

row_cluster_col <- setNames(
  colorspace::qualitative_hcl(
    length(cluster_levels),
    palette = "Pastel 1"
  ),
  cluster_levels
)

disease_col <- setNames(
  c("#4DBBD5", "#E64B35")[seq_along(disease_levels)],
  disease_levels
)

expr_col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#4575B4", "white", "#D73027")
)


## 4. 整理 GO 结果
## 假设 GO_res$FDR 已经是 -log10(qvalue)
GO_plot_df <- GO_res %>%
  filter(new_clutser_new %in% cluster_levels) %>%
  group_by(new_clutser_new) %>%
  arrange(desc(FDR), .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup()

GO_plot_df$new_clutser_new <- factor(
  GO_plot_df$new_clutser_new,
  levels = cluster_levels
)

GO_plot_df <- GO_plot_df %>%
  arrange(new_clutser_new, desc(FDR))

go_terms <- GO_plot_df$Description
go_score <- GO_plot_df$FDR
go_cluster <- GO_plot_df$new_clutser_new


## 5. 上方列注释
ha_col <- HeatmapAnnotation(
  disease = col_group,
  col = list(disease = disease_col),
  annotation_name_side = "left"
)


## 6. 左侧行 cluster 颜色
ha_row <- rowAnnotation(
  Cluster = row_cluster,
  col = list(Cluster = row_cluster_col),
  show_annotation_name = FALSE,
  width = unit(4, "mm")
)


## 7. 主热图
ht <- Heatmap(
  pd_mtx_scale_sub,
  name = "Relative\nexpression",
  col = expr_col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  
  row_split = row_cluster,
  column_split = col_group,
  
  show_row_names = FALSE,
  show_column_names = FALSE,
  
  top_annotation = ha_col,
  left_annotation = ha_row,
  
  row_title = NULL,
  column_title = NULL,
  
  heatmap_legend_param = list(
    title = "Relative\nexpression",
    at = c(-2, 0, 2),
    labels = c("Min", "0", "Max")
  )
)



library(dplyr)
library(ggplot2)
library(stringr)

GO_plot_df <- GO_res %>%
  filter(new_clutser_new %in% cluster_levels) %>%
  group_by(new_clutser_new) %>%
  arrange(desc(FDR), .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup() %>%
  mutate(
    new_clutser_new = factor(new_clutser_new, levels = cluster_levels),
    Description = str_wrap(Description, width = 40)
  ) %>%
  arrange(new_clutser_new, desc(FDR)) %>%
  mutate(
    Description = factor(Description, levels = rev(unique(Description)))
  )

max_val <- max(GO_plot_df$FDR)

p_go <- ggplot(GO_plot_df, aes(y = Description)) +
  
  ## bar（核心）
  geom_col(
    aes(x = FDR, fill = new_clutser_new),
    width = 0.7
  ) +
  
  ## label（放在 bar 起点附近，但在右侧区域）
  geom_text(
    aes(
      x = 0.5,   # ⭐ 关键：始终在正坐标
      label = Description
    ),
    hjust = 0,
    size = 3.2
  ) +
  
  scale_fill_manual(values = row_cluster_col) +
  
  scale_x_continuous(
    limits = c(0, max_val * 1.1),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "GO enriched pathway",
    x = expression(-log[10]("qvalue")),
    y = NULL
  ) +
  
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    
    ## 去掉左侧轴文字（因为已经手动画了）
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.line.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.4),
    
    plot.title = element_text(face = "bold", size = 12)
  )

library(ComplexHeatmap)
library(grid)
library(cowplot)

## 1. 抓取 heatmap
ht_grob <- grid.grabExpr({
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE
  )
})

## 2. 组合 heatmap + GO plot
p_final <- plot_grid(
  ht_grob,
  p_go,
  nrow = 1,
  align = "h",
  rel_widths = c(1.15, 0.85)
)

saveplot(p_final,'figure2_main_final',width = 6,height = 4)


row_split_vec <- factor(
  row_cluster[rownames(pd_mtx_scale_sub)],
  levels = cluster_levels
)

## 每个 cluster 对应 2 个 GO term
GO_text_df <- GO_res %>%
  filter(new_clutser_new %in% cluster_levels) %>%
  group_by(new_clutser_new) %>%
  arrange(desc(FDR), .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  summarise(
    text = list(str_wrap(Description, width = 35)),
    .groups = "drop"
  )

GO_text <- GO_text_df$text
names(GO_text) <- GO_text_df$new_clutser_new

## 按 cluster_levels 重新排序
GO_text <- GO_text[as.character(cluster_levels)]


go_anno <- rowAnnotation(
  GO = anno_textbox(
    align_to = row_split_vec,
    text = GO_text,
    which = "row",
    by = "anno_block",
    
    
    ## 文字样式
    gp = gpar(
      col = "black",
      fontsize = 8,
      fontface = "bold"
    ),
    
    ## 框的样式
    background_gp = gpar(
      fill = "grey90",
      col = NA
    ),
    
    ## 让多个 GO term 分行显示
    add_new_line = TRUE,
    word_wrap = TRUE
  )
)
ht2 <- Heatmap(
  pd_mtx_scale_sub,
  name = "Relative\nexpression",
  col = expr_col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  
  row_split = row_cluster,
  column_split = col_group,
  
  show_row_names = FALSE,
  show_column_names = FALSE,
  
  top_annotation = ha_col,
  left_annotation = ha_row,
  right_annotation = go_anno,
  row_title = NULL,
  column_title = NULL,
  
  heatmap_legend_param = list(
    title = "Relative\nexpression",
    at = c(-2, 0, 2),
    labels = c("Min", "0", "Max")
  )
)


saveplot(ht2,'figure2_temp',width = 3,height = 1.5)
