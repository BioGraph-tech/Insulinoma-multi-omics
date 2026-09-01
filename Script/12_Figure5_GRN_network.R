###############################################################################
# KK Final -- all labels, TF anti-overlap, TG cleanup, polished
###############################################################################

library(dplyr); library(tidyr); library(tibble)
library(igraph)
library(ggplot2); library(scales); library(ggnewscale); library(ggrepel)

setwd("/mnt/g/Work_NJU/PanNET/Figure4/GRN")
dir.create("GRN_kk_final", recursive=TRUE, showWarnings=FALSE)

# ============================================================================
# 1. DATA PIPELINE
# ============================================================================
cat("Loading data...\n")
grn <- readRDS("GRN_final.Rds")
go  <- readRDS("tf_go_top10_table.Rds")

grn_agg <- grn %>%
  group_by(TF,targetgene) %>%
  summarise(n_regions=n(), log2FoldChange=mean(log2FoldChange), .groups="drop")

main_TFs <- c("ZNF580","RREB1","PLAGL2","KLF2","KLF4","KLF7",
              "EGR1","EGR3","EGR4","RFX3","RFX6","FOS","FOSL2","BHLHE41","ZNF22")

modules_v3 <- list(
  "Stress-adaptive plasticity"=list(
    short="Stress-adaptive\nplasticity",
    tfs=c("EGR1","EGR3","EGR4"),
    candidates=c("ATF3","XBP1","SESN1","PDK4","GABARAPL1","ROBO1","SLIT2",
                 "IGF1R","BMPR2","MAPT","BCL2","STK24","PTPRO","METRN","ATP8A2")
  ),
  "Neuroendocrine vesicle trafficking"=list(
    short="Neuroendocrine\nvesicle trafficking",
    tfs=c("KLF2","KLF4","KLF7"),
    candidates=c("SNAP25","STX1A","VAMP2","SYT1","DNM1","PCLO","RIMS2",
                 "BLOC1S6","SLC18A2","NLGN1","SV2C","RAB27B","AP2S1","BRAF")
  ),
  "Ciliary signal sensing"=list(
    short="Ciliary\nsignal sensing",
    tfs=c("RFX3","RFX6"),
    candidates=c("WDR35","IFT22","IFT27","BBS9","TTC8","KIAA0586",
                 "ODF2","ASAP1","NEURL1","CEP126","CFAP44","CFAP53","ERICH3","BBOF1")
  ),
  "Beta-cell secretory memory"=list(
    short="Beta-cell\nsecretory memory",
    tfs=c("ZNF580","RREB1","PLAGL2"),
    candidates=c("INS","GCK","ABCC8","SLC30A8","PCSK1","TCF7L2","FOXO1",
                 "CHGA","TIAM1","PPP3CA","SYBU","NKX6-1","RFX3","TNS2","BCAR1","HDAC9")
  ),
  "Adhesion and migratory remodeling"=list(
    short="Adhesion &\nmigratory remodeling",
    tfs=c("FOS","FOSL2"),
    candidates=c("DDR1","FERMT2","TLN2","PTK2","TIAM2","TRIO","MYO10",
                 "GIT1","SORBS2","DST","DPYSL3","RUFY3","COBL")
  ),
  "Growth and information-processing rewiring"=list(
    short="Growth &\ninfo-processing\nrewiring",
    tfs=c("BHLHE41","ZNF22"),
    candidates=c("RPTOR","AKT1S1","LAMTOR4","SMARCA2","SMARCC1","PBRM1","UPF1",
                 "SIK3","PIP4P1","XBP1","TENT4B","SLF1","CDK5RAP2","NUP205","NUP214")
  )
)

tf_to_mod <- c()
for (m in names(modules_v3)) for (tf in modules_v3[[m]]$tfs) tf_to_mod[tf] <- m

selected_targets <- list()
for (m in names(modules_v3)) {
  mod <- modules_v3[[m]]
  for (g in mod$candidates) {
    e <- grn_agg %>% filter(TF %in% mod$tfs, targetgene==g)
    if (nrow(e)>0) {
      selected_targets[[paste(m,g)]] <- data.frame(
        module=m, targetgene=g, TF=paste(unique(e$TF),collapse="/"),
        log2FoldChange=mean(e$log2FoldChange), n_regions=max(e$n_regions),
        n_tfs=n_distinct(e$TF), stringsAsFactors=FALSE)
    }
  }
}
all_targets <- bind_rows(selected_targets) %>% filter(!targetgene %in% main_TFs)

tf_metrics <- grn_agg %>% filter(TF %in% main_TFs) %>%
  group_by(TF) %>%
  summarise(weighted_out_degree=sum(n_regions), full_out_degree=n(), .groups="drop") %>%
  mutate(node_size_raw=sqrt(weighted_out_degree), node_size=rescale(node_size_raw,to=c(8,12)))

edges_tf_target <- all_targets %>%
  tidyr::separate_rows(TF,sep="/") %>% select(from=TF, to=targetgene, n_regions)

tt <- grn_agg %>% filter(TF %in% main_TFs & targetgene %in% main_TFs & TF!=targetgene) %>%
  rename(from=TF, to=targetgene)
tt$from_module <- tf_to_mod[tt$from]; tt$to_module <- tf_to_mod[tt$to]
tt$is_cross <- tt$from_module != tt$to_module

edges_tf_tf_display <- tt %>% filter(is_cross) %>%
  left_join(tf_metrics%>%select(TF,weighted_out_degree),by=c("from"="TF")) %>%
  arrange(desc(weighted_out_degree)) %>% head(15) %>%
  distinct(from,to,.keep_all=TRUE)

tg_module <- setNames(all_targets$module, all_targets$targetgene)

module_colors <- c(
  "Stress-adaptive plasticity"="#4DAF4A",
  "Neuroendocrine vesicle trafficking"="#377EB8",
  "Ciliary signal sensing"="#984EA3",
  "Beta-cell secretory memory"="#E41A1C",
  "Adhesion and migratory remodeling"="#FF7F00",
  "Growth and information-processing rewiring"="#A65628"
)
fc_colors <- c("#67A9CF","#F7F7F7","#EF8A62")

tf_radius <- 0.28; tg_radius <- 0.15

# ============================================================================
# 2. KK LAYOUT + POST-PROCESSING
# ============================================================================
all_edges <- bind_rows(
  edges_tf_target %>% select(from,to),
  edges_tf_tf_display %>% select(from,to)
) %>% distinct(from,to)

g <- graph_from_data_frame(all_edges, directed=TRUE)

set.seed(1234)
kk_lay <- layout_with_kk(g)
kk_lay[,1] <- rescale(kk_lay[,1], to=c(-6,6))
kk_lay[,2] <- rescale(kk_lay[,2], to=c(-6,6))

nodes_kk <- bind_rows(
  tf_metrics %>% filter(TF %in% main_TFs) %>%
    mutate(node_id=TF, node_type="TF", node_label=TF,
           module=tf_to_mod[TF], log2FoldChange=NA_real_),
  all_targets %>%
    mutate(node_id=targetgene, node_type="target_gene", node_label=targetgene,
           log2FoldChange=log2FoldChange) %>%
    group_by(node_id) %>% slice(1) %>% ungroup()
) %>% mutate(
  module = ifelse(is.na(module), tg_module[node_id], module),
  plot_log2FC = pmax(pmin(log2FoldChange,2),-2),
  node_size   = ifelse(node_type=="TF", node_size, 5),
  x = kk_lay[match(node_id, V(g)$name), 1],
  y = kk_lay[match(node_id, V(g)$name), 2]
)

# ============================================================================
# 3. REMOVE ISOLATED TG NODES (far from their module center)
# ============================================================================

cat("Checking for isolated TG nodes...\n")
nodes_clean <- nodes_kk
removed <- c()

for (m in names(modules_v3)) {
  mod_tfs <- modules_v3[[m]]$tfs
  tf_pos <- nodes_clean %>% filter(node_id %in% mod_tfs)
  cx <- mean(tf_pos$x, na.rm=TRUE); cy <- mean(tf_pos$y, na.rm=TRUE)

  mod_tgs <- all_targets$targetgene[all_targets$module==m]
  for (tg in mod_tgs) {
    idx <- which(nodes_clean$node_id == tg)
    if (length(idx)!=1) next
    dist_to_center <- sqrt((nodes_clean$x[idx]-cx)^2 + (nodes_clean$y[idx]-cy)^2)
    # Compute median distance of TGs in this module
    all_tg_dists <- sapply(mod_tgs, function(g) {
      i <- which(nodes_clean$node_id==g)
      if(length(i)!=1) return(NA)
      sqrt((nodes_clean$x[i]-cx)^2 + (nodes_clean$y[i]-cy)^2)
    })
    med_dist <- median(all_tg_dists, na.rm=TRUE)
    mad_dist <- mad(all_tg_dists, na.rm=TRUE, constant=1)
    if (is.na(mad_dist) || mad_dist<0.5) mad_dist <- 0.5
    # Flag if >3 MAD from median
    if (!is.na(dist_to_center) && dist_to_center > med_dist + 3*mad_dist) {
      cat(sprintf("  Removing %s (dist=%.2f, median=%.2f, mad=%.2f)\n",
                  tg, dist_to_center, med_dist, mad_dist))
      removed <- c(removed, idx)
    }
  }
}

if (length(removed)>0) {
  nodes_clean <- nodes_clean[-unique(removed), ]
  cat(sprintf("  Removed %d isolated TG nodes\n", length(unique(removed))))
} else {
  cat("  No isolated nodes found\n")
}

# Also remove any nodes with NA coordinates
nodes_clean <- nodes_clean %>% filter(!is.na(x) & !is.na(y))

# ============================================================================
# 4. TF ANTI-OVERLAP: per-module centering + local TF spacing
# ============================================================================

cat("Adjusting TF positions to prevent overlap...\n")

# First: per-module TF centering (medium strength)
for (m in names(modules_v3)) {
  mod_tfs <- modules_v3[[m]]$tfs
  tf_pos <- nodes_clean %>% filter(node_id %in% mod_tfs)
  if (nrow(tf_pos)==0) next
  cx <- mean(tf_pos$x, na.rm=TRUE); cy <- mean(tf_pos$y, na.rm=TRUE)

  # Pull TFs toward center
  for (tf in mod_tfs) {
    idx <- which(nodes_clean$node_id==tf)
    if (length(idx)!=1) next
    nodes_clean$x[idx] <- cx + 0.55*(nodes_clean$x[idx]-cx)
    nodes_clean$y[idx] <- cy + 0.55*(nodes_clean$y[idx]-cy)
  }

  # Push TGs outward
  mod_tgs <- all_targets$targetgene[all_targets$module==m]
  for (tg in mod_tgs) {
    idx <- which(nodes_clean$node_id==tg)
    if (length(idx)!=1) next
    nodes_clean$x[idx] <- cx + 1.20*(nodes_clean$x[idx]-cx)
    nodes_clean$y[idx] <- cy + 1.20*(nodes_clean$y[idx]-cy)
  }
}

# Second: gentle local TF repulsion within same module to avoid label overlap
# TF labels are ~3.2pt bold text. In coordinate space, this is roughly 0.2-0.3 units.
# Minimum TF-to-TF distance should be ~0.45 units.

for (iter in 1:5) {
  for (m in names(modules_v3)) {
    mod_tfs <- modules_v3[[m]]$tfs
    if (length(mod_tfs)<2) next
    tf_idxs <- which(nodes_clean$node_id %in% mod_tfs)
    for (i in tf_idxs) {
      for (j in tf_idxs) {
        if (i>=j) next
        dx <- nodes_clean$x[i]-nodes_clean$x[j]
        dy <- nodes_clean$y[i]-nodes_clean$y[j]
        d <- sqrt(dx^2+dy^2)
        if (d<1e-6) { d<-0.01; dx<-0.01; dy<-0 }
        if (d < 0.5) {
          push <- (0.5-d)/2 * 0.8
          ux <- dx/d; uy <- dy/d
          nodes_clean$x[i] <- nodes_clean$x[i] + ux*push
          nodes_clean$y[i] <- nodes_clean$y[i] + uy*push
          nodes_clean$x[j] <- nodes_clean$x[j] - ux*push
          nodes_clean$y[j] <- nodes_clean$y[j] - uy*push
        }
      }
    }
  }
}

cat("TF overlap reduced\n")

# ============================================================================
# 5. EDGE SHORTENING
# ============================================================================

cat("Computing shortened edges...\n")

compute_edges <- function(nodes_df) {
  edges_tg <- edges_tf_target %>%
    left_join(nodes_df%>%select(node_id,x_from=x,y_from=y), by=c("from"="node_id")) %>%
    left_join(nodes_df%>%select(node_id,x_to=x,y_to=y),   by=c("to"="node_id")) %>%
    filter(!is.na(x_from) & !is.na(x_to)) %>%
    mutate(
      dx=x_to-x_from, dy=y_to-y_from, d=sqrt(dx^2+dy^2),
      ux=ifelse(d>1e-6,dx/d,0), uy=ifelse(d>1e-6,dy/d,0),
      x_start=x_from+ux*tf_radius, y_start=y_from+uy*tf_radius,
      x_end=x_to-ux*tg_radius,     y_end=y_to-uy*tg_radius
    ) %>% select(-dx,-dy,-d,-ux,-uy)

  edges_tt <- edges_tf_tf_display %>% select(from,to,n_regions) %>%
    left_join(nodes_df%>%select(node_id,x_from=x,y_from=y), by=c("from"="node_id")) %>%
    left_join(nodes_df%>%select(node_id,x_to=x,y_to=y),   by=c("to"="node_id")) %>%
    filter(!is.na(x_from) & !is.na(x_to)) %>%
    mutate(
      dx=x_to-x_from, dy=y_to-y_from, d=sqrt(dx^2+dy^2),
      ux=ifelse(d>1e-6,dx/d,0), uy=ifelse(d>1e-6,dy/d,0),
      x_start=x_from+ux*tf_radius, y_start=y_from+uy*tf_radius,
      x_end=x_to-ux*tf_radius,     y_end=y_to-uy*tf_radius
    ) %>% select(-dx,-dy,-d,-ux,-uy)

  list(tg=edges_tg, tt=edges_tt)
}

edges_list <- compute_edges(nodes_clean)
etg <- edges_list$tg; ett <- edges_list$tt

# ============================================================================
# 6. MODULE ELLIPSES
# ============================================================================

hull_df <- nodes_clean %>% filter(!is.na(module)) %>%
  group_by(module) %>%
  summarise(cx=mean(x,na.rm=TRUE), cy=mean(y,na.rm=TRUE),
            sx=sd(x,na.rm=TRUE)*1.6+0.5, sy=sd(y,na.rm=TRUE)*1.6+0.5, .groups="drop")

hull_ellipse <- do.call(rbind, lapply(1:nrow(hull_df), function(i) {
  theta <- seq(0,2*pi,length.out=60)
  data.frame(module=hull_df$module[i],
             x=hull_df$cx[i]+hull_df$sx[i]*cos(theta),
             y=hull_df$cy[i]+hull_df$sy[i]*sin(theta))
}))

# Module titles at top of each ellipse
mod_titles <- hull_df %>% mutate(x=cx, y=cy+sy+0.35) %>%
  left_join(
    bind_rows(lapply(names(modules_v3), function(m)
      data.frame(module=m, short=modules_v3[[m]]$short, stringsAsFactors=FALSE))),
    by="module"
  )

# ============================================================================
# 7. ALL-LABEL PLOT (every node name visible)
# ============================================================================

cat("Building final plot (ALL labels)...\n")

# Split TG labels: those inside module ellipse get ggrepel, rest use geom_text
# Actually ggrepel with larger canvas handles this

p_final <- ggplot() +
  geom_polygon(data=hull_ellipse, aes(x=x,y=y,fill=module,group=module),
               alpha=0.06, color=NA) +
  scale_fill_manual(values=module_colors, guide="none") +

  ggnewscale::new_scale_fill() +
  geom_segment(data=etg, aes(x=x_start,y=y_start,xend=x_end,yend=y_end),
               color="grey45", linewidth=0.35, alpha=0.5,
               arrow=arrow(length=unit(2,"pt"),type="closed")) +

  ggnewscale::new_scale_fill() +
  geom_curve(data=ett, aes(x=x_start,y=y_start,xend=x_end,yend=y_end),
             curvature=0.2, linewidth=0.7, color="grey15", alpha=0.7,
             arrow=arrow(length=unit(3,"pt"),type="closed")) +

  ggnewscale::new_scale_fill() +
  geom_point(data=nodes_clean%>%filter(node_type=="target_gene"),
             aes(x=x,y=y,fill=plot_log2FC), size=4.5, shape=21, stroke=0.6, color="grey60") +
  scale_fill_gradient2(
    name=expression(log[2]~"FC\n(clipped ±2)"),
    low=fc_colors[1], mid=fc_colors[2], high=fc_colors[3],
    midpoint=0, limits=c(-2,2), oob=squish, breaks=c(-2,-1,0,1,2),
    guide=guide_colorbar(order=2, barwidth=unit(10,"pt"), barheight=unit(70,"pt))
  ) +

  ggnewscale::new_scale_fill() +
  geom_point(data=nodes_clean%>%filter(node_type=="TF"),
             aes(x=x,y=y,size=node_size,fill=module),
             shape=21, stroke=0.4, color="grey55") +
  scale_fill_manual(values=module_colors, guide="none") +
  scale_size_continuous(range=c(8,12), guide="none") +

  # TF labels: bold, black, centered on node
  geom_text(data=nodes_clean%>%filter(node_type=="TF"),
            aes(x=x,y=y,label=node_label), size=3.2, fontface="bold", color="grey5") +

  # ALL target gene labels: repel to avoid overlap
  ggrepel::geom_text_repel(
    data=nodes_clean %>% filter(node_type=="target_gene"),
    aes(x=x,y=y,label=node_label),
    size=2.0, color="grey25", alpha=0.85,
    force=1.5, max.overlaps=100,
    segment.size=0.15, segment.color="grey70",
    box.padding=0.15, min.segment.length=0.05, fontface="italic"
  ) +

  geom_text(data=mod_titles, aes(x=x,y=y,label=short),
            size=3.0, fontface="bold", color="grey15", lineheight=0.9, hjust=0.5) +

  coord_fixed(clip="off") + theme_void() +
  theme(
    legend.position="right", legend.box="vertical",
    legend.spacing.y=unit(10,"pt"),
    legend.title=element_text(size=9,face="bold"),
    legend.text=element_text(size=8),
    legend.key.size=unit(14,"pt"),
    plot.margin=margin(15,30,15,25),
    plot.title=element_text(hjust=0.5,size=14,face="bold",margin=margin(b=6)),
    plot.subtitle=element_text(hjust=0.5,size=9.5,color="grey35")
  ) +
  labs(
    title="TF-family regulatory modules underlying PanNET functional reprogramming",
    subtitle="Kamada–Kawai layout -- TF-centered per module, all target genes labeled"
  )

# ============================================================================
# 8. SAVE
# ============================================================================

od <- "GRN_kk_final"
ggsave(file.path(od,"GRN_kk_final_all_labels.pdf"), p_final,
       width=18, height=17, device=cairo_pdf)
ggsave(file.path(od,"GRN_kk_final_all_labels.png"), p_final,
       width=18, height=17, dpi=300, bg="white")
ggsave(file.path(od,"GRN_kk_final_all_labels.svg"), p_final,
       width=18, height=17, device=svglite::svglite)

# Node/edge tables
nodes_out <- nodes_clean %>% select(
  node_id, node_type, node_label, module, log2FoldChange, plot_log2FC,
  node_size, weighted_out_degree, x, y
)
write.csv(nodes_out, file.path(od,"GRN_kk_final_nodes.csv"), row.names=FALSE)
write.csv(etg %>% select(from,to,n_regions), file.path(od,"GRN_kk_final_edges.csv"), row.names=FALSE)
write.csv(ett %>% select(from,to,n_regions), file.path(od,"GRN_kk_final_TF_TF_edges.csv"), row.names=FALSE)

# QC
cat("\n══════════ QC ══════════\n")
cat(sprintf("TFs:  %d\n", sum(nodes_clean$node_type=="TF")))
cat(sprintf("TGs:  %d\n", sum(nodes_clean$node_type=="target_gene")))
cat(sprintf("TF→TG edges: %d\n", nrow(etg)))
cat(sprintf("TF→TF edges: %d\n", nrow(ett)))
cat(sprintf("All TG labels: %d\n", sum(nodes_clean$node_type=="target_gene")))

# Check for TF overlap
tf_df <- nodes_clean %>% filter(node_type=="TF")
min_tf_dist <- Inf
for (i in 1:(nrow(tf_df)-1)) {
  for (j in (i+1):nrow(tf_df)) {
    d <- sqrt((tf_df$x[i]-tf_df$x[j])^2 + (tf_df$y[i]-tf_df$y[j])^2)
    if (d < min_tf_dist) min_tf_dist <- d
  }
}
cat(sprintf("Min TF-TF distance: %.3f\n", min_tf_dist))

cat(sprintf("\nFiles saved to %s/\n", od))
