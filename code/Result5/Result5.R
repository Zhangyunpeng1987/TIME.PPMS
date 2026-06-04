# Result 5. Multicellular interactions within cellular programs ----
## create CellChat object
# Step1  object creation and inference
library(Seurat)
library(ggsignif)
library(ggpubr)
library(cowplot)
library(dplyr)
library(NMF)
library(ggalluvial)
library(CellChat)
library(pheatmap)
library(qs)
library(patchwork)
library(future)
library(ComplexHeatmap)

pbmc <- subset(TIME, NMF == "NMF3")
table(pbmc$minorcell)
for (group in unique(pbmc$group)) {
  pbmc_subset <- pbmc[,pbmc$group == group]
  pbmc_cellchat <- createCellChat(pbmc_subset[["RNA"]]$data)
  pbmc_meta <- data.frame(celltype = pbmc_subset$celltype, row.names =  Cells(pbmc_subset))
  pbmc_cellchat <- addMeta(pbmc_cellchat, meta = pbmc_meta, meta.name = "celltype")
  pbmc_cellchat <- setIdent(pbmc_cellchat, ident.use = "celltype")
  groupSize <- as.numeric(table(pbmc_cellchat@idents)) 
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- CellChatDB
  pbmc_cellchat@DB <- CellChatDB.use 
  pbmc_cellchat <- subsetData(pbmc_cellchat) 
  pbmc_cellchat <- identifyOverExpressedGenes(pbmc_cellchat) 
  pbmc_cellchat <- identifyOverExpressedInteractions(pbmc_cellchat)
  pbmc_cellchat <- computeCommunProb(pbmc_cellchat)
  pbmc_cellchat <- filterCommunication(pbmc_cellchat, min.cells = 10)
  pbmc_cellchat <- computeCommunProbPathway(pbmc_cellchat)
  pbmc_cellchat <- aggregateNet(pbmc_cellchat)
  pbmc_cellchat <- netAnalysis_computeCentrality(pbmc_cellchat, slot.name = "netP") 
  saveRDS(pbmc_cellchat, paste0(group, '_Cellchat.rds'))
}

# step2 comparison analysis
cellchat.NT <- readRDS("NMF3_NT_Cellchat.rds")
cellchat.PT <- readRDS("NMF3_PT_Cellchat.rds")
object.list <- list(NT = cellchat.NT, PT = cellchat.PT)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

select_colors <- minorcell_colors[minorcell_colors$cell %in% unique(cellchat@meta$celltype), ]
select_colors_1 <- select_colors$colors
names(select_colors_1) <- select_colors$cell
gg1 <- netVisual_heatmap(cellchat, color.use = select_colors_1)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)- diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # 控制不同数据集中的气泡大小
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Plasma_c29_IGHM", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Fb_c20_DIO2", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "functional") 
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = 'uwot')
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
rankSimilarity(cellchat, type = "functional")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

i = 1
# outgoing
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i + 1], width = 10, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i + 1], width = 10, height = 20, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# all
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i + 1], width = 10, height = 20, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in PT", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in PT", angle.x = 45, remove.isolate = T)
gg1 + gg2

# DEGs
pos.dataset = "PT"
features.name <- paste0("diff_genes_", pos.dataset)
options(future.globals.maxSize = 10000000000000000)
plan(multicore, workers = 40)   
cellchat <- identifyOverExpressedGenes(
  cellchat, 
  group.dataset = "datasets", 
  pos.dataset = pos.dataset, 
  features.name = features.name, 
  only.pos = FALSE, 
  thresh.pc = 0.1, 
  thresh.fc = 0.1, 
  thresh.p = 0.05
)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PT", ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "NT", ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NT", "PT")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)

# save cellchat object and DEGs
save(cellchat.NT, cellchat.PT, cellchat, net.down, net.up, object.list, file = "cellchat_comparisonAnalysis_NMF3_NT_vs_PT_all.RData")


## 1. Figure 5A. Number of interactions and interaction strength between NT and PT ----
dt <- data.frame(NMF = c("NMF1", "NMF1", "NMF2", "NMF2", "NMF3", "NMF3", "NMF4", "NMF4"),
                 group = rep(c("NT", "PT"), 4),
                 number = c(7440, 8436, 10855, 11014, 2917, 3153, 4091, 6127),
                 strength = c(203.959, 231.858, 330.714, 276.603, 117.766, 75.953, 167.195, 221.033))
group_colors <- c("NT" = "#8dcdd5", "PT" = "#e6846d")

library(ggplot2)
library(tidyr)
dt_long <- pivot_longer(dt, 
                        cols = c(number, strength),
                        names_to = "variable",
                        values_to = "value")
ggplot(dt_long, aes(x = NMF, y = value, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), 
           width = 0.7,
           color = "black",
           linewidth = 0.3) +
  geom_text(aes(label = ifelse(variable == "number", 
                               format(round(value, 0), nsmall = 0),  
                               sprintf("%.2f", value))), 
            position = position_dodge(width = 0.8),
            vjust = -0.5,           
            size = 3,               
            color = "black",        
            fontface = "bold") +    
  facet_wrap(~variable, 
             scales = "free_y",
             ncol = 1,
             strip.position = "left",
             labeller = as_labeller(c(number = "Number", strength = "Strength"))) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "NMF Cluster", 
       y = NULL,
       fill = "Group") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
    axis.line = element_line(linewidth = 0.5),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  )
ggsave("Figure5A_Num_Str.pdf", height = 5, width = 7)


## 2. Figure 5B/Figure S6A-B. Comparison of differential numbers of interactions between NT and PT ----
minorcell_colors <- read.csv("minorcell_colors.csv")[, -1]
select_colors <- minorcell_colors[minorcell_colors$cell %in% unique(cellchat@meta$celltype), ]
select_colors <- select_colors[order(select_colors$cell), ]
select_colors_1 <- select_colors$colors

# NMF1
load("cellchat_comparisonAnalysis_NMF1_NT_vs_PT.RData")
pdf("Figure5B_NMF1_heatmap_Num.pdf", height = 7, width = 7.5)
netVisual_heatmap(cellchat, color.use = select_colors_1)
dev.off()

pdf("FigureS5A_NMF1_heatmap_Str.pdf", height = 7, width = 7.5)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

# NMF2
load("cellchat_comparisonAnalysis_NMF2_NT_vs_PT.RData")
pdf("Figure5B_NMF2_heatmap_Num.pdf", height = 7, width = 7.5)
netVisual_heatmap(cellchat, color.use = select_colors_1)
dev.off()

pdf("FigureS5A_NMF2_heatmap_Str.pdf", height = 7, width = 7.5)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

# NMF3
load("cellchat_comparisonAnalysis_NMF3_NT_vs_PT_all.RData")
load("cellchat_comparisonAnalysis_NMF3_NT_vs_PT.RData")
pdf("Figure5B_NMF3_heatmap_Num.pdf", height = 7, width = 7.5)
netVisual_heatmap(cellchat, color.use = select_colors_1)
netVisual_heatmap(cellchat.NT, color.use = select_colors_1, color.heatmap = c("white","#2166ac"))
dev.off()

pdf("FigureS5A_NMF3_heatmap_Str.pdf", height = 7, width = 7.5)
netVisual_heatmap(cellchat, measure = "weight")
netVisual_heatmap(cellchat.NT, color.use = select_colors_1, color.heatmap = c("white","#2166ac"), measure = "weight")
dev.off()

# NMF4
load("cellchat_comparisonAnalysis_NMF4_NT_vs_PT.RData")
pdf("Figure5B_NMF4_heatmap_Num.pdf", height = 7, width = 7.5)
netVisual_heatmap(cellchat, color.use = select_colors_1)
dev.off()

pdf("FigureS5A_NMF4_heatmap_Str.pdf", height = 7, width = 7.5)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()


## 3. Figure 5C. Shared pathways in each group ----
group_colors <- c("NT" = "#8dcdd5",
                  "PT" = "#e6846d")

#### NMF1
load("cellchat_comparisonAnalysis_NMF1_NT_vs_PT.RData")
pdf("Figure5C_NMF1_pathway.pdf", height = 10, width = 6)
p1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = group_colors)
print(p1)
dev.off()

#### NMF2
load("cellchat_comparisonAnalysis_NMF2_NT_vs_PT.RData")
pdf("Figure5C_NMF2_pathway.pdf", height = 10, width = 6)
p2 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = group_colors)
print(p2)
dev.off()

#### NMF3
load("cellchat_comparisonAnalysis_NMF3_NT_vs_PT.RData")
pdf("Figure5C_NMF3_pathway.pdf", height = 10, width = 6)
p3 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = group_colors)
print(p3)
dev.off()

#### NMF4
load("cellchat_comparisonAnalysis_NMF4_NT_vs_PT.RData")
pdf("Figure5C_NMF4_pathway.pdf", height = 10, width = 6)
p4 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = group_colors)
print(p4)
dev.off()


#### identify specific pathways in each group
library(tidyverse)

## NMF1
dt1 <- p1[["data"]] 
dt1 <- dt1[dt1$pvalues < 0.05, ] 
# dt1 <- dt1 %>% group_by(name) %>% top_n(1, wt = contribution.scaled)
dt1 <- dt1 %>%
  group_by(name) %>%
  slice(which.max(contribution.scaled)) %>%
  ungroup()
table(dt1$group)
dt1 <- dt1[dt1$name != "PTPRM", ]
dt1$asterisk <- case_when(
  dt1$pvalues < 0.0001 ~ "****",
  dt1$pvalues < 0.001 ~ "***",
  dt1$pvalues < 0.01 ~ "**",
  dt1$pvalues < 0.05 ~ "*",
  TRUE ~ ""
)
dt1$new_name <- paste0(dt1$asterisk, " ", dt1$name)
dt1 <- dt1 %>% select(name, group, asterisk, new_name, pvalues) %>% mutate(NMF = "NMF1")

## NMF2
dt2 <- p2[["data"]] 
dt2 <- dt2[dt2$pvalues < 0.05, ] 
dt2 <- dt2 %>%
  group_by(name) %>%
  slice(which.max(contribution.scaled)) %>%
  ungroup()
table(dt2$group)
dt2$asterisk <- case_when(
  dt2$pvalues < 0.0001 ~ "****",
  dt2$pvalues < 0.001 ~ "***",
  dt2$pvalues < 0.01 ~ "**",
  dt2$pvalues < 0.05 ~ "*",
  TRUE ~ ""
)
dt2$new_name <- paste0(dt2$asterisk, " ", dt2$name)
dt2 <- dt2 %>% select(name, group, asterisk, new_name, pvalues) %>% mutate(NMF = "NMF2")

## NMF3
dt3 <- p3[["data"]] 
dt3 <- dt3[dt3$pvalues < 0.05, ] 
dt3 <- dt3 %>%
  group_by(name) %>%
  slice(which.max(contribution.scaled)) %>%
  ungroup()
table(dt3$group)
dt3$asterisk <- case_when(
  dt3$pvalues < 0.0001 ~ "****",
  dt3$pvalues < 0.001 ~ "***",
  dt3$pvalues < 0.01 ~ "**",
  dt3$pvalues < 0.05 ~ "*",
  TRUE ~ ""
)
dt3$new_name <- paste0(dt3$asterisk, " ", dt3$name)
dt3 <- dt3 %>% select(name, group, asterisk, new_name, pvalues) %>% mutate(NMF = "NMF3")

## NMF4
dt4 <- p4[["data"]] 
dt4 <- dt4[dt4$pvalues < 0.05, ] 
dt4 <- dt4 %>%
  group_by(name) %>%
  slice(which.max(contribution.scaled)) %>%
  ungroup()
table(dt4$group)
dt4$asterisk <- case_when(
  dt4$pvalues < 0.0001 ~ "****",
  dt4$pvalues < 0.001 ~ "***",
  dt4$pvalues < 0.01 ~ "**",
  dt4$pvalues < 0.05 ~ "*",
  TRUE ~ ""
)
dt4$new_name <- paste0(dt4$asterisk, " ", dt4$name)
dt4 <- dt4 %>% select(name, group, asterisk, new_name, pvalues) %>% mutate(NMF = "NMF4")

dt <- rbind(dt1, dt2, dt3, dt4)

library(UpSetR)
library(RColorBrewer)
dt$group_nmf <- paste0(dt$group, "_", dt$NMF)
dt_all <- dt
dt_all$name <- as.character(dt_all$name)

library(ComplexHeatmap)
tmp <- dt_all
lists <- lapply(split(tmp, tmp$group_nmf), function(x)x$name)
m <- make_comb_mat(lists)
cs = comb_size(m)
labels = c("NT_NMF1", "PT_NMF1", "NT_NMF2", "PT_NMF2","NT_NMF3", "PT_NMF3", "NT_NMF4", "PT_NMF4")

pdf("Figure5C_upset_all.pdf", width = 8, height = 4)
colors <- c("#BEBADA", "#FB8072", "#CECB4A", "#B3DE69", "#FDB462", "#80B1D3", "#FCCDE5", "#8DD3C7")
ht = draw(UpSet(m,  
                top_annotation = upset_top_annotation(m, 
                                                      ylim = c(0, max(cs)*1.1),
                                                      annotation_name_rot = 90,
                                                      annotation_name_side = "left",height = unit(4, "cm")),
                bg_col = rev(colors),
                set_order = labels))
od = column_order(ht)

decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2., "pt"), default.units = "native", just = "left", gp = gpar(fontsize = 10),rot = 90)
})
dev.off()
write.csv(tmp, file = "Figure5C_tmp.csv")

#
df <- dt_all %>%
  count(group_nmf, name) %>%
  pivot_wider(
    names_from = name,
    values_from = n,
    values_fill = 0
  )
df <- as.data.frame(df)
rownames(df) <- df$group_nmf
df <- df[, -1]
df <- df[, colSums(df) == 1]

df1 <- dt_all[dt_all$name %in% colnames(df), ]
write.csv(df1, file = "Figure5C_df1.csv")


## 4. Figure 5D/Figure S7A-B. Comparison of signaling roles of NMFs between NT and PT ----
minorcell_colors <- read.csv("minorcell_colors.csv")[, -1]
select_colors <- minorcell_colors[minorcell_colors$cell %in% unique(cellchat@meta$celltype), ]
select_colors <- select_colors[order(select_colors$cell), ]
select_colors_1 <- select_colors$colors

df <- read.csv("Figure5C_df1.csv")[, -1] 

#### NMF1
load("cellchat_comparisonAnalysis_NMF1_NT_vs_PT.RData")
df1 <- subset(df, group_nmf %in% c("NT_NMF1", "PT_NMF1"))

i = 1
# outgoing
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
pathway.union_1 <- unique(df1$name)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union_1, title = names(object.list)[i], width = 10.5, height = 6, color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing", signaling = pathway.union_1, title = names(object.list)[i + 1], width = 10.5, height = 6, color.use = select_colors_1)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_outgoing_NMF1.pdf
# incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union_1, title = names(object.list)[i], width = 10.5, height = 6, color.heatmap = "GnBu", color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming", signaling = pathway.union_1, title = names(object.list)[i + 1], width = 10.5, height = 6, color.heatmap = "GnBu", color.use = select_colors_1)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_incoming_NMF1.pdf
# all
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union_1, title = names(object.list)[i], width = 10.5, height = 6, color.heatmap = "OrRd", color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "all", signaling = pathway.union_1, title = names(object.list)[i + 1], width = 10.5, height = 6, color.heatmap = "OrRd", color.use = select_colors_1)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_all_NMF1.pdf


#### NMF2
load("cellchat_comparisonAnalysis_NMF2_NT_vs_PT.RData")
df2 <- subset(df, group_nmf %in% c("NT_NMF2", "PT_NMF2"))

i = 1
# outgoing
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
pathway.union_2 <- unique(df2$name)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union_2, title = names(object.list)[i], width = 13.5, height = 6, color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing", signaling = pathway.union_2, title = names(object.list)[i + 1], width = 13.5, height = 6, color.use = select_colors_1)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_outgoing_NMF2.pdf
# incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union_2, title = names(object.list)[i], width = 13.5, height = 6, color.heatmap = "GnBu", color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming", signaling = pathway.union_2, title = names(object.list)[i + 1], width = 13.5, height = 6, color.heatmap = "GnBu", color.use = select_colors_1)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_incoming_NMF2.pdf
# all
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union_2, title = names(object.list)[i], width = 13.5, height = 6, color.heatmap = "OrRd", color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "all", signaling = pathway.union_2, title = names(object.list)[i + 1], width = 13.5, height = 6, color.heatmap = "OrRd", color.use = select_colors_1)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_all_NMF2.pdf


#### NMF3
load("cellchat_comparisonAnalysis_NMF3_NT_vs_PT_all.RData")
df3 <- subset(df, group_nmf %in% c("NT_NMF3", "PT_NMF3"))

i = 1
# outgoing
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
pathway.union_3 <- unique(df3$name)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union_3, title = names(object.list)[i], width = 10, height = 4, color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing", signaling = pathway.union_3, title = names(object.list)[i + 1], width = 10, height = 4, color.use = select_colors_1[-14])
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_outgoing_NMF3.pdf
# incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union_3, title = names(object.list)[i], width = 10, height = 4, color.heatmap = "GnBu", color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming", signaling = pathway.union_3, title = names(object.list)[i + 1], width = 10, height = 4, color.heatmap = "GnBu", color.use = select_colors_1[-14])
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_incoming_NMF3.pdf
# all
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union_3, title = names(object.list)[i], width = 10, height = 4, color.heatmap = "OrRd", color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "all", signaling = pathway.union_3, title = names(object.list)[i + 1], width = 10, height = 4, color.heatmap = "OrRd", color.use = select_colors_1[-14])
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_all_NMF3.pdf


#### NMF4
load("cellchat_comparisonAnalysis_NMF4_NT_vs_PT.RData")
df4 <- subset(df, group_nmf %in% c("NT_NMF4", "PT_NMF4"))

i = 1
# outgoing
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
pathway.union_4 <- unique(df4$name)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union_4, title = names(object.list)[i], width = 11, height = 5, color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing", signaling = pathway.union_4, title = names(object.list)[i + 1], width = 11, height = 5, color.use = select_colors_1)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_outgoing_NMF4.pdf
# incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union_4, title = names(object.list)[i], width = 11, height = 5, color.heatmap = "GnBu", color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming", signaling = pathway.union_4, title = names(object.list)[i + 1], width = 11, height = 5, color.heatmap = "GnBu", color.use = select_colors_1)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_incoming_NMF4.pdf
# all
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union_4, title = names(object.list)[i], width = 11, height = 5, color.heatmap = "OrRd", color.use = select_colors_1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "all", signaling = pathway.union_4, title = names(object.list)[i + 1], width = 11, height = 5, color.heatmap = "OrRd", color.use = select_colors_1)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # Figure5D_all_NMF4.pdf


## 5. Figure 5E. specific ligand-receptors between NT and PT ----
df <- read.csv("Figure5C_df1.csv")[, -1] 
df1_nt <- subset(df, group_nmf %in% c("NT_NMF1"))
df1_pt <- subset(df, group_nmf %in% c("PT_NMF1"))
df2_nt <- subset(df, group_nmf %in% c("NT_NMF2"))
df2_pt <- subset(df, group_nmf %in% c("PT_NMF2"))
df3_nt <- subset(df, group_nmf %in% c("NT_NMF3"))
df3_pt <- subset(df, group_nmf %in% c("PT_NMF3"))
df4_nt <- subset(df, group_nmf %in% c("NT_NMF4"))
df4_pt <- subset(df, group_nmf %in% c("PT_NMF4"))

# outgoing
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
pathway.union_1_nt <- unique(df1_nt$name)
pathway.union_1_pt <- unique(df1_pt$name)
pathway.union_2_nt <- unique(df2_nt$name)
pathway.union_2_pt <- unique(df2_pt$name)
pathway.union_3_nt <- unique(df3_nt$name)
pathway.union_3_pt <- unique(df3_pt$name)
pathway.union_4_nt <- unique(df4_nt$name)
pathway.union_4_pt <- unique(df4_pt$name)

#### NMF1 
load("cellchat_comparisonAnalysis_NMF1_NT_vs_PT.RData")
net.up <- net.up[net.up$pathway_name %in% pathway.union_1_pt, ]
net.down <- net.down[net.down$pathway_name %in% pathway.union_1_nt, ]
net_nmf1 <- rbind(net.down, net.up)
net_nmf1$NMF <- "NMF1"
net_nmf1$LR <- paste0(net_nmf1$ligand, " -> ", net_nmf1$receptor)
net_nmf1$LR_cell <- paste0(net_nmf1$source, " -> ", net_nmf1$target)
net_nmf1$LR_cell_group <- paste0(net_nmf1$LR_cell, " (", net_nmf1$datasets, ")")

ggplot(net_nmf1, aes(x = LR_cell_group, y = LR, fill = prob)) +
  geom_tile() +
  theme_test() +
  scale_fill_continuous_tableau() +
  RotatedAxis()


#### NMF2 
load("cellchat_comparisonAnalysis_NMF2_NT_vs_PT.RData")
net.up <- net.up[net.up$pathway_name %in% pathway.union_2_pt, ]
net.down <- net.down[net.down$pathway_name %in% pathway.union_2_nt, ]
net_nmf2 <- rbind(net.down, net.up)
net_nmf2$NMF <- "NMF2"
net_nmf2$LR <- paste0(net_nmf2$ligand, " -> ", net_nmf2$receptor)
net_nmf2$LR_cell <- paste0(net_nmf2$source, " -> ", net_nmf2$target)
net_nmf2$LR_cell_group <- paste0(net_nmf2$LR_cell, " (", net_nmf2$datasets, ")")
net_nmf2 <- net_nmf2[order(net_nmf2$datasets), ]
net_nmf2$LR_cell_group <- factor(net_nmf2$LR_cell_group, levels = unique(net_nmf2$LR_cell_group))

ggplot(net_nmf2, aes(x = LR_cell_group, y = LR, fill = prob)) +
  geom_tile() +
  theme_test() +
  scale_fill_continuous_tableau() +
  RotatedAxis()


#### NMF3 
load("cellchat_comparisonAnalysis_NMF3_NT_vs_PT_all.RData")
net.up <- net.up[net.up$pathway_name %in% pathway.union_3_pt, ]
net.down <- net.down[net.down$pathway_name %in% pathway.union_3_nt, ]
net_nmf3 <- rbind(net.down, net.up)
net_nmf3$NMF <- "NMF3"
net_nmf3$LR <- paste0(net_nmf3$ligand, " -> ", net_nmf3$receptor)
net_nmf3$LR_cell <- paste0(net_nmf3$source, " -> ", net_nmf3$target)
net_nmf3$LR_cell_group <- paste0(net_nmf3$LR_cell, " (", net_nmf3$datasets, ")")
net_nmf3 <- net_nmf3[order(net_nmf3$datasets), ]
net_nmf3$LR_cell_group <- factor(net_nmf3$LR_cell_group, levels = unique(net_nmf3$LR_cell_group))

ggplot(net_nmf3, aes(x = LR_cell_group, y = LR, fill = prob)) +
  geom_tile() +
  theme_test() +
  scale_fill_continuous_tableau() +
  RotatedAxis()


#### NMF4 
load("cellchat_comparisonAnalysis_NMF4_NT_vs_PT.RData")
net.up <- net.up[net.up$pathway_name %in% pathway.union_4_pt, ]
net.down <- net.down[net.down$pathway_name %in% pathway.union_4_nt, ]
net_nmf4 <- rbind(net.down, net.up)
net_nmf4$NMF <- "NMF4"
net_nmf4$LR <- paste0(net_nmf4$ligand, " -> ", net_nmf4$receptor)
net_nmf4$LR_cell <- paste0(net_nmf4$source, " -> ", net_nmf4$target)
net_nmf4$LR_cell_group <- paste0(net_nmf4$LR_cell, " (", net_nmf4$datasets, ")")
net_nmf4 <- net_nmf4[order(net_nmf4$datasets), ]
net_nmf4$LR_cell_group <- factor(net_nmf4$LR_cell_group, levels = unique(net_nmf4$LR_cell_group))

ggplot(net_nmf4, aes(x = LR_cell_group, y = LR, fill = prob)) +
  geom_tile() +
  theme_test() +
  scale_fill_continuous_tableau() +
  RotatedAxis()

###
net_nmf <- rbind(net_nmf1, net_nmf2, net_nmf3, net_nmf4)
write.csv(net_nmf, file = "Figure5E_net_nmf.csv", row.names = F)

a <- net_nmf %>% select(pathway_name, NMF, datasets) %>% unique()
write.csv(a, file = "Figure5E_net_nmf_1.csv", row.names = F)


### core pathways associated with immunity
core_pathways <- c("PD-L1", "BTLA", "IL2", "SEMA4", "SELL", "CD46", "OCLN", "CDH1", "NOTCH", "CHEMERIN", "VTN", "CD70", "APRIL", "JAM", "CLDN")
aa <- net_nmf[net_nmf$pathway_name %in% core_pathways, ]
aa$LR <- factor(aa$LR, levels = unique(aa$LR))
aa$LR_cell_group <- factor(aa$LR_cell_group, levels = unique(aa$LR_cell_group))

ggplot(aa, aes(x = LR, y = LR_cell_group, fill = prob)) +
  geom_tile() +
  theme_test() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Figure5E_core_LR.pdf", height = 20, width = 10)


## NMF-LRs
dt <- read.csv("Figure5E_net_nmf.csv")
dt <- dt[order(dt$NMF), ]
dt$LR <- factor(dt$LR, levels = unique(dt$LR))
dt$NMF <- factor(dt$NMF, levels = c("NMF4", "NMF3", "NMF2", "NMF1"))
df <- dt %>% select(LR, NMF, prob) %>% group_by(NMF, LR) %>% summarise(prob = mean(prob))
df$LR <- factor(df$LR, levels = unique(df$LR))
df$NMF <- factor(df$NMF, levels = c("NMF1", "NMF2", "NMF3", "NMF4"))

ggplot(df, aes(x = LR, y = NMF, fill = prob)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#ADD8E6", "#4682B4", "#191970"), name = "Role") +
  theme_minimal() +
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  Seurat::RotatedAxis()


## 6. Figure 5G/Figure S7C. ligand-receptors co-expression by ST ----
rm(list=ls())
library(Seurat)
library(hdf5r)
library(ggplot2)
library(data.table)
library(dplyr)

# ST Seurat
ccRCC_H46_T <-Load10X_Spatial(
  data.dir ="ST_visium/ST_visium/H46_T/outs",   
  filename = "filtered_feature_bc_matrix.h5", 
  slice ="ccRCC_H46_T")   
ccRCC_H46_T$orig.ident <-"ccRCC_H46_T" 
ccRCC_H46_T
SpatialFeaturePlot(ccRCC_H46_T, features = "nCount_Spatial")

ccRCC_H46_T <- SCTransform(ccRCC_H46_T, assay = "Spatial", verbose = FALSE)
ccRCC_H46_T <- RunPCA(ccRCC_H46_T, assay = "SCT", verbose = FALSE) 

ccRCC_H46_T<- FindNeighbors(ccRCC_H46_T, reduction = "pca", dims = 1:10)
ccRCC_H46_T <- FindClusters(ccRCC_H46_T, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(ccRCC_H46_T, label = TRUE, label.size = 5)

ccRCC_H46_T <- RunUMAP(ccRCC_H46_T, reduction = "pca", dims = 1:10)
p2 <- DimPlot(ccRCC_H46_T, reduction = "umap", label = TRUE)
p1+p2

SpatialFeaturePlot(ccRCC_H46_T, features =c("SOX10","SOD2"))
saveRDS(ccRCC_H46_T,"ccRCC_H46_T.rds")

# ST co-expression
library(SpaGene)
library(tidyverse)

load("ST/myST.rdata")
# count <- SpaCET_obj@assays$Spatial$counts
count <- ccRCC_H46_T@assays$Spatial$counts
location <- ccRCC_H46_T@images[["ccRCC_H46_T"]]@boundaries[["centroids"]]@coords
location <- as.data.frame(location)
location <- location %>% mutate(imagerow = y, imagecol = x) %>% select(imagecol, imagerow)
rownames(location) <- ccRCC_H46_T@images[["ccRCC_H46_T"]]@boundaries[["centroids"]]@cells
new_df <- data.frame(
  ligand_gene_symbol = c("DLL4", "DLL4", "CD274", "CD70", "JAG1", "JAG1", "CDH1"),
  receptor_gene_symbol = c("NOTCH3", "NOTCH4", "PDCD1", "CD27", "NOTCH3", "NOTCH4", "KLRG1")
)
obj_lr <- SpaGene_LR(count, location, LRpair = new_df)

## Plot Ligand-receptor pair 
plotLR(count, location, LRpair = c(new_df$ligand_gene_symbol[1], new_df$receptor_gene_symbol[1]), alpha.min = 0.2, pt.size = 1)
plotLR(count, location, LRpair = c(new_df$ligand_gene_symbol[2], new_df$receptor_gene_symbol[2]), alpha.min = 0.2, pt.size = 1)
plotLR(count, location, LRpair = c(new_df$ligand_gene_symbol[3], new_df$receptor_gene_symbol[3]), alpha.min = 0.2, pt.size = 1)
plotLR(count, location, LRpair = c(new_df$ligand_gene_symbol[4], new_df$receptor_gene_symbol[4]), alpha.min = 0.2, pt.size = 1)
plotLR(count, location, LRpair = c(new_df$ligand_gene_symbol[5], new_df$receptor_gene_symbol[5]), alpha.min = 0.2, pt.size = 1)
plotLR(count, location, LRpair = c(new_df$ligand_gene_symbol[6], new_df$receptor_gene_symbol[6]), alpha.min = 0.2, pt.size = 1)
plotLR(count, location, LRpair = c(new_df$ligand_gene_symbol[7], new_df$receptor_gene_symbol[7]), alpha.min = 0.2, pt.size = 1)

pdf("CDH1_KLRG1_LR.pdf", width = 10.5, height = 5)
plotLR(count, location, LRpair = c(new_df$ligand_gene_symbol[7], new_df$receptor_gene_symbol[7]), alpha.min = 0.2, pt.size = 1) +
  theme(
    panel.background = element_rect(fill = "black", color = "black"),
    panel.grid = element_blank(),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white")
  )
dev.off()

pdf("SEMA3C_PLXNA2_LR.pdf", width = 10.5, height = 5)
plotLR(count, location, LRpair = c("SEMA3C", "PLXNA2"), alpha.min = 0.2, pt.size = 1) +
  theme(
    panel.background = element_rect(fill = "black", color = "black"),
    panel.grid = element_blank(),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white")
  )
dev.off()

pdf("SEMA3C_NRP1_LR.pdf", width = 10.5, height = 5)
plotLR(count, location, LRpair = c("SEMA3C", "NRP1"), alpha.min = 0.2, pt.size = 1) +
  theme(
    panel.background = element_rect(fill = "black", color = "black"),
    panel.grid = element_blank(),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white")
  )
dev.off()