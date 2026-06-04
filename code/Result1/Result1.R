# Result 1. A large-scale pan-cancer single-cell atlas ----
setwd("/data1/yang/test")
library(tidyverse)
library(ggplot2)
library(Seurat)


## 1. Figure1A. Design process ----
# The number of cells, patients and samples for each type of cancer
obj <- read_rds("all.rds")

# patient
dt <- as.data.frame(table(obj$cancer, obj$patient))
dt <- dt[dt$Freq > 0,]
dt <- as.data.frame(table(dt$Var1))
ggplot(dt, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_text(aes(label = Freq), vjust = -0.5) + 
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) +
  ylab("The number of patient")
ggsave("Figure1A_Num.patient.pdf", width = 20, height = 12, units = "cm")  

# sample
dt <- as.data.frame(table(obj$cancer, obj$sample))
dt <- dt[dt$Freq > 0,]
dt <- as.data.frame(table(dt$Var1))
ggplot(dt, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "#d8d0e7") +
  geom_text(aes(label = Freq), vjust = -0.5) + 
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) +
  ylab("The number of sample")
ggsave("Figure1A_Num.sample.pdf", width = 20, height = 12, units = "cm")  

# cell
dt <- as.data.frame(table(obj$cancer))
dt$Freq <- round(dt$Freq / 10000, digits = 2)
ggplot(dt, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "#eaeaa2") +
  geom_text(aes(label = Freq), vjust = -0.5) + 
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) +
  ylab("The number of cell (*10000)")
ggsave("Figure1A_Num.cell.pdf", width = 20, height = 12, units = "cm")  


## 2. Figure1B. umap plotting ----
# ALL
DimPlot(
  obj,
  reduction = "umap.harmony",
  label = F,
  cols = cell_colors,
  group.by = "celltype",
  raster = F
) + theme_void() + NoLegend() + labs(title = "")
ggsave("Figure1B_celltype_umap_void.png", width = 20, height = 20, units = "cm")

# NT
obj1 <- subset(obj, group == "NT")
DimPlot(
  obj1,
  reduction = "umap.harmony",
  label = F,
  cols = cell_colors,
  group.by = "celltype",
  raster = F
) + theme_void() + NoLegend() + labs(title = "")
ggsave("Figure1B_NT_celltype_umap_void.png", width = 20, height = 20, units = "cm")

# PT
obj2 <- subset(obj, group == "PT")
DimPlot(
  obj2,
  reduction = "umap.harmony",
  label = F,
  cols = cell_colors,
  group.by = "celltype",
  raster = F
) + theme_void() + NoLegend() + labs(title = "")
ggsave("Figure1B_PT_celltype_umap_void.png", width = 20, height = 20, units = "cm")


## 3. Figure1C. The proportion of cells between NT and PT groups ----
ggplot(data = obj@meta.data) +
  geom_bar(mapping = aes(x = group, fill = celltype, ), position = "fill", width = 0.75) +
  scale_fill_manual(values = cell_colors) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  labs(y = "Total cell proportion (%) \n", x="") +
  theme(axis.title.x = element_text(size = 15, colour = "black")) +
  theme(axis.text.y = element_text(size = 15, colour = "black")) +
  theme(axis.text.x = element_text(size = 10, colour = "black")) +
  scale_y_continuous(expand = c(0,0)) 
ggsave("Figure1C_NT-PT.pdf", width = 10, height = 10, units = "cm")


## 4. Figure1D. The proportion of cells across different cancer types in NT and PT groups ----
meta <- obj@meta.data
meta_NT <- meta[meta$group == "NT", ]
meta_PT <- meta[meta$group == "PT", ]

# NT
ggplot(data = meta_NT) +
  geom_bar(mapping = aes(x = cancer, fill = celltype, ), position = "fill", width = 0.75) +
  scale_fill_manual(values = cell_colors) + 
  coord_flip() +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  labs(y = "Total cell proportion (%) \n", x="") +
  theme(axis.title.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black")) +
  theme(axis.text.x = element_text(colour = "black")) +
  scale_y_continuous(expand = c(0,0)) 
ggsave("Figure1D_NT_cancer.pdf", width = 13, height = 10, units = "cm")

# PT
ggplot(data = meta_PT) +
  geom_bar(mapping = aes(x = cancer, fill = celltype, ), position = "fill", width = 0.75) +
  scale_fill_manual(values = cell_colors) + 
  coord_flip() +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  labs(y = "Total cell proportion (%) \n", x="") +
  theme(axis.title.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black")) +
  theme(axis.text.x = element_text(colour = "black")) +
  scale_y_continuous(expand = c(0,0)) 
ggsave("Figure1D_PT_cancer.pdf", width = 13, height = 10, units = "cm")


## 5. Figure1E. cell marker ----
mainmarkers <- c("CD3D","NKG7","CD79A","MS4A1","CD68","CD163","CD14","CD1C","LAMP3","CPA3","VWF","COL1A1","RERGL","ALB","EPCAM")
DotPlot(obj, features = mainmarkers, group.by = "celltype") +
  coord_flip() +
  scale_color_gradientn(colors = c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = NULL, y = "", x = "") +
  guides(colour = guide_colourbar(title = "avg.exp\n(scaled)"),
         size = guide_legend(title = "pct.exp"))
ggsave("Figure1E_cell_marker.pdf", width = 15, height = 12, units = "cm")


## 6. Figure1F. umap of cell marker ----
library(scCustomize)

all.marker <- read.csv("all.marker.csv")
df <- all.marker %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%
  arrange(cluster)

marker_umap <- function(obj, gene) {
  FeaturePlot_scCustom(seurat_object = obj, 
                       colors_use = viridis_plasma_dark_high, na_color = "lightgray",
                       features = gene) + theme_void() + NoLegend() + labs(title = "")
}
marker_umap(obj, "CD3D"); ggsave("Figure1F_CD3D.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "NKG7"); ggsave("Figure1F_NKG7.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "CD79A"); ggsave("Figure1F_CD79A.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "MS4A1"); ggsave("Figure1F_MS4A1.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "CD68"); ggsave("Figure1F_CD68.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "CD163"); ggsave("Figure1F_CD163.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "CD14"); ggsave("Figure1F_CD14.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "CD1C"); ggsave("Figure1F_CD1C.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "LAMP3"); ggsave("Figure1F_LAMP3.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "CPA3"); ggsave("Figure1F_CPA3.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "VWF"); ggsave("Figure1F_VWF.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "COL1A1"); ggsave("Figure1F_COL1A1.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "RERGL"); ggsave("Figure1F_RERGL.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "ALB"); ggsave("Figure1F_ALB.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "EPCAM"); ggsave("Figure1F_EPCAM.png", width = 10, height = 10, units = "cm")
marker_umap(obj, "ACAT2"); ggsave("Figure1F_ACAT2.png", width = 10, height = 10, units = "cm")


## 7. Figure1G. subclustering umap ----
# T/NK cell
load("T_5.0.RData")
DimPlot(
  obj,
  reduction = "umap.harmony",
  label = T,
  cols = tnk_colors,
  group.by = "RNA_snn_res.0.9",
  raster = F
)
ggsave("Figure1G_T_NK_clusters_umap.png", width = 21, height = 20, units = "cm")

# B cell
load("B_5.0.RData")
DimPlot(
  obj,
  reduction = "umap.harmony",
  label = T,
  cols = B_colorS,
  group.by = "RNA_snn_res.1.3",
  raster = F
)
ggsave("Figure1G_B_clusters_umap.png", width = 21, height = 20, units = "cm")

# Myeloid cell
load("Myeloid_5.0.RData")
DimPlot(
  obj,
  reduction = "umap.harmony",
  label = T,
  cols = Myeloid_colors,
  group.by = "RNA_snn_res.0.3",
  raster = F
)
ggsave("Figure1G_Myeloid_clusters_umap.png", width = 21, height = 20, units = "cm")

# Stromal cell
load("Stromal_5.0.RData")
DimPlot(
  obj,
  reduction = "umap.harmony",
  label = T,
  cols = Stromal_colors,
  group.by = "RNA_snn_res.0.5",
  raster = F
)
ggsave("Figure1G_Stromal_clusters_umap.png", width = 21, height = 20, units = "cm")


## 8. FigureS2A. umap plots ----
# patient
DimPlot(
  all,
  reduction = "umap.harmony",
  label = F,
  group.by = "patient",
  raster = F
) + theme_void() + NoLegend() + labs(title = "")
ggsave("FigureS2A_patient_umap.png", width = 20, height = 20, units = "cm")

# study
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
pie(rep(1,n), col = sample(col_vector, 29))
col_1 = sample(col_vector, 29)

DimPlot(
  all,
  reduction = "umap.harmony",
  label = F,
  cols = col_1,
  group.by = "study",
  raster = F
) + theme_void() + NoLegend() + labs(title = "")
ggsave("FigureS2A_study_umap.png", width = 20, height = 20, units = "cm")

# cancer
DimPlot(
  obj,
  reduction = "umap.harmony",
  label = F,
  cols = cancer.color,
  group.by = "cancer",
  raster = F
) + theme_void() + NoLegend() + labs(title = "")
ggsave("FigureS2A_cancer_umap.png", width = 20, height = 20, units = "cm")


## 9. FigureS2B. The leave-one-dataset/cancer-type-out analyses of major cell types ----
meta <- obj@meta.data
group_col  <- "group"
sample_col <- "sample"
study_col  <- "study"
cancer_col <- "cancer"
cell_col   <- "celltype"
normal_label <- "NT"
tumor_label  <- "PT"
major_cells <- c(
  "T cell", "B cell", "Macrophage", "Monocyte", "DC",
  "Mast cell", "Endothelial cell", "Fibroblast", "SMC", "non-TIME"
)
meta2 <- meta %>%
  filter(.data[[group_col]] %in% c(normal_label, tumor_label)) %>%
  filter(.data[[cell_col]] %in% major_cells)
sample_info <- meta2 %>%
  distinct(
    sample = .data[[sample_col]],
    group  = .data[[group_col]],
    study  = .data[[study_col]],
    cancer = .data[[cancer_col]]
  )

# Overall cell-level proportion
overall_prop_wide <- meta2 %>%
  count(
    group = .data[[group_col]],
    celltype = .data[[cell_col]],
    name = "cell_number"
  ) %>%
  group_by(group) %>%
  mutate(
    total_cells = sum(cell_number),
    proportion = cell_number / total_cells
  ) %>%
  ungroup() %>%
  select(group, celltype, cell_number, proportion) %>%
  pivot_wider(
    names_from = group,
    values_from = c(cell_number, proportion),
    values_fill = 0
  ) %>%
  mutate(
    log2FC_PT_vs_NT = log2((.data[[paste0("proportion_", tumor_label)]] + 1e-6) /
                             (.data[[paste0("proportion_", normal_label)]] + 1e-6))
  ) %>%
  arrange(desc(log2FC_PT_vs_NT))

# sample-level proportion
sample_prop <- meta2 %>%
  count(
    sample = .data[[sample_col]],
    celltype = .data[[cell_col]],
    name = "cell_number"
  ) %>%
  right_join(
    expand_grid(sample = sample_info$sample, celltype = major_cells),
    by = c("sample", "celltype")
  ) %>%
  mutate(cell_number = ifelse(is.na(cell_number), 0, cell_number)) %>%
  left_join(sample_info, by = "sample") %>%
  group_by(sample) %>%
  mutate(
    total_cells = sum(cell_number),
    proportion = cell_number / total_cells
  ) %>%
  ungroup()

sample_mean_prop_wide <- sample_prop %>%
  group_by(group, celltype) %>%
  summarise(
    n_sample = n_distinct(sample),
    mean_prop = mean(proportion, na.rm = T),
    median_prop = median(proportion, na.rm = T),
    sd_prop = sd(proportion, na.rm = T),
    se_prop = sd_prop / sqrt(n_sample),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(n_sample, mean_prop, median_prop, sd_prop, se_prop),
    values_fill = 0
  ) %>%
  mutate(
    log2FC_PT_vs_NT = log2((.data[[paste0("mean_prop_", tumor_label)]] + 1e-6) /
                             (.data[[paste0("mean_prop_", normal_label)]] + 1e-6))
  ) %>%
  arrange(desc(log2FC_PT_vs_NT))

# leave-one-out senesitivity analysis
calc_loo_effect <- function(remove_col) {
  remove_values <- sort(unique(meta2[[remove_col]]))
  map_dfr(remove_values, function(x) {
    dat <- meta2 %>%
      filter(.data[[remove_col]] != x)
    dat_sample_info <- dat %>%
      distinct(
        sample = .data[[sample_col]],
        group  = .data[[group_col]],
        study  = .data[[study_col]],
        cancer = .data[[cancer_col]]
      )
    dat_sample_prop <- dat %>%
      count(
        sample = .data[[sample_col]],
        celltype = .data[[cell_col]],
        name = "cell_number"
      ) %>%
      right_join(
        expand_grid(sample = dat_sample_info$sample, celltype = major_cells),
        by = c("sample", "celltype")
      ) %>%
      mutate(cell_number = ifelse(is.na(cell_number), 0, cell_number)) %>%
      left_join(dat_sample_info, by = "sample") %>%
      group_by(sample) %>%
      mutate(
        total_cells = sum(cell_number),
        proportion = cell_number / total_cells
      ) %>%
      ungroup()
    dat_sample_prop %>%
      group_by(group, celltype) %>%
      summarise(
        n_sample = n_distinct(sample),
        mean_prop = mean(proportion, na.rm = T),
        median_prop = median(proportion, na.rm = T),
        sd_prop = sd(proportion, na.rm = T),
        se_prop = sd_prop / sqrt(n_sample),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = group,
        values_from = c(n_sample, mean_prop, median_prop, sd_prop, se_prop),
        values_fill = 0
      ) %>%
      mutate(
        removed_type = remove_col,
        removed_value = x,
        log2FC_PT_vs_NT = log2((.data[[paste0("mean_prop_", tumor_label)]] + 1e-6) /
                                 (.data[[paste0("mean_prop_", normal_label)]] + 1e-6))
      )
  })
}
loo_study_effect  <- calc_loo_effect(study_col)
loo_cancer_effect <- calc_loo_effect(cancer_col)
loo_effect <- bind_rows(loo_study_effect, loo_cancer_effect) %>%
  left_join(
    sample_mean_prop_wide %>%
      select(celltype, full_log2FC_PT_vs_NT = log2FC_PT_vs_NT),
    by = "celltype"
  ) %>%
  mutate(
    same_direction = sign(log2FC_PT_vs_NT) == sign(full_log2FC_PT_vs_NT),
    delta_log2FC = log2FC_PT_vs_NT - full_log2FC_PT_vs_NT
  )
loo_summary <- loo_effect %>%
  group_by(removed_type, celltype) %>%
  summarise(
    full_log2FC_PT_vs_NT = unique(full_log2FC_PT_vs_NT),
    mean_loo_log2FC = mean(log2FC_PT_vs_NT, na.rm = T),
    median_loo_log2FC = median(log2FC_PT_vs_NT, na.rm = T),
    min_loo_log2FC = min(log2FC_PT_vs_NT, na.rm = T),
    max_loo_log2FC = max(log2FC_PT_vs_NT, na.rm = T),
    sd_loo_log2FC = sd(log2FC_PT_vs_NT, na.rm = T),
    max_abs_delta_log2FC = max(abs(delta_log2FC), na.rm = T),
    same_direction_rate = mean(same_direction, na.rm = T),
    n_leave_one = n(),
    robust_direction = ifelse(same_direction_rate >= 0.9, "Yes", "No"),
    .groups = "drop"
  ) %>%
  arrange(removed_type, desc(abs(full_log2FC_PT_vs_NT)))
write.csv(loo_effect, "leave_one_each_result_major_cell_with_nonTIME.csv", row.names = F)

# Stability summary plot：same-direction rate
p1 <- ggplot(loo_summary, aes(x = celltype, y = same_direction_rate)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0.9, linetype = 2, color = "grey50") +
  facet_wrap(~ removed_type, nrow = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 11),
  ) +
  labs(
    x = NULL,
    y = "Same-direction rate"
  )
p1
ggsave("FigureS2B_Cell_Same_direction.pdf", p1, width = 5.27, height = 3.78)


## 10. FigureS2C. The main cell proportions of each patient across cancer types ----
load("FigureS2C_all.rds_meta.data.RData")
dt <- a %>% select(cancer, patient, celltype)

ggplot(data = dt) +
  geom_bar(mapping = aes(x = patient, fill = celltype, ), position = "fill", width = 0.75) +
  scale_fill_manual(values = cell_colors) + 
  theme_classic() +
  labs(y = "Total cell proportion (%) \n", x="") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) 
ggsave("FigureS2C_patient.pdf", width = 20, height = 3.5)


## 11. FigureS2E. Similarity of transcriptomes of cell subpopulations ----
library(qs)
library(ape)
library(factoextra)

TMEHluster <- qread("TIME.qs")
Idents(TMEHluster) <- "minorcell"; 
PCA <- Embeddings(TMEHluster, reduction = "pca"); 
TME.PCA <- matrix(0,nrow=length(SelectClus),ncol=50); 
rownames(TME.PCA) <- SelectClus; 
colnames(TME.PCA) <- colnames(PCA);
for(i in 1:length(SelectClus)) { SelectCell <- WhichCells(TMEHluster, idents=SelectClus[i]); 
SelectMatrix <- PCA[match(SelectCell,rownames(PCA)),]; 
for(j in 1:ncol(SelectMatrix)) { TME.PCA[i,j] <- mean(SelectMatrix[,j]); } }
Dist.PCA <- dist(TME.PCA[match(SelectClus,rownames(TME.PCA)),c(1:30)], method = "euclidean"); 
HClus.PCA <- hclust(Dist.PCA, method = "ward.D"); 
Cluster <- integer(length(as.phylo(HClus.PCA)$tip.label));
names(Cluster) <- as.phylo(HClus.PCA)$tip.label
Cluster[which(is.na(match(names(Cluster),TNKclus))==F)] <- 1; 
Cluster[which(is.na(match(names(Cluster),Bclus))==F)] <- 2; 
Cluster[which(is.na(match(names(Cluster),Mclus))==F)] <- 3;
Cluster[which(is.na(match(names(Cluster),Sclus))==F)] <- 4;
Cluster[which(is.na(match(names(Cluster),Tclus))==F)] <- 5;

library(ggtree)
ggtree(HClus.PCA) + 
  geom_tiplab(hjust = 0, size = 3)
ggsave("FigureS2E.ggtree.pdf", height = 20, width = 10)


## 12. FigureS2F. The proportion of non-TIME cells in NT and PT groups ----
# cell proportion
library(ggalluvial)
obj <- readRDS("FigureS2_nonTIME_1.0.rds")
df_pie <- obj@meta.data %>%
  group_by(tissue) %>%
  mutate(percentage = number / sum(number))
ggplot(df_pie, aes(x = tissue, 
                   stratum = cell, 
                   alluvium = cell,
                   y = percentage,
                   fill = cell,
                   label = cell)) +
  scale_y_continuous(label = scales::percent_format(),
                     expand=c(0,0)) +
  scale_fill_manual(values = as.character(major_color),breaks = names(major_color)) +
  geom_flow() +
  geom_stratum(width = 0.6,colour = "grey40") +
  theme(axis.title=element_blank(),
        axis.text.x.bottom = element_text(angle = 0,hjust = 1,vjust = 1),
        legend.title=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(),
        axis.ticks=element_line(),
        legend.position = "none",
        axis.text = element_text(size = 12))
ggsave("FigureS2F_non-TIME_tissue_proportion.pdf", width = 3, height = 5)

# umap
DimPlot(
  obj,
  reduction = "umap.harmony",
  label = T,
  cols = Stromal_colors,
  group.by = "RNA_snn_res.0.5",
  raster = F
)
ggsave("Figure2F_non-TIME_umap.png", width = 21, height = 20, units = "cm")


## 13. FigureS1H. umap of copycat results ----
load("non-TIME_copycat.RData")
# cell types
DimPlot(
  obj,
  reduction = "umap.harmony",
  label = T,
  group.by = "copycat.predicted.id",
  raster = F
)
ggsave("Figure2H_copycat_cell_umap.png", width = 21, height = 20, units = "cm")

# cancer types
DimPlot(
  obj,
  reduction = "umap.harmony",
  label = T,
  group.by = "cancer",
  cols = cancer.color,
  raster = F
)
ggsave("Figure2H_copycat_cancer_umap.png", width = 21, height = 20, units = "cm")


## 14. FigureS3A. umap of cell marker ----
library(Nebulosa)
library(scCustomize)
library(Seurat)
library(readr)

# T/NK cell
load("T_5.0.RData")
FeaturePlot_scCustom(seurat_object = obj, 
                     colors_use = viridis_magma_dark_high, 
                     features = "CD3D") + NoAxes() + theme_void() + ggtitle(label = "") + NoLegend()
ggsave("FigureS3A_T_NK_CD3D.png", height = 5, width = 5)
FeaturePlot_scCustom(seurat_object = obj, 
                     colors_use = viridis_magma_dark_high, 
                     features = "NKG7") + NoAxes() + theme_void() + ggtitle(label = "") + NoLegend()
ggsave("FigureS3A_T_NK_NKG7.png", height = 5, width = 5)

# B cell
load("B_5.0.RData")
FeaturePlot_scCustom(seurat_object = obj, 
                     colors_use = viridis_magma_dark_high, 
                     features = "MS4A1") + NoAxes() + theme_void() + ggtitle(label = "") + NoLegend()
ggsave("FigureS3A_B_MS4A1.png", height = 5, width = 5)
FeaturePlot_scCustom(seurat_object = obj, 
                     colors_use = viridis_magma_dark_high, 
                     features = "JCHAIN") + NoAxes() + theme_void() + ggtitle(label = "") + NoLegend()
ggsave("FigureS3A_B_JCHAIN.png", height = 5, width = 5)

# Myeloid cell
load("Myeloid_5.0.RData")
FeaturePlot_scCustom(seurat_object = obj, 
                     colors_use = viridis_magma_dark_high, 
                     features = "CD68") + NoAxes() + theme_void() + ggtitle(label = "") + NoLegend()
ggsave("FigureS3A_Myeloid_CD68.png", height = 5, width = 5)

# Stromal cell
load("Stromal_5.0.RData")
FeaturePlot_scCustom(seurat_object = obj, 
                     colors_use = viridis_magma_dark_high, 
                     features = "VWF") + theme_void() + ggtitle(label = "") + NoLegend()
ggsave("FigureS3A_Stromal_VWF.png", height = 5, width = 5)
FeaturePlot_scCustom(seurat_object = obj,
                     colors_use = viridis_magma_dark_high, 
                     features = "COL1A1") + theme_void() + ggtitle(label = "") + NoLegend()
ggsave("FigureS3A_Stromal_COL1A1.png", height = 5, width = 5)


## 15. FigureS3B. dotplot of T/NK cell marker ----
mainmarkers <- c("IL7R","CCR7","LTB","NPC2","FOXP3","TNFRSF4","HSPA1B","TMEM66","NBEAL1","FYB",
                 "H3-3B","CD79A","MZB1","CD8A","GZMH","GZMB","IGHA1","TYMS","NBEAL1","GAS5","HSPH1",
                 "MT1E","CMC1","SELK","SNHG29","C12ORF57","C4ORF3","FGFBP2","TYROBP","GNLY","MTRNR2L12","FCGR3A","FCGR2A")

DotPlot(obj, features = unique(mainmarkers), group.by = "minorcell") +
  scale_color_gradientn(colors = c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right") +
  labs(title = NULL, y = "", x = "") +
  guides(colour = guide_colourbar(title = "avg.exp\n(scaled)"),
         size = guide_legend(title = "pct.exp"))
ggsave("FigureS3B_T_NK_marker.pdf", width = 11, height = 7)


## 16. FigureS3C. dotplot of B cell marker ----
mainmarkers <- c("YBX3","CD83","CD52","CAPG","GNB2L1","HSP90AB1","PTPRCAP","GLTSCR2","CD97","EEF1G",
                 "H4C3","NIBAN3","MT1E","IGHD","KIAA0226L","SDF2L1","HSPA1B","IGHG4","MZB1","TNFRSF17","IGKV1-5",
                 "XIST","IGHG4","IGHM","MT1E","RGS13","PCLAF","CDC20","FGB","MKI67","CD79A","HSPA1A","CD3D")

DotPlot(obj, features = unique(mainmarkers), group.by = "minorcell") +
  scale_color_gradientn(colors = c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right") +
  labs(title = NULL, y = "", x = "") +
  guides(colour = guide_colourbar(title = "avg.exp\n(scaled)"),
         size = guide_legend(title = "pct.exp"))
ggsave("FigureS3C_B_marker.pdf", width = 12.5, height = 7.5)


## 17. FigureS3D. dotplot of Stromal cell marker ----
mainmarkers <- c("CAVIN2","TM4SF1","VWA1","PLVAP","SEMA3G","CXCR4","ALB","CRHBP","TFF3","H2AZ1",
                 "FBLN1","MFAP5","COL1A1","CENPF","CFD","DIO2","KCNMA1","HIGD1B","COL4A1","H2AJ","MYH11","RAMP1","APOC1","ACTA2")

DotPlot(obj, features = unique(mainmarkers), group.by = "minorcell") +
  scale_color_gradientn(colors = c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right") +
  labs(title = NULL, y = "", x = "") +
  guides(colour = guide_colourbar(title = "avg.exp\n(scaled)"),
         size = guide_legend(title = "pct.exp"))
ggsave("FigureS3D_Stromal_marker.pdf", width = 8.75, height = 6.3)


## 18. FigureS3E. dotplot of Myeloid cell marker ----
mainmarkers <- c("FCN1","VCAN","CD14","ATP5MJ","STAB1","C1QC","FBP1","APOA2","APOC3","SLC40A1",
                 "FOLR2","IDO1","LAMP3","TOP2A","MKI67","LILRA4","GZMB","IRF4","FBP1","TPSAB1","TYMS","G0S2")

DotPlot(obj, features = unique(mainmarkers), group.by = "minorcell") +
  scale_color_gradientn(colors = c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right") +
  labs(title = NULL, y = "", x = "") +
  guides(colour = guide_colourbar(title = "avg.exp\n(scaled)"),
         size = guide_legend(title = "pct.exp"))
ggsave("FigureS3E_Myeloid_marker.pdf", width = 8.75, height = 4.5)
