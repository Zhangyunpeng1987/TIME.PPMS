# Result 3. Program-specific metabolic reprogramming in TIME ----
## collecting human metabolic pathways from KEGG 
library(KEGGREST)
organism.info <- keggList('organism') 
head(organism.info) 
listDatabases()
pathways <- keggList("pathway","hsa")
tryCatch({
  hsa_pathways <- keggLink("pathway", "hsa")
  print(hsa_pathways)
}, error = function(e) {
  message("keggLink error：", e$message)
})

meta_pathways <- unique(hsa_pathways)[grepl("hsa00", unique(hsa_pathways))]
meta_info <- lapply(meta_pathways, keggGet)

meta_info <- lapply(meta_pathways, function(pathway) {
  tryCatch({
    keggGet(pathway)
  }, error = function(e) {
    message("Error retrieving pathway: ", pathway)
    return(NULL)
  })
})
meta_info <- Filter(Negate(is.null), meta_info)

names <- unlist(lapply(meta_info, function(x) x[[1]]$NAME))
genes <- unlist(lapply(meta_info, function(x) {
  g <- x[[1]]$GENE
  paste(str_split(g[seq(2, length(g), by = 2)], ";", simplify = TRUE)[, 1], collapse = ";")
}))

class.id <- unlist(lapply(meta_info, function(x) x[[1]]$CLASS))

df <- data.frame(
  Pathway_ID = meta_pathways,
  Pathway_Name = names,
  Genes = genes,
  Class = class.id
)
df$Pathway_ID[1]
df$Pathway_ID <- gsub("path:", "", df$Pathway_ID)
df$Pathway_Name[1]
df$Pathway_Name <- gsub(" - Homo sapiens \\(human\\)", "", df$Pathway_Name)
df$Class[1]
df$Class <- gsub("Metabolism; ", "", df$Class)
df <- df %>% separate_rows(Genes, sep = ";", convert = TRUE)
write.csv(df, file = "Figure3A_KEGG_metabolism.csv")


## 1. Figure 3A. metabolic activity between NT and PT ----
# metabolic score
TIME <- read_rds("TIME.rds")
dt <- read.csv("Figure3A_KEGG_metabolism.csv")
dt <- dt[order(dt$Class), ]
meta_list <- split(dt$Genes, dt$Pathway_Name); rm(dt)
TIME <- AddModuleScore(TIME,
                       features = meta_list,
                       ctrl = 100,
                       name = "features")
colnames(TIME@meta.data)[26:109] <- names(meta_list); rm(meta_list)

## pheatmap
pt_TIME <- subset(TIME, group == "PT")
nt_TIME <- subset(TIME, group == "NT")

# NT
sce_Metal_exp <- nt_TIME
sce_Metal_exp$celltype <- sce_Metal_exp$NMF
mscore_data <- nt_TIME@meta.data %>% select(25:109)
avg_sM <- aggregate(mscore_data[, 2:ncol(mscore_data)], list(mscore_data$NMF), mean)
rownames(avg_sM) <- avg_sM$Group.1
avg_sM <- data.frame(t(avg_sM[,-1]))
avg_sM$KEGG <- rownames(avg_sM)
dt <- read.csv("Figure3A_KEGG_metabolism.csv")
dt <- dt %>% select(Pathway_Name, Class) %>% unique()
colnames(dt)[1] <- "KEGG"
avg_sM_dt <- inner_join(avg_sM, dt, by = "KEGG")
avg_sM_dt <- avg_sM_dt[order(avg_sM_dt$Class), ]
rownames(avg_sM_dt) <- avg_sM_dt$KEGG
avg_sM_dt[, 5] <- NULL

avg_sM_class <- aggregate(avg_sM_dt[, 1:ncol(avg_sM_dt)-1], list(avg_sM_dt$Class), mean)
rownames(avg_sM_class) <- avg_sM_class$Group.1
avg_sM_class[, 1] <- NULL

pdf("Figure3A_scMetabolism_nmf_nt.pdf", height = 15, width = 8)
pheatmap::pheatmap(avg_sM_dt[, 1:4],
                   scale = 'row',
                   border_color = "white",
                   cluster_cols = F,
                   cluster_rows = F,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()

pdf("Figure3A_class_nmf_nt.pdf", height = 5, width = 5)
pheatmap::pheatmap(avg_sM_class,
                   scale = 'row',
                   border_color = "white",
                   cluster_cols = F,
                   cluster_rows = F,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()

# PT
sce_Metal_exp <- pt_TIME
sce_Metal_exp$celltype <- sce_Metal_exp$NMF
mscore_data <- pt_TIME@meta.data %>% select(25:109)
avg_sM <- aggregate(mscore_data[, 2:ncol(mscore_data)], list(mscore_data$NMF), mean)
rownames(avg_sM) <- avg_sM$Group.1
avg_sM <- data.frame(t(avg_sM[,-1]))
avg_sM$KEGG <- rownames(avg_sM)
dt <- read.csv("Figure3A_KEGG_metabolism.csv")
dt <- dt %>% select(Pathway_Name, Class) %>% unique()
colnames(dt)[1] <- "KEGG"
avg_sM_dt <- inner_join(avg_sM, dt, by = "KEGG")
avg_sM_dt <- avg_sM_dt[order(avg_sM_dt$Class), ]
rownames(avg_sM_dt) <- avg_sM_dt$KEGG
avg_sM_dt[, 5] <- NULL

avg_sM_class <- aggregate(avg_sM_dt[, 1:ncol(avg_sM_dt)-1], list(avg_sM_dt$Class), mean)
rownames(avg_sM_class) <- avg_sM_class$Group.1
avg_sM_class[, 1] <- NULL

pdf("Figure3A_scMetabolism_nmf_pt.pdf", height = 15, width = 8)
pheatmap::pheatmap(avg_sM_dt[, 1:4],
                   scale = 'row',
                   border_color = "white",
                   cluster_cols = F,
                   cluster_rows = F,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()

pdf("Figure3A_class_nmf_pt.pdf", height = 5, width = 5)
pheatmap::pheatmap(avg_sM_class,
                   scale = 'row',
                   border_color = "white",
                   cluster_cols = F,
                   cluster_rows = F,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()

## metabolic score of more than nine cancer types
# NT
meta <- nt_TIME@meta.data %>% dplyr::select(10, 26:109)
meta <- meta %>%
  dplyr::group_by(cancer) %>%
  summarise(across(1:84, mean, na.rm = TRUE)) %>%
  mutate(across(2:85, ~ifelse(.x > mean(.x, na.rm = TRUE), 1, 0)))

nonzero_counts <- meta %>%
  select(2:85) %>%  
  summarise(across(everything(), ~sum(.x == 1))) %>%  
  t() %>%  
  as.data.frame() %>%
  rownames_to_column("metabolism")
nonzero_counts$metabolism <- factor(nonzero_counts$metabolism, levels = rev(rownames(avg_sM_dt)))

ggplot(nonzero_counts, aes(x = V1, y = metabolism)) +
  geom_bar(stat = "identity", fill = "#8DCDD5") +
  theme_classic() +
  ylab("") + xlab("Number of cancer") +
  geom_vline(xintercept = c(3, 6, 9, 12), colour = 'black', lwd = 0.36, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12))
ggsave("Figure3A_Num.count_NT.pdf", width = 8, height = 15)

# PT
meta <- pt_TIME@meta.data %>% dplyr::select(10, 26:109)
meta <- meta %>%
  dplyr::group_by(cancer) %>%
  summarise(across(1:84, mean, na.rm = TRUE)) %>%
  mutate(across(2:85, ~ifelse(.x > mean(.x, na.rm = TRUE), 1, 0)))

nonzero_counts <- meta %>%
  select(2:85) %>%  
  summarise(across(everything(), ~sum(.x == 1))) %>%  
  t() %>%  
  as.data.frame() %>%
  rownames_to_column("metabolism")
nonzero_counts$metabolism <- factor(nonzero_counts$metabolism, levels = rev(rownames(avg_sM_dt)))

ggplot(nonzero_counts, aes(x = V1, y = metabolism)) +
  geom_bar(stat = "identity", fill = "#E6846d") +
  theme_classic() +
  ylab("") + xlab("Number of cancer") +
  geom_vline(xintercept = c(3, 6, 9, 12), colour = 'black', lwd = 0.36, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12))
ggsave("Figure3A_Num.count_PT.pdf", width = 8, height = 15)


## 2. Figure 3B. metabolic genes ----
dt <- read.csv("Figure3A_KEGG_metabolism.csv")
dt1 <- read.csv("Figure2E_markers_nt_pt.csv")
markers_nmf <- read.csv("Figure2E_markers_nmf.csv")
markers_nmf <- markers_nmf[markers_nmf$avg_log2FC > 0, ]
meta_gene <- unique(dt$Genes)
gene_sets <- list(
  KEGG_MG = meta_gene,
  NMF_M = unique(markers_nmf$gene)
  # Group_M = unique(dt1$gene)
)
common_genes <- Reduce(intersect, gene_sets)

library(UpSetR)
library(RColorBrewer)
pdf('Figure3A_upset_plot.pdf', height = 6, width = 8, onefile = F)
upset(fromList(gene_sets),
      nsets = 2, 
      order.by = "freq", 
      point.size = 5, 
      line.size = 1.3, 
      mainbar.y.label = "IntersectionSize", 
      sets.x.label = "", 
      mb.ratio = c(0.60, 0.40), 
      text.scale = c(2, 2, 2, 2, 2, 2),
      main.bar.color = "gray30", 
      queries = list(list(query = intersects,params = list("KEGG_MG", "NMF_M"), color = "red", active = T)) #给自己想要展示的特定组的交集设置颜色，可设置多个
)
dev.off()

## Venn diagram
options(stringsAsFactors = F)
library(pacman)
p_load(ggvenn,
       ggVennDiagram,
       UpSetR,data.table,
       stringr)
df <- markers_nmf[markers_nmf$gene %in% common_genes, ]
GE_list <- split(df$gene, 
                 factor(df$cluster, levels = paste0("NMF", 1:4)))
ggvenn(GE_list,
       fill_color = nmf_colors)
ggsave("Figure3B_ggvenn.pdf", width = 5, height = 5)


## 3. Figure 3C. core metabolic genes ----
common_genes <- intersect(gene_sets[[1]], gene_sets[[2]])
expr_all <- GetAssayData(TIME, layer = "counts")[common_genes, ]
expr <- matr.filter(as.matrix(expr_all))
expr <- t(expr) %>% as.data.frame()
sce_Metal_exp <- TIME[, colnames(TIME) %in% rownames(expr)]
expr$celltype <- sce_Metal_exp$NMF
avg_sM <- aggregate(expr[, 1:ncol(expr) - 1], list(expr$celltype), mean)
rownames(avg_sM) <- avg_sM$Group.1
avg_sM <- data.frame(t(avg_sM[,-1]))

pdf("Figure3C_scMetabolism_NMF_gene.pdf", height = 3, width = 6)
pheatmap::pheatmap(t(avg_sM),
                   scale = 'column',
                   border_color = "white",
                   cluster_cols = T,
                   cluster_rows = F,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()

## group
mainmarkers <- c("GGT5", "PTGS2", "CD38", "HPGD", "CHPF", "PLPP5", "SDS", "PLD4", "LIPA", "GATM", "GPX1",
                 "AKR1C3", "COX7A1", "PLCG2", "MGLL", "HYAL2", "ENTPD1", "FUCA1", "EZH2", "TYMP" , "HMOX1")
## Dotplot
DotPlot(TIME, features = mainmarkers, group.by = "group", cols = c("#424da7","#dd2b19")) +
  # coord_flip() +
  scale_color_gradientn(colors = c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2),
        panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="right") +
  labs(title = NULL, y = "", x = "") +
  guides(colour = guide_colourbar(title = "avg.exp\n(scaled)"),
         size = guide_legend(title = "pct.exp"))
ggsave("Figure3E_commongene_group.pdf", width = 15, height = 12, units = "cm")


## 4. Figure 3D. metabolic gene pheatmap ----
dt <- dt[dt$Genes %in% common_genes, ]
dt <- dt[, 3:5]
dt <- dt[order(dt$Class), ]
dt1 <- dt[, 1:2]
dt1$Genes <- factor(dt1$Genes, levels = c("GGT5", "PTGS2", "CD38", "HPGD",  "CHPF", "PLPP5", "SDS", "PLD4", "LIPA", "GATM", "GPX1",
                                          "AKR1C3", "COX7A1", "PLCG2", "MGLL", "HYAL2", "ENTPD1", "FUCA1", "EZH2", "TYMP" , "HMOX1"))
dt1$Pathway_Name <- factor(dt1$Pathway_Name, levels = unique(dt$Pathway_Name))                    

dt2 <- dt[, -1]
dt3 <- dt[, -2]

df1 <- dt1 %>%
  count(Pathway_Name, Genes) %>%
  pivot_wider(
    names_from = Genes,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("Pathway_Name")
df1$Pathway_Name <- rownames(df1)
df1_dt3 <- inner_join(df1, dt3, by = "Pathway_Name") %>% unique()
rownames(df1_dt3) <- df1_dt3$Pathway_Name
df1_dt3 <- df1_dt3[, -c(22:23)]
df1_dt3 <- df1_dt3 %>% select("GGT5", "PTGS2", "CD38", "HPGD",  "CHPF", "PLPP5", "SDS", "PLD4", "LIPA", "GATM", "GPX1",
                              "AKR1C3", "COX7A1", "PLCG2", "MGLL", "HYAL2", "ENTPD1", "FUCA1", "EZH2", "TYMP" , "HMOX1")
pdf("Figure3D_GENE_META.pdf", width = 8, height = 6)
pheatmap(
  df1_dt3,
  color = colorRampPalette(c("white", "#f8bda1"))(2),  
  number_color = "black",
  cluster_rows = F,
  cluster_cols = F
)
dev.off()

## PT vs. NT
# NT
nt_time <- subset(obj, group == "NT")
sce_Metal_exp <- nt_time
sce_Metal_exp$celltype <- sce_Metal_exp$cancer
mscore_data <- nt_time@meta.data %>% select(10, 27:109)
avg_sM <- aggregate(mscore_data[, 2:ncol(mscore_data)], list(mscore_data$cancer), mean)
rownames(avg_sM) <- avg_sM$Group.1
avg_sM <- data.frame(t(avg_sM[,-1]))
avg_sM$KEGG <- rownames(avg_sM)
avg_sM <- avg_sM[avg_sM$KEGG %in% unique(dt3$Pathway_Name), ]
avg_sM$KEGG <- factor(avg_sM$KEGG, levels = unique(dt3$Pathway_Name))
avg_sM <- avg_sM[order(avg_sM$KEGG), ]
avg_sM_nt <- avg_sM[, -20]

# PT
sce_Metal_exp <- pt_time
sce_Metal_exp$celltype <- sce_Metal_exp$cancer
mscore_data <- pt_time@meta.data %>% select(10, 27:109)
avg_sM <- aggregate(mscore_data[, 2:ncol(mscore_data)], list(mscore_data$cancer), mean)
rownames(avg_sM) <- avg_sM$Group.1
avg_sM <- data.frame(t(avg_sM[,-1]))
avg_sM$KEGG <- rownames(avg_sM)
avg_sM <- avg_sM[avg_sM$KEGG %in% unique(dt3$Pathway_Name), ]
avg_sM$KEGG <- factor(avg_sM$KEGG, levels = unique(dt3$Pathway_Name))
avg_sM <- avg_sM[order(avg_sM$KEGG), ]
avg_sM_pt <- avg_sM[, -20]

avg <- avg_sM_pt - avg_sM_nt
write.csv(avg, file = "Figure3D_cancer_df.exp.csv")

pdf("Figure3D_scMetabolism_cancer.pdf", height = 7, width = 10)
pheatmap::pheatmap(avg,
                   scale = 'row',
                   border_color = "white",
                   cluster_cols = F,
                   cluster_rows = F,
                   color = c(colorRampPalette(colors = c("#f0cfac", "white"))(21),  
                             colorRampPalette(colors = c("white", "#425a04"))(21)), 
                   # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
)
dev.off()


## 5. Figure 3E. risk metabolic gene ----
library(GSVA)
library(TCGAplot)
library(tidyverse)
library(survival)
library(meta)
library(doParallel) 

registerDoParallel(cores = 20)

tpm <- get_all_tpm()
meta <- get_all_meta()
cancers <- unique(meta$Cancer)

kegg_meta <- read_csv("Figure3A_KEGG_metabolism.csv")
genelist <- split(kegg_meta$Genes, kegg_meta$Pathway_ID)

final_pan <- foreach(geneset_name = names(genelist), 
                     .combine = bind_rows,
                     .packages = c("GSVA", "tidyverse", "survival", "meta", "TCGAplot"),
                     .export = c("tpm", "meta", "cancers")) %dopar% {
                       current_geneset <- genelist[[geneset_name]]
                       current_genelist <- list(current_geneset)
                       names(current_genelist) <- geneset_name
                       
                       cox_results <- list()
                       
                       for (cancer in cancers) {
                         exprSet <- subset(tpm, Group == "Tumor" & Cancer == cancer) %>% 
                           tibble::add_column(ID = stringr::str_sub(rownames(.), 1, 12), .before = "Cancer") %>% 
                           dplyr::filter(!duplicated(ID)) %>% 
                           tibble::remove_rownames() %>% 
                           tibble::column_to_rownames("ID") %>% 
                           dplyr::filter(rownames(.) %in% rownames(subset(meta, Cancer == cancer)))
                         
                         exprSet <- exprSet[, -(1:2)] %>% as.matrix() %>% t()
                         
                         gsvapar <- gsvaParam(exprData = exprSet, 
                                              geneSets = current_genelist, 
                                              kcdf = "Gaussian")
                         exprSet_gsva <- gsva(gsvapar)
                         
                         cl <- meta[colnames(exprSet_gsva), ]
                         cl$symbol <- exprSet_gsva[geneset_name, ]
                         
                         if (sum(!is.na(cl$time) & !is.na(cl$event)) > 0) {
                           m <- tryCatch(
                             coxph(Surv(time, event) ~ symbol + age, data = cl),
                             error = function(e) NULL
                           )
                           
                           if (!is.null(m)) {
                             beta <- coef(m)
                             se <- sqrt(diag(vcov(m)))
                             tmp <- round(cbind(
                               HR = exp(beta),
                               se = se,
                               lower = exp(beta - 1.96 * se),
                               upper = exp(beta + 1.96 * se),
                               p = 1 - pchisq((beta/se)^2, 1)
                             ), 3)
                             cox_results[[cancer]] <- tmp["symbol", ]
                           }
                         }
                       }
                       
                       if (length(cox_results) > 0) {
                         a <- do.call(rbind, cox_results) %>% 
                           as.data.frame() %>% 
                           rownames_to_column("cancer") %>% 
                           select(cancer, HR, se, lower, upper, p)
                         
                         # Pan-cancer meta-analysis
                         meta_res <- metagen(log(a$HR), a$se, sm = "HR")
                         pan_row <- data.frame(
                           cancer = "Pan cancer",
                           HR = exp(meta_res$TE.random),
                           lower = exp(meta_res$lower.random),
                           upper = exp(meta_res$upper.random),
                           p = meta_res$pval.random,
                           stringsAsFactors = FALSE
                         )
                         
                         a <- a %>% select(cancer, HR, lower, upper, p)
                         
                         pan <- bind_rows(a, pan_row) %>% 
                           mutate(
                             celltype = geneset_name,
                             significance = case_when(
                               p <= 0.0001 ~ "****",
                               p <= 0.001 ~ "***",
                               p <= 0.01 ~ "**",
                               p <= 0.05 ~ "*",
                               TRUE ~ ""
                             ),
                             color = case_when(
                               HR > 1 ~ "Worse survival",
                               HR < 1 ~ "Better survival",
                               TRUE ~ "Neutral"
                             )
                           )
                         pan
                       } else {
                         NULL
                       }
                     }

stopImplicitCluster()
write.csv(final_pan, file = "TCGAplot_84_metabolism_HR.csv")

final_result1 <- final_pan[final_pan$cancer == "Pan cancer", ]
final_result1 <- final_result1 %>% filter(p < 0.05)
table(final_result1$color) 
# Better survival  Worse survival 
# 22               18 
a <- kegg_meta[kegg_meta$Pathway_ID %in% unique(final_result1$celltype), ]
write.csv(a, file = "pancancer_Hightriskmetagenes.csv")
