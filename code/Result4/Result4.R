# Result 4. Transcriptional regulatory programs orchestrating immune remodeling ----
## 1. Figure 4A. TF activity in each NMF ----
library(tidyverse)
library(Seurat)
library(pheatmap)
library(viper)
library(decoupleR)
library(patchwork)
library(OmnipathR)

# NMF-group
TIME <- qread("TIME.qs")
nt_time <- subset(TIME, group == "NT")
pt_time <- subset(TIME, group == "PT")

## NT
net <- get_collectri(organism = 'human', split_complexes = FALSE)
write_csv(net, "Figure4_decoupleR_net.csv")
net <- read.csv("Figure4_decoupleR_net.csv")
# mat <- as.matrix(TIME@assays[["RNA"]]@layers[["data"]])
mat <- GetAssayData(nt_time, layer = "data") %>% as.matrix()
plan("multisession", workers = 40) 
acts <- run_ulm(mat, net, minsize = 1, .source = 'source', .target = 'target', .mor = 'mor')

nt_time[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', 
              names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(object = nt_time) <- "tfsulm"
nt_time <- ScaleData(nt_time)
nt_time@assays$tfsulm@data <- nt_time@assays$tfsulm@scale.data

n_tfs <- 25
df <- t(as.matrix(nt_time@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(nt_time)) %>% # Levels: NMF4 NMF2 NMF1 NMF3
  pivot_longer(cols = -cluster,
               names_to = "source",
               values_to = "score") %>% 
  group_by(cluster, source) %>% 
  summarise(mean = mean(score))

df1 <- df %>% group_by(cluster) %>% top_n(10, wt = mean)
df1$cluster <- factor(df1$cluster, levels = c("NMF1", "NMF2", "NMF3", "NMF4"))
df1 <- df1[order(df1$cluster), ]
df1$source <- factor(df1$source, levels = unique(df1$source))
top_acts_mat <- df %>%
  filter(source %in% unique(df1$source)) %>%
  pivot_wider(id_cols = 'cluster', 
              names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster')
top_acts_mat$NMF <- rownames(top_acts_mat)
top_acts_mat <- top_acts_mat[order(top_acts_mat$NMF), ]
top_acts_mat <- top_acts_mat %>% select(unique(as.character(df1$source)))

range(top_acts_mat[, -ncol(top_acts_mat)])
green2purple.less.white <- colorRampPalette(c("darkgreen","#5F9E5F","white","#A074B6","darkorchid4"))

pdf("Figure4A_NT_NMF_top10TF.pdf", height = 3, width = 8)
Heatmap(top_acts_mat[, -ncol(top_acts_mat)], 
        # col = circlize::colorRamp2(c(-0.6, 0, 1), c("#F9A80A", "white", "#03018D")),
        col = circlize::colorRamp2(seq(-1,1,length.out=21), green2purple.less.white(21)),
        cluster_rows = F,
        cluster_columns = F
)
dev.off()
qsave(nt_time, file = "Figure4_nt_time.qs")
write.csv(df, file = "Figure4A_tf_NT.csv")

## PT
net <- read.csv("Figure4_decoupleR_net.csv")
# mat <- as.matrix(TIME@assays[["RNA"]]@layers[["data"]])
mat <- GetAssayData(pt_time, layer = "data") %>% as.matrix()
plan("multisession", workers = 40) 
acts <- run_ulm(mat, net, minsize = 1, .source = 'source', .target = 'target', .mor = 'mor')

pt_time[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', 
              names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(object = pt_time) <- "tfsulm"
pt_time <- ScaleData(pt_time)
pt_time@assays$tfsulm@data <- pt_time@assays$tfsulm@scale.data

df <- t(as.matrix(pt_time@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(pt_time)) %>% # Levels: NMF4 NMF2 NMF1 NMF3
  pivot_longer(cols = -cluster, 
               names_to = "source", 
               values_to = "score") %>% 
  group_by(cluster, source) %>% 
  summarise(mean = mean(score))

df1 <- df %>% group_by(cluster) %>% top_n(10, wt = mean)
df1$cluster <- factor(df1$cluster, levels = c("NMF1", "NMF2", "NMF3", "NMF4"))
df1 <- df1[order(df1$cluster), ]
df1$source <- factor(df1$source, levels = unique(df1$source))

top_acts_mat <- df %>%
  filter(source %in% unique(df1$source)) %>%
  pivot_wider(id_cols = 'cluster', 
              names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster')
top_acts_mat$NMF <- rownames(top_acts_mat)
top_acts_mat <- top_acts_mat[order(top_acts_mat$NMF), ]
top_acts_mat <- top_acts_mat %>% select(unique(as.character(df1$source)))

pdf("Figure4A_PT_NMF_top10TF.pdf", height = 3, width = 8)
Heatmap(top_acts_mat[, -ncol(top_acts_mat)], 
        # col = circlize::colorRamp2(c(-0.6, 0, 1), c("#F9A80A", "white", "#03018D")),
        col = circlize::colorRamp2(seq(-1, 1, length.out = 21), green2purple.less.white(21)),
        cluster_rows = F,
        cluster_columns = F
)
dev.off()

qsave(pt_time, file = "Figure4_pt_time.qs")
write.csv(df, file = "Figure4A_tf_PT.csv")


## 2. Figure 4B. TFs in each NMF ----
dt <- read.csv("Figure4_decoupleR_net.csv")
exp_tf <- intersect(rownames(TIME), unique(dt$source))
exp_gene <- intersect(rownames(TIME), unique(dt$target))
dt <- dt[dt$source %in% exp_tf & dt$target %in% exp_gene, ]

#### NT
nt_tf <- read.csv("Figure4A_tf_NT.csv")[, -1]
# nt_tf <- nt_tf[nt_tf$source %in% exp_tf, ]
nt_tf <- nt_tf %>% group_by(source) %>% 
  filter(mean == max(mean) & mean > 0)
nt_tf <- nt_tf[order(nt_tf$cluster, -nt_tf$mean), ]
nt_tf_1 <- split(nt_tf$source, factor(nt_tf$cluster, levels = c("NMF1", "NMF2", "NMF3", "NMF4")))

#### PT
pt_tf <- read.csv("Figure4A_tf_PT.csv")[, -1]
# pt_tf <- pt_tf[pt_tf$source %in% exp_tf, ]
pt_tf <- pt_tf %>% group_by(source) %>% 
  filter(mean == max(mean) & mean > 0)
pt_tf <- pt_tf[order(pt_tf$cluster, -pt_tf$mean), ]
pt_tf_1 <- split(pt_tf$source, factor(pt_tf$cluster, levels = c("NMF1", "NMF2", "NMF3", "NMF4")))

####
nmf1_j <- intersect(nt_tf_1[["NMF1"]], pt_tf_1[["NMF1"]])
nmf2_j <- intersect(nt_tf_1[["NMF2"]], pt_tf_1[["NMF2"]])
nmf3_j <- intersect(nt_tf_1[["NMF3"]], pt_tf_1[["NMF3"]])
nmf4_j <- intersect(nt_tf_1[["NMF4"]], pt_tf_1[["NMF4"]])

####
nmf1_nt_s <- setdiff(nt_tf_1[["NMF1"]], nmf1_j)
nmf2_nt_s <- setdiff(nt_tf_1[["NMF2"]], nmf2_j)
nmf3_nt_s <- setdiff(nt_tf_1[["NMF3"]], nmf3_j)
nmf4_nt_s <- setdiff(nt_tf_1[["NMF4"]], nmf4_j)

####
nmf1_pt_s <- setdiff(pt_tf_1[["NMF1"]], nmf1_j)
nmf2_pt_s <- setdiff(pt_tf_1[["NMF2"]], nmf2_j)
nmf3_pt_s <- setdiff(pt_tf_1[["NMF3"]], nmf3_j)
nmf4_pt_s <- setdiff(pt_tf_1[["NMF4"]], nmf4_j)

tf <- data.frame(tf = c(nmf1_j, nmf2_j, nmf3_j, nmf4_j, nmf1_nt_s, nmf2_nt_s, nmf3_nt_s, nmf4_nt_s, nmf1_pt_s, nmf2_pt_s, nmf3_pt_s, nmf4_pt_s),
                 type = c(rep("NMF1_NPT", length(nmf1_j)), rep("NMF2_NPT", length(nmf2_j)), rep("NMF3_NPT", length(nmf3_j)), rep("NMF4_NPT", length(nmf4_j)),
                          rep("NMF1_NT", length(nmf1_nt_s)), rep("NMF2_NT", length(nmf2_nt_s)), rep("NMF3_NT", length(nmf3_nt_s)), rep("NMF4_NT", length(nmf4_nt_s)),
                          rep("NMF1_PT", length(nmf1_pt_s)), rep("NMF2_PT", length(nmf2_pt_s)), rep("NMF3_PT", length(nmf3_pt_s)), rep("NMF4_PT", length(nmf4_pt_s))
                 ))
write.csv(tf, file = "Figure4B_tf_type.csv")
dt_tf <- split(tf$tf, tf$type)

df <- tf %>%
  count(type, tf) %>%
  pivot_wider(
    names_from = tf,
    values_from = n,
    values_fill = 0
  )
df <- as.data.frame(df)
rownames(df) <- df$type
df <- df[, -1]
df <- df[, colSums(df) == 1]
df$sum <- rowSums(df)

####
library(ggvenn)
venn_list_1 <- list(NT_NMF1 = nt_tf_1[["NMF1"]], PT_NMF1 = pt_tf_1[["NMF1"]])
ggvenn(venn_list_1, 
       fill_color = c("#BEBADA","#FB8072")
)
ggsave("Figure4B_ggvenn_NMF1.pdf", width = 6, height = 4)

venn_list_2 <- list(NT_NMF2 = nt_tf_1[["NMF2"]], PT_NMF2 = pt_tf_1[["NMF2"]])
ggvenn(venn_list_2, 
       fill_color = c("#CECB4A","#B3DE69")
)
ggsave("Figure4B_ggvenn_NMF2.pdf", width = 6, height = 4)

venn_list_3 <- list(NT_NMF3 = nt_tf_1[["NMF3"]], PT_NMF3 = pt_tf_1[["NMF3"]])
ggvenn(venn_list_3, 
       fill_color = c("#FDB462","#80B1D3")
)
ggsave("Figure4B_ggvenn_NMF3.pdf", width = 6, height = 4)

venn_list_4 <- list(NT_NMF4 = nt_tf_1[["NMF4"]], PT_NMF4 = pt_tf_1[["NMF4"]])
ggvenn(venn_list_4, 
       fill_color = c("#FCCDE5","#8DD3C7")
)
ggsave("Figure4B_ggvenn_NMF4.pdf", width = 6, height = 4)


## 3. Figure 4C. Identification of core regulons ----
TIME <- qread("TIME.qs")
dt <- read.csv("Figure4_decoupleR_net.csv")
exp_tf <- intersect(rownames(TIME), unique(dt$source))
exp_gene <- intersect(rownames(TIME), unique(dt$target))
dt <- dt[dt$source %in% exp_tf, ]
dt <- dt[dt$target %in% exp_gene, ]
write.csv(dt, file = "Figure4C_decoupleR_net_exp.csv")

# ④
data <- read.csv("Figure4C_decoupleR_net_exp.csv")
data <- data[data$source %in% common, ]
length(unique(data$target))

# ①ImmReg database
dt <- fread("TF_pathway_all.txt", data.table = F)
dt <- dt[dt$`P Adjust` < 0.05 & abs(dt$Score) > 0.995, ]
length(unique(dt$`TF Symbol`))
dt <- dt[dt$`TF Symbol` %in% common, ]
dt <- dt %>% separate_rows(`Marker Gene`, sep = ",") %>%
  mutate(ENSG = trimws(`Marker Gene`))

library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = unlist(strsplit(dt$ENSG, ",")),
  mart = mart
)
colnames(mapping)[1] <- "ENSG"
mapping <- unique(mapping)
df_final <- inner_join(dt, mapping, by = "ENSG") 

# ②NMF-spoecific TFs
df <- read.csv("Figure4B_tf_type.csv")
df <- df[df$type %in% c("NMF1_NT", "NMF1_PT", "NMF2_NT", "NMF2_PT", "NMF3_NT", "NMF3_PT", "NMF4_NT", "NMF4_PT"), ]
length(unique(df$tf))

# ③intersection
common <- intersect(dt$`TF Symbol`, df$tf)
common_target <- intersect(unique(data$target), unique(df_final$external_gene_name))

# ⑤DEGs of NMFs
nmf_markers <- read.csv("Figure2E_markers_nmf.csv")
nmf_markers <- unique(nmf_markers$gene)
nmf_markers_target <- intersect(common_target, nmf_markers)

# ⑥core TFs
data1 <- data[data$target %in% nmf_markers_target, ]
length(unique(data1$source))
df_final1 <- df_final[df_final$external_gene_name %in% nmf_markers_target, ]
length(unique(df_final1$`TF Symbol`))
core_TFs <- intersect(unique(df_final1$`TF Symbol`), unique(data1$source))
df_final2 <- df_final1[df_final1$`TF Symbol` %in% rownames(avg_sM), ]
length(unique(df_final2$`TF Symbol`))
unique(df_final2$Cancer)

# ⑦75%
nt <- read.csv("Figure4A_tf_NT.csv")[, -c(1, 2)]
nt$type <- paste0("NT", "_", nt$cluster)
pt <- read.csv("Figure4A_tf_PT.csv")[, -c(1, 2)]
pt$type <- paste0("PT", "_", pt$cluster)

nt1 <- nt[nt$source %in% core_TFs, ]
colnames(nt1)[2:4] <- paste0("NT", "_", colnames(nt1)[2:4])
pt1 <- pt[pt$source %in% core_TFs,]
colnames(pt1)[2:4] <- paste0("PT", "_", colnames(pt1)[2:4])
identical(unique(nt1$source), unique(pt1$source))

npt <- cbind(nt1, pt1)
npt$df <- npt$PT_mean - npt$NT_mean
npt <- npt[, -5]

npt1 <- npt %>% filter(!(NT_mean < 0 & PT_mean < 0))
length(unique(npt1$NT_source))
length(unique(npt1$PT_source))

npt2 <- npt1[abs(npt1$df) > as.numeric(quantile(abs(npt1$df))[4]), ]
length(unique(npt2$NT_source))
length(unique(npt2$PT_source))

# ⑧core Regulon
data2 <- data1[data1$source %in% unique(npt2$NT_source), ]
write.csv(data2[, -1], file = "Figure4C_core_tf.csv")

data3 <- data2[, c(2, 3)]
data4 <- data.frame(source = unique(data3$source),
                    target = unique(data3$source))
data5 <- rbind(data3, data4)
regulon <- split(data5$target, data5$source)

# ⑨Regulon score
library(Seurat)
library(future)
plan(multisession, workers = 20) 
score <- AddModuleScore(TIME,
                        features = regulon,
                        ctrl = 100,
                        name = "Regulon")
colnames(score@meta.data)[112:196] <- names(regulon)
meta <- score@meta.data

nt_meta <- subset(meta, group == "NT")
pt_meta <- subset(meta, group == "PT")

mscore_data <- pt_meta %>% select(10, 112:196)
avg_sM <- aggregate(mscore_data[, 2:ncol(mscore_data)], list(mscore_data$cancer), mean)
rownames(avg_sM) <- avg_sM$Group.1
avg_sM <- data.frame(t(avg_sM[,-1]))
avg_sM[avg_sM < 0] <- 0
# write.csv(avg_sM, file = "Figure4C_regulon_score_cancer_PT.csv")

avg_sM_nt <- avg_sM
avg_sM_pt <- avg_sM
avg_sM_npt <- avg_sM_pt - avg_sM_nt
write.csv(avg_sM_npt, file = "Figure4C_regulon_score_cancer_NPT.csv")

#⑩HR analysis of Regulon
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

genelist <- regulon
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
                         
                         # GSVA
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
                         
                         # Pan-cancer meta
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
write.csv(final_pan, file = "Figure4C_regulon_HR.csv")


## 4. Figure 4D. Regulon survival analysis ----
HR <- read.csv("Figure4C_regulon_HR.csv")[, -1]
HR$celltype <- factor(HR$celltype, levels = order_TFs)
HR <- HR[order(HR$celltype), ]
HR <- HR[HR$cancer == "Pan cancer", ]
HR <- HR %>%
  mutate(
    significance = ifelse(significance == "" | is.na(significance), "none", significance)
  )

HR11 <- HR[HR$cancer == "Pan cancer", ]
HR12 <- HR11[HR11$significance != "none", ]
HR1 <- HR[HR$celltype %in% unique(HR12$celltype), ]
HR1$color[which(HR1$significance == "none")] <- "No"
HR1$significance <- ifelse(HR1$significance == "none", "", HR1$significance)
HR1$cancer <- factor(HR1$cancer, levels = unique(HR1$cancer))
HR1 <- HR1[order(HR1$HR, decreasing = T),]
HR1$celltype <- factor(HR1$celltype, levels = rev(colnames(score1)))

color_color <- c("Better survival" = "#0F7B9F",
                 "Worse survival" = "#C3423F",
                 "No" = "white")

ggplot(HR1, aes(x = cancer, y = celltype, fill = color)) +
  geom_tile(colour = "grey50") +
  scale_fill_manual(values = color_color) +
  geom_text(aes(label = significance), size = 3) +
  theme_test() +
  labs(
    x = "",
    y = ""
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = "none")
ggsave("Figure4E_Regulon_HR1.pdf", width = 9.39, height = 6.67)

score_nt <- read_csv("Figure4C_regulon_score_cancer_NT.csv")
score_nt <- as.data.frame(score_nt)
rownames(score_nt) <- score_nt$...1
score_nt <- score_nt[, -1]
score_nt <- t(score_nt)
score_nt <- as.data.frame(score_nt)
rownames(score_nt) <- paste0("NT_", rownames(score_nt))

score_pt <- read_csv("Figure4C_regulon_score_cancer_PT.csv")
score_pt <- as.data.frame(score_pt)
rownames(score_pt) <- score_pt$...1
score_pt <- score_pt[, -1]
score_pt <- t(score_pt)
score_pt <- as.data.frame(score_pt)
rownames(score_pt) <- paste0("PT_", rownames(score_pt))

score <- rbind(score_nt, score_pt)
score$new <- rownames(score)
score <- score %>% separate(col = new, into = c("Group", "cancer"), sep = "_")
score <- score %>% select(Group, cancer, everything())
score <- score %>% arrange(cancer)
score <- score[, -c(1, 2)]

score1 <- score
score1 <- score1[, colnames(score1) %in% unique(as.character(HR1$celltype))]
score1 <- score1[, colnames(score1) %in% unique(as.character(HR1$celltype))]
score1 <- score1 %>% select(rev(HR1$celltype))

# define color palette
mycolors <- c("#3b374c", "#44598e", "#64a0c0", "#7ec4b7", "#deebcd") 
mycolors <- c("#073f82", "#1b71b4", "#58a4cf", "#a2cbe3", "#f2f9fe") 
mycolors <- c("#eeecdf", "#becdd2", "#6f9ad1", "#44679f", "#3f4f71") 
mycolors <- c("#492952", "#82677e", "white", "#59829e", "#1e4668") 
mycolors <- c("#57121d", "#d56e5e", "#eaebea", "#5390b5", "#1f294e") 
mycol1 <- colorRampPalette(c("#06a7cd", "white", "#e74a32"), alpha = TRUE)(21)

# scale_to_minusone_to_one <- function(x) {
#   (x - min(x)) / (max(x) - min(x)) * 2 - 1
# }
# 
# df_scaled <- as.data.frame(lapply(df, scale_to_minusone_to_one))

bk = c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
pdf("Figure4E_Regulon_HR_1.pdf", height = 5.61, width = 9.39)
pheatmap(t(score1),
         scale = "row",
         border_color = "black",
         cluster_cols = F,
         angle_col = 45,
         # color = mycol1,
         color = colorRampPalette(rev(mycolors))(21),
         # color = c(colorRampPalette(colors = c("#f0cfac", "white"))(21),  
         #           colorRampPalette(colors = c("white", "#425a04"))(21)),  
         cluster_rows = F
)
dev.off()

annotation_row <- data.frame(
  Group = rep(c("NT", "PT"), 19)
)
row.names(annotation_row) <- rownames(score)
annotation_col <- data.frame(
  Survival_association = HR$color,
  significance = HR$significance
)
row.names(annotation_col) <- order_TFs
annotation_colors <- list(
  Group = c("NT" = "#8dcdd5", "PT" = "#e6846d"),
  Survival_association = c("Better survival" = "#0F7B9F",
                           "Worse survival" = "#C3423F"),
  significance = c("none" = "#C6DBEF",
                   "*" = "#9ECAE1",
                   "**" = "#6BAED6",
                   "***" = "#2171B5",
                   "****" = "#08519C")
)

library(RColorBrewer)
p <- pheatmap::pheatmap(score,
                        scale = "column",
                        border_color = "white",
                        cluster_cols = T,
                        cluster_rows = F,
                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pdf("Figure4D_Regulon.pdf", height = 7, width = 20)
pheatmap::pheatmap(score,
                   scale = "column",
                   border_color = "white",
                   cluster_cols = T,
                   cluster_rows = F,
                   annotation_row = annotation_row,
                   annotation_col = annotation_col,
                   annotation_colors = annotation_colors,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()

order_TFs <- p[["tree_col"]][["labels"]][p[["tree_col"]][["order"]]]


## TF-targets
dt <- read.csv("Figure4C_core_tf.csv")
df <- data.frame(table(dt$source))
df$Var1 <- factor(df$Var1, levels = order_TFs)

library(ggplot2)
ggplot(df, aes(x = Var1, y = Freq)) +
  geom_col() +
  cowplot::theme_cowplot() +
  theme(
    text = element_text(size = 7, family = "ArialMT"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(linetype = 1, color = "black", linewidth = 0.3),
    axis.ticks = element_line(linetype = 1, color = "black", linewidth = 0.3)
  ) 
ggsave("Figure4D_targets.pdf", height = 4, width = 17)


## 5. Figure 4E. TFs, genes, and immune pathways ----
HR <- HR[HR$significance != "none", ]

dt <- read.csv("Figure4C_core_tf.csv")
dt <- dt[dt$source %in% unique(HR$celltype), ]

df <- read_csv("Figure4B_tf_type.csv")
df <- df[df$tf %in% unique(HR$celltype), ]
df <- df %>% separate(type, into = c("NMF", "Group"), sep = "_")
df <- df[, -1]
colnames(df)[1] <- "source"
df1 <- df[, 1:2]

dtf <- inner_join(dt, df, by = "source")

data <- as.data.frame.matrix(table(df1$source, df1$NMF))
data$source <- rownames(data)

datadt <- inner_join(dtf, data, by = "source")

library(PieGlyph)
ggplot(data = datadt, aes(y = source, x = target, color = Group))+
  geom_pie_glyph(slices = c('NMF1', 'NMF2', 'NMF3', 'NMF4'))+
  cowplot::theme_cowplot() +
  theme(
    text = element_text(size = 7, family = "ArialMT"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(linetype = 1, color = "black", linewidth = 0.3),
    axis.ticks = element_line(linetype = 1, color = "black", linewidth = 0.3)
  ) +
  scale_fill_manual(values = c("#F9F4D3","#B9D9F2","#B7938F","#CFEAE4")) +
  scale_color_manual(values = c("red", "blue"))
ggsave("Figure4E_1_tf-target_pieglyph.pdf", height = 6, width = 8.5)


#②TF-pathway
dt <- fread("TF_pathway_all.txt", data.table = F)
dt <- dt[dt$`P Adjust` < 0.05 & abs(dt$Score) > 0.995, ]
dt <- dt[dt$`TF Symbol` %in% HR$celltype, ]
dt <- dt %>% select(`TF Symbol`, `Immune Pathway`) %>% distinct()
dt$Freq <- 1

library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)

df <- dt %>%
  make_long(`Immune Pathway`, `TF Symbol`)

ggplot(df, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = node,
               label = node)) +
  geom_sankey(type = "sankey") +
  geom_sankey_text(size = 3.2, 
                   color = "black")+ 
  theme(legend.position = 'none')

c4a_gui()
mycol <- c4a('rainbow_wh_rd', 40)
mycol2 <- sample(mycol, length(mycol))
p1 <- ggplot(df, aes(x = x,
                     next_x = next_x,
                     node = node,
                     next_node = next_node,
                     fill = node,
                     label = node)) +
  scale_fill_manual(values = mycol2) + 
  geom_sankey(flow.alpha = 0.5, 
              smooth = 8, 
              width = 0.12) + 
  geom_sankey_text(size = 3.2,
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')
p1
ggsave("Figure4E_2_sankey.pdf", height = 6, width = 6)









