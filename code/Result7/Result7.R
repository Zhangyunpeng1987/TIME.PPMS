# Result 7. PPMS captures immune characteristics ----
## 1. Figure 7A. paired tumor-normal analysis of PPMS genes ----
library(TCGAplot)
gs_pan_paired_boxplot(sig.genes,"PPMS")


## 2. Figure 7B. PPMS and NMF ----
library(Seurat)
dt <- read.csv("sig.genes.csv")
dt[which(dt$x == "H1.3"), ] <- "H1-3"
dt[which(dt$x == "H1.4"), ] <- "H1-4"
ppms <- list(ppms = dt$x)
TIME <- AddModuleScore(obj,
                       features = ppms,
                       ctrl = 100,
                       name = "PPMS.score")

meta <- TIME@meta.data %>% dplyr::select(22, 25, 12, 112)

library(ggplot2)
library(ggpubr) 
P <- ggplot(meta, aes(x = NMF, y = PPMS.score1, fill = NMF)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_classic() + 
  scale_fill_manual(values = c("#F9F4D3","#B9D9F2","#B7938F","#CFEAE4")) +
  geom_hline(yintercept = mean(meta$PPMS.score1), linetype = "dashed", color = "red") +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10, angle = 45, hjust = 1, vjust = 1)) +
  xlab("") 
P + stat_compare_means(method = "anova", 
                       aes(label = "p.format"), size = 5) 
ggsave("Figure7_PPMS.score-NMF.pdf", width = 5, height = 2.5)


## 3. Figure 7C. PPMS risk score across cancer types ----
library(tidyverse)
library(ggradar)

load("Figure6_model.TCGA-validation.RData")
dt <- exp_meta %>% select(Cancer) %>% rownames_to_column()
colnames(dt)[1] <- "ID"
dt <- inner_join(all.tmp, dt, by = "ID")
df <- dt %>% group_by(Cancer) %>% summarise(mean = mean(RS))
df1 <- as.data.frame(t(df))
df1 <- df1 %>% rownames_to_column()
colnames(df1) <- df1[1, ]
df1 <- df1[-1, ]
df1 <- df1 %>% rownames_to_column() %>% as.data.frame()
df1 <- df1[, -1]
df1[, -1] <- as.numeric(df1[, -1] )

range(as.numeric(df1[1, -1]))
ggradar(as_tibble(df1[, -1]),
        grid.min = 11.19393,
        grid.mid = 72.60243,
        grid.max = 134.01093, 
        group.point.size = 2,
        group.line.width = 1, 
        background.circle.colour = 'grey', 
        background.circle.transparency = 0, 
        legend.position = 'right', 
        legend.text.size = 12, 
        fill = TRUE, 
        fill.alpha = 0.3 
)
ggsave("riskscore_TCGA.pdf", width = 6, height = 6)


## 4. Figure 7D. PPMS gene KEGG ----
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

entrez_ids <- mget(sigs, org.Hs.egSYMBOL2EG, ifnotfound = NA)  
entrez_ids <- as.character(entrez_ids)           
gene_group_diff_matrix <- data.frame(sigs, entrez_id = entrez_ids)  
gene <- entrez_ids[entrez_ids != "NA"]           
kegg_enrich <- enrichKEGG(gene = gene,
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

KEGG <- as.data.frame(kegg_enrich)  
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x) paste(gene_group_diff_matrix$sigs[match(strsplit(x, "/")[[1]], as.character(gene_group_diff_matrix$entrez_id))], collapse = "/")))  # 将geneID转换为基因名称，并加入到KEGG结果中
write.table(KEGG, file = "genes_kegg_enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)  # 将KEGG结果写入名为kegg.tsv的文件中，以制表符分隔，不包含行名
barplot(kegg_enrich, drop = TRUE, showCategory = 15, label_format = 60) 


## ploting
library(dplyr)         
library(ggplot2)        
library(ggrepel)       
library(stringr)        
library(ggforce)        
library(tidyverse)
library(data.table)
dt <- fread("Sig.genes_kegg_enrichment.tsv")

dat <- dt %>%
  arrange(pvalue) %>%       
  head(15) %>%             
  arrange(category) %>%     
  mutate(Description = factor(Description, levels = rev(unique(Description)))) 

dat$GeneRatio <- sapply(strsplit(dat$GeneRatio, "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])  
})

dat <- dat %>% mutate(Description = factor(Description, levels = rev(unique(Description))))

custom_colors <- colorRampPalette(c("#eeeeee", "#ff5743"))(100)  
levelcolor <- c("#AEC7E8", "#9c93e5", "#008bd0","#00a08f","#f2b603") 

dat$Wrapped_Description <- str_wrap(dat$Description, width = 40)

ggplot(dat) +
  ggforce::geom_link(
    aes(
      x = 0,  
      y = Wrapped_Description,  
      xend = GeneRatio,  
      yend = Wrapped_Description,  
      color = category,  
      alpha = after_stat(index),  
      size = after_stat(index)  
    ),
    n = 500,  
    show.legend = c(color = TRUE, linewidth = FALSE)  
  ) +
  geom_point(
    aes(
      x = GeneRatio, 
      y = Wrapped_Description,
      fill = -log10(pvalue)  
    ),
    color = "black",  
    size = 6,         
    shape = 21,        
    show.legend = TRUE
  ) +
  scale_color_manual(values = levelcolor, name = "Category") +  
  scale_fill_gradientn(
    colors = custom_colors,  
    values = scales::rescale(seq(0, 6, length.out = 100)),  
    name = "-Log10(pvalue)"
  ) +
  guides(alpha = "none", size = "none") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = NA), 
    panel.grid = element_blank(),                
    panel.border = element_blank(),              
    axis.line = element_line(color = "black", linewidth = 0.75),  
    axis.text.y = element_text(color = "black", size = 12, lineheight = 0.8),  
    axis.text.x = element_text(color = "black", size = 12),  
    axis.title = element_text(size = 12),       
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 12),      
    plot.margin = margin(1, 1, 1, 2, "cm")     
  ) +
  ylab("") +
  xlab("GeneRatio") + 
  scale_y_discrete(limits = rev)
ggsave("keggresult.pdf", width = 10, height = 6.5)


## 5. Figure 7E-F. PPMS and immratio ----
gs_immucell_heatmap(geneset = sig.genes)
gs_immunescore_heatmap(geneset = sig.genes)


## 6. Figure 7G. High/Low-PPMS ----
library(tidyverse)
library(TCGAplot)

load("Figure6_all_mime1.RData")
res[["riskscore"]][["RSF"]][["Dataset1"]] -> testing
testing$RS_type <- ifelse(testing$RS >= median(testing$RS), 1, -1) # 1-High; -1-Low
testing$group <- "Testing"
res[["riskscore"]][["RSF"]][["Dataset2"]] -> validation
validation$RS_type <- ifelse(validation$RS >= median(validation$RS), 1, -1)
validation$group <- "Validation"
tcga <- rbind(testing, validation)
rm(list = ls()[which(ls() != "tcga")])

dt <- get_immu_ratio()
dt <- dt %>% rownames_to_column() %>% as.data.frame()
colnames(dt)[1] <- "ID"
dt$ID <- sub("-01A$", "", dt$ID)
j_names <- intersect(dt$ID, tcga$ID)
tcga <- tcga[tcga$ID %in% j_names, ]
dt <- dt[dt$ID %in% j_names, ]
dt_tmp <- inner_join(tcga, dt, by = "ID")
df <- dt_tmp[, -c(1,2,3,4,6)]
df1 <- df %>% pivot_longer(cols = -RS_type, 
                           names_to = "cell", 
                           values_to = "value" )
df1$RS_type <- ifelse(df1$RS_type == 1, "High", "Low")

p <- ggplot(df1, aes(x = cell, y = value, fill = RS_type)) +
  geom_boxplot(aes(fill = RS_type), alpha = 0.6, outliers = F) +
  scale_fill_manual(values = c("lightgreen", "skyblue"), name = NULL) +
  labs(x = NULL, y = "immune cell ratio",
       title = "Low PPMS vs High PPMS") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p + stat_compare_means(aes(group = RS_type), label = "p.signif", 
                       method = "t.test", label.y.npc = "top", hide.ns = T)


## 7. Figure 7H. PPMS risk score and expression of immune-related genes ----
# chemokine
Chemokine <- c("CCL1","CCL2","CCL3","CCL4","CCL5","CCL7","CCL8","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25","CCL26","CCL28","CX3CL1","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CXCL9","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL17")
dt <- exp_meta %>% 
  select(all_of(Chemokine)) %>%
  mutate(new = rowMeans(.)) %>%
  select(new) %>%
  rownames_to_column() %>%
  as.data.frame()
colnames(dt)[1] <- "ID"
dt_tmp <- inner_join(all.tmp, dt)

my_theme <- theme_bw() +
  theme(
    aspect.ratio        = 1,                                        
    legend.position     = "none",                                   
    axis.text           = element_text(size = 18, color = "black"), 
    axis.ticks          = element_line(color = "black"),            
    axis.ticks.length   = unit(0.2, "cm"),                          
    panel.spacing       = unit(0, "cm"),                            
    axis.title          = element_blank()                           
  )

blue_range <- function(n) {
  colorRampPalette(c("#00204d", "#00326f", "#31446b", "#4e576c",
                     "#666970", "#958f78", "#b0a473", "#e7d159", "#ffea46"))(n)
}

cor_coef  <- cor(dt_tmp$new, dt_tmp$RS, method = "pearson")  
cor_test  <- cor.test(dt_tmp$new, dt_tmp$RS, method = "pearson")

ggplot(dt_tmp, aes(x = RS, y = new)) +
  geom_density2d_filled(bins = 9) +                                                  
  scale_fill_manual(values = blue_range(9)) +                                         
  geom_smooth(method = "lm", se = FALSE, color = "#bebebe", linetype = "dashed") +    
  annotate("text",                                                                    
           x = 0.85, y = 0.35,   
           label = sprintf("R = %.3f\nP = %.2g", cor_coef, cor_test$p.value),
           size = 7,            
           colour = "white",
           hjust = 1, vjust = 1) +
  labs(x = "tran", y = "meta") +
  my_theme
ggsave("ICGs.pdf", width = 5, height = 5)

# immune-related genes in PT and NT
library(GSVA)
immuinhibitor <- c("ADORA2A","BTLA","CD160","CD244","CD274","CD96","CSF1R","CTLA4","HAVCR2","IDO1","IL10","IL10RB","KDR","KIR2DL1","KIR2DL3","LAG3","LGALS9","PDCD1","PDCD1LG2","TGFB1","TGFBR1","TIGIT","VTCN1")
immustimulator <- c("CD27","CD276","CD28","CD40","CD40LG","CD48","CD70","CD80","CD86","CXCL12","CXCR4","ENTPD1","HHLA2","ICOS","ICOSLG","IL2RA","IL6","IL6R","KLRC1","KLRK1","LTA","MICB","NT5E","PVR","RAET1E","TMIGD2","TNFRSF13B","TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF13","TNFSF13B","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF9","ULBP1")
ICGs <- c("CD274","CTLA4","HAVCR2","LAG3","PDCD1","PDCD1LG2","SIGLEC15","TIGIT")
Chemokine <- c("CCL1","CCL2","CCL3","CCL4","CCL5","CCL7","CCL8","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25","CCL26","CCL28","CX3CL1","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CXCL9","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL17")
receptor <- c("CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10", "CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","XCR1")
geneset <- list(immuinhibitor = immuinhibitor,
                immustimulator = immustimulator,
                ICGs = ICGs,
                Chemokine = Chemokine,
                receptor = receptor)

# TCGA 
tcga <- get_all_paired_tpm()
df = tcga %>% select(-c(Cancer, Group, ID)) %>% t()
gsvapar <- gsvaParam(exprData = df, geneSets = geneset, 
                     kcdf = "Gaussian")
df <- gsva(gsvapar)
df = t(df)
df1 = tcga %>% select(Group, ID)
identical(rownames(df), rownames(df1))
dat = cbind(df1, df)
ggplot(dat, aes(x = Group, y = immuinhibitor, fill = Group)) +
  geom_violin() +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outliers = F) +
  stat_compare_means(method = 'wilcox.test', label = 'p.format',
                     paired = F, label.x.npc = 'centre', label.y.npc = 'top') +
  labs(y = 'immuinhibitor', x = "Group") +
  guides(fill = 'none', color = 'none') +
  theme_test() +
  scale_fill_manual(values = c("#8dcdd5", "#e6846d")) 
ggsave("FigureS6_Group-immuinhibitor.pdf", width = 3.98, height = 4.55)


## 8. Figure 7I. PPMS and TMB/MSI ----
## TMB
tmb <- get_tmb()
tmb <- tmb %>% tibble::rownames_to_column() %>% as.data.frame()
colnames(tmb)[1] <- "ID"
tmb$ID <- sub("-01A$", "", tmb$ID)

dt <- all.tmp[all.tmp$ID %in% tmb$ID, ] %>% select(ID, RS)
tmb <- tmb[tmb$ID %in% dt$ID, ]
dt_tmb <- inner_join(dt, tmb, by = "ID")
p1 <- ggscatter(dt_tmb, x = "RS", y = "TMB",
                add = "reg.line", conf.int = TRUE, color = "#dbdba7", fill = "#dbdba7",
                add.params = list(             
                  color = "#c0c05a",              
                  fill = "lightgray"           
                ),
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Risk score", ylab = "TMB")
p1

dt_tmb$Group <- ifelse(dt_tmb$RS > median(dt_tmb$RS), "High", "Low")
dt_tmb$tmb_group <- ifelse(dt_tmb$TMB > median(dt_tmb$TMB), "High", "Low")

corT = cor.test(dt_tmb$RS, dt_tmb$TMB, method = 'pearson')
cor = corT$estimate 
pvalue = corT$p.value
ggplot(dt_tmb, aes(x = RS, y = TMB))+
  geom_point()+
  geom_smooth(method = 'lm') +
  labs(y = 'TMB')+
  theme_test() +
  ggtitle(label = "R = 0.0211, p = 0.054")
ggsave("FigureS6_RS-TMB_cor.pdf", width = 5, height = 5)

ggplot(dt_tmb, aes(x = tmb_group, y = RS)) +
  stat_boxplot(geom = 'errorbar', width = 0.3,
               position = position_dodge(0.75)
  ) +
  geom_boxplot(mapping = aes(fill = tmb_group), position = position_dodge(0.75), size = 0.8, outliers = F #outlier.shape = NA
  ) +
  stat_compare_means(method = 't.test', label = 'p.format',
                     paired = F,#vjust=0,hjust=0,
                     label.x.npc = 'centre', label.y.npc = 'top') +
  geom_jitter(mapping = aes(color = tmb_group, alpha = 0.5), width = 0.1, size = 2, height = 0) +
  labs(x = 'TMB', y = "RS") +
  guides(fill = 'none', color = 'none') +
  theme_test() +
  scale_fill_manual(values = c("#b1bfe2", "#e3ddee")) +
  scale_color_manual(values = c("#b1bfe2", "#e3ddee")) +
  ggtitle(label = "High(n=4181), Low(n=4193)")
ggsave("FigureS6_RS-TMB_t.test.pdf", width = 5, height = 5)


## MSI
msi <- get_msi()
msi <- msi %>% tibble::rownames_to_column() %>% as.data.frame()
colnames(msi)[1] <- "ID"
msi$ID <- sub("-01A$", "", msi$ID)

dt <- all.tmp[all.tmp$ID %in% msi$ID, ] %>% select(ID, RS)
msi <- msi[msi$ID %in% dt$ID, ]
dt_msi <- inner_join(dt, msi, by = "ID")
p2 <- ggscatter(dt_msi, x = "RS", y = "MSI",
                add = "reg.line", conf.int = TRUE, color = "#dbdba7", fill = "#dbdba7",
                add.params = list(             
                  color = "#c0c05a",              
                  fill = "lightgray"           
                ),
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Risk score", ylab = "MSI")
p2
p1 + p2
ggsave("FigureS6_RS-MSI-TMB_cor.pdf", width = 10, height = 5)

dt_msi$Group <- ifelse(dt_msi$RS > median(dt_msi$RS), "High", "Low")
dt_msi$msi_group <- ifelse(dt_msi$MSI > median(dt_msi$MSI), "High", "Low")

corT = cor.test(dt_msi$RS, dt_msi$MSI, method = 'pearson')
cor = corT$estimate 
pvalue = corT$p.value
ggplot(dt_msi, aes(x = RS, y = MSI)) +
  geom_point()+
  geom_smooth(method = 'lm') +
  labs(y = 'MSI')+
  theme_test() +
  ggtitle(label = "R = -0.00173, p = 0.873")
ggsave("FigureS6_RS-MSI_cor.pdf", width = 5, height = 5)

ggplot(dt_msi, aes(x = msi_group, y = RS))+
  stat_boxplot(geom = 'errorbar', width = 0.3,
               position = position_dodge(0.75)
  ) +
  geom_boxplot(mapping = aes(fill = msi_group), position = position_dodge(0.75), size = 0.8, outliers = F #outlier.shape = NA
  ) +
  stat_compare_means(method = 't.test', label = 'p.format',
                     paired = F,#vjust=0,hjust=0,
                     label.x.npc = 'centre', label.y.npc = 'top') +
  geom_jitter(mapping = aes(color = msi_group, alpha = 0.5), width = 0.1, size = 2, height = 0) +
  labs(x = 'MSI', y = "RS") +
  guides(fill = 'none', color = 'none') +
  theme_test() +
  scale_fill_manual(values = c("#b1bfe2", "#e3ddee")) +
  scale_color_manual(values = c("#b1bfe2", "#e3ddee")) +
  ggtitle(label = "High(n=4181), Low(n=4193)")
ggsave("FigureS6_RS-MSI_t.test.pdf", width = 5, height = 5)


## 9. Figure 7J. PPMS and immunotherapy response ----
library(Seurat)
library(qs)
library(tidyverse)
library(ggplot2)
library(ggpubr)
obj <- qread("ICI/ICI.qs")
meta <- obj@meta.data
meta$group <- ifelse(meta$Sig.genes_UCell > median(meta$Sig.genes_UCell), "High", "Low")
obj@meta.data$group <- paste0(obj$study, "_", obj$cancer_type)
ggplot(obj@meta.data, aes(x = Response, y = Sig.genes_UCell, fill = Response)) +
  geom_boxplot(
    alpha = 0.7, 
    outlier.shape = NA,  
    width = 0.6
  ) +
  facet_wrap(~ group, nrow = 2, scales = "free_y") +
  stat_compare_means(
    method = "t.test", 
    label = "p.signif", 
    label.x = 1.5,
    comparisons = list(c("NR", "R"))
  ) +
  scale_fill_manual(values = c("NR" = "#F8766D", "R" = "#00BFC4")) +
  scale_color_manual(values = c("NR" = "#F8766D", "R" = "#00BFC4")) +
  labs(
    y = "Ucell_score"
  ) +
  theme_test() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top",
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank()
  )
ggsave("Ucell_score_by_Cancer_and_Response.png", p, width = 12, height = 8, dpi = 300)


## 10. Figure 7K. PPMS and ICI based on bulk ----
coef <- read.csv("Sig.genes_coefficients.csv")[, -1]
dt <- read.csv("ICI/ICI.tpm.clinical.csv")
dt1 <- read.csv("ICI/ICI.tpm.count.csv")
rownames(dt1) <- dt1$X
dt1$X <- NULL

common_genes <- intersect(coef$Sig.genes, rownames(dt1))
expr_common <- dt1[common_genes, ]
coef <- coef[coef$Sig.genes %in% common_genes, ]

ppms_coef <- coef$Coefficients
names(ppms_coef) <- coef$Sig.genes

risk_scores <- data.frame(
  Sample = colnames(expr_common),
  Risk_Score = apply(expr_common, 2, function(x) sum(x * ppms_coef))
)

data1 <- dt[, 2:4]
colnames(risk_scores)[1] <- "patient"
risk_scores$patient <- gsub("-", ".", risk_scores$patient)
data1$patient <- gsub("-", ".", data1$patient)
data2 <- inner_join(risk_scores, data1, by = "patient")
data2$group <- ifelse(data2$Risk_Score > median(data2$Risk_Score), "High", "Low")
result <- data2 %>%
  group_by(cancer, group, Response) %>%
  tally() %>%
  pivot_wider(names_from = Response, values_from = n, values_fill = 0)

cancer_types <- unique(result$cancer)
p_values <- data.frame(cancer = character(), p_value = numeric())

for (cancer in cancer_types) {
  current_data <- result[result$cancer == cancer, ]
  current_data <- as.data.frame(current_data)
  contingency_table <- matrix(c(
    current_data[current_data$group == "High", "NR"],
    current_data[current_data$group == "High", "R"],
    current_data[current_data$group == "Low", "NR"],
    current_data[current_data$group == "Low", "R"]
  ), nrow = 2, byrow = TRUE)
  test_result <- chisq.test(contingency_table)
  p_values <- rbind(p_values, data.frame(
    cancer_type = cancer,
    p_value = test_result$p.value
  ))
}


library(reshape2)
# long_dt <- result %>% pivot_longer(cols = -c("group", "cancer"),
#                                    names_to = "Response",
#                                    values_to = "Num")
long_dt <- result %>%
  group_by(cancer) %>%
  mutate(total = NR + R) %>%
  mutate(NR_ratio = NR / total, R_ratio = R / total) %>%
  ungroup()
long_dt <- melt(long_dt, id.vars = c("cancer", "group", "NR", "R", "total"), 
                measure.vars = c("NR_ratio", "R_ratio"), 
                variable.name = "response_ratio", 
                value.name = "ratio")
long_dt$level <- paste0(long_dt$group, "_", long_dt$response_ratio)

ggplot(long_dt, aes(x = ratio, y = cancer, fill = level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  theme_test() +
  scale_fill_manual(values = c("High_NR_ratio" = "#BE8061",
                               "High_R_ratio" = "#A84C26",
                               "Low_NR_ratio" = "#C5C7C4",
                               "Low_R_ratio" = "#909191"))
ggsave("Figure7_bulk_NRR.pdf", width = 4, height = 6)


## 11. Figure S12A. PPMS genes in ST ----
library(SpaCET)
library(Seurat)
library(SingleCellExperiment)
visiumPath = file.path("ST_visium/LUAD/")
SpaCET_LUAD = create.SpaCET.object.10X(visiumPath = visiumPath)
SpaCET_LUAD = SpaCET.quality.control(SpaCET_LUAD)
SpaCET.visualize.spatialFeature(
  SpaCET_LUAD,
  spatialType = "QualityControl",
  spatialFeatures = c("UMI", "Gene"),
  imageBg = T
)
library(MUDAN)
SpaCET_LUAD = SpaCET.deconvolution(SpaCET_LUAD, cancerType = "LUAD", coreNo = 8)
myST = SpaCET.GeneSetScore(SpaCET_obj = SpaCET_LUAD, GeneSets = pt_NMF_gene.list)
myST@results$GeneSetScore[1:4, 1:6]
SpaCET.visualize.spatialFeature(
  myST,
  spatialType = "GeneSetScore",
  spatialFeatures = sig.genes
)
save(SpaCET_HCC_HCC5NR,file = "SpaCET_LUAD.RData")

pdf("LUAD.pdf",width = 12,height = 2.5)
SpaCET.visualize.spatialFeature(
  myST,
  spatialType = "GeneSetScore",
  spatialFeatures = sig.genes
)
dev.off()




