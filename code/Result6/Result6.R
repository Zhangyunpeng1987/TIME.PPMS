# Result 6. Survival prognostication based on PPMS ----
## step 1. TCGA risk genes
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
genes <- colnames(tpm)[-c(1, 2)]

final_result <- foreach(gene = genes, 
                        .combine = bind_rows,
                        .packages = c("survival", "meta", "dplyr", "tibble", "stringr")) %dopar% {
                          
                          cox_results <- list()
                          for (cancer in cancers) {
                            exprSet <- subset(tpm, Group == "Tumor" & Cancer == cancer) %>%
                              tibble::add_column(ID = stringr::str_sub(rownames(.), 1, 12), .before = "Cancer") %>%
                              dplyr::filter(!duplicated(ID)) %>%
                              tibble::remove_rownames() %>%
                              tibble::column_to_rownames("ID") %>%
                              dplyr::filter(rownames(.) %in% rownames(subset(meta, Cancer == cancer)))
                            
                            exprSet <- exprSet[, -(1:2)]
                            exprSet <- as.matrix(t(exprSet))
                            cl_data <- meta[colnames(exprSet), ]
                            
                            if (!gene %in% rownames(exprSet)) next
                            cl_data$symbol <- exprSet[gene, ]
                            
                            if (sum(!is.na(cl_data$time) & !is.na(cl_data$event)) == 0) next
                            
                            m <- tryCatch(
                              coxph(Surv(time, event) ~ symbol + age, data = cl_data),
                              error = function(e) NULL
                            )
                            
                            if (is.null(m)) next
                            
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
                          
                          if (length(cox_results) > 0) {
                            a <- do.call(rbind, cox_results) %>% 
                              as.data.frame() %>% 
                              rownames_to_column("cancer") %>% 
                              select(cancer, HR, se, lower, upper, p) %>%
                              filter(
                                se > 0,
                                !is.na(HR),
                                HR > 0,              
                                is.finite(log(HR))
                              )
                            
                            if (nrow(a) >= 1) {
                              # Pan-cancer meta分析
                              meta_res <- tryCatch(
                                {
                                  metagen(
                                    log(a$HR),
                                    a$se,
                                    sm = "HR",
                                    method = "DL",       
                                    control = list(maxiter = 1000)
                                  )
                                },
                                error = function(e) {
                                  message("Meta-analysis failed: ", e$message)
                                  return(NULL)
                                }
                              )
                              if (!is.null(meta_res)) {
                                pan_row <- data.frame(
                                  cancer = "Pan cancer",
                                  HR = exp(meta_res$TE.random),
                                  lower = exp(meta_res$lower.random),
                                  upper = exp(meta_res$upper.random),
                                  p = meta_res$pval.random,
                                  stringsAsFactors = FALSE
                                )
                              }
                            } else {
                              message("No valid studies for meta-analysis.")
                              pan_row <- NULL
                            }
                            a <- a %>% select(cancer, HR, lower, upper, p)
                            bind_rows(a, pan_row) %>%
                              mutate(
                                gene = gene,
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
                          } else {
                            NULL
                          }
                        }
stopCluster(cl)
write.csv(final_result, file = "all_Hightriskgenes.csv")

final_result1 <- final_result[final_result$cancer == "Pan cancer", ]
final_result1 <- final_result1 %>% filter(p < 0.05)
table(final_result1$color) 
# Better survival  Worse survival 
# 3785            4626 
write.csv(final_result1, file = "pancancer_Hightriskgenes.csv")

## step 2. scRNA-seq risk genes
## step 3. candidate risk genes
library(tidyverse)
library(readr)

dir_for_data <- "Figure6/TCGA/"

# Figure2
fig2 <- read_csv(file.path(dir_for_data, "TCGAplot_103_HR.csv")) %>%
  filter(cancer == "Pan cancer", color == "Worse survival", p < 0.05)

dir_for_signatures <- file.path(dir_for_data, "01.all_subset_signature_in_TIME")
result_df <- list()

for (i_file in list.files(dir_for_signatures)) {
  markers_1 <- read_tsv(file.path(dir_for_signatures, i_file)) %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 20)
  markers_1$celltype <- str_extract(i_file, "^.*(?=_wilcoxon)")
  result_df[[i_file]] <- markers_1
}

genelist <- bind_rows(result_df) %>%
  filter(celltype %in% fig2$celltype, avg_log2FC > 1)
fig2_genes <- unique(genelist$gene)
dt1 <- genelist %>% select(gene, celltype)
colnames(dt1)[2] <- "class"
dt1$group <- "Minorcell"

# Figure3
fig3 <- read_csv(file.path(dir_for_data, "TCGAplot_84_metabolism_HR.csv")) %>%
  filter(cancer == "Pan cancer", color == "Worse survival", p < 0.05)

fig3_1 <- read_csv(file.path(dir_for_data, "pancancer_Hightriskmetagenes.csv")) %>%
  filter(Pathway_ID %in% fig3$celltype)
fig3_genes <- unique(fig3_1$Genes)
dt2 <- fig3_1 %>% select(Genes, Class)
colnames(dt2) <- c("gene", "class")
dt2$group <- "Metabolism"

# Figure4
fig4 <- read_csv(file.path(dir_for_data, "Figure4C_regulon_HR.csv")) %>%
  filter(cancer == "Pan cancer", color == "Worse survival", p < 0.05)

fig4_1 <- read_csv(file.path(dir_for_data, "Figure4C_core_tf.csv")) %>%
  filter(source %in% fig4$celltype)
fig4_genes <- union(unique(fig4_1$target), unique(fig4_1$source))
dt3 <- fig4_1 %>% select(source, target)
dt3$class <- dt3$source
a <- dt3[,c(1,3)] %>% as.data.frame()
colnames(a)[1] <- "gene"
aa <- dt3[,c(2,3)] %>% as.data.frame()
colnames(aa)[1] <- "gene"
dt3 <- rbind(a, aa)
dt3$group <- "Regulon"

# Figure5
fig5 <- read_csv(file.path(dir_for_data, "TCGAplot_60_LRs_HR.csv")) %>%
  filter(cancer == "Pan cancer", color == "Worse survival", p < 0.05)

fig5_1 <- read_csv(file.path(dir_for_data, "pancancer_HightriskLRsgenes.csv")) %>%
  filter(interaction_name_2 %in% fig5$celltype) %>%
  dplyr::select(interaction_name_2, interaction_name) %>%
  separate_rows(interaction_name, sep = "_") %>%
  distinct()
fig5_genes <- unique(fig5_1$interaction_name)
colnames(fig5_1) <- c("class", "gene")
dt4 <- fig5_1 %>% select(gene, class) %>% mutate(group = "LRs")

scRNA_highrisk_genes <- c(fig2_genes, fig3_genes, fig4_genes, fig5_genes) %>% unique()
write.csv(scRNA_highrisk_genes, file = file.path(dir_for_data, "scRNA_highrisk_genes_precise.csv"))

bulk_highrisk_genes <- read_csv(file.path(dir_for_data, "pancancer_Hightriskgenes.csv")) %>%
  filter(cancer == "Pan cancer", color == "Worse survival", p < 0.05) %>%
  pull(gene)
write.csv(bulk_highrisk_genes, file = file.path(dir_for_data, "bulk_highrisk_genes.csv"))
dt5 <- data.frame(gene = bulk_highrisk_genes,
                  class = "TCGA",
                  group = "Bulk")

sc_bulk <- intersect(bulk_highrisk_genes, scRNA_highrisk_genes)
write.csv(sc_bulk, file = file.path(dir_for_data, "sc_bulk_genes_precise.csv"))


## step 4. data split
library(TCGAplot)
library(tidyverse)
library(caret)

tpm <- get_all_tpm()
exprSet <- tpm %>%
  filter(Group == "Tumor") %>%
  mutate(ID = stringr::str_sub(row.names(.), 1, 12)) %>%
  distinct(ID, .keep_all = TRUE)
exprSet <- exprSet[, -c(1, 2, 3)]

meta <- get_all_meta() %>% mutate(ID = row.names(.))

samsample <- intersect(row.names(exprSet), row.names(meta))
exprSet <- exprSet[samsample, ]
exprSet$ID <- row.names(exprSet)
exp_meta <- inner_join(meta, exprSet, by = "ID") %>%
  select(ID, time, event, everything())
genelist <- read_csv("Figure6/TCGA/sc_bulk_genes_precise.csv")$x

set.seed(1234)
sample_indices <- sample(c(TRUE, FALSE), nrow(exprSet), replace = TRUE, prob = c(0.7, 0.3))
list_train_vali_Data <- list(
  Dataset1 = exp_meta[exp_meta$ID %in% rownames(exprSet[sample_indices, ]), ],
  Dataset2 = exp_meta[exp_meta$ID %in% rownames(exprSet[!sample_indices, ]), ]
)

list_train_vali_Data$Dataset1 <- list_train_vali_Data$Dataset1 %>% rename(OS.time = time, OS = event)
list_train_vali_Data$Dataset2 <- list_train_vali_Data$Dataset2 %>% rename(OS.time = time, OS = event)
save(genelist, list_train_vali_Data, file = "mime1.RData")

## step 5. model construction and evaluation
library(Mime1)
load("Figure6_all_mime1.RData")
source("/data3/home/yang/GI/ML.Dev.Prog.Sig_all.R")
res <- ML.Dev.Prog.Sig_all(train_data = list_train_vali_Data$Dataset1,
                           list_train_vali_Data = list_train_vali_Data,
                           unicox.filter.for.candi = T,
                           unicox_p_cutoff = 0.05,
                           candidate_genes = genelist,
                           mode = 'all', nodesize = 5, seed = 5201314)


## 1. Figure 6B-C. C-index and AUC ----
df("Figure6B_mime1.pdf", height = 6, width = 6)
cindex_dis_all(res, validate_set = names(list_train_vali_Data)[-1], 
               order = names(list_train_vali_Data), width = 0.35)
dev.off()

all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res, train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data, mode = 'all', AUC_time = 1,
                             auc_cal_method = "KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res, train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = 'all', AUC_time = 3,
                             auc_cal_method = "KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res, train_data = list_train_vali_Data[["Dataset1"]],
                             inputmatrix.list = list_train_vali_Data, mode = 'all', AUC_time = 5,
                             auc_cal_method = "KM")

pdf("Figure6C_mime1.pdf", height = 6, width = 6)
auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set = names(list_train_vali_Data)[-1],
            order = names(list_train_vali_Data),
            width = 0.35, year = 1)
dev.off()


## 2. Figure 6D. the importance of top PPMS genes ----
dt <- res[["ml.res"]][["RSF"]][["importance"]] %>% as.data.frame()
dt$sig <- rownames(dt)
dt <- dt[order(dt$., decreasing = T), ]
top20_dt <- dt[1:20, ]
colnames(top20_dt)[1] <- "importance"
top20_dt$sig <- factor(top20_dt$sig, levels = top20_dt$sig)

library(ggplot2)
ggplot(top20_dt, aes(x = sig, y = importance)) +
  geom_bar(stat = "identity", fill = "#7AA8C7") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))
ggsave("top20.sig.pdf", width = 8.49, height = 4.16)


## 3. Figure 6E. Kaplan-Meier curves ----
survplot <- vector("list", 2) 
for (i in c(1:2)) {
  print(survplot[[i]] <- rs_sur(res, model_name = "RSF",
                                dataset = names(list_train_vali_Data)[i],
                                median.line = "hv",
                                cutoff = 0.5, conf.int = T,
                                xlab = "Day", pval.coord = c(1000, 0.9)))
}
pdf("Figure6E_survplot.pdf", width = 10, height = 5)


## 4. Figure 6F/Figure S10A. model comparison ----
library(Mime1)
library(tidyverse)

# pan-cancer comparison
load("Figure6_all_mime1.RData"); rm(all.auc.1y, all.auc.3y, genelist); gc()
set.seed(123)
my_dt1 <- Mime1::pre.prog.sig
my_dt11 <- my_dt1$Glioma[my_dt1$Glioma$model %in% sample(my_dt1[["Glioma"]]$model, 10), ]
my_dt12 <- my_dt1$GBM[my_dt1$GBM$model %in% sample(my_dt1[["GBM"]]$model, 10), ]
my_dt13 <- my_dt1$LGG[my_dt1$LGG$model %in% sample(my_dt1[["LGG"]]$model, 10), ]
dt <- rbind(my_dt11, my_dt12, my_dt13)
rs.cancer <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = T, 
                                     collected_sig_table = my_dt,
                                     list_input_data = list_train_vali_Data)
pdf("model_comparison.pdf", height = 25, width = 10)
cindex_comp(rs.cancer,
            res,
            model_name = "RSF",
            dataset = names(list_train_vali_Data))
dev.off()

my_dt <- read.csv("model.comparison.csv")
c_names <- unique(my_dt$Cancer)
my_dt$Cancer[which(my_dt$Cancer == "Glioma")] <- "GBM"
my_dt$Cancer[which(my_dt$Cancer == "LUNG")] <- "LUAD"

load("Figure6_TCGA_exp_meta.RData")
result <- list()
for (cancer in c_names) {
  my_dt2 <- my_dt[my_dt$Cancer == cancer, ]
  dt <- exp_meta[exp_meta$Cancer == cancer, -c(1, 4, 5, 6)] %>% 
    rownames_to_column(var = "ID") %>% 
    as.data.frame() %>% 
    mutate(time = time * 30) %>%
    select(ID, time, event, everything()) %>%
    rename(OS.time = time, OS = event)
  list_train_vali_Data1 <- list(Dataset1 = dt, Dataset2 = dt)
  re.cancer <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = T, 
                                       collected_sig_table = my_dt2,
                                       list_input_data = list_train_vali_Data1)
  result <- c(re.cancer, result)
}
all_result <- bind_rows(result)
all_result <- all_result[all_result$ID == "Dataset2", ]
all_result$model <- names(result)
a <- my_dt %>% select(model, Cancer) %>% distinct()
aa <- inner_join(a, all_result)
aa[140, ] <- c("Our model", "pancancer", "Dataset2", "0.750")
aa <- aa[order(aa$Cindex, decreasing = T), ]
aa$model <- factor(aa$model, levels = rev(unique(aa$model)))
aaa <- as.data.frame(table(aa$Cancer))
aaa$new <- paste0(aaa$Var1, " (n = ", aaa$Freq, ")")
colnames(aaa)[1] <- "Cancer"
aaaa <- inner_join(aa, aaa)

ggplot(aaaa, aes(x = ID, y = model, colour = new)) +
  geom_point() +
  theme_test() +
  scale_color_manual(values = sample(col_vector, 17))
ggsave("model_comparison_color-num.pdf", height = 25, width = 10)

pdf("model_comparison.pdf", height = 25, width = 10)
cindex_comp(result,
            res,
            model_name = "RSF",
            dataset = names(list_train_vali_Data1))
dev.off()


## 5. Figure 6G. C-index of PPMS and clinical factors ----
## all
surv_obj <- Surv(dt$OS.time, dt$OS)
cox_model1 <- coxph(surv_obj ~ RS, data = dt)
cox_model2 <- coxph(surv_obj ~ age, data = dt)
cox_model3 <- coxph(surv_obj ~ gender, data = dt)
cox_model4 <- coxph(surv_obj ~ stage, data = dt)
summary(cox_model1)
summary(cox_model2)
summary(cox_model3)
summary(cox_model4)

data <- data.frame(
  Variable = c("RS", "Age", "Gender", "Stage"),
  C_index = c(0.904, 0.626, 0.53, 0.648),
  SE = c(0.003, 0.007, 0.006, 0.007)
)

comparison <- data %>%
  rowwise() %>%
  mutate(
    z_score = ifelse(Variable == "RS", NA,
                     (C_index - data$C_index[data$Variable == "RS"]) / 
                       sqrt(SE^2 + data$SE[data$Variable == "RS"]^2)),
    p_value = ifelse(Variable == "RS", NA,
                     2 * pnorm(abs(z_score), lower.tail = FALSE))
  )
data <- bind_cols(data, comparison[, c("z_score", "p_value")])

data$Significance <- case_when(
  data$p_value < 0.0001 ~ "****",
  data$p_value < 0.001 ~ "***",
  data$p_value < 0.01 ~ "**",
  data$p_value < 0.05 ~ "*",
  TRUE ~ ""
)
data$Variable <- factor(data$Variable, levels = c("RS", "Age", "Gender", "Stage"))

ggplot(data, aes(x = Variable, y = C_index, fill = Variable)) +
  geom_col(color = "black", width = 0.7) +
  # geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) +
  labs(title = "C-index (Compared with RS)",
       y = "C-index",
       x = "") +
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("#015493", "#019092", "#999999", "#f4a99b"))
ggsave("Figure6G_RS-two-sided_z-score_test_all.pdf", width = 4.77, height = 4.18)


## test
surv_obj <- Surv(test$OS.time, test$OS)
cox_model1 <- coxph(surv_obj ~ RS, data = test)
cox_model2 <- coxph(surv_obj ~ age, data = test)
cox_model3 <- coxph(surv_obj ~ gender, data = test)
cox_model4 <- coxph(surv_obj ~ stage, data = test)
summary(cox_model1)
summary(cox_model2)
summary(cox_model3)
summary(cox_model4)

data <- data.frame(
  Variable = c("RS", "Age", "Gender", "Stage"),
  C_index = c(0.953, 0.616, 0.531, 0.646),
  SE = c(0.002, 0.008, 0.007, 0.009)
)

comparison <- data %>%
  rowwise() %>%
  mutate(
    z_score = ifelse(Variable == "RS", NA,
                     (C_index - data$C_index[data$Variable == "RS"]) / 
                       sqrt(SE^2 + data$SE[data$Variable == "RS"]^2)),
    p_value = ifelse(Variable == "RS", NA,
                     2 * pnorm(abs(z_score), lower.tail = FALSE))
  )

data <- bind_cols(data, comparison[, c("z_score", "p_value")])

data$Significance <- case_when(
  data$p_value < 0.0001 ~ "****",
  data$p_value < 0.001 ~ "***",
  data$p_value < 0.01 ~ "**",
  data$p_value < 0.05 ~ "*",
  TRUE ~ ""
)
data$Variable <- factor(data$Variable, levels = c("RS", "Age", "Gender", "Stage"))

ggplot(data, aes(x = Variable, y = C_index, fill = Variable)) +
  geom_col(color = "black", width = 0.7) +
  # geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) +
  labs(title = "C-index (Compared with RS)",
       y = "C-index",
       x = "") +
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("#015493", "#019092", "#999999", "#f4a99b"))
ggsave("Figure6G_RS-two-sided_z-score_test_test.pdf", width = 4.77, height = 4.18)


surv_obj <- Surv(val$OS.time, val$OS)
cox_model1 <- coxph(surv_obj ~ RS, data = val)
cox_model2 <- coxph(surv_obj ~ age, data = val)
cox_model3 <- coxph(surv_obj ~ gender, data = val)
cox_model4 <- coxph(surv_obj ~ stage, data = val)
summary(cox_model1)
summary(cox_model2)
summary(cox_model3)
summary(cox_model4)

data <- data.frame(
  Variable = c("RS", "Age", "Gender", "Stage"),
  C_index = c(0.75, 0.649, 0.528, 0.654),
  SE = c(0.01, 0.012, 0.012, 0.014)
)

comparison <- data %>%
  rowwise() %>%
  mutate(
    z_score = ifelse(Variable == "RS", NA,
                     (C_index - data$C_index[data$Variable == "RS"]) / 
                       sqrt(SE^2 + data$SE[data$Variable == "RS"]^2)),
    p_value = ifelse(Variable == "RS", NA,
                     2 * pnorm(abs(z_score), lower.tail = FALSE))
  )

data <- bind_cols(data, comparison[, c("z_score", "p_value")])

data$Significance <- case_when(
  data$p_value < 0.0001 ~ "****",
  data$p_value < 0.001 ~ "***",
  data$p_value < 0.01 ~ "**",
  data$p_value < 0.05 ~ "*",
  TRUE ~ ""
)
data$Variable <- factor(data$Variable, levels = c("RS", "Age", "Gender", "Stage"))

ggplot(data, aes(x = Variable, y = C_index, fill = Variable)) +
  geom_col(color = "black", width = 0.7) +
  # geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) +
  labs(title = "C-index (Compared with RS)",
       y = "C-index",
       x = "") +
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("#015493", "#019092", "#999999", "#f4a99b"))
ggsave("Figure6G_RS-two-sided_z-score_test_validation.pdf", width = 4.77, height = 4.18)


## 5. Figure S8A. random-split validation of PPMS ----
library(randomForestSRC)
library(caret)
library(purrr)

tcga_data <- rbind(list_train_vali_Data$Dataset1, list_train_vali_Data$Dataset2)
sig_genes <- names(res$ml.res$RSF$importance)
sig_genes <- gsub("-", ".", sig_genes)
sig_genes <- gsub("_", ".", sig_genes)
sig_genes <- make.names(sig_genes)
colnames(tcga_data) <- gsub("-", ".", colnames(tcga_data))
colnames(tcga_data) <- gsub("_", ".", colnames(tcga_data))
colnames(tcga_data) <- make.names(colnames(tcga_data))
sig_genes <- intersect(sig_genes, colnames(tcga_data))

tcga_data2 <- tcga_data %>%
  dplyr::select(ID, OS.time, OS, all_of(sig_genes)) %>%
  filter(!is.na(OS.time), !is.na(OS), OS.time > 0) %>%
  mutate(
    ID = as.character(ID),
    OS.time = as.numeric(OS.time),
    OS = as.numeric(OS)
  )
tcga_data2[, sig_genes] <- apply(
  tcga_data2[, sig_genes, drop = FALSE],
  2,
  as.numeric
) %>% as.data.frame()
tcga_data2[, sig_genes] <- apply(
  tcga_data2[, sig_genes, drop = FALSE],
  2,
  function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  }
) %>% as.data.frame()

# time-dependent AUC
calc_auc <- function(dat, times = c(365, 1095, 1825)) {
  auc <- tryCatch({
    roc <- timeROC(
      T = dat$OS.time,
      delta = dat$OS,
      marker = dat$RS,
      cause = 1,
      weighting = "marginal",
      times = times,
      iid = FALSE
    )
    as.numeric(roc$AUC)
  }, error = function(e) {
    rep(NA_real_, length(times))
  })
  auc
}

# repeated random split validation
run_one_repeat_fast <- function(i, data, genes, train_prop = 0.7) {
  set.seed(1234 + i)
  inTrain <- caret::createDataPartition(
    y = data$OS,
    p = train_prop,
    list = FALSE
  )
  train_dat <- data[inTrain, ]
  valid_dat <- data[-inTrain, ]
  train_fit <- train_dat %>%
    dplyr::select(OS.time, OS, all_of(genes))
  valid_fit <- valid_dat %>%
    dplyr::select(ID, OS.time, OS, all_of(genes))
  fit <- tryCatch({
    rfsrc(
      Surv(OS.time, OS) ~ .,
      data = train_fit,
      ntree = 500,
      nodesize = 5,
      splitrule = "logrank",
      importance = FALSE,
      proximity = FALSE,
      forest = TRUE,
      seed = 5201314 + i
    )
  }, error = function(e) NULL)
  if (is.null(fit)) {
    return(data.frame(
      repeat_id = i,
      n_train = nrow(train_dat),
      n_valid = nrow(valid_dat),
      n_event_train = sum(train_dat$OS == 1),
      n_event_valid = sum(valid_dat$OS == 1),
      C_index = NA_real_,
      AUC_1year = NA_real_,
      AUC_3year = NA_real_,
      AUC_5year = NA_real_
    ))
  }
  pred <- predict(
    fit,
    newdata = valid_fit[, genes, drop = FALSE]
  )
  valid_pred <- valid_fit %>%
    mutate(RS = as.numeric(pred$predicted))
  cindex <- tryCatch({
    as.numeric(
      summary(
        coxph(Surv(OS.time, OS) ~ RS, data = valid_pred)
      )$concordance[1]
    )
  }, error = function(e) NA_real_)
  auc <- calc_auc(valid_pred)
  data.frame(
    repeat_id = i,
    n_train = nrow(train_dat),
    n_valid = nrow(valid_dat),
    n_event_train = sum(train_dat$OS == 1),
    n_event_valid = sum(valid_dat$OS == 1),
    C_index = cindex,
    AUC_1year = auc[1],
    AUC_3year = auc[2],
    AUC_5year = auc[3]
  )
}

n_repeat <- 100
repeated_validation_results <- map_dfr(
  1:n_repeat,
  ~run_one_repeat_fast(
    i = .x,
    data = tcga_data2,
    genes = sig_genes,
    train_prop = 0.7
  )
)
write.csv(
  repeated_validation_results,
  "PPMS_RSF_repeated_random_split_validation_test2_each_repeat.csv",
  row.names = FALSE
)
repeated_validation_results

summary_metric <- function(x) {
  x <- as.numeric(x)
  data.frame(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    lower_95CI = as.numeric(quantile(x, 0.025, na.rm = TRUE)),
    upper_95CI = as.numeric(quantile(x, 0.975, na.rm = TRUE))
  )
}
repeated_validation_summary <- bind_rows(
  C_index  = summary_metric(repeated_validation_results$C_index),
  AUC_1year = summary_metric(repeated_validation_results$AUC_1year),
  AUC_3year = summary_metric(repeated_validation_results$AUC_3year),
  AUC_5year = summary_metric(repeated_validation_results$AUC_5year),
  .id = "Metric"
)
write.csv(
  repeated_validation_summary,
  "PPMS_RSF_repeated_random_split_validation_summary.csv",
  row.names = FALSE
)
repeated_validation_summary

# plot
plot_df <- repeated_validation_results %>%
  select(repeat_id, C_index, AUC_1year, AUC_3year, AUC_5year) %>%
  pivot_longer(
    cols = -repeat_id,
    names_to = "Metric",
    values_to = "Value"
  )
p <- ggplot(plot_df, aes(x = Metric, y = Value)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Repeated random-split validation of the PPMS-RSF model",
    x = NULL,
    y = "Performance"
  )
p
pdf("PPMS_RSF_repeated_random_split_validation.pdf", width = 4.58, height = 3.7)
print(p)
dev.off()


## 6. Figure S8B. ROC ----
## Time-dependent ROC
library(timeROC)
library(survival)
library(ggplot2)
library(tidyverse)
library(patchwork)

# train and test
train_rs <- res$riskscore$RSF$Dataset1
valid_rs <- res$riskscore$RSF$Dataset2
train_rs$Cohort <- "Training"
valid_rs$Cohort <- "Validation"
all_rs <- bind_rows(train_rs, valid_rs)
all_rs$Cohort <- "All TCGA"
time_points <- c(365, 1095, 1825)
time_labels <- c("1-year", "3-year", "5-year")
calc_timeROC <- function(dat, cohort_name) {
  dat <- dat %>%
    filter(!is.na(OS.time), !is.na(OS), !is.na(RS)) %>%
    mutate(
      OS.time = as.numeric(OS.time),
      OS = as.numeric(OS),
      RS = as.numeric(RS)
    )
  roc <- timeROC(
    T = dat$OS.time,
    delta = dat$OS,
    marker = dat$RS,
    cause = 1,
    weighting = "marginal",
    times = time_points,
    iid = TRUE
  )
  auc_df <- data.frame(
    Cohort = cohort_name,
    Time = time_labels,
    AUC = as.numeric(roc$AUC)
  )
  list(
    roc = roc,
    auc = auc_df
  )
}
roc_train <- calc_timeROC(train_rs, "Training")
roc_valid <- calc_timeROC(valid_rs, "Validation")
roc_all   <- calc_timeROC(all_rs, "All TCGA")
auc_summary <- bind_rows(
  roc_train$auc,
  roc_valid$auc,
  roc_all$auc
)
write.csv(auc_summary, "PPMS_RSF_time_dependent_AUC.csv", row.names = FALSE)
auc_summary

# plot
plot_timeROC <- function(roc_obj, cohort_name) {
  plot_df <- bind_rows(
    data.frame(
      FPR = roc_obj$roc$FP[, 1],
      TPR = roc_obj$roc$TP[, 1],
      Time = paste0("1-year AUC = ", round(roc_obj$roc$AUC[1], 3))
    ),
    data.frame(
      FPR = roc_obj$roc$FP[, 2],
      TPR = roc_obj$roc$TP[, 2],
      Time = paste0("3-year AUC = ", round(roc_obj$roc$AUC[2], 3))
    ),
    data.frame(
      FPR = roc_obj$roc$FP[, 3],
      TPR = roc_obj$roc$TP[, 3],
      Time = paste0("5-year AUC = ", round(roc_obj$roc$AUC[3], 3))
    )
  )
  ggplot(plot_df, aes(x = FPR, y = TPR, color = Time)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      title = cohort_name,
      x = "False positive rate",
      y = "True positive rate",
      color = NULL
    )
}
p_train <- plot_timeROC(roc_train, "Training set")
p_valid <- plot_timeROC(roc_valid, "Validation set")
p_all   <- plot_timeROC(roc_all, "All TCGA")

pdf("PPMS_RSF_time_dependent_ROC.pdf", width = 17, height = 4)
print(p_train | p_valid | p_all)
dev.off()


## 7. Figure S8C. calibration curve ----
# Cox calibration model
train_rs <- train_rs %>%
  filter(!is.na(OS.time), !is.na(OS), !is.na(RS)) %>%
  mutate(
    OS.time = as.numeric(OS.time),
    OS = as.numeric(OS),
    RS = as.numeric(RS)
  )
cox_cal <- coxph(
  Surv(OS.time, OS) ~ RS,
  data = train_rs,
  x = TRUE
)

# Cox model
predict_event_prob <- function(cox_model, newdata, times) {
  newdata <- newdata %>%
    filter(!is.na(OS.time), !is.na(OS), !is.na(RS)) %>%
    mutate(
      OS.time = as.numeric(OS.time),
      OS = as.numeric(OS),
      RS = as.numeric(RS)
    )
  base_haz <- basehaz(cox_model, centered = FALSE)
  beta <- coef(cox_model)["RS"]
  lp <- beta * newdata$RS
  pred_list <- lapply(times, function(t0) {
    bh_sub <- base_haz %>% filter(time <= t0)
    if (nrow(bh_sub) == 0) {
      H0_t <- 0
    } else {
      H0_t <- max(bh_sub$hazard, na.rm = TRUE)
    }
    pred_surv <- exp(-H0_t * exp(lp))
    pred_event <- 1 - pred_surv
    data.frame(
      ID = newdata$ID,
      OS.time = newdata$OS.time,
      OS = newdata$OS,
      RS = newdata$RS,
      time = t0,
      pred_event = pred_event
    )
  })
  bind_rows(pred_list)
}

# calibration
calc_calibration <- function(dat, cohort_name, n_group = 5) {
  dat <- dat %>%
    filter(!is.na(OS.time), !is.na(OS), !is.na(RS)) %>%
    mutate(
      OS.time = as.numeric(OS.time),
      OS = as.numeric(OS),
      RS = as.numeric(RS)
    )
  pred_df <- predict_event_prob(cox_cal, dat, time_points)
  cal_df <- pred_df %>%
    group_by(time) %>%
    mutate(
      risk_group = ntile(pred_event, n_group)
    ) %>%
    ungroup()
  obs_df <- cal_df %>%
    group_by(time, risk_group) %>%
    group_modify(~{
      t0 <- .y$time[1]
      fit <- survfit(Surv(OS.time, OS) ~ 1, data = .x)
      s <- summary(fit, times = t0, extend = TRUE)
      data.frame(
        n = nrow(.x),
        mean_pred_event = mean(.x$pred_event, na.rm = TRUE),
        observed_event = 1 - s$surv,
        lower_event = 1 - s$upper,
        upper_event = 1 - s$lower
      )
    }) %>%
    ungroup() %>%
    mutate(
      Cohort = cohort_name,
      Time = case_when(
        time == 365 ~ "1-year",
        time == 1095 ~ "3-year",
        time == 1825 ~ "5-year",
        TRUE ~ as.character(time)
      )
    )
  obs_df
}
cal_train <- calc_calibration(train_rs, "Training")
cal_valid <- calc_calibration(valid_rs, "Validation")
cal_all   <- calc_calibration(all_rs, "All TCGA")
cal_summary <- bind_rows(cal_train, cal_valid, cal_all)
write.csv(cal_summary, "PPMS_RSF_calibration_summary.csv", row.names = FALSE)
cal_summary

plot_calibration <- function(cal_df, cohort_name) {
  ggplot(cal_df %>% filter(Cohort == cohort_name),
         aes(x = mean_pred_event, y = observed_event)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    geom_errorbar(
      aes(ymin = lower_event, ymax = upper_event),
      width = 0.015,
      linewidth = 0.4
    ) +
    geom_point(size = 2.2) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ Time, nrow = 1) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 11),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    ) +
    labs(
      title = cohort_name,
      x = "Predicted event probability",
      y = "Observed event probability"
    )
}
p_train <- plot_calibration(cal_summary, "Training")
p_valid <- plot_calibration(cal_summary, "Validation")
p_all   <- plot_calibration(cal_summary, "All TCGA")
pdf("PPMS_RSF_calibration_curve.pdf", width = 9.32, height = 11.1)
print(p_train / p_valid / p_all)
dev.off()


## 8. Figure S8D-E. Related to Figure 6E ----
## step function
# Test
rs_data <- res[["riskscore"]][["RSF"]][["Dataset1"]]
rs_data$RS_group <- ifelse(
  rs_data$RS >= median(rs_data$RS, na.rm = TRUE),
  "High",
  "Low"
)
fit <- survfit(Surv(OS.time, OS) ~ RS_group, data = rs_data)
pdf("Test_stepfunction.pdf", width = 5.5, height = 5)
plot(fit, conf.int = TRUE, col = c("#868686", "#B24745"),
     lwd = 2, xlab = "Time", ylab = "Survival probability")
legend("topright", 
       legend = c("Low", "High"), 
       col = c("#868686", "#B24745"), 
       lwd = 2)
dev.off()

# Validation
rs_data <- res[["riskscore"]][["RSF"]][["Dataset2"]]
rs_data$RS_group <- ifelse(
  rs_data$RS >= median(rs_data$RS, na.rm = TRUE),
  "High",
  "Low"
)
fit <- survfit(Surv(OS.time, OS) ~ RS_group, data = rs_data)
pdf("Validation_stepfunction.pdf", width = 5.5, height = 5)
plot(fit, conf.int = TRUE, col = c("#868686", "#B24745"),
     lwd = 2, xlab = "Time", ylab = "Survival probability")
legend("topright", 
       legend = c("Low", "High"), 
       col = c("#868686", "#B24745"), 
       lwd = 2)
dev.off()

## cumulative hazard plot
# Test
pdf("Test_cumulativeHazard.pdf", width = 5.5, height = 6.8)
ggsurvplot(
  fit,
  data = rs_data,
  fun = "cumhaz",
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#B24745", "#868686")
)
dev.off()

# Vlidation
pdf("Validation_cumulativeHazard.pdf", width = 5.5, height = 6.8)
ggsurvplot(
  fit,
  data = rs_data,
  fun = "cumhaz",
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#B24745", "#868686")
)
dev.off()


## 9. Figure S8F. Additional validation ----
load("/data3/home/yang/GI/test/Figure6_all_mime1.RData")
fit_rsf <- res$ml.res$RSF
out_prefix <- "External_GEO_RSF_validation"

library(randomForestSRC)

# 1. 
sig_genes <- fit_rsf$xvar.names
if (is.null(sig_genes)) {
  sig_genes <- names(fit_rsf$importance)
}
sig_genes <- gsub("-", ".", sig_genes)
sig_genes <- gsub("_", ".", sig_genes)
sig_genes <- make.names(sig_genes)
length(sig_genes)

# 2. 
tcga_train <- list_train_vali_Data$Dataset1 %>%
  as.data.frame()

colnames(tcga_train) <- gsub("-", ".", colnames(tcga_train))
colnames(tcga_train) <- gsub("_", ".", colnames(tcga_train))
colnames(tcga_train) <- make.names(colnames(tcga_train))
missing_tcga <- setdiff(sig_genes, colnames(tcga_train))

length(missing_tcga)
missing_tcga

# 3. 
tcga_train <- tcga_train %>%
  mutate(
    OS.time = as.numeric(OS.time),
    OS = as.numeric(OS)
  ) %>%
  filter(
    !is.na(OS.time),
    !is.na(OS),
    OS.time > 0
  )

tcga_train[, sig_genes] <- lapply(
  tcga_train[, sig_genes, drop = FALSE],
  as.numeric
)

# 4. Z-score
zscore_fun <- function(x) {
  
  mu <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- mu
  
  sigma <- sd(x, na.rm = TRUE)
  
  if (!is.finite(sigma) || sigma == 0) {
    return(rep(0, length(x)))
  }
  
  (x - mu) / sigma
}

tcga_train_z <- tcga_train

tcga_train_z[, sig_genes] <- lapply(
  tcga_train_z[, sig_genes, drop = FALSE],
  zscore_fun
)

# 5. 
train_fit <- tcga_train_z %>%
  select(
    OS.time,
    OS,
    all_of(sig_genes)
  )

set.seed(5201314)

fit_rsf_z <- rfsrc(
  Surv(OS.time, OS) ~ .,
  data = train_fit,
  ntree = 1000,
  nodesize = 5,
  mtry = 18,
  splitrule = "logrank",
  nsplit = 10,
  importance = TRUE,
  proximity = FALSE,
  forest = TRUE,
  seed = 5201314
)

fit_rsf_z

save(
  fit_rsf_z,
  sig_genes,
  file = "PPMS_RSF_Zscore.RData"
)

length(sig_genes)
length(missing_tcga)
fit_rsf_z

##
load("/data3/home/yang/GI/test/SecondV/PPMS_RSF_Zscore.RData")
load("/data3/home/yang/GI/test/FirstV/Survival_data/PDAC_survival.RData")
load("/data3/home/yang/GI/test/FirstV/LUAD_survival.RData")
dt <- LUAD_survival
external_dat <- as.data.frame(dt)
external_dat$ID <- rownames(external_dat)

colnames(external_dat) <- gsub("-", ".", colnames(external_dat))
colnames(external_dat) <- gsub("_", ".", colnames(external_dat))
colnames(external_dat) <- make.names(colnames(external_dat))

external_dat <- external_dat %>%
  mutate(
    ID = as.character(ID),
    OS.time = as.numeric(OS.time),
    OS = as.numeric(OS)
  ) %>%
  filter(
    !is.na(OS.time),
    !is.na(OS),
    OS.time > 0
  )

sig_genes <- fit_rsf_z$xvar.names

length(sig_genes)

available_genes <- intersect(
  sig_genes,
  colnames(external_dat)
)

missing_genes <- setdiff(
  sig_genes,
  colnames(external_dat)
)

cat("Available PPMS genes:",
    length(available_genes), "/", length(sig_genes), "\n")

cat("Missing PPMS genes:",
    length(missing_genes), "\n")

cat("Gene coverage:",
    round(length(available_genes) / length(sig_genes) * 100, 2),
    "%\n")

external_dat[, available_genes] <- lapply(
  external_dat[, available_genes, drop = FALSE],
  as.numeric
)

zscore_external <- function(x) {
  
  mu <- mean(x, na.rm = TRUE)
  
  if (!is.finite(mu)) {
    return(rep(0, length(x)))
  }
  
  x[is.na(x)] <- mu
  
  sigma <- sd(x, na.rm = TRUE)
  
  if (!is.finite(sigma) || sigma == 0) {
    return(rep(0, length(x)))
  }
  
  (x - mu) / sigma
}

external_dat[, available_genes] <- lapply(
  external_dat[, available_genes, drop = FALSE],
  zscore_external
)

if (length(missing_genes) > 0) {
  external_dat[, missing_genes] <- 0
}

external_x <- external_dat[
  ,
  sig_genes,
  drop = FALSE
]

dim(external_x)

stopifnot(
  identical(
    colnames(external_x),
    fit_rsf_z$xvar.names
  )
)

stopifnot(!anyNA(external_x))

pred <- predict(
  fit_rsf_z,
  newdata = external_x
)

external_dat$RS <- as.numeric(pred$predicted)

summary(external_dat$RS)

length(available_genes)
length(missing_genes)
dim(external_x)
summary(external_dat$RS)

library(survival)

cox_fit <- coxph(
  Surv(OS.time, OS) ~ RS,
  data = external_dat
)

cox_sum <- summary(cox_fit)

cox_result <- data.frame(
  Cohort = "CRC_GSE161158",
  n_sample = nrow(external_dat),
  n_event = sum(external_dat$OS == 1, na.rm = TRUE),
  n_PPMS_genes = length(sig_genes),
  n_available_genes = length(available_genes),
  gene_coverage = length(available_genes) / length(sig_genes),
  HR = as.numeric(cox_sum$coefficients[1, "exp(coef)"]),
  lower_95CI = as.numeric(cox_sum$conf.int[1, "lower .95"]),
  upper_95CI = as.numeric(cox_sum$conf.int[1, "upper .95"]),
  cox_p = as.numeric(cox_sum$coefficients[1, "Pr(>|z|)"]),
  C_index = as.numeric(cox_sum$concordance[1]),
  C_index_se = as.numeric(cox_sum$concordance[2])
)

cox_result

library(timeROC)
library(dplyr)
library(ggplot2)

# 1/3/5 years
time_points <- c(365, 1095, 1825)
time_labels <- c("1-year", "3-year", "5-year")

roc_obj <- timeROC(
  T = external_dat$OS.time,
  delta = external_dat$OS,
  marker = external_dat$RS,
  cause = 1,
  weighting = "marginal",
  times = time_points,
  iid = TRUE
)

# AUC
auc_result <- data.frame(
  Cohort = "CRC_GSE161158",
  Time = time_labels,
  Time_days = time_points,
  AUC = as.numeric(roc_obj$AUC)
)

auc_result

roc_df <- bind_rows(
  data.frame(
    FPR = roc_obj$FP[, 1],
    TPR = roc_obj$TP[, 1],
    Time = paste0("1-year AUC = ", round(roc_obj$AUC[1], 3))
  ),
  data.frame(
    FPR = roc_obj$FP[, 2],
    TPR = roc_obj$TP[, 2],
    Time = paste0("3-year AUC = ", round(roc_obj$AUC[2], 3))
  ),
  data.frame(
    FPR = roc_obj$FP[, 3],
    TPR = roc_obj$TP[, 3],
    Time = paste0("5-year AUC = ", round(roc_obj$AUC[3], 3))
  )
)

p_roc <- ggplot(
  roc_df,
  aes(x = FPR, y = TPR, color = Time)
) +
  geom_line(linewidth = 1) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    color = "grey50"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  labs(
    title = "CRC_GSE161158",
    x = "False positive rate",
    y = "True positive rate",
    color = NULL
  )
p_roc


## 10. Figure S9A-B. PPMS in TCGA pan-cancer ----
load("/data3/home/yang/GI/test/Figure6_all_mime1.RData")
load("/data3/home/yang/GI/test/Figure6_TCGA_exp_meta.RData")
train_rs <- res$riskscore$RSF$Dataset1
valid_rs <- res$riskscore$RSF$Dataset2
train_rs$Cohort <- "Training"
valid_rs$Cohort <- "Validation"
all_rs <- bind_rows(train_rs, valid_rs)
all_rs <- bind_rows(train_rs, valid_rs)
exp_meta <- exp_meta[all_rs$ID, ]
all(rownames(exp_meta) == all_rs$ID)
all_rs$Cancer <- exp_meta$Cancer

# data processing for KM
plot_dat <- all_rs %>%
  filter(!is.na(OS.time), !is.na(OS), !is.na(RS), !is.na(Cancer)) %>%
  mutate(
    OS.time = as.numeric(OS.time),
    OS = as.numeric(OS),
    RS = as.numeric(RS),
    Cancer = as.character(Cancer)
  ) %>%
  group_by(Cancer) %>%
  mutate(
    RiskGroup = ifelse(RS >= median(RS, na.rm = TRUE), "High", "Low"),
    RiskGroup = factor(RiskGroup, levels = c("Low", "High"))
  ) %>%
  ungroup()

# KM
km_df <- plot_dat %>%
  group_by(Cancer) %>%
  group_modify(~{
    
    fit <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = .x)
    ss <- surv_summary(fit, data = .x)
    
    ss %>%
      mutate(
        Cancer = unique(.x$Cancer),
        RiskGroup = gsub("RiskGroup=", "", strata),
        time_year = time / 365.25
      )
  }) %>%
  ungroup()

# log-rank p
stat_df <- plot_dat %>%
  group_by(Cancer) %>%
  group_modify(~{
    n_sample <- nrow(.x)
    n_event <- sum(.x$OS == 1, na.rm = TRUE)
    p_val <- tryCatch({
      sd <- survdiff(Surv(OS.time, OS) ~ RiskGroup, data = .x)
      1 - pchisq(sd$chisq, length(sd$n) - 1)
    }, error = function(e) NA_real_)
    data.frame(
      n_sample = n_sample,
      n_event = n_event,
      logrank_p = p_val,
      x_pos = max(.x$OS.time, na.rm = TRUE) / 365.25 * 0.05,
      y_pos = 0.15,
      label = paste0(
        "n=", n_sample,
        ", events=", n_event,
        "\nP=", ifelse(is.na(p_val), "NA", format.pval(p_val, digits = 2, eps = 1e-4))
      )
    )
  }) %>%
  ungroup()
write.csv(
  stat_df,
  "PPMS_KM_by_cancer_logrank_statistics.csv",
  row.names = FALSE
)

# KM
p_all_km <- ggplot(km_df, aes(x = time_year, y = surv, color = RiskGroup)) +
  geom_step(linewidth = 0.6) +
  geom_text(
    data = stat_df,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 4,
    hjust = 0
  ) +
  scale_color_manual(values = c("#B24745", "#868686")) +
  facet_wrap(~ Cancer, ncol = 5, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 10),
    legend.position = "top"
  ) +
  labs(
    title = "Cancer-specific Kaplan–Meier curves for the PPMS-RSF risk score",
    x = "Time (years)",
    y = "Overall survival probability",
    color = NULL
  )

pdf("PPMS_KM_curves_by_cancer_all.pdf", width = 15, height = 21)
print(p_all_km)
dev.off()

## TCGAplot pan-cancer ROC of PPMS genes
setwd("pancancerROC")
sig <- read.csv("sig.genes.csv")
sig.genes <- sig$x
sig.genes[which(sig.genes == "H1.3")] <- "H1-3"
sig.genes[which(sig.genes == "H1.4")] <- "H1-4"

library(TCGAplot)
a <- get_cancers() %>% as.data.frame()
cancers <- as.character(unique(a$Var1))
survplot <- vector("list", 33) 
for (i in c(1:33)) {
  survplot[[i]] <- tryCatch({
    gs_roc(cancers[i], sig.genes, "RSF_Sig.genes")
  }, error = function(e) {
    message(paste("Error occurred for cancer", cancers[i], ":", e$message))
    NULL
  })
}
survplot1 <- Filter(function(x) !is.null(x), survplot)
pdf("ROC-Sig.genes-TCGA.pdf", width = 25, height = 25)
aplot::plot_list(gglist = survplot1, ncol = 5)
dev.off()

dt <- data.frame(cancer = c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
                            "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PCPG", "PRAD", "READ", "STAD", "THCA", "UCEC"),
                 AUC = c(0.871, 0.933, 0.875, 0.746, 0.975, 0.965, 0.99, 0.958, 0.768, 0.914,
                         0.862, 0.797, 0.881, 0.965, 0.742, 0.547, 0.736, 0.946, 0.845, 0.856, 0.976))
dt <- dt[order(dt$AUC), ]
dt$cancer <- factor(dt$cancer, levels = dt$cancer)
ggplot(dt, aes(y = cancer, x = AUC)) +
  geom_col(fill = "#d0e0ef") +
  geom_text(aes(label = AUC), position = position_dodge(0.9), vjust = 0) +
  xlab(label = "AUC") + ylab("") +
  theme_test()
ggsave("ROC-Sig.genes-TCGA_1.pdf", width = 5, height = 6.5)


## 11. Figure S11A-B PPMS with other clinical features ----
load("Figure6_all_mime1.RData")
load("Figure6_TCGA_exp_meta.RData")
exp_meta <- exp_meta[, c(1, 4, 5, 6)]
exp_meta$ID <- rownames(exp_meta)

res[["riskscore"]][["RSF"]][["Dataset1"]] -> testing
testing$RS_type <- ifelse(testing$RS >= median(testing$RS), 1, -1) # 1-High; -1-Low
testing$group <- "Testing"
res[["riskscore"]][["RSF"]][["Dataset2"]] -> validation
validation$RS_type <- ifelse(validation$RS >= median(validation$RS), 1, -1)
validation$group <- "Validation"

tcga <- rbind(testing, validation)

tcga_all <- inner_join(tcga, exp_meta, by = "ID")
tcga_all <- tcga_all %>%
  mutate(age_type = case_when(
    age <= 40 ~ 1,
    age > 40 & age <= 60 ~ 2,
    .default = 3
  ))
tcga_all <- tcga_all %>%
  mutate(stage_type = case_when(
    stage == "I" ~ 1,
    stage == "II" ~ 2,
    stage == "III" ~ 3,
    stage == "IV" ~ 4,
    .default = 0
  ))
tcga_all <- tcga_all %>%
  mutate(gender_type = case_when(
    gender == "M" ~ 1,
    .default = -1
  ))
write.csv(tcga_all, file = "FigureS11_tcga-RS-clinical.csv")

# ploting
gs_gender("COAD",sig.genes,"PPMS")
gs_3age("COAD",Sig.genes,"PPMS")
gs_age("COAD",Sig.genes,"PPMS")

## forest
ibrary(survival)
library(survminer)
dt <- read.csv("tcga-RS-clinical.csv")[, -1]

### HR森林图
## all
surv_obj <- Surv(dt$OS.time, dt$OS)
cox_model <- coxph(surv_obj ~ RS + age + gender + stage, data = dt)
summary(cox_model)
pdf("FigureS11B_RS-clinical_ggforest_all.pdf", width = 6.05, height = 4.18)
ggforest(cox_model, data = dt)
dev.off()

## testing
test <- dt[dt$group == "Testing", ]
surv_obj <- Surv(test$OS.time, test$OS)
cox_model <- coxph(surv_obj ~ RS + age + gender + stage, data = test)
pdf("FigureS11B_RS-clinical_ggforest_test.pdf", width = 6.05, height = 4.18)
ggforest(cox_model, data = test)
dev.off()

## validation
val <- dt[dt$group == "Validation", ]
surv_obj <- Surv(val$OS.time, val$OS)
cox_model <- coxph(surv_obj ~ RS + age + gender + stage, data = val)
pdf("FigureS11B_RS-clinical_ggforest_val.pdf", width = 6.05, height = 4.18)
ggforest(cox_model, data = val)
dev.off()

