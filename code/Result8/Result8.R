# Figure 8. CRISPR/Cas9 Screening prioritizes candidate therapeutic targets ----
# CRISPR, Chronos score
Chronos <- fread("CRISPRGeneEffect.csv", data.table = F)
Chronos <- Chronos %>%
  column_to_rownames(var = "V1")
Chronos <- as.data.frame(t(Chronos))
Chronos <- Chronos %>%
  rownames_to_column(var = "gene") %>%
  mutate(gene = sub(" \\(.*\\)", "", gene))

# information of models
meta <- read_csv("Model.csv")
meta <- meta[meta$ModelID %in% colnames(Chronos)[-1], ]
my_Chronos <- Chronos %>%
  select(gene, all_of(my_cells)) 
rownames(my_Chronos) <- my_Chronos$gene
my_Chronos <- my_Chronos[ , -1]

# mean Chronos score
my_Chronos$mean <- rowMeans(my_Chronos[,-1], na.rm = T)
my_Chronos <- my_Chronos %>%
  rownames_to_column(var = "gene") %>%
  arrange(mean)
my_Chronos1 <- my_Chronos %>%
  select(gene, mean)
my_Chronos2 <- my_Chronos
write.csv(my_Chronos1, "Pancancer_Chronos_Mean_CRISPRGeneEffect.csv", row.names = F)


## 1. Figure 8A-B. CRISPR/Cas9 Screening prioritizes candidate therapeutic targets ----
dt <- read.csv("Pancancer_Chronos_Mean_CRISPRGeneEffect.csv")
genes <- dt[dt$mean < -0.5, ]$gene
sig <- read.csv("sig.genes.csv")
common_genes <- intersect(genes, sig$x)

## 
dt1 <- readxl::read_xlsx("CTRP2.xlsx")
dt1$Correlation <- paste0("CTRP2_" ,dt1$Correlation)
dt1 <- split(dt1$Gene_name, dt1$Correlation)

dt2 <- readxl::read_xlsx("GDSC.xlsx")
dt2$Correlation <- paste0("GDSC_" ,dt2$Correlation)
dt2<- split(dt2$Gene_name, dt2$Correlation)

intersect(dt1$CTRP2_sensitivity, dt2$GDSC_sensitivity)
intersect(dt1$CTRP2_resistance, dt2$GDSC_resistance)
dt <- data.frame(
  Gene = c(intersect(dt1$CTRP2_sensitivity, dt2$GDSC_sensitivity),
           intersect(dt1$CTRP2_resistance, dt2$GDSC_resistance)),
  Drug_correlation = c(rep("sensitivity", length(intersect(dt1$CTRP2_sensitivity, dt2$GDSC_sensitivity))),
                       rep("resistance", length(intersect(dt1$CTRP2_resistance, dt2$GDSC_resistance))))
)
write.csv(dt, "Pancancer_Sig_Drug_Genes.csv", row.names = F)
dt <- read.csv("Pancancer_Sig_Drug_Genes.csv")
dt <- dt[dt$Drug_correlation == "sensitivity", ]
pdrug_genes <- dt$Gene
##

final_genes <- intersect(genes, pdrug_genes)
dt1 <- dt[dt$gene %in% final_genes, ]
colnames(dt1) <- c("Gene", "Pancancer_Chronos_Mean")
write.csv(dt1, "Pancancer_Chronos_Mean_CRISPRGeneEffect_Sig_Drug_Genes.csv", row.names = F)

ggplot(dt1, aes(x = reorder(Gene, -Pancancer_Chronos_Mean), y = Pancancer_Chronos_Mean)) +
  geom_bar(stat = "identity", fill = "#4c98a3") +
  coord_flip() +
  theme_test() +
  labs(x = "Gene", y = "Pancancer Chronos Mean Score") +
  theme(axis.text = element_text(color = "black"))
ggsave("Pancancer_Chronos_Mean_CRISPRGeneEffect_Sig_Drug_Genes_Barplot.pdf", width = 3.15, height = 6.75)

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    Set1 = genes,
    Set2 = pdrug_genes
  ),
  filename = NULL,
  fill = c("#66C2A5", "#FC8D62"),  
  alpha = 0.6,
  cex = 2,
  cat.cex = 2,
  lwd = 2
)
griE::grid.newpage()
griE::grid.draw(venn.plot)


## 2. Figure S14A-D. Cancer-lineage-specific dependency and tumor-normal expression profiles of PPMS-associated candidate targets ----
## Cancer-lineage-specific dependency analysis
# CRISPR, Chronos score
Chronos <- fread("CRISPRGeneEffect.csv", data.table = F)
Chronos <- Chronos %>%
  column_to_rownames(var = "V1")
Chronos <- as.data.frame(t(Chronos))
Chronos <- Chronos %>%
  rownames_to_column(var = "gene") %>%
  mutate(gene = sub(" \\(.*\\)", "", gene))

meta <- read_csv("Model.csv")
meta <- meta[meta$ModelID %in% colnames(Chronos)[-1], ]

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(readr)

outdir <- "PPMS_DepMap_lineage_dependency"
dir.create(outdir, showWarnings = FALSE)
lineage_col <- "OncotreeLineage"
min_cell_lines <- 5

# Chronos gene effect threshold
dependency_cutoff <- -0.5
strong_dependency_cutoff <- -1.0
top_genes <- c("RRM1", "PLK1", "CDC45", "RRM2", "SPC24")
top_genes <- intersect(top_genes, final_genes)

# DepMap metadata
meta_clean <- meta %>%
  as.data.frame() %>%
  dplyr::select(
    ModelID,
    CellLineName,
    DepmapModelType,
    OncotreeLineage,
    OncotreePrimaryDisease,
    OncotreeSubtype,
    OncotreeCode,
    PrimaryOrMetastasis,
    SampleCollectionSite
  ) %>%
  distinct(ModelID, .keep_all = TRUE) %>%
  filter(
    !is.na(ModelID),
    !is.na(.data[[lineage_col]]),
    .data[[lineage_col]] != ""
  ) %>%
  rename(lineage = all_of(lineage_col))

model_ids <- intersect(colnames(Chronos), meta_clean$ModelID)
message("Matched DepMap models: ", length(model_ids))
meta_clean <- meta_clean %>%
  filter(ModelID %in% model_ids)

# Chronos gene effect matrix
Chronos_clean <- Chronos %>%
  as.data.frame()
Chronos_clean$gene <- sub(" \\(.*\\)", "", Chronos_clean$gene)
Chronos_clean$gene <- gsub("-", ".", Chronos_clean$gene)
Chronos_clean$gene <- gsub("_", ".", Chronos_clean$gene)
final_genes_clean <- final_genes
final_genes_clean <- gsub("-", ".", final_genes_clean)
final_genes_clean <- gsub("_", ".", final_genes_clean)

Chronos_clean[, model_ids] <- lapply(
  Chronos_clean[, model_ids, drop = FALSE],
  function(x) as.numeric(as.character(x))
)

Chronos_gene <- Chronos_clean %>%
  dplyr::select(gene, all_of(model_ids)) %>%
  group_by(gene) %>%
  summarise(
    across(all_of(model_ids), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
available_genes <- intersect(final_genes_clean, Chronos_gene$gene)
missing_genes <- setdiff(final_genes_clean, Chronos_gene$gene)
message("Available PPMS candidate genes in DepMap: ",
        length(available_genes), "/", length(final_genes_clean))
if (length(missing_genes) > 0) {
  write.csv(
    data.frame(missing_gene = missing_genes),
    file.path(outdir, "missing_final_genes_in_DepMap.csv"),
    row.names = FALSE
  )
}

# gene × model × lineage
dep_long <- Chronos_gene %>%
  filter(gene %in% available_genes) %>%
  pivot_longer(
    cols = all_of(model_ids),
    names_to = "ModelID",
    values_to = "gene_effect"
  ) %>%
  left_join(meta_clean, by = "ModelID") %>%
  filter(
    !is.na(gene_effect),
    !is.na(lineage)
  )
write.csv(
  dep_long,
  file.path(outdir, "PPMS_46_DepMap_gene_effect_long.csv"),
  row.names = FALSE
)

# dependency summary for gene × cancer lineage
lineage_summary <- dep_long %>%
  group_by(gene, lineage) %>%
  summarise(
    n_cell_lines = n(),
    mean_gene_effect = mean(gene_effect, na.rm = TRUE),
    median_gene_effect = median(gene_effect, na.rm = TRUE),
    q25_gene_effect = quantile(gene_effect, 0.25, na.rm = TRUE),
    q75_gene_effect = quantile(gene_effect, 0.75, na.rm = TRUE),
    dependency_rate_0.5 = mean(gene_effect < dependency_cutoff, na.rm = TRUE),
    dependency_rate_1.0 = mean(gene_effect < strong_dependency_cutoff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cell_lines >= min_cell_lines)
write.csv(
  lineage_summary,
  file.path(outdir, "PPMS_46_DepMap_lineage_dependency_summary.csv"),
  row.names = FALSE
)

# Cancer-lineage-specific dependency test
test_one_gene_lineage <- function(g, lin, dat, min_n = 5) {
  tmp <- dat %>%
    filter(gene == g) %>%
    mutate(group = ifelse(lineage == lin, lin, "Other"))
  n_lineage <- sum(tmp$group == lin, na.rm = TRUE)
  n_other <- sum(tmp$group == "Other", na.rm = TRUE)
  if (n_lineage < min_n || n_other < min_n) {
    return(NULL)
  }
  x <- tmp$gene_effect[tmp$group == lin]
  y <- tmp$gene_effect[tmp$group == "Other"]
  wt <- wilcox.test(
    x = x,
    y = y,
    alternative = "two.sided",
    exact = FALSE
  )
  data.frame(
    gene = g,
    lineage = lin,
    n_lineage = n_lineage,
    n_other = n_other,
    median_lineage = median(x, na.rm = TRUE),
    median_other = median(y, na.rm = TRUE),
    delta_median = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
    p_value = wt$p.value
  )
}
gene_lineage_grid <- expand.grid(
  gene = available_genes,
  lineage = sort(unique(dep_long$lineage)),
  stringsAsFactors = FALSE
)
lineage_test <- pmap_dfr(
  gene_lineage_grid,
  function(gene, lineage) {
    test_one_gene_lineage(
      g = gene,
      lin = lineage,
      dat = dep_long,
      min_n = min_cell_lines
    )
  }
) %>%
  mutate(
    FDR_global = p.adjust(p_value, method = "BH"),
    stronger_dependency_in_lineage = delta_median < 0,
    significant_lineage_dependency = FDR_global < 0.05 & delta_median < 0
  ) %>%
  group_by(gene) %>%
  mutate(
    FDR_within_gene = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup()
write.csv(
  lineage_test,
  file.path(outdir, "PPMS_46_DepMap_lineage_specific_dependency_tests.csv"),
  row.names = FALSE
)

lineage_test_sig <- lineage_test %>%
  filter(significant_lineage_dependency) %>%
  arrange(gene, FDR_global, delta_median)
write.csv(
  lineage_test_sig,
  file.path(outdir, "PPMS_46_DepMap_significant_lineage_biased_dependencies.csv"),
  row.names = FALSE
)

top_dependent_lineages <- lineage_summary %>%
  group_by(gene) %>%
  arrange(median_gene_effect, .by_group = TRUE) %>%
  mutate(rank_dependency = row_number()) %>%
  filter(rank_dependency <= 5) %>%
  ungroup()
write.csv(
  top_dependent_lineages,
  file.path(outdir, "PPMS_46_top5_dependent_lineages_per_gene.csv"),
  row.names = FALSE
)

# Heatmap: 46 genes × cancer lineages
heatmap_df <- lineage_summary %>%
  dplyr::select(gene, lineage, median_gene_effect) %>%
  pivot_wider(
    names_from = lineage,
    values_from = median_gene_effect
  ) %>%
  as.data.frame()

rownames(heatmap_df) <- heatmap_df$gene
heatmap_df$gene <- NULL
heatmap_mat <- as.matrix(heatmap_df)

mat <- heatmap_mat
min_val <- min(mat, na.rm = TRUE)
max_val <- max(mat, na.rm = TRUE)
cutoff  <- -0.5
n1 <- 60  # < -0.5
n2 <- 60  # > -0.5
cols_low <- colorRampPalette(c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF"))(n1)
cols_high <- colorRampPalette(c("#FFF5EB", "#FDBE85", "#FD8D3C", "#D94701"))(n2)
my_cols <- c(cols_low, cols_high)
breaks_low <- seq(min_val, cutoff, length.out = n1 + 1)
breaks_high <- seq(cutoff, max_val, length.out = n2 + 1)
my_breaks <- c(breaks_low, breaks_high[-1])
pdf(
  file.path(outdir, "PPMS_46_DepMap_lineage_dependency_heatmap_cutoff_neg0.5.pdf"),
  width = 5.76,
  height = 8.5
)
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  na_col = "grey90",
  color = my_cols,
  breaks = my_breaks,
  main = "Cancer-lineage-specific dependency of PPMS candidate targets",
  fontsize_row = 7,
  fontsize_col = 8,
  border_color = NA
)
dev.off()

# Boxplot: selected top genes across lineages
if (length(top_genes) > 0) {
  major_lineages <- meta_clean %>%
    count(lineage, name = "n_cell_lines") %>%
    filter(n_cell_lines >= 10) %>%
    arrange(desc(n_cell_lines)) %>%
    pull(lineage)
  plot_df <- dep_long %>%
    filter(
      gene %in% top_genes,
      lineage %in% major_lineages
    ) %>%
    group_by(gene) %>%
    mutate(
      lineage = factor(
        lineage,
        levels = lineage_summary %>%
          filter(gene == unique(gene)[1], lineage %in% major_lineages) %>%
          arrange(median_gene_effect) %>%
          pull(lineage)
      )
    ) %>%
    ungroup()
  p_box <- ggplot(
    plot_df,
    aes(x = lineage, y = gene_effect)
  ) +
    geom_boxplot(outlier.shape = NA, width = 0.65) +
    geom_jitter(width = 0.15, size = 0.4, alpha = 0.35) +
    facet_wrap(~ gene, ncol = 1, scales = "free_y") +
    geom_hline(yintercept = dependency_cutoff, linetype = 2, linewidth = 0.3) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
      strip.text = element_text(face = "bold", size = 10)
    ) +
    labs(
      x = "Cancer lineage",
      y = "Chronos gene effect score"
    )
  ggsave(
    file.path(outdir, "PPMS_top5_DepMap_lineage_dependency_boxplot.pdf"),
    p_box,
    width = 14,
    height = 12
  )
}

# Summary table for selected top genes
top_gene_summary <- lineage_summary %>%
  filter(gene %in% top_genes) %>%
  arrange(gene, median_gene_effect)
write.csv(
  top_gene_summary,
  file.path(outdir, "PPMS_top5_DepMap_lineage_dependency_summary.csv"),
  row.names = FALSE
)

## compare target expression in tumor versus normal tissues
setwd("PPMS_DepMap_lineage_dependency/")
library(TCGAplot)
library(ggplot2)
library(pheatmap)
library(tidyverse)

dt <- get_all_paired_tpm()
df <- read.csv("PPMS_46_top5_dependent_lineages_per_gene.csv", header = T, stringsAsFactors = F)
genes <- unique(df$gene)

# input
outdir <- "PPMS_46_target_expression_PT_vs_NT"
dir.create(outdir, showWarnings = FALSE)
target_genes <- genes
available_genes <- intersect(target_genes, colnames(dt))
missing_genes <- setdiff(target_genes, colnames(dt))
message("Available genes: ", length(available_genes), "/", length(target_genes))
message("Missing genes: ", paste(missing_genes, collapse = ", "))

# long format
expr_long <- dt %>%
  as.data.frame() %>%
  mutate(
    Group = as.character(Group),
    Group = recode(Group, "Normal" = "NT", "Tumor" = "PT"),
    Group = factor(Group, levels = c("NT", "PT"))
  ) %>%
  filter(Group %in% c("NT", "PT")) %>%
  select(Cancer, Group, ID, all_of(available_genes)) %>%
  pivot_longer(
    cols = all_of(available_genes),
    names_to = "gene",
    values_to = "expr"
  ) %>%
  mutate(expr = as.numeric(expr))

# paired wide format
paired_wide <- expr_long %>%
  group_by(Cancer, gene, ID, Group) %>%
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = Group,
    values_from = expr
  ) %>%
  filter(!is.na(NT), !is.na(PT)) %>%
  mutate(
    diff_PT_minus_NT = PT - NT
  )
write.csv(
  paired_wide,
  file.path(outdir, "PPMS_46_target_expression_paired_wide.csv"),
  row.names = FALSE
)

# paired Wilcoxon test per cancer and gene
paired_stat_fun <- function(x) {
  x <- x %>%
    filter(!is.na(NT), !is.na(PT))
  n_pair <- nrow(x)
  if (n_pair < 3) {
    return(data.frame(
      n_pairs = n_pair,
      median_NT = median(x$NT, na.rm = TRUE),
      median_PT = median(x$PT, na.rm = TRUE),
      mean_NT = mean(x$NT, na.rm = TRUE),
      mean_PT = mean(x$PT, na.rm = TRUE),
      median_diff_PT_minus_NT = median(x$PT - x$NT, na.rm = TRUE),
      mean_diff_PT_minus_NT = mean(x$PT - x$NT, na.rm = TRUE),
      p_value = NA_real_
    ))
  }
  wt <- wilcox.test(
    x$PT,
    x$NT,
    paired = TRUE,
    alternative = "two.sided",
    exact = FALSE
  )
  data.frame(
    n_pairs = n_pair,
    median_NT = median(x$NT, na.rm = TRUE),
    median_PT = median(x$PT, na.rm = TRUE),
    mean_NT = mean(x$NT, na.rm = TRUE),
    mean_PT = mean(x$PT, na.rm = TRUE),
    median_diff_PT_minus_NT = median(x$PT - x$NT, na.rm = TRUE),
    mean_diff_PT_minus_NT = mean(x$PT - x$NT, na.rm = TRUE),
    p_value = wt$p.value
  )
}

expr_stats <- paired_wide %>%
  group_by(Cancer, gene) %>%
  group_modify(~ paired_stat_fun(.x)) %>%
  ungroup() %>%
  mutate(
    FDR_global = p.adjust(p_value, method = "BH"),
    direction = case_when(
      median_diff_PT_minus_NT > 0 ~ "Higher in PT",
      median_diff_PT_minus_NT < 0 ~ "Higher in NT",
      TRUE ~ "No difference"
    ),
    significance = case_when(
      FDR_global < 0.001 ~ "***",
      FDR_global < 0.01  ~ "**",
      FDR_global < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  group_by(gene) %>%
  mutate(FDR_within_gene = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write.csv(
  expr_stats,
  file.path(outdir, "PPMS_46_target_expression_PT_vs_NT_statistics.csv"),
  row.names = FALSE
)

# pan-cancer paired comparison per gene
pan_stat <- paired_wide %>%
  group_by(gene) %>%
  group_modify(~ paired_stat_fun(.x)) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    direction = case_when(
      median_diff_PT_minus_NT > 0 ~ "Higher in PT",
      median_diff_PT_minus_NT < 0 ~ "Higher in NT",
      TRUE ~ "No difference"
    ),
    significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(FDR)

write.csv(
  pan_stat,
  file.path(outdir, "PPMS_46_target_expression_pan_cancer_paired_statistics.csv"),
  row.names = FALSE
)

# heatmap of median paired expression difference
heatmap_df <- expr_stats %>%
  select(gene, Cancer, median_diff_PT_minus_NT) %>%
  pivot_wider(
    names_from = Cancer,
    values_from = median_diff_PT_minus_NT
  ) %>%
  as.data.frame()
rownames(heatmap_df) <- heatmap_df$gene
heatmap_df$gene <- NULL
heatmap_mat <- as.matrix(heatmap_df)

# significance matrix
sig_df <- expr_stats %>%
  select(gene, Cancer, significance) %>%
  pivot_wider(
    names_from = Cancer,
    values_from = significance
  ) %>%
  as.data.frame()

rownames(sig_df) <- sig_df$gene
sig_df$gene <- NULL

sig_mat <- as.matrix(sig_df)
sig_mat <- sig_mat[rownames(heatmap_mat), colnames(heatmap_mat)]

# gene order: put selected top genes first if present
top_targets <- intersect(c("RRM1", "PLK1", "CDC45", "RRM2", "SPC24"), rownames(heatmap_mat))
gene_order <- c(top_targets, setdiff(rownames(heatmap_mat), top_targets))
heatmap_mat <- heatmap_mat[gene_order, , drop = FALSE]
sig_mat <- sig_mat[gene_order, , drop = FALSE]

# color centered at 0
max_abs <- max(abs(heatmap_mat), na.rm = TRUE)

hm_cols <- colorRampPalette(c(
  "#2166AC", "white", "#B2182B"
))(100)

hm_breaks <- seq(-max_abs, max_abs, length.out = 101)

pdf(
  file.path(outdir, "PPMS_46_target_expression_PT_vs_NT_heatmap.pdf"),
  width = 5,
  height = 9
)
pheatmap(
  heatmap_mat,
  color = hm_cols,
  breaks = hm_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = sig_mat,
  number_color = "black",
  fontsize_number = 9,
  fontsize_row = 7,
  fontsize_col = 9,
  border_color = NA,
  na_col = "grey90",
  main = "PT-versus-NT expression differences of PPMS candidate targets"
)
dev.off()

# paired boxplots for selected top targets
plot_genes <- intersect(c("RRM1", "PLK1", "CDC45", "RRM2", "SPC24"), available_genes)

plot_df <- expr_long %>%
  filter(gene %in% plot_genes)

label_df <- expr_stats %>%
  filter(gene %in% plot_genes) %>%
  left_join(
    plot_df %>%
      group_by(Cancer, gene) %>%
      summarise(
        y_pos = max(expr, na.rm = TRUE) + 0.08 * diff(range(expr, na.rm = TRUE)),
        .groups = "drop"
      ),
    by = c("Cancer", "gene")
  ) %>%
  mutate(
    label = case_when(
      FDR_global < 0.001 ~ "***",
      FDR_global < 0.01  ~ "**",
      FDR_global < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

p_top <- ggplot(
  plot_df,
  aes(x = Group, y = expr, fill = Group)
) +
  geom_violin() +
  # geom_line(
  #   aes(group = interaction(Cancer, ID)),
  #   color = "grey75",
  #   linewidth = 0.25,
  #   alpha = 0.5
  # ) +
  # geom_boxplot(
  #   width = 0.55,
  #   outlier.shape = NA,
  #   alpha = 0.85
  # ) +
  # geom_jitter(
  #   width = 0.08,
  #   size = 0.45,
  #   alpha = 0.45
  # ) +
  geom_text(
    data = label_df,
    aes(x = 1.5, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  facet_grid(gene ~ Cancer, scales = "free_y") +
  scale_fill_manual(values = c("NT" = "#8DCDD5", "PT" = "#E6846D")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = NULL,
    y = "Expression"
  )

ggsave(
  file.path(outdir, "PPMS_top5_target_expression_PT_vs_NT_paired_boxplot.pdf"),
  p_top,
  width = 10.5,
  height = 7.5
)
p_top