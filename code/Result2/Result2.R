# Result 2. Coordinated cellular programs across the tumor-normal ecosystem ----
setwd("/data1/yang/test")
# obj <- read_rds("sce.all_int.rds")
# obj <- subset(obj, celltype != "non-TIME"); all_seu <- obj; rm(obj); gc()
# load("T_5.0.RData"); all_t <- obj; rm(obj); gc()
# load("B_5.0.RData"); all_b <- obj; rm(obj); gc()
# load("Myeloid_5.0.RData"); all_m <- obj; rm(obj); gc()
# load("Stromal_5.0.RData"); all_s <- obj; rm(obj); gc()
# load("nonTIME_1.0.RData"); nonTIME <- obj; rm(obj); gc()
# 
# meta <- all_seu@meta.data %>%
#   mutate(cellid = rownames(.))
# meta_t <- all_t@meta.data %>%
#   mutate(cellid = rownames(.)) %>%
#   select(cellid, majorcell, minorcell)
# meta_b <- all_b@meta.data %>%
#   mutate(cellid = rownames(.)) %>%
#   select(cellid, majorcell, minorcell)
# meta_m <- all_m@meta.data %>%
#   mutate(cellid = rownames(.)) %>%
#   select(cellid, majorcell, minorcell)
# meta_s <- all_s@meta.data %>%
#   mutate(cellid = rownames(.)) %>%
#   select(cellid, majorcell, minorcell)
# 
# meta_all <- rbind(meta_t, meta_b, meta_m, meta_s)
# meta_int <- inner_join(meta, meta_all, by = "cellid")
# rownames(meta_int) <- meta_int$cellid
# TIME_seu <- AddMetaData(all_seu, metadata = meta_int)
# saveRDS(TIME_seu, file = "TIME.rds")


## 1. Figure2A/FigureS4B. NMF ----
library(qs)
TIME <- qread("TIME.qs")
cell_counts <- TIME@meta.data %>%
  group_by(sample, majorcell, minorcell) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  na.omit()
total_majorcell_counts <- TIME@meta.data %>%
  group_by(sample, majorcell) %>%
  summarise(total_majorcell_count = n(), .groups = 'drop') %>%
  na.omit()
abundance_matrix <- cell_counts %>%
  left_join(total_majorcell_counts, by = c("sample", "majorcell")) %>%
  mutate(proportion = cell_count / total_majorcell_count) %>%
  select(sample, majorcell, minorcell, proportion)
V <- abundance_matrix %>%
  tidyr::pivot_wider(names_from = minorcell, values_from = proportion, values_fill = list(proportion = 0))
V <- as.data.frame(V)
V <- V %>% dplyr::select(-majorcell)
dt <- matrix(rep(0, 104*391), nrow = 104, ncol = 391)
dt <- as.data.frame(dt)
rownames(dt) <- colnames(V)[-1]
colnames(dt) <- unique(V$sample)
update_data2 <- function(data1, data2) {
  for (i in 1:nrow(data1)) {
    sample_value <- data1$sample[i]
    for (j in 2:ncol(data1)) {
      value <- data1[i, j]
      if (value > 0 && sample_value %in% colnames(data2) && colnames(data1)[j] %in% rownames(data2)) {
        data2[rownames(data2) == colnames(data1)[j], colnames(data2) == sample_value] <- value
      }
    }
  }
  return(data2)
}
data3 <- update_data2(V, dt)
mtx <- as.matrix(data3)
ranks <- 3:8  
nrun <- 20   
nmf.res = nmf(mtx, rank = ranks, nrun = 20, seed = 123456, .options = "v")
saveRDS(
  nmf.res,
  "nmf.res.rds"
)

# plot
cophenetic_coeffs <- numeric(length(ranks))
for (i in seq_along(ranks)) {
  rank <- ranks[i]
  nmf_results <- nmf(mtx, rank = rank, nrun = nrun, .options = "v", seed = 123456)
  cophenetic_coeffs[i] <- cophcor(nmf_results)
}
cophenetic_df <- data.frame(FactorizationRank = ranks, CopheneticCoefficient = cophenetic_coeffs)
ggplot(cophenetic_df, aes(x = FactorizationRank, y = CopheneticCoefficient)) +
  geom_point(color = "purple", size = 3) +
  geom_line(color = "purple") +
  labs(title = "Cophenetic correlation survey", x = "Factorization rank", y = "Coefficient") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16))
ggsave("FigureS4B_Coefficient.pdf", width = 20, height = 15, units = "cm")

# k = 4
nmf.res = nmf(mtx, rank = 4, nrun = 20, seed = 123456, .options = "v")
w <- basis(nmf.res)
colnames(w) <- paste0('NMF', seq(1, 4))
w_ <- w %>% t %>% scale
p <- ComplexHeatmap::Heatmap(w_, width = 11, height = 4, column_km = 4, 
                             col = colorRamp2(c(-1,0,2), c("#415288", "white", "#BD4342")),
                             cluster_rows = F, clustering_method_columns = 'single') 
heatmap <- draw(p)
pdf("Figure2A_top.pdf", width = 19, height = 6)
heatmap
dev.off()
save(cophenetic_df, w_, file = "Figure2A_NMF.RData")


## 2. Figure2B. NMF score ----
# marker genes of NMF
TIME$new_NMF <- factor(paste(TIME$group, TIME$NMF, sep = "_"))
Idents(TIME) <- TIME$new_NMF
table(TIME$new_NMF)
marker_condition <- data.frame()
for ( ci in sort(as.character(unique(TIME$NMF))) ) {
  tmp.marker <- FindMarkers(TIME, 
                            logfc.threshold = 0.25, min.pct = 0.1,
                            only.pos = F, test.use = "wilcox",
                            ident.1 = paste0("PT_", ci), ident.2 = paste0("NT_", ci)
  )
  tmp.marker$gene = rownames(tmp.marker)
  tmp.marker$condition = ifelse(tmp.marker$avg_log2FC > 0,'PT','NT')
  tmp.marker$cluster = ci
  tmp.marker = tmp.marker %>% filter(p_val_adj < 0.05)
  tmp.marker = as.data.frame(tmp.marker)
  tmp.marker = tmp.marker %>% arrange(desc(avg_log2FC))
  marker_condition = marker_condition %>% rbind(tmp.marker)
}
write.csv(marker_condition, file = "Figure2C_marker_condition.csv")
NMF_gene.list <- split(marker_condition$gene, marker_condition$cluster)
save(NMF_gene.list, file = "NMF_gene.list.RData")

# ST
library(SpaCET)
visiumPath <- '~/data1/yang/test/ST'
# Create SpaCET object
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)
str(SpaCET_obj)
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes=1)
SpaCET.visualize.spatialFeature( 
  SpaCET_obj,  
  spatialType = "QualityControl",  
  spatialFeatures = c("UMI","Gene"), 
  imageBg = TRUE)
# To perform deconvolution on the ST data, the type of tumor needs to be specified.
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj,           
                                   cancerType = "LGG",           
                                   coreNo = 1)
myST = SpaCET.GeneSetScore(SpaCET_obj = SpaCET_obj, GeneSets = NMF_gene.list)
myST@results$GeneSetScore[1:4, 1:6]
pdf("Figure2B_LGG.pdf", width = 20, height = 15)
SpaCET.visualize.spatialFeature(
  myST,
  spatialType = "GeneSetScore",
  spatialFeatures = c("NMF1", "NMF2", "NMF3", "NMF4")
)
dev.off()


## 3. Figure2C. NMF marker genes between PT and NT ----
# 1.4 paired pseudobulk DESeq2：PT vs NT ----
suppressPackageStartupMessages({
  library(DESeq2)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(readr)
  library(BiocParallel)
})
register(MulticoreParam(workers = 20))

run_deseq2_one_nmf <- function(pb_item, nmf_i,
                               min_count = 10,
                               min_pb_samples = 2,
                               lfc_cutoff = 0.25) {
  
  message("Running paired pseudobulk DESeq2 for ", nmf_i, " ...")
  pb_counts <- pb_item$counts
  coldata <- as.data.frame(pb_item$coldata)
  
  coldata <- coldata[match(colnames(pb_counts), coldata$pb_id), , drop = FALSE]
  if (any(is.na(coldata$pb_id))) {
    stop(nmf_i, ": coldata matching failed. Some pb_id values are NA.")
  }
  rownames(coldata) <- coldata$pb_id
  stopifnot(identical(colnames(pb_counts), rownames(coldata)))
  
  coldata$patient_pb <- factor(coldata$patient_pb)
  coldata$group_raw <- factor(coldata$group_raw, levels = c("NT", "PT"))
  
  message(
    nmf_i, ": ",
    ncol(pb_counts), " pseudobulk samples; ",
    dplyr::n_distinct(coldata$patient_pb), " paired patients before gene filtering."
  )
  
  keep_gene <- Matrix::rowSums(pb_counts >= min_count) >= min_pb_samples &
    Matrix::rowSums(pb_counts) > 0
  
  pb_counts_use <- pb_counts[keep_gene, , drop = FALSE]
  
  message(
    nmf_i, ": ",
    nrow(pb_counts_use), " genes retained; ",
    ncol(pb_counts_use), " pseudobulk samples; ",
    dplyr::n_distinct(coldata$patient_pb), " paired patients."
  )
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(pb_counts_use)),
    colData = coldata,
    design = ~ patient_pb + group_raw
  )
  
  dds <- DESeq(dds, quiet = FALSE, parallel = TRUE)
  
  res <- results(
    dds,
    contrast = c("group_raw", "PT", "NT"),
    alpha = 0.05,
    parallel = TRUE
  )
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      NMF_program = nmf_i,
      avg_log2FC = log2FoldChange,
      p_val = pvalue,
      p_val_adj = padj,
      condition = dplyr::case_when(
        !is.na(p_val_adj) & p_val_adj < 0.05 & avg_log2FC > 0 ~ "PT",
        !is.na(p_val_adj) & p_val_adj < 0.05 & avg_log2FC < 0 ~ "NT",
        TRUE ~ "NS"
      ),
      direction = dplyr::case_when(
        condition == "PT" ~ "PT-up",
        condition == "NT" ~ "NT-up",
        TRUE ~ "NS"
      ),
      significant_FDR_0.05 = !is.na(p_val_adj) & p_val_adj < 0.05,
      significant_FDR_0.05_logFC_0.25 =
        !is.na(p_val_adj) & p_val_adj < 0.05 & abs(avg_log2FC) > lfc_cutoff,
      n_paired_patients = dplyr::n_distinct(coldata$patient_pb),
      n_pseudobulk_samples = nrow(coldata),
      method = "patient-level pseudobulk paired DESeq2; design = ~ patient_pb + group_raw; contrast = PT vs NT"
    ) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))
  
  summary_df <- tibble::tibble(
    NMF_program = nmf_i,
    n_paired_patients = dplyr::n_distinct(coldata$patient_pb),
    n_pseudobulk_samples = nrow(coldata),
    n_tested_genes = nrow(res_df),
    n_sig_FDR_0.05 = sum(res_df$significant_FDR_0.05, na.rm = TRUE),
    n_sig_FDR_0.05_logFC_0.25 = sum(res_df$significant_FDR_0.05_logFC_0.25, na.rm = TRUE),
    n_PT_up_FDR_0.05_logFC_0.25 = sum(
      res_df$p_val_adj < 0.05 & res_df$avg_log2FC > lfc_cutoff,
      na.rm = TRUE
    ),
    n_NT_up_FDR_0.05_logFC_0.25 = sum(
      res_df$p_val_adj < 0.05 & res_df$avg_log2FC < -lfc_cutoff,
      na.rm = TRUE
    )
  )
  
  list(
    res = res_df,
    summary = summary_df,
    dds = dds
  )
}

nmf_vec <- c("NMF1", "NMF2", "NMF3", "NMF4")
for (nmf_i in nmf_vec) {
  test_i <- run_deseq2_one_nmf(
    pb_count_list[[nmf_i]],
    nmf_i
  )
  result_i <- test_i$res
  write.csv(
    result_i,
    file = paste0(
      "/data3/home/yang/GI/test/SecondV/",
      tolower(nmf_i),
      "_pseudobulk_DEGs.csv"
    ),
    row.names = FALSE
  )
  save(
    test_i,
    file = paste0(
      "/data3/home/yang/GI/test/SecondV/",
      tolower(nmf_i),
      "_pseudobulk_DEGs.RData"
    )
  )
  message(nmf_i, " finished and saved.")
}

# plot
library(tidyverse)
library(ggrepel)
library(HGNChelper)

path <- "/data3/home/yang/GI/test/SecondV"

dt <- list.files(path, pattern = "nmf[1-4]_pseudobulk_DEGs\\.csv$", full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  mutate(
    NMF_program = factor(NMF_program, levels = paste0("NMF", 1:4)),
    change = ifelse(log2FoldChange > 0, "sigUp", "sigDown")
  )

dt <- dt %>%
  filter(!grepl("^Tissue:", gene))

table(dt$NMF_program, dt$direction)

set.seed(123)

dt <- dt %>%
  mutate(
    x = as.numeric(NMF_program),
    x_jitter = x + runif(n(), -0.28, 0.28)
  )

hgnc <- checkGeneSymbols(unique(dt$gene))
approved_genes <- hgnc$x[hgnc$Approved]
lab <- dt %>%
  filter(gene %in% approved_genes) %>%
  group_by(NMF_program, change) %>%
  slice_max(abs(log2FoldChange), n = 5, with_ties = FALSE) %>%
  ungroup()

p <- ggplot(dt, aes(x_jitter, log2FoldChange, color = change)) +
  geom_point(size = 1.3, alpha = 1) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text_repel(
    data = lab,
    aes(x = x_jitter, y = log2FoldChange, label = gene),
    size = 3, color = "black",
    box.padding = 0.3,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_x_continuous(
    breaks = 1:4,
    labels = paste0("NMF", 1:4)
  ) +
  scale_color_manual(values = c(sigUp = "#e6846d", sigDown = "#8dcdd5")) +
  labs(x = NULL, y = "log2 fold change (PT vs. NT)", color = "PT vs. NT") +
  theme_classic() +
  theme(
    axis.text.x = element_text(face = "bold"),
    legend.position = "right"
  )
p
ggsave(file.path(path, "Fig2C_pseudobulk_DEGs.pdf"), p, width = 8, height = 5)

## 4. Figure2D. KEGG pathway analysis of NMF marker genes ----
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

path <- "/data3/home/yang/GI/test/SecondV"

dt <- list.files(path, pattern = "nmf[1-4]_pseudobulk_DEGs\\.csv$", full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  filter(!grepl("^Tissue:", gene))

enrich_nmf <- function(dat, nmf, logfc.cut = 0.5, top_n = 5) {
  
  x <- dat %>% filter(NMF_program == nmf)
  bg <- bitr(unique(x$gene), "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  
  sig <- x %>%
    filter(padj < 0.05, abs(log2FoldChange) > logfc.cut) %>%
    mutate(Group = ifelse(log2FoldChange > 0, "PT", "NT"))
  
  res <- map_dfr(c("PT", "NT"), function(g) {
    genes <- sig %>% filter(Group == g) %>% pull(gene)
    ids <- bitr(genes, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
    
    enrichGO(
      gene = ids, universe = bg, OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID", ont = "BP",
      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
      readable = TRUE
    )@result %>%
      filter(p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      # slice_head(n = top_n) %>%
      mutate(Group = g)
  })
  
  res <- res %>%
    mutate(
      score = -log10(p.adjust) * ifelse(Group == "PT", 1, -1),
      term = paste0(Description, "_", Group),
      Group = factor(Group, levels = c("NT", "PT"))
    ) %>%
    arrange(Group, score) %>%
    mutate(term = factor(term, levels = unique(term)))
  
  p <- ggplot(res, aes(score, term, fill = Group)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.25) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_y_discrete(labels = setNames(res$Description, res$term)) +
    scale_fill_manual(values = c(NT = "#8dcdd5", PT = "#e6846d")) +
    labs(
      title = nmf,
      x = "Signed -log10(FDR)",
      y = NULL,
      fill = "PT vs. NT"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9),
      legend.position = "top"
    )
  
  list(enrichment = res, plot = p)
}

res_NMF1 <- enrich_nmf(dt, "NMF1")
res_NMF2 <- enrich_nmf(dt, "NMF2")
res_NMF3 <- enrich_nmf(dt, "NMF3")
res_NMF4 <- enrich_nmf(dt, "NMF4")

nmf <- res_NMF4$enrichment
write.csv(nmf, file = "NMF4_GO.csv")

# p_all <- (res_NMF1$plot | res_NMF2$plot) /
#   (res_NMF3$plot | res_NMF4$plot)
# 
# p_all
# 
# ggsave(
#   file.path(path, "Fig2D_pseudobulk_GO.pdf"),
#   p_all, width = 14, height = 10
# )

dt <- res_NMF3$enrichment 
dt <- dt %>%
  mutate(
    Cluster = factor(as.character(Group), levels = c("NT", "PT")),
    Description = str_wrap(Description, width = 35)
  ) %>%
  arrange(Cluster, p.adjust)

dt$term <- factor(dt$term, levels = rev(dt$term))

dt$lable <- ifelse(dt$Cluster == "PT", 1, -1)
dt$pvalue_loc <- -log10(dt$p.adjust) * dt$lable
dt$lable_hjust <- ifelse(dt$Cluster == "PT", 1, 0)
dt$lable_xloc <- ifelse(dt$Cluster == "PT", -0.5, 0.5)

dt$label_color <- "black"
dt$label_color[dt$Cluster == "PT"][1:min(3, sum(dt$Cluster == "PT"))] <- "red"
dt$label_color[dt$Cluster == "NT"][1:min(3, sum(dt$Cluster == "NT"))] <- "red"

color <- c(NT = "#8dcdd5", PT = "#e6846d")
xmax <- ceiling(max(abs(dt$pvalue_loc))) + 1

p <- ggplot(dt, aes(x = pvalue_loc, y = term, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_manual(values = color) +
  scale_x_continuous(limits = c(-xmax, xmax)) +
  scale_y_discrete(labels = setNames(dt$Description, dt$term)) +
  geom_text(
    aes(x = lable_xloc, y = term, label = Description, hjust = lable_hjust),
    color = dt$label_color, size = 4.5, lineheight = 0.9
  ) +
  labs(x = expression(-log[10]~"(FDR)"), y = NULL, title = "") +
  # annotate("text", x = xmax * 0.75, y = length(levels(dt$term)) - 1, label = "PT", size = 6, fontface = "bold") +
  # annotate("text", x = -xmax * 0.75, y = 2, label = "NT", size = 6, fontface = "bold") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

p
ggsave("NMF4_GO.pdf", plot = p, width = 4.5, height = 4.5)


## 5. Figure2E. NMF survial analysis ----
TIME <- qread("TIME.qs")
pt_time <- subset(TIME, group == "PT")
celltypes <- levels(Idents(pt_time))

## calculate DEGs
dir_for_data <- "/data1/yang/test"
dir.create(file.path(dir_for_data, "01.NMF_subset_signature_in_TIME"))
cl <- makeCluster(20)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths(c("/data3/home/yang/R/x86_64-pc-linux-gnu-library/4.4", "/data3/software/R/4.4.1/lib64/R/library")))
foreach(
  ident.1 = celltypes,
  .export = c("pt_time", "dir_for_data"),
  .packages = c("Seurat", "readr", "tidyverse")
) %dopar% {
  markers <- Seurat::FindMarkers(pt_time, ident.1 = ident.1, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox") %>%
    rownames_to_column("gene") %>%
    arrange(desc(avg_log2FC))
  write_tsv(markers, file.path(dir_for_data, "01.NMF_subset_signature_in_TIME", paste0(ident.1, "_wilcoxon_DE_genes.tsv")))
}
stopCluster(cl)


library(GSVA)
library(TCGAplot)
library(tidyverse)
library(survival)
library(meta)
library(doParallel)  
## 
registerDoParallel(cores = 20)
tpm <- get_all_tpm()
meta <- get_all_meta()
meta <- meta[meta$Cancer != "TGCT", ]
cancers <- unique(meta$Cancer)
dir_for_signatures <- file.path(dir_for_data, "Figure2F.NMF_subset_signature_in_TIME")
remove_cell_cycle_genes <- FALSE
top_n_DE_genes <- 20
result_df <- list()
for (i_file in list.files(dir_for_signatures)) {
  # get signature with cell cycle genes removed
  markers_1 <- read_tsv(file.path(dir_for_signatures, i_file))
  if (remove_cell_cycle_genes) {
    markers_1 <- markers_1 %>%
      dplyr::filter(
        p_val_adj < 0.05,
        avg_log2FC > 1,
        gene %in% colnames(tpm)[-c(1, 2)],
        !gene %in% cell_cycle_genes
      ) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = top_n_DE_genes)
  } else {
    markers_1 <- markers_1 %>%
      dplyr::filter(
        p_val_adj < 0.05,
        avg_log2FC > 1,
        gene %in% colnames(tpm)[-c(1, 2)]
      ) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = top_n_DE_genes)
  }
  markers_1$celltype <- str_extract(i_file, "^.*(?=_wilcoxon)")
  result_df[[i_file]] <- markers_1
}
genelist <- bind_rows(result_df)
genelist <- split(genelist$gene, genelist$celltype)
# foreach
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
write.csv(final_pan, file = "TCGAplot_NMF_HR.csv")

# plot pan-cancer survival for each B subset
plot_df <- final_pan[final_pan$celltype == "NMF3", ]
# plot_df$HR <- round(plot_df$HR, digits = 2)
plot_df <- plot_df %>%
  arrange(HR)
plot_df$cancer <- factor(plot_df$cancer, levels = c(setdiff(plot_df$cancer, "Pan cancer"), "Pan cancer"))

## color
plot_df$color <- as.character(plot_df$color)
plot_df$color[plot_df$p >= 0.05 & plot_df$HR < 1] <- "Better survival (P > 0.05)"
plot_df$color[plot_df$p >= 0.05 & plot_df$HR > 1] <- "Worse survival (P > 0.05)"
plot_df$color <- factor(plot_df$color, levels = c(
  "Better survival", "Better survival (P > 0.05)",
  "Worse survival (P > 0.05)", "Worse survival"
))
mycolor <- c(
  "Better survival" = "#0F7B9F",
  "Better survival (P > 0.05)" = "#E0F3F8FF",
  "Worse survival (P > 0.05)" = "#FDDBC7FF",
  "Worse survival" = "#C3423F"
)
## plot
## plot
ggplot(data = plot_df, aes(x = cancer, y = HR, ymin = lower, ymax = upper, color = color)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, colour = "grey40", linetype = "dashed", size = 0.2) + # add a dotted line at x=1 after flip
  geom_text(aes(x = cancer, y = max(upper), label = significance), show.legend = FALSE) +
  scale_color_manual(values = mycolor, name = "Survival association") +
  xlab("") +
  ylab("Hazard ratio (95% CI)") +
  cowplot::theme_cowplot() +
  theme(
    text = element_text(size = 10, family = "ArialMT"),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 8),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(linetype = 1, color = "black", size = 0.3),
    axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
  ) +
  scale_y_continuous(trans = "log2")
ggsave(
  filename = "Figure2E.survival_NMF3.pdf",
  height = 3.5,
  width = 10
)


## 6. FigureS4A/FigureS4G. Ro/e analysis ----
library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)

meta <- TIME@meta.data
R_oe <- calTissueDist(meta,
                      byPatient = F,
                      colname.cluster = "group",
                      colname.patient = "patient",
                      colname.tissue = "minorcell",
                      method = "chisq", 
                      min.rowSum = 0) 

# R_oe <- calTissueDist(meta,
#                       byPatient = F,
#                       colname.cluster = "group",
#                       colname.patient = "patient",
#                       colname.tissue = "NMF",
#                       method = "chisq", 
#                       min.rowSum = 0) 

R_oe
data <- as.matrix(t(R_oe))
annotation_colors <- c(
  "+++" = "#E4540E",  
  "++" = "#F5904C",    
  "+" = "#FCBE8C",     
  "+/-" = "#FFE7CF",   
  "-" = "#0f86a9"      
)

# cell_fun function
cell_fun <- function(j, i, x, y, width, height, fill) {
  value <- data[i, j]
  if (value > 3) {
    label <- "+++"
  } else if (value > 1.5) {
    label <- "++"
  } else if (value > 1) {
    label <- "+"
  } else if (value >= 0) {
    label <- "+/-"
  } else {
    label <- "-"
  }
  color <- annotation_colors[label]
  grid.rect(x, y, width, height, gp = gpar(fill = color, col = NA))
  grid.text(label, x, y, gp = gpar(fontsize = 13, col = "gray10"))
}

pdf(file="FigureS3A.pdf", width = 4, height = 20)
ht <- Heatmap(
  data,
  show_heatmap_legend = FALSE,  
  cluster_rows = T,
  cluster_columns = T,
  row_names_side = 'right',
  show_column_names = T,
  show_row_names = T,
  col = colorRamp2(c(-2, 3), c("white", "white")),  
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 0,
  cell_fun = cell_fun
)
lgd <- Legend(
  title = "R_o/e",
  labels = names(annotation_colors),
  legend_gp = gpar(fill = annotation_colors),
  labels_gp = gpar(fontsize = 10)
)
draw(ht, heatmap_legend_list = list(lgd))
dev.off()

col_fun <- colorRamp2(
  c(0, 1, 1.5, 3),
  c("white", "#FCBE8C", "#F5904C", "#E4540E")
)
# cell_fun
cell_fun <- function(j, i, x, y, width, height, fill) {
  value <- data[i, j]
  grid.text(
    sprintf("%.2f", value),
    x,
    y,
    gp = gpar(
      fontsize = 10,
      col = ifelse(value > 2, "white", "gray10")
    )
  )
}
pdf(file = "FigureS3A_numeric_gradient.pdf", width = 4, height = 20)
ht <- Heatmap(
  data,
  name = "Ro/e",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_side = "right",
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 0,
  cell_fun = cell_fun,
  heatmap_legend_param = list(
    title = "Ro/e",
    at = c(0, 1, 1.5, 3),
    labels = c("0", "1", "1.5", "3")
  )
)
draw(ht)
dev.off()

# TCGA pan-cancer DEGs
run_tcga_de <- function(cancer_type, 
                        data_dir = "TCGA/33-pan-cancer-TCGA/",
                        logFC_t = 1, 
                        pvalue_t = 0.05) {
  require(data.table)
  require(edgeR)
  require(stringr)
  require(dplyr)
  file_path <- file.path(data_dir, paste0("TCGA-", cancer_type, "_data.csv"))
  exp <- fread(file_path, data.table = FALSE)
  rownames(exp) <- exp$V1
  exp <- dplyr::select(exp, -V1) 
  samples <- colnames(exp)
  sample_info <- data.frame(
    sample = samples,
    type = str_sub(samples, 14, 15)
  )
  valid_samples <- sample_info %>%
    filter(type %in% c("01", "11")) %>%
    pull(sample)
  exp <- exp[, colnames(exp) %in% valid_samples]
  if (ncol(exp) == 0) stop("No valid tumor/normal samples found")
  table(str_sub(colnames(exp), 14, 15)) %>% print()
  Group <- case_when(
    str_sub(colnames(exp), 14, 15) == "01" ~ "tumor",
    str_sub(colnames(exp), 14, 15) == "11" ~ "normal"
  ) %>% factor(levels = c("normal", "tumor"))
  dge <- DGEList(counts = exp, group = Group)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~Group)
  dge <- estimateDisp(dge, design) 
  fit <- glmQLFit(dge, design)     
  qlf <- glmQLFTest(fit)
  DEG <- topTags(qlf, n = Inf, adjust.method = "BH") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(
      change = case_when(
        logFC > logFC_t & FDR < pvalue_t ~ "UP",
        logFC < -logFC_t & FDR < pvalue_t ~ "DOWN",
        TRUE ~ "NOT"
      ),
      cancer = cancer_type
    )
  attr(DEG, "cancer_type") <- cancer_type
  attr(DEG, "sample_counts") <- table(Group)
  return(DEG)
}
tcga_cancers <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
                  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS",
                  "UVM")
combined_degs <- lapply(tcga_cancers, function(ct) {
  tryCatch({
    message("Analyzing ", ct)
    run_tcga_de(ct)
  }, error = function(e) {
    message(paste("Failed for", ct, ":", e$message))
    return(NULL)
  })
}) %>% 
  bind_rows()


## 7. FigureS4C. top two marker genes of NMF ----
DotPlot(TIME, features = genes <- c(
  "SPP1", "COL4A1", 
  "CCR7", "RGS11", 
  "TYROBP", "MS4A6A", 
  "BCL2A1", "TPSAB1" 
)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


## 8. FigureS4D. Leave-one-out analysis of NMF ----
meta <- TIME@meta.data
group_col  <- "group"
sample_col <- "sample"
study_col  <- "study"
cancer_col <- "cancer"
nmf_col    <- "NMF"
normal_label <- "NT"
tumor_label  <- "PT"
nmf_programs <- c("NMF1", "NMF2", "NMF3", "NMF4")
meta2 <- meta %>%
  filter(.data[[group_col]] %in% c(normal_label, tumor_label)) %>%
  filter(.data[[nmf_col]] %in% nmf_programs)
sample_info <- meta2 %>%
  distinct(
    sample = .data[[sample_col]],
    group  = .data[[group_col]],
    study  = .data[[study_col]],
    cancer = .data[[cancer_col]]
  )

# cell-level NMF proportion
overall_nmf_prop_wide <- meta2 %>%
  count(
    group = .data[[group_col]],
    NMF = .data[[nmf_col]],
    name = "cell_number"
  ) %>%
  group_by(group) %>%
  mutate(
    total_cells = sum(cell_number),
    proportion = cell_number / total_cells
  ) %>%
  ungroup() %>%
  select(group, NMF, cell_number, proportion) %>%
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

# sample-level NMF proportion
sample_nmf_prop <- meta2 %>%
  count(
    sample = .data[[sample_col]],
    NMF = .data[[nmf_col]],
    name = "cell_number"
  ) %>%
  right_join(
    expand_grid(sample = sample_info$sample, NMF = nmf_programs),
    by = c("sample", "NMF")
  ) %>%
  mutate(cell_number = ifelse(is.na(cell_number), 0, cell_number)) %>%
  left_join(sample_info, by = "sample") %>%
  group_by(sample) %>%
  mutate(
    total_cells = sum(cell_number),
    proportion = cell_number / total_cells
  ) %>%
  ungroup()

sample_nmf_mean_prop_wide <- sample_nmf_prop %>%
  group_by(group, NMF) %>%
  summarise(
    n_sample = n_distinct(sample),
    mean_prop = mean(proportion, na.rm = TRUE),
    median_prop = median(proportion, na.rm = TRUE),
    sd_prop = sd(proportion, na.rm = TRUE),
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

# leave-one-study-out / leave-one-cancer-out NMF 
calc_nmf_loo_effect <- function(remove_col) {
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
        NMF = .data[[nmf_col]],
        name = "cell_number"
      ) %>%
      right_join(
        expand_grid(sample = dat_sample_info$sample, NMF = nmf_programs),
        by = c("sample", "NMF")
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
      group_by(group, NMF) %>%
      summarise(
        n_sample = n_distinct(sample),
        mean_prop = mean(proportion, na.rm = TRUE),
        median_prop = median(proportion, na.rm = TRUE),
        sd_prop = sd(proportion, na.rm = TRUE),
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
loo_study_nmf_effect  <- calc_nmf_loo_effect(study_col)
loo_cancer_nmf_effect <- calc_nmf_loo_effect(cancer_col)
loo_nmf_effect <- bind_rows(loo_study_nmf_effect, loo_cancer_nmf_effect) %>%
  left_join(
    sample_nmf_mean_prop_wide %>%
      select(NMF, full_log2FC_PT_vs_NT = log2FC_PT_vs_NT),
    by = "NMF"
  ) %>%
  mutate(
    same_direction = sign(log2FC_PT_vs_NT) == sign(full_log2FC_PT_vs_NT),
    delta_log2FC = log2FC_PT_vs_NT - full_log2FC_PT_vs_NT
  )
loo_nmf_summary <- loo_nmf_effect %>%
  group_by(removed_type, NMF) %>%
  summarise(
    full_log2FC_PT_vs_NT = unique(full_log2FC_PT_vs_NT),
    mean_loo_log2FC = mean(log2FC_PT_vs_NT, na.rm = TRUE),
    median_loo_log2FC = median(log2FC_PT_vs_NT, na.rm = TRUE),
    min_loo_log2FC = min(log2FC_PT_vs_NT, na.rm = TRUE),
    max_loo_log2FC = max(log2FC_PT_vs_NT, na.rm = TRUE),
    sd_loo_log2FC = sd(log2FC_PT_vs_NT, na.rm = TRUE),
    max_abs_delta_log2FC = max(abs(delta_log2FC), na.rm = TRUE),
    same_direction_rate = mean(same_direction, na.rm = TRUE),
    n_leave_one = n(),
    robust_direction = ifelse(same_direction_rate >= 0.9, "Yes", "No"),
    .groups = "drop"
  ) %>%
  arrange(removed_type, desc(abs(full_log2FC_PT_vs_NT)))
write.csv(overall_nmf_prop_wide, "01_overall_cell_level_PT_NT_NMF_proportion.csv", row.names = FALSE)
write.csv(sample_nmf_mean_prop_wide, "02_sample_level_PT_NT_NMF_mean_proportion.csv", row.names = FALSE)
write.csv(loo_nmf_effect, "03_leave_one_each_result_NMF.csv", row.names = FALSE)
write.csv(loo_nmf_summary, "04_leave_one_summary_NMF.csv", row.names = FALSE)

# plot
plot_nmf_direction <- loo_nmf_summary %>%
  mutate(
    NMF = factor(NMF, levels = nmf_levels),
    removed_type = recode(
      removed_type,
      cancer = "Leave-one-cancer-type-out",
      study = "Leave-one-study-out"
    )
  )
p1 <- ggplot(plot_nmf_direction, aes(x = NMF, y = same_direction_rate)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0.9, linetype = 2, color = "grey50") +
  facet_wrap(~ removed_type, nrow = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 11)
  ) +
  labs(
    x = NULL,
    y = "Same-direction rate"
  )
p1
ggsave("NMF_Same_direction.pdf", p1, width = 5.27, height = 3.78)


## 9. FigureS4E. Cosine similarity of NMF programs ----
cosine_sim <- function(x, y) {
  sum(x * y) / sqrt(sum(x^2) * sum(y^2))
}
permute_vec <- function(x) {
  if (length(x) == 1) return(list(x))
  unlist(
    lapply(seq_along(x), function(i) {
      lapply(permute_vec(x[-i]), function(y) c(x[i], y))
    }),
    recursive = FALSE
  )
}

match_programs <- function(W_ref, W_boot) {
  cor_mat <- cor(W_ref, W_boot, method = "pearson")
  cos_mat <- matrix(
    NA,
    nrow = ncol(W_ref),
    ncol = ncol(W_boot),
    dimnames = list(colnames(W_ref), colnames(W_boot))
  )
  for (i in seq_len(ncol(W_ref))) {
    for (j in seq_len(ncol(W_boot))) {
      cos_mat[i, j] <- cosine_sim(W_ref[, i], W_boot[, j])
    }
  }
  perms <- permute_vec(seq_len(ncol(W_boot)))
  scores <- sapply(perms, function(p) {
    sum(cor_mat[cbind(seq_len(ncol(W_ref)), p)])
  })
  best_perm <- perms[[which.max(scores)]]
  data.frame(
    ref_program = colnames(W_ref),
    boot_program = colnames(W_boot)[best_perm],
    pearson_cor = cor_mat[cbind(seq_len(ncol(W_ref)), best_perm)],
    cosine_similarity = cos_mat[cbind(seq_len(ncol(W_ref)), best_perm)]
  )
}

set.seed(123456)
nmf.ref <- readRDS("nmf.res.rds")
W_ref <- basis(nmf.ref)
colnames(W_ref) <- paste0("NMF", 1:4)
n_boot <- 100
sample_fraction <- 0.8
boot_results <- map_dfr(seq_len(n_boot), function(b) {
  set.seed(123456 + b)
  boot_samples <- sample(
    colnames(mtx),
    size = floor(ncol(mtx) * sample_fraction),
    replace = TRUE
  )
  mtx_boot <- mtx[, boot_samples, drop = FALSE]
  nmf.boot <- nmf(
    mtx_boot,
    rank = 4,
    nrun = 20,
    method = "brunet",
    seed = 123456 + b,
    .options = "v"
  )
  W_boot <- basis(nmf.boot)
  colnames(W_boot) <- paste0("Boot_NMF", 1:4)
  matched <- match_programs(W_ref, W_boot)
  matched$bootstrap_id <- b
  matched
})
boot_summary <- boot_results %>%
  group_by(ref_program) %>%
  summarise(
    n_boot = n(),
    mean_pearson = mean(pearson_cor),
    median_pearson = median(pearson_cor),
    min_pearson = min(pearson_cor),
    sd_pearson = sd(pearson_cor),
    mean_cosine = mean(cosine_similarity),
    median_cosine = median(cosine_similarity),
    min_cosine = min(cosine_similarity),
    sd_cosine = sd(cosine_similarity),
    success_rate_pearson_0.7 = mean(pearson_cor >= 0.7),
    success_rate_pearson_0.8 = mean(pearson_cor >= 0.8),
    success_rate_cosine_0.8 = mean(cosine_similarity >= 0.8),
    .groups = "drop"
  )
write.csv(boot_results, "NMF_bootstrap_each_iteration.csv", row.names = FALSE)
write.csv(boot_summary, "NMF_bootstrap_summary.csv", row.names = FALSE)

# plot
plot_boot <- boot_results %>%
  select(bootstrap_id, ref_program, pearson_cor, cosine_similarity) %>%
  pivot_longer(
    cols = c(pearson_cor, cosine_similarity),
    names_to = "metric",
    values_to = "similarity"
  ) %>%
  mutate(
    ref_program = factor(ref_program, levels = c("NMF1", "NMF2", "NMF3", "NMF4")),
    metric = recode(
      metric,
      pearson_cor = "Pearson correlation",
      cosine_similarity = "Cosine similarity"
    )
  )
p_boot <- ggplot(plot_boot, aes(x = ref_program, y = similarity)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5) +
  geom_hline(yintercept = 0.7, linetype = 2, color = "grey50") +
  facet_wrap(~ metric, nrow = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 11)
  ) +
  labs(
    x = NULL,
    y = "Similarity to reference program"
  )
pdf("NMF_bootstrap_reproducibility_boxplot.pdf", width = 5.27, height = 3.78)
print(p_boot)
dev.off()


## 10. FigureS4F. NMF marker genes score in ST ----
library(SpaCET)
visiumPath <- 'HNSCC_GSM8633895/'
# Create SpaCET object
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)
str(SpaCET_obj)
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes=1)
SpaCET.visualize.spatialFeature( 
  SpaCET_obj,  
  spatialType = "QualityControl",  
  spatialFeatures=c("UMI","Gene"), 
  imageBg = TRUE)
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj,           
                                   cancerType="HNSC",           
                                   coreNo=1) 
genes <- c(
  "SPP1", "COL4A1", 
  "CCR7", "RGS11", 
  "TYROBP", "MS4A6A", 
  "BCL2A1", "TPSAB1" 
)
feature_cols <- c(
  "#91BED9", 
  "#E7D5A0",  
  "#F5B574",  
  "#D8663E"   
)
pdf("ST/HNSCC_GSM8633895_SpatialFeature.pdf", width=10.81, height=4.03)
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "GeneExpression",
  spatialFeatures = genes,
  colors = feature_cols,
  nrow = 2,
  imageBg = T
)
dev.off()


## 11. FigureS4H. NMF programs in each cancer type ----
TIME$NMF <- factor(TIME$NMF, levels = c("NMF1","NMF2","NMF3","NMF4"))
meta <- TIME@meta.data
meta_NT <- meta[meta$group == "NT", ]
meta_PT <- meta[meta$group == "PT", ]

# NT
ggplot(data = meta_NT) +
  geom_bar(mapping = aes(x = cancer, fill = NMF, ), position = "fill", width = 0.75) +
  scale_fill_manual(values = nmf_colors) + 
  coord_flip() +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  labs(y = "Total NMF proportion (%) \n", x="") +
  theme(axis.title.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black")) +
  theme(axis.text.x = element_text(colour = "black")) +
  scale_y_continuous(expand = c(0,0)) 
ggsave("FigureS4H_NT_cancer_NMF.pdf", width = 13, height = 10, units = "cm")

# PT
ggplot(data = meta_PT) +
  geom_bar(mapping = aes(x = cancer, fill = NMF, ), position = "fill", width = 0.75) +
  scale_fill_manual(values = nmf_colors) + 
  coord_flip() +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  labs(y = "Total NMF proportion (%) \n", x="") +
  theme(axis.title.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black")) +
  theme(axis.text.x = element_text(colour = "black")) +
  scale_y_continuous(expand = c(0,0)) 
ggsave("FigureS4H_PT_cancer_NMF.pdf", width = 13, height = 10, units = "cm")


## 12. FigureS4I. NMF programs and GEs ----
# Jaccard correlation
load("GE_list.Rdata")
pt_TIME <- subset(TIME, group == "PT")
Idents(pt_TIME) <- pt_TIME$NMF
all.markers <- FindAllMarkers(pt_TIME, only.pos = T, logfc.threshold = 0.25, min.pct = 0.1)
pt_nmf_markers <- all.markers[all.markers$avg_log2FC > 1 & all.markers$p_val_adj < 0.05, ]
pt_nmf <- split(pt_nmf_markers$gene, 
                factor(pt_nmf_markers$cluster, levels = paste0("NMF", 1:4)))
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) return(0) 
  return(round(intersection / union, 4))
}
jaccard_matrix <- matrix(nrow = length(GE_list),
                         ncol = length(pt_nmf),
                         dimnames = list(names(GE_list), names(pt_nmf)))
for(i in seq_along(GE_list)){
  for(j in seq_along(pt_nmf)){
    jaccard_matrix[i, j] <- jaccard(GE_list[[i]], pt_nmf[[j]])
  }
}
pdf("Jaccard_NMF_GE.pdf", width = 6, height = 6)
pheatmap(jaccard_matrix,
         col = colorRampPalette(c("white", "darkgreen"))(256),
         cluster_rows = F,
         main = "Jaccard Index Matrix",
         cluster_cols = F)
dev.off()

# boxplot of GEs
library(UCell)
library(ggcorrplot)
library(ggthemes)
library(ggplot2)
GE.color <- c("#9fd4ca","#f8f6bf","#c9c7de","#f1958a","#95bbd4","#f6bf7f","#bbd985","#f8d5e5","#b3ddf6","#d2e8cd")

GE.gene <- readxl::read_xlsx("cancer_GE.xlsx")
GE.gene$GE <- paste0("GE", GE.gene$GE)
GE_list <- split(GE.gene$gene, 
                 factor(GE.gene$GE, levels = paste0("GE", 1:10)))
pt_TIME <- subset(TIME, group == "PT")
score <- AddModuleScore_UCell(pt_TIME, features = GE_list)
raw.data <- as.matrix(score@assays$RNA$counts)

colnames(raw.data) <- score$NMF
GEs <- NULL
for (var in paste0("GE", 1:10)) {
  ucell <- paste0(var, "_UCell")
  b = data.frame(expression = score@meta.data %>% dplyr::pull(ucell),
                 group = colnames(raw.data))
  data = as.data.frame(b)
  data$group = factor(data$group, levels = names(table(pt_TIME$NMF)))
  data = na.omit(data)
  data$GE <- rep(var, dim(data)[1])
  GEs <- rbind(GEs, data)
}
df1 <- GEs[GEs$GE %in% paste0("GE", 1:5), ]
df2 <- GEs[GEs$GE %in% paste0("GE", 6:10), ]
df1$GE <- factor(df1$GE, levels = paste0("GE", 1:5))
df2$GE <- factor(df2$GE, levels = paste0("GE", 6:10))

p1 <- ggplot(df1, aes(GE, expression, fill = group))+
  geom_boxplot(linewidth = 0.6, outlier.shape = NA) +
  annotate("rect", xmin = 0.4, xmax = 1.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill="#9fd4ca") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill="#f8f6bf") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill="#c9c7de")+
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill="#f1958a") +
  annotate("rect", xmin = 4.5, xmax = 5.6, ymin = -Inf, ymax = Inf, alpha = 0.2, fill="#95bbd4")+
  # geom_dotplot(dotsize = 0.8, binaxis = "y", stackdir = "center", position = position_dodge(0.8))+
  geom_vline(xintercept = 1.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 2.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 3.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 4.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  theme_bw()+
  theme(axis.text.y = element_text(size=10, color = "#204056"),
        axis.text.x = element_text(size=10, angle = 45, hjust = 1, vjust = 1, color = "#204056"),
        axis.title = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values = nmf_colors)
p1
ggsave("GE1_5_Score_NMF.pdf", width = 30, height = 20, units = "cm")

p2 <- ggplot(df, aes(GE, expression, fill = group))+
  geom_boxplot(linewidth = 0.6) +
  annotate("rect", xmin = 0.4, xmax = 1.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill="#f6bf7f") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill="#bbd985") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill="#f8d5e5")+
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill="#b3ddf6") +
  annotate("rect", xmin = 4.5, xmax = 5.6, ymin = -Inf, ymax = Inf, alpha = 0.2, fill="#d2e8cd")+
  geom_vline(xintercept = 1.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 2.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 3.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  geom_vline(xintercept = 4.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  theme_bw()+
  theme(axis.text.y = element_text(size=10, color = "#204056"),
        axis.text.x = element_text(size=10, angle = 45, hjust = 1, vjust = 1, color = "#204056"),
        axis.title = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values = nmf_colors)
p2
ggsave("GE6_10_Score_NMF.pdf", width = 30, height = 20, units = "cm")


## FigureS5 GEs
# Analyze the gene elements composition of tumor cells
library(BiocManager)
library(GEOquery) 
library(plyr)
library(dplyr) 
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot) 
library(multtest)
library(msigdbr)
library(fgsea)
library(loomR)
library(clustree)
library(tibble)
library(SeuratData)
library(matrixStats)
library(sparseMatrixStats)
library(DESeq2)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(viridis)
library(gridExtra)
library(ggplotify)
library(multtest)
library(metap)
library(writexl)
library(Rcpp)
library(RcppZiggurat)
library(Rfast)
library(ggh4x)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(cola)
library(msigdbr)
library(UCell)
library(RColorBrewer)
# DEG_remove_mito 
DEG_Remove_mito <- function (df){
  df_rm_mito <- df[!grepl("^MT-|^MT.",df$gene),]
  return(df_rm_mito)
}

# DEG_remove_heat 
DEG_Remove_heat <- function (df){
  df_rm_hsp <- df[!grepl("^HSP",rownames(df)),]
  return(df_rm_hsp)
}

# plotSimilarityMatrix 
# create a similarity matrix plot from a dataframe
# adapted from klic package (adjusted heatmap parameters)
plotSimilarityMatrix = function(X, y = NULL, 
                                min.val = 0, 
                                max.val = 1,
                                clusLabels = NULL, 
                                colX = NULL, colY = NULL, 
                                clr = FALSE, clc = FALSE, 
                                annotation_col = NULL, 
                                annotation_row = NULL, 
                                annotation_colors = NULL, 
                                myLegend = NULL, 
                                fileName = "posteriorSimilarityMatrix", 
                                savePNG = FALSE, 
                                semiSupervised = FALSE, 
                                showObsNames = FALSE) {
  
  if (!is.null(y)) {
    # Check if the rownames correspond to the ones in the similarity matrix
    check <- sum(1 - rownames(X) %in% row.names(y))
    if (check == 1)
      stop("X and y must have the same row names.")
  }
  
  if (!is.null(clusLabels)) {
    if (!is.integer(clusLabels))
      stop("Cluster labels must be integers.")
    
    n_clusters <- length(table(clusLabels))
    riordina <- NULL
    for (i in 1:n_clusters) {
      riordina <- c(riordina, which(clusLabels == i))
    }
    
    X <- X[riordina, riordina]
    y <- y[riordina, ]
    y <- as.data.frame(y)
  }
  
  if (savePNG)
    grDevices::png(paste(fileName, ".png", sep = ""))
  
  if (!is.null(y)) {
    ht <- ComplexHeatmap::pheatmap(X, legend = TRUE,  
                                   color = rev(brewer.pal(11, "RdBu")), 
                                   breaks = seq(min.val, max.val, length.out = 11),
                                   cluster_rows = clr, 
                                   cluster_cols = clc, 
                                   #annotation_col = y,
                                   show_rownames = showObsNames, 
                                   show_colnames = showObsNames, 
                                   drop_levels = TRUE, 
                                   treeheight_row = -1,
                                   treeheight_col = -1,
                                   #annotation_row = annotation_row, 
                                   #annotation_col = annotation_col, 
                                   annotation_colors = annotation_colors)
  } else {
    ht <- ComplexHeatmap::pheatmap(X, legend = TRUE,
                                   color = rev(brewer.pal(11, "RdBu")),  
                                   breaks = seq(min.val, max.val, length.out = 11),
                                   cluster_rows = clr, 
                                   cluster_cols = clc,
                                   show_rownames = showObsNames, 
                                   show_colnames = showObsNames,
                                   treeheight_row = -1,
                                   treeheight_col = -1,
                                   drop_levels = TRUE, 
                                   #annotation_row = annotation_row, 
                                   #annotation_col = annotation_col, 
                                   annotation_colors = annotation_colors)
  }
  
  if (savePNG)
    grDevices::dev.off()
  
  return(ht)
}

# simil function 
# calculate similarity matrix 
# param df is dataframe of all cells for similarity matrix computation
# param drop is list of genes to drop from analysis
# param file is file name for output
# param method ("jaccard" or "corr") for method used
simil <- function(df, drop, file, method) {
  # drop the input gene list from analysis
  if (length(drop) > 0) {
    jc <- df[-drop, , drop = TRUE]
  }
  else {
    jc <- df
  }
  
  #jc[jc > 0] <- 1
  #jc[jc <= 0] <- -1
  jc <- as.matrix(jc)
  
  # calculate jaccard similarity matrix if method == "jaccard"
  if (method == "jaccard") {
    jc <- prabclus::jaccard(jc)
    jc <- 1 - jc
  }
  
  # calculate correlation matrix if method == "corr"
  else if (method == "corr") {
    jc <- Rfast::cora(jc)
  }
  
  # save file
  saveRDS(jc, file = file) 
  
  # return quantiles and mean similarity value
  v <- c(quantile(as.vector(jc), na.rm = TRUE), mean(as.vector(jc), na.rm = TRUE))
  names(v) <- c('0%','25%', '50%', '75%', '100%', 'mean')
  print(v)
  return(v)
}

# simil_plot 
# plot interactive similarity heatmap 
# param a is an .rds file generated by simil_calc
# param annot is a vector used for heatmap annotation (metadata column from Seurat object)

simil_plot <- function(a, min.val, max.val, annot) {
  jc <- readRDS(a) #read in similarity matrix
  
  # add annotation
  if (length(annot) > 0) {
    row_annot <- annot[rownames(jc), , drop = FALSE] # select rows to use for heatmap annotation
    col_annot <- annot[colnames(jc), , drop = FALSE] # select columns to use for heatmap annotation
    colors <- mako(n_distinct(annot)) # create color vector to use for annotation
    names(colors) <- base::unique(annot)[[1]]
    colors <- list(colors, colors)
    names(colors) <- c(as.name(names(annot)), as.name(names(annot)))
  } 
  else {
    row_annot <- NULL
    col_annot <- NULL
    colors <- NULL
  }
  
  # generate simialrity matrix plot (using plotSimilarity Matrix)
  hm <- plotSimilarityMatrix(jc, clr = TRUE, clc = TRUE, 
                             min.val = min.val, max.val = max.val,
                             annotation_row = row_annot, 
                             annotation_col = col_annot, 
                             annotation_colors = colors,
                             showObsNames = T) # plot full matrix
  return(hm)
  #hm <- draw(hm)
  #htShiny(hm) # generate interactive heatmap 
}

# simil_GE function 
setwd("simil_cancerepi/091322/110322_rerun")

GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- which(rownames(cancer.epi) %in% GElist$gene) # get indices of GE genes

# calculate Jaccard similarity matrix 
# param df is dataframe for comparison
# param file is file name for output
# param method ("jaccard" or "corr") for similarity method used
simil_GE <- function(df, file, method) {
  jc <- df[GElist, , drop = FALSE] # subset to sc50 genes
  # jc[jc > 0.5] <- 1
  # jc[jc < 0.5 & jc > -0.5] <- 0
  # jc[jc < -0.5] <- -1
  jc <- as.matrix(jc)
  
  # calculate jaccard similarity matrix if method == "jaccard"
  if (method == "jaccard") {
    jc <- prabclus::jaccard(jc)
    jc <- 1 - jc
  }
  
  # calculate correlation matrix if method == "corr"
  else if (method == "corr") {
    jc <- Rfast::cora(jc)
  }
  
  # save file
  saveRDS(jc, file = file) 
  
  # return quantiles and mean similarity value
  v <- c(quantile(as.vector(jc), na.rm = TRUE), mean(as.vector(jc), na.rm = TRUE))
  names(v) <- c('0%','25%', '50%', '75%', '100%', 'mean')
  print(v)
  return(v)
}

# Load in Seurat object 
setwd("/home/data/t100553/wss/genemodule_epicancer")
load("~/wss/genemodule_epicancer/epi_pt.RData")
cancer.epi <- subset(obj,copykat.pred=="aneuploid")
DefaultAssay(cancer.epi) <- "RNA"

# Scale data 

DefaultAssay(cancer.epi) <- "RNA"
cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")
cancer.epi <- FindVariableFeatures(cancer.epi, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
all.genes <- rownames(cancer.epi)
cancer.epi <- ScaleData(cancer.epi, features = all.genes)

# Drop samples with too few cells 

# create dataframe of Seurat object metadata
sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("sample",
                                "patient", 
                                "cancer", 
                                "group"))

# determine number of cancer epi cells per sample
drop_samples <- sobjlists %>% dplyr::count(sobjlists$sample) # generate dataframe with counts of # cells per sample
drop_samples <- drop_samples[drop_samples$n < 50, ] # get list of samples with n < 10 cells
nrow(drop_samples)

# drop samples with too few cells 
sobjlists <- sobjlists[!(sobjlists$sample %in% drop_samples$`sobjlists$sample`), ] 

# create dataframe for bar graph analysis
sobjlists <- sobjlists %>% dplyr::group_by(sample, 
                                           patient, 
                                           cancer, 
                                           group) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C*100) # add column with % of cells by sc50.Pred out of all cells in sample

# reorder sobjlists
sobjlists <- sobjlists[order(sobjlists$cancer),] # order by BC Subtype

# get list of tumor samples with enough cells for analysis
samples <- unique(sobjlists$sample)

# All cancer epi clustering 

# prepare cancer epithelial cell object
Epi.all.combo <- cancer.epi
DefaultAssay(Epi.all.combo) <- "RNA"

# cluster cells at various resolutions 
resolution.range <- c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)

# generate cancer epithelial cell UMAP
calculate_pcs <- function(scobj) {
  pct <- scobj[["pca"]]@stdev / sum(scobj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- if (length(pct) > 1) {
    sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  } else {
    NA
  }
  pcs <- ifelse(is.na(co2), co1, min(co1, co2, na.rm = TRUE))
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  list(
    pcs = pcs,
    co1 = co1,
    co2 = co2,
    plot_df = plot_df
  )
}
pca <- calculate_pcs(Epi.all.combo)$pcs
pca
Epi.all.combo <- FindNeighbors(Epi.all.combo, reduction = "pca", dims = 1:pca)
Epi.all.combo <- FindClusters(Epi.all.combo, resolution = resolution.range)
Epi.all.combo <- RunUMAP(Epi.all.combo, reduction = "pca", dims = 1:pca, verbose = TRUE, seed.use = 123)

# plot UMAP by clustering
p <- DimPlot(Epi.all.combo, reduction = "umap",group.by = "RNA_snn_res.0.1", 
             label = F, 
             repel = TRUE, 
             raster = FALSE) + SeuratAxes()
ggsave("cancerepi_UMAP.pdf", plot = as.ggplot(p), width = 5.8, height = 5.5)

# plot UMAP by original dataset
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "study") + SeuratAxes() #+ NoLegend()
ggsave("cancerepi_UMAP_origdataset.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# plot UMAP by clinical subtype
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "sample") + SeuratAxes()+ NoLegend()
ggsave("cancerepi_UMAP_sample.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# plot UMAP by sc50 subtype
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "cancer") + SeuratAxes()
ggsave("cancerepi_UMAP_cancer.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# update cancer epithelial cell object
cancer.epi <- Epi.all.combo

# list of clinically actionable targets
targets <- list("ESR1",
                "ERBB2", #HER2
                "PIK3CA",
                c("NTRK1", "NTRK2", "NTRK3"),
                "CD274", #PD-L1
                "ERBB3", #HER3
                "EGFR",
                c("FGFR1", "FGFR2", "FGFR3", "FGFR4"),
                "TACSTD2", #TROP2
                c("CDK4", "CDK6"), 
                "AR",
                "NECTIN2", 
                "LAG3")

# Add UCell score for clinically actionable targets
cancer.epi <- AddModuleScore_UCell(cancer.epi,
                                   features = targets,
                                   name = names(targets),
                                   assay = "RNA")

cancer.epi@meta.data <- cancer.epi@meta.data[,-c(99:4942,4944:15186)]

colnames(cancer.epi@meta.data)[102:114] <- c("ESR1",
                                             "ERBB2", 
                                             "PIK3CA",
                                             "NTRK", 
                                             "CD274", 
                                             "ERBB3", 
                                             "EGFR",
                                             "FGFR",
                                             "TACSTD2", 
                                             "CDK",
                                             "AR",
                                             "NECTIN2", 
                                             "LAG3")

getwd()
saveRDS(cancer.epi, "cancerepi_umap.rds")

# Unsupervised DGE generation 
setwd("/home/data/t100553/wss/genemodule_epicancer")
j <- readRDS("cancerepi_umap.rds")

# initialize objects for unsupervised DGE gene signature generation
all_DGE <- data.frame()

# perform DGE analysis at each of the various clustering resolutions
for (i in colnames(j@meta.data)[c(15,17, 21:33)]) {
  Idents(j) <- j@meta.data[,i]
  
  # DGE analysis (cluster biomarkers)
  j.markers_DGE <- FindAllMarkers(j, only.pos = T, 
                                  min.cells.group = 50, 
                                  min.diff.pct = 0.25, 
                                  # test.use = "MAST",
                                  logfc.threshold = 0.25)
  j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
  j.markers_DGE$res <- strsplit(i, "res.")[[1]][2]
  j.markers_DGE$cluster_res <- paste0(j.markers_DGE$cluster, "_", j.markers_DGE$res)
  
  all_DGE <- rbind(all_DGE, j.markers_DGE)
}

# save full output of all unsupervised DGE signatures
setwd("/home/data/t100553/wss/genemodule_epicancer")
write_xlsx(all_DGE,
           path = "cancerepi_DGE_unsupervised_3.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)


# Unsupervised sample-level DGE generation

setwd("D:/DI/test/genemodule_epicancer")
j <- readRDS("cancerepi_withtargets_110922.rds")
cancer.epi <- j

# initialize objects for unsupervised DGE gene signature generation
all_DGE <- data.frame()

for (i in samples[20:34]) {
  subset <- subset(j, subset = sample == i)
  DefaultAssay(subset) <- "RNA"
  subset <- NormalizeData(subset, assay = "RNA")
  subset <- FindVariableFeatures(subset, 
                                 selection.method = "vst", 
                                 nfeaøtures = 2000)
  all.genes <- rownames(subset)
  subset <- ScaleData(subset, features = all.genes)
  
  subset <- FindNeighbors(subset, reduction = "pca", dims = 1:30)
  resolution.range <- c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
  subset <- FindClusters(subset, 
                         graph.name = "RNA_snn", 
                         resolution = resolution.range)
  subset <- RunUMAP(subset, reduction = "pca", dims = 1:30, verbose = TRUE, seed.use=123)
  
  # perform DGE analysis at each of the various clustering resolutions
  for (k in colnames(subset@meta.data)[(grep("RNA_snn", colnames(subset@meta.data)[85:115]) + 84)]) {
    Idents(subset) <- subset@meta.data[,k]
    
    # DGE analysis (cluster biomarkers)
    j.markers_DGE <- FindAllMarkers(subset, only.pos = T, 
                                    min.cells.group = 5, 
                                    # test.use = "MAST",
                                    #logfc.threshold = 0.25,
                                    min.diff.pct = 0.1)
    if(dim(j.markers_DGE)[1] > 0) {
      j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
      j.markers_DGE$res <- strsplit(k, "res.")[[1]][2]
      j.markers_DGE$cluster_res <- paste0(i, "_",j.markers_DGE$cluster, "_", j.markers_DGE$res)
    }
    
    all_DGE <- rbind(all_DGE, j.markers_DGE)
  }
}

# save full output of all unsupervised DGE signatures
setwd("simil_cancerepi/091322/110322_rerun")
write_xlsx(all_DGE,
           path = "cancerepi_DGE_unsupervised_samplelevel_4.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)



# Calculate Jaccard similarity 

setwd("genemodule_epicancer")

# create Jaccard similarity matrix for supervised and unsupervised DGE lists
files <- list.files(pattern = ".xlsx")

cancer_DGEs <- data.frame()
for (i in files) {
  cancer_DGEs <- rbind(cancer_DGEs, readxl::read_xlsx(i))
}

# filter DGE signatures
cancer_DGEs <- readxl::read_xlsx(files)
cancer_DGEs <- DEG_Remove_mito(cancer_DGEs) # remove mitochondrial genes
#cancer_DGEs <- DEG_Remove_heat(cancer_DGEs) # remove HSP gene
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$p_val_adj < 0.05),] # filter with p-value threshold
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$avg_log2FC > 0),] # filter with log2FC threshold
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$pct.1 > 0.2),] # filter with percent threshold

# select top n genes per signature by adjusted p-value (includes ties)
cancer_DGEs <- cancer_DGEs %>%
  group_by(cluster_res) %>%
  top_n(-200, p_val_adj)

# select top n genes per signature by avg_log2FC
cancer_DGEs <- cancer_DGEs %>% 
  group_by(cluster_res) %>%
  top_n(200, avg_log2FC) 

# get number of genes in each DGE signature
DGE_counts <- cancer_DGEs %>% dplyr::count(cancer_DGEs$cluster_res)
View(DGE_counts)

# prepare signature matrix to use for calculating pairwise Jaccard indices between signatures
cancer_DGEs_forjaccard <- as.data.frame(unique(cancer_DGEs$gene))
colnames(cancer_DGEs_forjaccard) <- c('gene')
rownames(cancer_DGEs_forjaccard) <- cancer_DGEs_forjaccard$gene
for (i in unique(cancer_DGEs$cluster_res)) {
  j <- as.data.frame(cancer_DGEs[which(cancer_DGEs$cluster_res == i), ]$gene)
  j$i <- 1
  colnames(j) <- c("gene", i)
  rownames(j) <- j$gene
  cancer_DGEs_forjaccard <- left_join(cancer_DGEs_forjaccard, j, by = "gene")
}

rownames(cancer_DGEs_forjaccard) <- cancer_DGEs_forjaccard$gene
cancer_DGEs_forjaccard <- cancer_DGEs_forjaccard[,-1]
cancer_DGEs_forjaccard[is.na(cancer_DGEs_forjaccard)] <- 0

# drop unsupervised DGE signatures with fewer than n genes
drop <- vector()
for (i in 1:(dim(cancer_DGEs_forjaccard)[2])) {
  if (DGE_counts[which(DGE_counts$cluster_res ==
                       colnames(cancer_DGEs_forjaccard)[i]),3] < 20) {
    drop <- append(drop, i)
  }
}
cancer_DGEs_forjaccard <- cancer_DGEs_forjaccard[, -drop]

# calculate Jaccard similarity between all DGE signatures
setwd("/home/data/t100553/wss/genemodule_epicancer")
simil(cancer_DGEs_forjaccard, NULL, "cancerepi_DGEs_forjaccard_unsupervised_only.rds", "jaccard")
cancerepi_jaccard <- as.data.frame(readRDS("cancerepi_DGEs_forjaccard_unsupervised_only.rds"))

# remove redundant unsupervised DGE signatures
drop <- vector()
dim(cancerepi_jaccard)

for (i in 1:(ncol(cancerepi_jaccard))) {
  for (j in (i + 1):ncol(cancerepi_jaccard)) {
    if (!is.na(cancerepi_jaccard[i, j]) && is.numeric(cancerepi_jaccard[i, j]) && cancerepi_jaccard[i, j] > 0.95) {
      drop <- append(drop, j)
    }
  }
}
drop <- unique(drop)
cancerepi_jaccard_noredundant <- cancerepi_jaccard[-drop, -drop]
saveRDS(cancerepi_jaccard_noredundant, "cancerepi_DGEs_forjaccard_noredundant_unsupervised_only.rds")

# plot unsupervised DGE signature Jaccard similarities
p <- simil_plot("cancerepi_DGEs_forjaccard_noredundant_unsupervised_only.rds", 0, 1, NULL)
p <- as.ggplot(p)
ggsave("GE_prelim_unsupervised_only.pdf", plot = p, width = 40, height = 40)

# cola package to define GEs 
set.seed(123)
rh = consensus_partition(cancerepi_jaccard_noredundant, 
                         top_value_method = "ATC",
                         partition_method = "skmeans",
                         top_n = c(dim(cancerepi_jaccard_noredundant)[1]), 
                         p_sampling = 0.8,
                         max_k = 15)

k <- suggest_best_k(rh) #k <- 10

pdf("GE_stats_hclust.pdf", width = 7.3, height = 4)
select_partition_number(rh, mark_best = F)
dev.off()

write.csv(get_stats(rh), "GE_stats_hclust.csv")

pdf("GE_consensusplot.pdf", width = 7, height = 4)
ht_opt$message = FALSE
consensus_heatmap(rh, k = k)
dev.off()

write.csv(get_classes(rh, k = k), "GE_classes_hclust.csv")

# Define GEs
GEs <- read.csv("GE_classes_hclust.csv", stringsAsFactors = FALSE, na.strings = "unknown")[,c(1,2)]
colnames(GEs) <- c("cluster_res", "GE")

cancer_DGEs_GE <- left_join(cancer_DGEs, GEs, by = "cluster_res")

cancer_DGEs_GE <- cancer_DGEs_GE[-which(is.na(cancer_DGEs_GE$GE)),] 
cancer_DGEs_GE$totavg_log2FC <- NA
cancer_DGEs_GE$totmin.pct.diff <- NA

for (i in unique(cancer_DGEs_GE$gene)) {
  mean <-  mean(cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$avg_log2FC)
  cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$totavg_log2FC <- mean
  
  mean <-  mean(cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$percent)
  cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$totmin.pct.diff <- mean
}

# Filter GEs 
cancer_GE <- cancer_DGEs_GE %>% dplyr::group_by(GE, gene) %>%
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 

#cancer_GE <- cancer_GE[which(cancer_GE$percent > 0.1),]

cancer_GE <- cancer_GE %>%
  group_by(GE) %>%
  top_n(350, Nb)

cancer_GE <- left_join(cancer_GE, unique(cancer_DGEs_GE[,c(7,12)]))
cancer_GE <- cancer_GE %>% 
  group_by(GE) %>%
  top_n(200, totavg_log2FC) 
GE_summary <- cancer_GE %>% dplyr::group_by(GE) %>%
  dplyr::summarise(Nb = n())
View(GE_summary)
cancer_GE$old_GE <- cancer_GE$GE
cancer_GE$GE[which(cancer_GE$old_GE == 1)] <- 1
cancer_GE$GE[which(cancer_GE$old_GE == 2)] <- 3
cancer_GE$GE[which(cancer_GE$old_GE == 3)] <- 6
cancer_GE$GE[which(cancer_GE$old_GE == 4)] <- 4
cancer_GE$GE[which(cancer_GE$old_GE == 5)] <- 5
cancer_GE$GE[which(cancer_GE$old_GE == 6)] <- 9
cancer_GE$GE[which(cancer_GE$old_GE == 7)] <- 2
cancer_GE$GE[which(cancer_GE$old_GE == 8)] <- 7
cancer_GE$GE[which(cancer_GE$old_GE == 9)] <- 10
cancer_GE$GE[which(cancer_GE$old_GE == 10)] <- 8
write_xlsx(cancer_GE, "cancer_GE.xlsx")

# Add GE scores for cancer cells 
setwd("/home/data/t100553/wss/genemodule_epicancer")
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- split(GElist$gene,GElist$GE)
length(GElist) # get number of GEs

cancer.epi <- AddModuleScore_UCell(cancer.epi,
                                   features = GElist,
                                   assay = "RNA")
saveRDS(cancer.epi, file = "cancerepi_withGEs.rds")

# GE heatmap
set.seed(123)
expdata <- t(cancer.epi@meta.data[,c(8:12,34:43)])
collapse_expdata <- as.data.frame(rownames(expdata))

for (i in samples) {
  subset <- expdata[,which(expdata[2,] == i)]
  subset <- subset[,sample.int(dim(subset)[2],min(dim(subset)[2],20000))]
  collapse_expdata <- cbind(collapse_expdata, subset)
}

labels <- rownames(collapse_expdata)[6:15]
labels <- labels %>% gsub("_UCell", "",.)
labels <- labels %>% paste0("GE",.)

collapse_expdata <- collapse_expdata[,-1]
collapse_zscore <- collapse_expdata[-c(1:5),]
collapse_zscore <- as.matrix(sapply(collapse_zscore, as.numeric))
collapse_zscore <- t(apply(collapse_zscore, 1, function(x) (x-mean(x))/sd(x)))
sort <- rbind(apply(collapse_zscore, 2, function(x) which.max(x)),
              apply(collapse_zscore, 2, function(x) max(x)), 
              collapse_zscore)
sort <- sort[,order(sort[2,],decreasing = T)] # sort by max z-score
sort <- sort[,order(sort[1,],decreasing = F)] #sort by GE
collapse_zscore <- sort[-c(1,2),]
collapse_expdata <- collapse_expdata[,colnames(collapse_zscore)]

orig_anno <- t(as.matrix(collapse_expdata[1,]))
colnames(orig_anno) <- c("origin")
sample_anno <- t(as.matrix(collapse_expdata[2,]))
colnames(sample_anno) <- c("sample")
BC_anno <- t(as.matrix(collapse_expdata[3,]))
colnames(BC_anno) <- c("BC subtype")
pam50_anno <- t(as.matrix(collapse_expdata[5,]))
colnames(pam50_anno) <- c("PAM50 subtype")
sc50_anno <- t(as.matrix(collapse_expdata[4,]))
colnames(sc50_anno) <- c("SC50 subtype") 
anno <- HeatmapAnnotation("orig" = orig_anno,
                          "sample" = sample_anno,
                          "BC" = BC_anno,
                          "pam50" = pam50_anno,
                          "sc50" = sc50_anno, 
                          show_legend = c("sample" = FALSE), 
                          col = list("orig" = c("Pal_Prim" = "#18A900",
                                                "Qian" = "#FFC300",
                                                "Wu" = "#C70039", 
                                                "Karaayvaz" = "#eb7d34",
                                                "Wu2021prim" = "#006CA9",
                                                "Xu" = "#A300DB"),
                                     "BC" = c("HER2+" = "#700639", "HR+" = "#397006", "TNBC" = "#063970"),
                                     "pam50" = c("Basal" = "#2596be", "Her2" = "#9925be", "LumA" = "#49be25", "LumB" = "#be4d25"),
                                     "sc50" = c("Basal_SC" = "#76b5c5", "Her2E_SC" = "#ad76c5", "LumA_SC" = "#c58676", "LumB_SC" = "#8dc576"))
)
pdf("heatmap.pdf",width = 8,height = 15)
collapse_zscore <- t(collapse_zscore)
pheatmap(collapse_zscore, cluster_rows = FALSE, cluster_cols = FALSE,show_rownames = FALSE, show_colnames = FALSE,
         color = colorRampPalette(c("#6c94c7","white","#f0b2aa","#dd7e74","#d66c61","#d46258","#bd2b2e","#bc272c","#B8222E"))(1000))
dev.off()
pheatmap(collapse_zscore, cluster_rows = FALSE, cluster_cols = FALSE)

# GE function analysis
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))  
  union <- length(union(set1, set2))            
  return(intersection / union)                 
}

library(msigdbr)
msigdbr_species() 

human_KEGG = msigdbr(species = "Homo sapiens",category = "H" 
                     
) %>% 
  dplyr::select(gs_name,gene_symbol)
#category = "H"
H.geneset = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)#list
num_sets_A <- length(GElist)
num_sets_B <- length(H.geneset)
similarity_matrix <- matrix(0, nrow = num_sets_A, ncol = num_sets_B)
for (i in 1:num_sets_A) {
  for (j in 1:num_sets_B) {
    similarity_matrix[i, j] <- jaccard_similarity(GElist[[i]], H.geneset[[j]])
  }
}
rownames(similarity_matrix) <- names(GElist)
colnames(similarity_matrix) <- names(H.geneset)
rownames(similarity_matrix) <- paste0("GE",rownames(similarity_matrix))
print(similarity_matrix)
similarity_matrix <- data.frame(t(similarity_matrix))
pheatmap(similarity_matrix,color = colorRampPalette(c("#6c94c7","white","#fdedea","#fbdedb","#f0b2aa","#dd7e74","#d66c61","#d46258","#bd2b2e","#bc272c"))(1000))

write_xlsx(similarity_matrix,
           path = "GE_H.geneset_jaccard.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)

# ClusterProfiler for each GE
library(org.Hs.eg.db)
library(enrichplot)
converted_gene_sets <- lapply(GElist, function(genes) {
  gene <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(gene$ENTREZID)
})

go_results <- lapply(converted_gene_sets, function(genes) {
  enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
})

kegg_results <- lapply(converted_gene_sets, function(genes) {
  enrichKEGG(gene = genes, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
})
getwd()
save(go_results,files = "./go_results.RData")
save(kegg_results,files = "kegg_results.RData")

g1 <- go_results[[1]]@result
k1 <- kegg_results[[1]]@result
g2 <- go_results[[2]]@result
k2 <- kegg_results[[2]]@result
g3 <- go_results[[3]]@result
k3 <- kegg_results[[3]]@result
g4 <- go_results[[4]]@result
k4 <- kegg_results[[4]]@result
g5 <- go_results[[5]]@result
k5 <- kegg_results[[5]]@result
g6 <- go_results[[6]]@result
k6 <- kegg_results[[6]]@result
g7 <- go_results[[7]]@result
k7 <- kegg_results[[7]]@result
g8 <- go_results[[8]]@result
k8 <- kegg_results[[8]]@result
g9 <- go_results[[9]]@result
k9 <- kegg_results[[9]]@result
g10 <- go_results[[10]]@result
k10 <- kegg_results[[10]]@result

#ploting
#GE1
go.GE1 <- as.data.frame(g1) %>%
  head(n = 10) %>% 
  arrange(desc(qvalue)) 
kegg.GE1 <- as.data.frame(k1) %>% 
  head(n = 10) %>% 
  arrange(desc(qvalue)) 

ggplot(data = go.GE1, 
       aes(x = Count, y = reorder(Description, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "go.GE1") +theme_bw()

ggplot(data = kegg.GE1, 
       aes(x = Count, y = reorder(subcategory, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "kegg.GE1") +theme_bw()

#GE2
go.GE2 <- as.data.frame(g2) %>%
  head(n = 10) %>% 
  arrange(desc(qvalue))
kegg.GE2 <- as.data.frame(k2) %>%
  head(n = 10) %>% 
  arrange(desc(qvalue))

ggplot(data = go.GE2, 
       aes(x = Count, y = reorder(Description, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "go.GE2") +theme_bw()

ggplot(data = kegg.GE2, 
       aes(x = Count, y = reorder(subcategory, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "kegg.GE2") +theme_bw()

#GE3
go.GE3 <- as.data.frame(g3) %>%
  head(n = 10) %>% # 
  arrange(desc(qvalue))
kegg.GE3 <- as.data.frame(k3) %>%
  head(n = 10) %>% # 
  arrange(desc(qvalue))

ggplot(data = go.GE3, 
       aes(x = Count, y = reorder(Description, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "go.GE3") +theme_bw()

ggplot(data = kegg.GE3, 
       aes(x = Count, y = reorder(subcategory, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "kegg.GE3") +theme_bw()

#GE4
go.GE4 <- as.data.frame(g4) %>%
  head(n = 10) %>% # 
  arrange(desc(qvalue))
kegg.GE4 <- as.data.frame(k4) %>%
  head(n = 10) %>% # 
  arrange(desc(qvalue))

ggplot(data = go.GE4, 
       aes(x = Count, y = reorder(Description, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "go.GE4") +theme_bw()

ggplot(data = kegg.GE4, 
       aes(x = Count, y = reorder(subcategory, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "kegg.GE4") +theme_bw()

#GE5
go.GE5 <- as.data.frame(g5) %>%
  head(n = 10) %>% # 
  arrange(desc(qvalue))
kegg.GE5 <- as.data.frame(k5) %>%
  head(n = 10) %>% 
  arrange(desc(qvalue))

ggplot(data = go.GE5, 
       aes(x = Count, y = reorder(Description, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "go.GE5") +theme_bw()

ggplot(data = kegg.GE5, 
       aes(x = Count, y = reorder(subcategory, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "kegg.GE5") +theme_bw()


#GE6
go.GE6 <- as.data.frame(g6) %>%
  head(n = 10) %>% # 
  arrange(desc(qvalue))
kegg.GE6 <- as.data.frame(k6) %>%
  head(n = 10) %>% # 
  arrange(desc(qvalue))

ggplot(data = go.GE6, 
       aes(x = Count, y = reorder(Description, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "go.GE6") +theme_bw()

ggplot(data = kegg.GE6, 
       aes(x = Count, y = reorder(subcategory, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "kegg.GE6") +theme_bw()

#GE7
go.GE7 <- as.data.frame(g7) %>%
  head(n = 10) %>% # 
  arrange(desc(qvalue))
kegg.GE7 <- as.data.frame(k7) %>%
  head(n = 10) %>% 
  arrange(desc(qvalue))

ggplot(data = go.GE7, 
       aes(x = Count, y = reorder(Description, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "go.GE7") +theme_bw()

ggplot(data = kegg.GE7, 
       aes(x = Count, y = reorder(subcategory, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "kegg.GE7") +theme_bw()

#GE8
go.GE8 <- as.data.frame(g8) %>%
  head(n = 10) %>% # 
  arrange(desc(qvalue))
kegg.GE8 <- as.data.frame(k8) %>%
  head(n = 10) %>% #
  arrange(desc(qvalue))

ggplot(data = go.GE8, 
       aes(x = Count, y = reorder(Description, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "go.GE8") +theme_bw()

ggplot(data = kegg.GE8, 
       aes(x = Count, y = reorder(subcategory, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "kegg.GE8") +theme_bw()

#GE9
go.GE9 <- as.data.frame(g9) %>%
  head(n = 10) %>% 
  arrange(desc(qvalue))
kegg.GE9 <- as.data.frame(k9) %>%
  head(n = 10) %>% 
  arrange(desc(qvalue))

ggplot(data = go.GE9, 
       aes(x = Count, y = reorder(Description, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "go.GE9") +theme_bw()
ggplot(data = kegg.GE9, 
       aes(x = Count, y = reorder(subcategory, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "kegg.GE9") +theme_bw()

#GE10
go.GE10 <- as.data.frame(g10) %>%
  head(n = 10) %>% 
  arrange(desc(qvalue))
kegg.GE10 <- as.data.frame(k10) %>%
  head(n = 10) %>% 
  arrange(desc(qvalue))
ggplot(data = go.GE10, 
       aes(x = Count, y = reorder(Description, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "go.GE10") +theme_bw()
ggplot(data = kegg.GE10, 
       aes(x = Count, y = reorder(subcategory, Count), fill = -log10(qvalue)))+
  geom_bar(stat = "identity", width = 0.8) + 
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "kegg.GE10") +theme_bw()

# ClusterProfiler
library(ggplot2)
library(clusterProfiler)
data(gcSample)
xx.go <- compareCluster(converted_gene_sets,fun="enrichGO", OrgDb="org.Hs.eg.db")
dotplot(xx.go, showCategory=5, includeAll=FALSE) 
xx.kegg <- compareCluster(
  geneClusters = converted_gene_sets,
  fun = "enrichKEGG",
  organism = "hsa"  
)
dotplot(xx.kegg, showCategory=5, includeAll=FALSE) 
save(xx.go,files = "./go.RData")

# Calculate the Jaccard correlation between the NMF modules of the GE module and the TIME cells
setwd("genemodule_epicancer")
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist$GE <- paste0("GE",GElist$GE)
GElist <- split(GElist$gene,GElist$GE)
TIME.NMF <- read.csv("./TIME_NMF_markers2.csv")
TIME.NMFlist <- split(TIME.NMF$gene,TIME.NMF$cluster)

jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))  
  union <- length(union(set1, set2))            
  return(intersection / union)                  
}
num_sets_A <- length(GElist)
num_sets_B <- length(TIME.NMFlist)
similarity_matrix <- matrix(0, nrow = num_sets_A, ncol = num_sets_B)
for (i in 1:num_sets_A) {
  for (j in 1:num_sets_B) {
    similarity_matrix[i, j] <- jaccard_similarity(GElist[[i]], TIME.NMFlist[[j]])
  }
}
rownames(similarity_matrix) <- names(GElist)
colnames(similarity_matrix) <- names(TIME.NMFlist)
print(similarity_matrix)
similarity_matrix <- data.frame(t(similarity_matrix))
pheatmap(similarity_matrix,color = colorRampPalette(c("#6c94c7","white","#fdedea","#fbdedb","#f0b2aa","#dd7e74","#d66c61","#d46258","#bd2b2e","#bc272c"))(1000))
write_xlsx(similarity_matrix,
           path = "GE_TIME.NMFlist_jaccard.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)
