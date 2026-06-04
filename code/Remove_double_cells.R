## Remove the double cells according to the samples after QC
library(Seurat)
library(tidyverse)
library(scuttle)
library(scDblFinder)
library(ggplot2)
library(writexl)
library(SeuratDisk)
library(celldex)
library(sctransform)
library(readxl)

sce <- readRDS("~/GI/cancerData_postQC/OPSCC.rds")
dim(sce)
seurat_list <- SplitObject(sce, split.by = "sample")

finddoublet <- function(seurat_obj){
  set.seed(34) 
  sce_obj <- as.SingleCellExperiment(seurat_obj)
  sce_obj <- scDblFinder(sce_obj)
  seu_new <- as.Seurat(sce_obj)
  return(seu_new)
}
nMAD <- function(x, nmads = 3){
  xm <- median(x)
  md <- median(abs(x - xm))
  mads <- xm + nmads * md
  return(mads)
}

nmad = 5
count = 1
sample.nfeaure.cut <- c()
sample.ncount.cut <- c()
sample.mt.cut <- c()
seurat_list_qc <- c()
for (obj in seurat_list){
  print(names(seurat_list)[count])
  # add percent.mt for qc
  obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = "^MT-")
  # visual nfeature ncount percent.mt
  print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  # add doublet info
  obj = finddoublet(obj)
  # nMAD cut
  nfeature.upcut <- ceiling(nMAD(obj$nFeature_RNA, nmad))
  ncount.upcut <- ceiling(nMAD(obj$nCount_RNA, nmad))
  permt.upcut <- ceiling(nMAD(obj$percent.mt, nmad))
  sample.nfeaure.cut <- c(sample.nfeaure.cut, nfeature.upcut)
  sample.ncount.cut <- c(sample.ncount.cut, ncount.upcut)
  sample.mt.cut <- c(sample.mt.cut, permt.upcut)
  # perform QC
  obj.filt <- subset(obj, subset = nFeature_RNA <= nfeature.upcut & nFeature_RNA > nfeature.upcut/20)
  obj.filt <- subset(obj.filt, subset = nCount_RNA <= ncount.upcut & nCount_RNA > ncount.upcut/20)
  obj.filt <- subset(obj.filt, subset = percent.mt <= min(permt.upcut, 25))
  obj.filt <- subset(obj.filt, subset = scDblFinder.class %in% c('singlet'))
  Idents(obj.filt) <- names(seurat_list)[count]
  # visualize again
  print(VlnPlot(obj.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  seurat_list_qc <- c(seurat_list_qc, obj.filt)
  count = count + 1
}
names(sample.mt.cut) <- names(seurat_list)
names(sample.ncount.cut) <- names(seurat_list)
names(sample.nfeaure.cut) <- names(seurat_list)
names(seurat_list_qc) <- names(seurat_list)

sce.big <- merge(x = seurat_list_qc[[1]],
                 y = seurat_list_qc[-1],
                 add.cell.ids = names(seurat_list_qc),
                 project = "OPSCC")
dim(sce.big)
dim(sce)
saveRDS(sce.big, file = "~/GI/cancerData_postQC_finddoublet/OPSCC_QC.rds")
rm(list = ls()); gc()
