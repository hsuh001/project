########################################
# transformation between different scRNA-seq data formats
# date: 2021.01.25 - 01.25
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-15
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# 不同R包数据的相互转换 -----------------------------------------------------------------
##### SeuratObject与SingleCellExperiment的相互转换
library(scater)
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
library(loomR)
remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
library(Seurat)
BiocManager::install("patchwork")
library(patchwork)

## SeuratObject转SingleCellExperiment
data("pbmc_small")
pbmc_small
# An object of class Seurat
# 230 features across 80 samples within 1 assay 
# Active assay: RNA (230 features, 20 variable features)
#  2 dimensional reductions calculated: pca, tsne

# 使用as.SingleCellExperiment()
sce_pbmc <- as.SingleCellExperiment(pbmc_small)
sce_pbmc
# class: SingleCellExperiment
# dim: 230 80
# metadata(0):
# assays(2): counts logcounts
# rownames(230): MS4A1 CD79B ... SPON2 S100B
# rowData names(5): vst.mean vst.variance vst.variance.expected
#   vst.variance.standardized vst.variable
# colnames(80): ATGCCAGAACGACT CATGGCCTGTGCAT ... GGAACACTTCAGAC
#   CTTGATTGATCTTC
# colData names(8): orig.ident nCount_RNA ... RNA_snn_res.1 ident
# reducedDimNames(2): PCA TSNE
# spikeNames(0):
# altExpNames(0):

# 接下来使用scater操作
p1 <- plotExpression(
    sce_pbmc,
    features = "MS4A1",
    x = "ident"
) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
)
p2 <- plotPCA(sce_pbmc, colour_by = "ident")
p1 + p2

## SingleCellExperiment转SeuratObject
# https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/manno_human.rds
sce_manno <- readRDS(file = "manno_human.rds")
sce_manno
# class: SingleCellExperiment
# dim: 20560 4029
# metadata(0):
# assays(2): counts logcounts
# rownames(20560): 'MARC1' 'MARC2' ... ZZEF1 ZZZ3
# rowData names(10): feature_symbol
#   is_feature_control ... total_counts
#   log10_total_counts
# colnames(4029): 1772122_301_C02 1772122_180_E05
#   ... 1772116-063_G02 1772099-259_H03
# colData names(34): Species cell_type1 ...
#   pct_counts_ERCC is_cell_control
# reducedDimNames(0):
# altExpNames(0):

sce_manno <- runPCA(sce_manno)
# 使用as.Seurat()
seurat_manno <- as.Seurat(
    sce_manno,
    counts = "counts",
    data = "logcounts"
)

seurat_manno
# An object of class Seurat
# 20560 features across 4029 samples within 1 assay 
# Active assay: RNA (20560 features)
#  1 dimensional reduction calculated: PCA

Idents(seurat_manno) <- "cell_type1"
p1 <- DimPlot(
    seurat_manno,
    reduction = "PCA",
    group.by = "Source"
) + NoLegend()
p2 <- RidgePlot(
    seurat_manno,
    features = "ACTB",
    group.by = "Source"
)
p1 + p2

##### SeuratObject与loom的相互转换
## SeuratObject转loom
loom_pbmc <- as.loom(pbmc, filename = "pbmc3k.loom", verbose = FALSE)
loom_pbmc
# Class: loom
# Filename: /__w/1/s/output/pbmc3k.loom
# Access type: H5F_ACC_RDWR
# Attributes: version, chunks, LOOM_SPEC_VERSION, assay, last_modified
# Listing:
#        name    obj_type dataset.dims dataset.type_class
#   col_attrs   H5I_GROUP         <NA>               <NA>
#  col_graphs   H5I_GROUP         <NA>               <NA>
#      layers   H5I_GROUP         <NA>               <NA>
#      matrix H5I_DATASET 2638 x 13714          H5T_FLOAT
#   row_attrs   H5I_GROUP         <NA>               <NA>
#  row_graphs   H5I_GROUP         <NA>               <NA>

# 最后使用完要记得关上loom对象
loom_pbmc$close_all()

## loom转SeuratObject
# 用loomR::connect()读取
loom_l6_immune <- connect(
    filename = "../data/l6_r1_immune_cells.loom",
    mode = "r"
)
loom_l6_immune
## Class: loom
## Filename: /__w/1/s/data/l6_r1_immune_cells.loom
## Access type: H5F_ACC_RDONLY
## Attributes: CreationDate, last_modified
## Listing:
##        name    obj_type  dataset.dims dataset.type_class
##   col_attrs   H5I_GROUP          <NA>               <NA>
##  col_graphs   H5I_GROUP          <NA>               <NA>
##      layers   H5I_GROUP          <NA>               <NA>
##      matrix H5I_DATASET 14908 x 27998          H5T_FLOAT
##   row_attrs   H5I_GROUP          <NA>               <NA>
##  row_graphs   H5I_GROUP          <NA>               <NA>

# as.Seurat()转换
seurat_l6_immune <- as.Seurat(loom_l6_immune)
VlnPlot(
    seurat_l6_immune,
    features = c("Sparc", "Ftl1", "Junb", "Ccl4"),
    ncol = 2,
    pt.size = 0.1
)

# 关闭loom对象
loom_l6_immune$close_all()

## 补充，Seurat V2，Convert()
data("pbmc_small")
pfile <- Convert(
    from = pbmc_small,
    to = "loom",
    filename = "pbmc_small.loom",
    display.progress = FALSE
)
pfile
# Class: loom
# Filename: /home/paul/Documents/Satija/pbmc_small.loom
# Access type: H5F_ACC_RDWR
# Attributes: version, chunks
# Listing:
#        name    obj_type dataset.dims dataset.type_class
#   col_attrs   H5I_GROUP         <NA>               <NA>
#  col_graphs   H5I_GROUP         <NA>               <NA>
#      layers   H5I_GROUP         <NA>               <NA>
#      matrix H5I_DATASET     80 x 230          H5T_FLOAT
#   row_attrs   H5I_GROUP         <NA>               <NA>
#  row_graphs   H5I_GROUP         <NA>               <NA>

##### AnnData转SeuratObject
pbmc_3k <- ReadH5AD(file = "pbmc3k.h5ad")
# 利用Seurat操作
Idents(pbmc_3k) <- "louvain"
p1 <- DimPlot(pbmc_3k, label = TRUE)
p2 <- VlnPlot(pbmc_3k, features = c("CST3", "NKG7", "PPBP"), combine = FALSE)
wrap_plots(c(list(p1), p2), ncol = 2) & NoLegend()
