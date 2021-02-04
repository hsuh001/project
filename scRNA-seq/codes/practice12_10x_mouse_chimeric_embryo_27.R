########################################
# practice 11, 10X Genomics, mouse chimeric embryo
# date: 2021.02.04 - 02.04
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-12
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
# 使用小鼠E8.5时期的嵌合体胚胎数据
# BiocManager::install("MouseGastrulationData")
library(MouseGastrulationData)
sce_chimera <- WTChimeraData(samples = 5:10)
sce_chimera
# class: SingleCellExperiment
# dim: 29453 20935
# metadata(0):
# assays(1): counts
# rownames(29453): ENSMUSG00000051951 ENSMUSG00000089699 ... ENSMUSG00000095742
#   tomato-td
# rowData names(2): ENSEMBL SYMBOL
# colnames(20935): cell_9769 cell_9770 ... cell_30702 cell_30703
# colData names(11): cell barcode ... doub.density sizeFactor
# reducedDimNames(2): pca.corrected.E7.5 pca.corrected.E8.5
# altExpNames(0):


# gene annotation -------------------------------------------------------------
names(colData(sce_chimera))
#  [1] "cell"            "barcode"         "sample"          "stage"
#  [5] "tomato"          "pool"            "stage.mapped"    "celltype.mapped"
#  [9] "closest.cell"    "doub.density"    "sizeFactor"

table(sce_chimera$sample)
#    5    6    7    8    9   10
# 2411 1047 3007 3097 4544 6829

library(scater)
rownames(sce_chimera) <- uniquifyFeatureNames(
    rowData(sce_chimera)$ENSEMBL,
    rowData(sce_chimera)$SYMBOL
)


# quality control -------------------------------------------------------------
drop <- sce_chimera$celltype.mapped %in% c("stripped", "Doublet")

table(drop)
# drop
# FALSE  TRUE
# 19426  1509

sce_chimera_filtered <- sce_chimera[, !drop]
sce_chimera_filtered
# class: SingleCellExperiment
# dim: 29453 19426
# metadata(0):
# assays(1): counts
# rownames(29453): Xkr4 Gm1992 ... CAAA01147332.1 tomato-td
# rowData names(2): ENSEMBL SYMBOL
# colnames(19426): cell_9769 cell_9770 ... cell_30701 cell_30702
# colData names(11): cell barcode ... doub.density sizeFactor
# reducedDimNames(2): pca.corrected.E7.5 pca.corrected.E8.5
# altExpNames(0):


# normalization ---------------------------------------------------------------
# 数据已经计算过size factors
sce_chimera_filtered <- logNormCounts(sce_chimera_filtered)


# variance modeling -----------------------------------------------------------
library(scran)
dec_chimera <- modelGeneVar(
    sce_chimera_filtered,
    block = sce_chimera_filtered$sample
)
chosen_hvgs_chimera <- dec_chimera$bio > 0
table(chosen_hvgs_chimera)
# chosen_hvgs_chimera
# FALSE  TRUE
# 14754 14699


# 数据整合 + 矫正批次效应 ---------------------------------------------------------------
library(batchelor)
set.seed(01001001)
sce_chimera_merged <- correctExperiments(
    sce_chimera_filtered,
    batch = sce_chimera_filtered$sample,
    subset.row = chosen_hvgs_chimera,
    PARAM = FastMnnParam(
        merge.order = list(
            list(1, 3, 5),  # WT (3 replicates)
            list(2, 4, 6)   # td-Tomato (3 replicates)
        )
    )
)

# 查看lost.var值
metadata(sce_chimera_merged)$merge.info$lost.var
#                 5            6            7            8           9          10
# [1,] 0.000000e+00 0.0204433327 0.000000e+00 0.0169566707 0.000000000 0.000000000
# [2,] 0.000000e+00 0.0007388531 0.000000e+00 0.0004408687 0.000000000 0.015473671
# [3,] 3.089886e-02 0.0000000000 2.011941e-02 0.0000000000 0.000000000 0.000000000
# [4,] 9.023683e-05 0.0000000000 8.272484e-05 0.0000000000 0.018046768 0.000000000
# [5,] 4.320774e-03 0.0072517998 4.123814e-03 0.0078279624 0.003831263 0.007786139

sce_chimera_merged
# class: SingleCellExperiment
# dim: 14699 19426
# metadata(2): merge.info pca.info
# assays(3): reconstructed counts logcounts
# rownames(14699): Xkr4 Rp1 ... Vmn2r122 CAAA01147332.1
# rowData names(3): rotation ENSEMBL SYMBOL
# colnames(19426): cell_9769 cell_9770 ... cell_30701 cell_30702
# colData names(12): batch cell ... doub.density sizeFactor
# reducedDimNames(3): corrected pca.corrected.E7.5 pca.corrected.E8.5
# altExpNames(0):


# clustering ------------------------------------------------------------------
snn_gr_chimera <- buildSNNGraph(sce_chimera_merged, use.dimred = "corrected")
clusters_chimera <- igraph::cluster_louvain(snn_gr_chimera)$membership
colLabels(sce_chimera_merged) <- factor(clusters_chimera)

# 查看分群和细胞类型之间的关系
tab <- table(
    Cluster = clusters_chimera,
    Sample = sce_chimera_merged$sample
)
pheatmap::pheatmap(
    log(tab + 10),
    color = viridis::viridis(100)
)


# dimension reduce ------------------------------------------------------------
sce_chimera_merged <- runTSNE(
    sce_chimera_merged,
    dimred = "corrected"
)
sce_chimera_merged <- runUMAP(
    sce_chimera_merged,
    dimred = "corrected"
)

gridExtra::grid.arrange(
    plotTSNE(
        sce_chimera_merged,
        colour_by = "label",
        text_by = "label",
        text_col = "red"
    ),
    plotTSNE(sce_chimera_merged, colour_by = "batch"),
    ncol = 2
)
