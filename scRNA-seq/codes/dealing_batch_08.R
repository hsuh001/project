########################################
# dealing with batch effect
# date: 2021.01.11 - 01.13
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-8
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
# 使用多个不同批次的10X PBMC数据
BiocManager::install("TENxPBMCData")
library(TENxPBMCData)
all_pbmc_sce <- list(
    pbmc3k = TENxPBMCData("pbmc3k"),
    pbmc4k = TENxPBMCData("pbmc4k"),
    pbmc8k = TENxPBMCData("pbmc8k")
)
all_pbmc_sce
# $pbmc3k
# class: SingleCellExperiment
# dim: 32738 2700
# metadata(0):
# assays(1): counts
# rownames(32738): ENSG00000243485 ENSG00000237613 ... ENSG00000215616
#   ENSG00000215611
# rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
# colnames: NULL
# colData names(11): Sample Barcode ... Individual Date_published
# reducedDimNames(0):
# altExpNames(0):
# 
# $pbmc4k
# class: SingleCellExperiment
# dim: 33694 4340
# metadata(0):
# assays(1): counts
# rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#   ENSG00000268674
# rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
# colnames: NULL
# colData names(11): Sample Barcode ... Individual Date_published
# reducedDimNames(0):
# altExpNames(0):
# 
# $pbmc8k
# class: SingleCellExperiment
# dim: 33694 8381
# metadata(0):
# assays(1): counts
# rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#   ENSG00000268674
# rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
# colnames: NULL
# colData names(11): Sample Barcode ... Individual Date_published
# reducedDimNames(0):
# altExpNames(0):


# quality control -------------------------------------------------------------
library(scater)
stats <- high_mito <- list()
for (n in names(all_pbmc_sce)) {
    current <- all_pbmc_sce[[n]]
    is_mito <- grep("MT", rowData(current)$Symbol_TENx)
    stats[[n]] <- perCellQCMetrics(current, subsets = list(Mito = is_mito))
    high_mito[[n]] <- isOutlier(stats[[n]]$subsets_Mito_percent, type = "higher")
    all_pbmc_sce[[n]] <- current[, !high_mito[[n]]]
}


# normalization ---------------------------------------------------------------
all_pbmc_sce <- lapply(all_pbmc_sce, logNormCounts)


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
library(scran)
all_pbmc_dec <- lapply(all_pbmc_sce, modelGeneVar)
all_pbmc_hvgs <- lapply(all_pbmc_dec, getTopHVGs, prop = 0.1)
length(all_pbmc_hvgs[[1]])
# [1] 1001

length(all_pbmc_hvgs[[2]])
# [1] 1352


# dimension reduce ------------------------------------------------------------
library(BiocSingular)
set.seed(10000)
all_pbmc_sce <- mapply(
    FUN = runPCA,
    x = all_pbmc_sce,
    subset_row = all_pbmc_hvgs,
    MoreArgs = list(ncomponents = 25, BSPARAM = RandomParam()),
    SIMPLIFY = FALSE
)

set.seed(100000)
all_pbmc_sce <- lapply(all_pbmc_sce, runTSNE, dimred = "PCA")

set.seed(1000000)
all_pbmc_sce <- lapply(all_pbmc_sce, runUMAP, dimred = "PCA")


# clustering ------------------------------------------------------------------
for (n in names(all_pbmc_sce)) {
    g <- buildSNNGraph(all_pbmc_sce[[n]], k = 10, use.dimred = "PCA")
    clust <- igraph::cluster_walktrap(g)$membership
    colLabels(all_pbmc_sce[[n]]) <- factor(clust)
}


# re-normalization ------------------------------------------------------------
# 提取数据pbmc_3k和pbmc_4k
pbmc_3k <- all_pbmc_sce$pbmc3k
dec_3k <- all_pbmc_dec$pbmc3k
pbmc_3k
# class: SingleCellExperiment
# dim: 32738 2609
# metadata(0):
# assays(2): counts logcounts
# rownames(32738): ENSG00000243485 ENSG00000237613 ... ENSG00000215616
#   ENSG00000215611
# rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
# colnames: NULL
# colData names(13): Sample Barcode ... sizeFactor label
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):

pbmc_4k <- all_pbmc_sce$pbmc4k
dec_4k <- all_pbmc_dec$pbmc4k
pbmc_4k
# class: SingleCellExperiment
# dim: 33694 4182
# metadata(0):
# assays(2): counts logcounts
# rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#   ENSG00000268674
# rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
# colnames: NULL
# colData names(13): Sample Barcode ... sizeFactor label
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):

# 选出两个数据的gene交集，并各自取数据子集
universe <- intersect(rownames(pbmc_3k), rownames(pbmc_4k))
length(universe)
# [1] 31232

pbmc_3k_sub <- pbmc_3k[universe, ]
dec_3k_sub <- dec_3k[universe, ]

pbmc_4k_sub <- pbmc_4k[universe, ]
dec_4k_sub <- dec_4k[universe, ]

# 矫正批次间测序深度的影响
BiocManager::install("batchelor")
library(batchelor)
re_normalized <- multiBatchNorm(pbmc_3k_sub, pbmc_4k_sub)
pbmc_3k_sub_renorm <- re_normalized[[1]]
pbmc_3k_sub_renorm
# class: SingleCellExperiment
# dim: 31232 2609
# metadata(0):
# assays(2): counts logcounts
# rownames(31232): ENSG00000243485 ENSG00000237613 ... ENSG00000198695
#   ENSG00000198727
# rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
# colnames: NULL
# colData names(13): Sample Barcode ... sizeFactor label
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):

pbmc_4k_sub_renorm <- re_normalized[[2]]
pbmc_4k_sub_renorm
# class: SingleCellExperiment
# dim: 31232 4182
# metadata(0):
# assays(2): counts logcounts
# rownames(31232): ENSG00000243485 ENSG00000237613 ... ENSG00000198695
#   ENSG00000198727
# rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
# colnames: NULL
# colData names(13): Sample Barcode ... sizeFactor label
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):


# re-modeling and re-selecting HVGs -------------------------------------------
library(scran)
combined_dec <- combineVar(dec_3k_sub, dec_4k_sub)
combined_chosen_hvgs <- combined_dec$bio > 0
sum(combined_chosen_hvgs)
# [1] 13431


# check out the batch effects -------------------------------------------------
##### 合并数据
# cbind()会报错
cbind(pbmc_3k_sub_renorm, pbmc_4k_sub_renorm)
# Error in FUN(X[[i]], ...) :
#   column(s) 'Symbol_TENx' in 'mcols' are duplicated and the data do not match

# 检查原因, Symbol_TENx。同一个的Enseml ID的Symbol就不同
mcols(pbmc_3k_sub_renorm)
# DataFrame with 31232 rows and 3 columns
#                      ENSEMBL_ID  Symbol_TENx       Symbol
#                     <character>  <character>  <character>
# ENSG00000243485 ENSG00000243485   MIR1302-10           NA
# ENSG00000237613 ENSG00000237613      FAM138A      FAM138A
# ENSG00000186092 ENSG00000186092        OR4F5        OR4F5
# ENSG00000238009 ENSG00000238009 RP11-34P13.7 LOC100996442
# ENSG00000239945 ENSG00000239945 RP11-34P13.8           NA
# ...                         ...          ...          ...
# ENSG00000212907 ENSG00000212907      MT-ND4L         ND4L
# ENSG00000198886 ENSG00000198886       MT-ND4          ND4
# ENSG00000198786 ENSG00000198786       MT-ND5          ND5
# ENSG00000198695 ENSG00000198695       MT-ND6          ND6
# ENSG00000198727 ENSG00000198727       MT-CYB         CYTB

mcols(pbmc_4k_sub_renorm)
# DataFrame with 31232 rows and 3 columns
#                      ENSEMBL_ID  Symbol_TENx       Symbol
#                     <character>  <character>  <character>
# ENSG00000243485 ENSG00000243485 RP11-34P13.3           NA
# ENSG00000237613 ENSG00000237613      FAM138A      FAM138A
# ENSG00000186092 ENSG00000186092        OR4F5        OR4F5
# ENSG00000238009 ENSG00000238009 RP11-34P13.7 LOC100996442
# ENSG00000239945 ENSG00000239945 RP11-34P13.8           NA
# ...                         ...          ...          ...
# ENSG00000212907 ENSG00000212907      MT-ND4L         ND4L
# ENSG00000198886 ENSG00000198886       MT-ND4          ND4
# ENSG00000198786 ENSG00000198786       MT-ND5          ND5
# ENSG00000198695 ENSG00000198695       MT-ND6          ND6
# ENSG00000198727 ENSG00000198727       MT-CYB         CYTB

# 暂时先不考虑Symbol的问题，检查Ensembl ID是否一致
identical(mcols(pbmc_3k_sub_renorm)[, 1], mcols(pbmc_4k_sub_renorm)[, 1])
# [1] TRUE

# Ensembl ID一致，变成同一个即可。暂时不考虑Symbol的问题
rowData(pbmc_3k_sub_renorm) <- rowData(pbmc_4k_sub_renorm)
pbmc_3k_sub_renorm$batch <- "3k"
pbmc_4k_sub_renorm$batch <- "4k"
uncorrected <- cbind(pbmc_3k_sub_renorm, pbmc_4k_sub_renorm)

##### PCA
# 使用合并两个数据集后得到的HVGs进行PCA
library(scater)
library(BiocSingular)
set.seed(0010101010)
uncorrected <- runPCA(
    uncorrected,
    subset_row = combined_chosen_hvgs,
    BSPARAM = RandomParam()
)

# 聚类分群
library(scran)
combined_snn_gr <- buildSNNGraph(uncorrected, use.dimred = "PCA")
combined_clusters <- igraph::cluster_walktrap(combined_snn_gr)$membership
combined_tab <- table(
    Cluster = combined_clusters,
    Batch = uncorrected$batch
)
combined_tab
#        Batch
# Cluster   3k   4k
#      1     1  781
#      2     0 1309
#      3     0  535
#      4    14   51
#      5     0  605
#      6   489    0
#      7     0  184
#      8  1272    0
#      9     0  414
#      10  151    0
#      11    0   50
#      12  155    0
#      13    0   65
#      14    0   61
#      15    0   88
#      16   30    0
#      17  339    0
#      18  145    0
#      19   11    3
#      20    2   36

# t-SNE
set.seed(1111001)
uncorrected <- runTSNE(uncorrected, dimred = "PCA")
plotTSNE(uncorrected, colour_by = "batch")


# correcting batch effects by linear regression -------------------------------
library(batchelor)
rescaled <- rescaleBatches(pbmc_3k_sub_renorm, pbmc_4k_sub_renorm)
rescaled
# class: SingleCellExperiment
# dim: 31232 6791
# metadata(0):
# assays(1): corrected
# rownames(31232): ENSG00000243485 ENSG00000237613 ... ENSG00000198695
#   ENSG00000198727
# rowData names(0):
# colnames: NULL
# colData names(1): batch
# reducedDimNames(0):
# altExpNames(0):

# PCA, clustering
rescaled <- runPCA(
    rescaled,
    subset_row = combined_chosen_hvgs,
    exprs_values = "corrected"
)
rescaled_snn_gr <- buildSNNGraph(rescaled, use.dimred = "PCA")
rescaled_clusters <- igraph::cluster_walktrap(rescaled_snn_gr)$membership
rescaled_tab <- table(
    Cluster = rescaled_clusters,
    Batch = rescaled$batch
)
rescaled_tab
#        Batch
# Cluster    1    2
#      1   147   93
#      2   336  605
#      3   415   28
#      4   148  185
#      5    65  763
#      6   271  517
#      7   471  481
#      8    13    3
#      9    17   78
#      10   17   50
#      11   13   40
#      12  566 1057
#      13    4   34
#      14  109  201
#      15    3   36
#      16    3    8
#      17   11    3

# t-SNE
rescaled <- runTSNE(rescaled, dimred = "PCA")
rescaled$batch <- factor(rescaled$batch)
plotTSNE(rescaled, colour_by = "batch")


# correcting batch effects by MNN ---------------------------------------------
set.seed(1000101001)
# d, 指定的主成分个数
# k, 近邻的数量
mnn_out <- fastMNN(
    pbmc_3k_sub_renorm,
    pbmc_4k_sub_renorm,
    d = 50,
    k = 20,
    subset.row = combined_chosen_hvgs,
    BSPARAM = BiocSingular::RandomParam(deferred = TRUE)
)
mnn_out
# class: SingleCellExperiment
# dim: 13431 6791
# metadata(2): merge.info pca.info
# assays(1): reconstructed
# rownames(13431): ENSG00000239945 ENSG00000228463 ... ENSG00000198695
#   ENSG00000198727
# rowData names(1): rotation
# colnames: NULL
# colData names(1): batch
# reducedDimNames(1): corrected
# altExpNames(0):

head(mnn_out$batch)
# [1] 1 1 1 1 1 1

dim(reducedDim(mnn_out, "corrected"))
# [1] 6791   50

reducedDim(mnn_out, "corrected")[1:3, 1:3]
#             [,1]       [,2]         [,3]
# [1,] -0.12783777  0.0409469 -0.000116194
# [2,] -0.03364614 -0.1466529  0.161690348
# [3,] -0.09895631  0.1219669  0.010321992

dim(assay(mnn_out, "reconstructed"))
# [1] 13431  6791

assay(mnn_out, "reconstructed")[1:3, 1:3]
# <3 x 3> matrix of class LowRankMatrix and type "double":
#                          [,1]          [,2]          [,3]
# ENSG00000239945 -2.522191e-06 -1.851424e-06 -1.198984e-05
# ENSG00000228463 -6.626821e-04 -6.724341e-04 -4.820230e-04
# ENSG00000237094 -8.077231e-05 -8.038006e-05 -9.630608e-05

# clustering
library(scran)
mnn_snn_gr <- buildSNNGraph(mnn_out, use.dimred = "corrected")
mnn_clusters <- igraph::cluster_walktrap(mnn_snn_gr)$membership
mnn_tab <- table(
    Cluster = mnn_clusters,
    Batch = mnn_out$batch
)
mnn_tab
#        Batch
# Cluster    1    2
#      1   337  606
#      2   289  542
#      3   152  181
#      4    12    4
#      5   517  467
#      6    17   19
#      7   313  661
#      8   162  118
#      9    11   56
#      10  547 1083
#      11   17   59
#      12   16   58
#      13  144   93
#      14   67  191
#      15    4   36
#      16    4    8

# t-SNE
library(scater)
mnn_out <- runTSNE(mnn_out, dimred = "corrected")
mnn_out$batch <- factor(mnn_out$batch)
plotTSNE(mnn_out, colour_by = "batch")

# batchelor::fastMNN()自带检测结果
metadata(mnn_out)$merge.info$lost.var
#             [,1]        [,2]
# [1,] 0.006617087 0.003315395


# 更为细致的检查 ---------------------------------------------------------------------
##### 同一批次内clusters的比较
library(pheatmap)
## 方法一，绘制两个批次混合后的cluster与各自之前的cluster的对比图
# batch 1, pbmc3k
after_3k_tab <- table(
    paste("after", mnn_clusters[rescaled$batch == 1]),
    paste("before", colLabels(pbmc_3k_sub_renorm))
)
heat_3k <- pheatmap(
    log10(after_3k_tab + 10),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    main = "PBMC 3K comparison",
    silent = TRUE
)

# batch 2, pbmc4k
after_4k_tab <- table(
    paste("after", mnn_clusters[rescaled$batch == 2]),
    paste("before", colLabels(pbmc_4k_sub_renorm))
)
heat_4k <- pheatmap(
    log10(after_4k_tab + 10),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    main = "PBMC 4K comparison",
    silent = TRUE
)

gridExtra::grid.arrange(heat_3k[[4]], heat_4k[[4]], ncol = 2)

## 方法二，计算rand index
install.packages("fossil")
library(fossil)
# batch 1, pbmc3k
rand_index_3k <- rand.index(
    as.integer(mnn_clusters[rescaled$batch == 1]),
    as.integer(colLabels(pbmc_3k_sub_renorm))
)
rand_index_3k
# [1] 0.9335226

# batch 2, pbmc4k
rand_index_4k <- rand.index(
    as.integer(mnn_clusters[rescaled$batch == 2]),
    as.integer(colLabels(pbmc_4k_sub_renorm))
)
rand_index_4k
# [1] 0.9575746

##### 利用marker gene检查
# 首先分别从原来的两个batch中寻找marker gene
stats_3k <- pairwiseWilcox(pbmc_3k_sub_renorm, direction = "up")
markers_3k <- getTopMarkers(stats_3k[[1]], stats_3k[[2]], n = 10)
length(markers_3k)
# [1] 10

stats_4k <- pairwiseWilcox(pbmc_4k_sub_renorm, direction = "up")
markers_4k <- getTopMarkers(stats_4k[[1]], stats_4k[[2]], n = 10)
length(markers_4k)
# [1] 13

# 取两个批次的marker合集，作为矫正的HVGs
markers_3k_4k <- unique(unlist(c(unlist(markers_3k), unlist(markers_4k))))
length(markers_3k_4k)
# [1] 314

set.seed(1000110)
mnn_out_2 <- fastMNN(
    pbmc_3k_sub_renorm,
    pbmc_4k_sub_renorm,
    d = 50,
    k = 20,
    subset.row = markers_3k_4k,
    BSPARAM = BiocSingular::RandomParam(deferred = TRUE)
)

# 绘图, t-SNE
library(scater)
mnn_out_2 <- runTSNE(mnn_out_2, dimred = "corrected")
gridExtra::grid.arrange(
    plotTSNE(
        mnn_out_2[, mnn_out_2$batch == 1],
        colour_by = I(colLabels(pbmc_3k_sub_renorm))
    ),
    plotTSNE(
        mnn_out_2[, mnn_out_2$batch == 2],
        colour_by = I(colLabels(pbmc_4k_sub_renorm))
    ),
    ncol = 2
)

##### 矫正后数据的应用
# 使用合并后、但未矫正的数据寻找marker gene
uncorrected <- cbind(pbmc_3k_sub_renorm, pbmc_4k_sub_renorm)

mnn_clusters <- igraph::cluster_walktrap(mnn_snn_gr)$membership

markers_out <- findMarkers(
    uncorrected,
    mnn_clusters,
    block = uncorrected$batch,
    direction = "up",
    lfc = 1,
    row.data = rowData(uncorrected)[, 3, drop = FALSE]
)
markers_out
# List of length 16
# names(16): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
