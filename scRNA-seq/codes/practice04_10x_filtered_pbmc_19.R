########################################
# practice 4, 10x, filtered PBMC
# date: 2021.01.31 - 01.31
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-4
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
####################10X Genomics PBMC 3K, 4K, 8K
BiocManager::install("TENxPBMCData")
library(TENxPBMCData)
sce_pbmc_all <- list(
    pbmc3k = TENxPBMCData("pbmc3k"),
    pbmc4k = TENxPBMCData("pbmc4k"),
    pbmc8k = TENxPBMCData("pbmc8k")
)
sce_pbmc_all
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


# qc --------------------------------------------------------------------------
# 备份数据用于质控探索
unfiltered <- sce_pbmc_all

# 通过线粒体含量计算质控结果，随后过滤
library(scater)
stats <- high_mito <- list()
for (n in names(sce_pbmc_all)) {
    current <- sce_pbmc_all[[n]]
    is_mito <- grep("MT", rowData(current)$Symbol_TENx)
    stats[[n]] <- perCellQCMetrics(current, subsets = list(Mito = is_mito))
    high_mito[[n]] <- isOutlier(
        stats[[n]]$subsets_Mito_percent,
        type = "higher"
    )
    sce_pbmc_all[[n]] <- current[, !high_mito[[n]]]
}

# 查看根据线粒体含量过滤的结果
lapply(high_mito, summary)
# $pbmc3k
#    Mode   FALSE    TRUE
# logical    2609      91
# 
# $pbmc4k
#    Mode   FALSE    TRUE
# logical    4182     158
# 
# $pbmc8k
#    Mode   FALSE    TRUE
# logical    8157     224

# 批量作图，将结果存在list，方便后期导出
qc_plots <- list()

for (n in names(sce_pbmc_all)) {
    current <- unfiltered[[n]]
    colData(current) <- cbind(colData(current), stats[[n]])
    current$discard <- high_mito[[n]]
    qc_plots[[n]] <- plotColData(
        current,
        x = "sum",
        y = "subsets_Mito_percent",
        colour_by = "discard"
    ) + scale_x_log10()
}
do.call(gridExtra::grid.arrange, c(qc_plots, ncol = 3))


# 归一化, normalization ----------------------------------------------------------
sce_pbmc_all <- lapply(sce_pbmc_all, logNormCounts)

# 查看size factor
lapply(sce_pbmc_all, function(x) summary(sizeFactors(x)))
# $pbmc3k
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.2338  0.7478  0.9262  1.0000  1.1571  6.6042
# 
# $pbmc4k
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.3155  0.7109  0.8903  1.0000  1.1272 11.0267
# 
# $pbmc8k
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.2963  0.7043  0.8772  1.0000  1.1177  6.7942


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
library(scran)
pbmc_dec_all <- lapply(sce_pbmc_all, modelGeneVar)
pbmc_hvgs_all <- lapply(pbmc_dec_all, getTopHVGs, prop = 0.1)

lapply(pbmc_hvgs_all, length)
# $pbmc3k
# [1] 1001
# 
# $pbmc4k
# [1] 1352
# 
# $pbmc8k
# [1] 1509

# 查看方差大小
par(mfrow = c(1, 3))
for (n in names(pbmc_dec_all)) {
    cur_dec <- pbmc_dec_all[[n]]
    plot(
        cur_dec$mean,
        cur_dec$total,
        pch = 16,
        cex = 0.5,
        xlab = "Mean of log-expression",
        ylab = "Variance of log-expression"
        )
        cur_fit <- metadata(cur_dec)
        # 蓝线指的是所有gene都会存在的一种偏差
        curve(cur_fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
}


# dimension reduce ------------------------------------------------------------
library(BiocSingular)
set.seed(10000)
sce_pbmc_all <- mapply(
    FUN = runPCA,
    x = sce_pbmc_all,
    subset_row = pbmc_hvgs_all,
    MoreArgs = list(ncomponents = 25, BSPARAM = RandomParam()),
    SIMPLIFY = FALSE
)

set.seed(100000)
sce_pbmc_all <- lapply(sce_pbmc_all, runTSNE, dimred = "PCA")

set.seed(1000000)
sce_pbmc_all <- lapply(sce_pbmc_all, runUMAP, dimred = "PCA")


# clustering ------------------------------------------------------------------
for (n in names(sce_pbmc_all)) {
    snn_gr <- buildSNNGraph(sce_pbmc_all[[n]], k = 10, use.dimred = "PCA")
    cluster <- igraph::cluster_walktrap(snn_gr)$membership
    colLabels(sce_pbmc_all[[n]]) <- factor(cluster)
}

# 查看各自分的cluster数目
lapply(sce_pbmc_all, function(x) table(colLabels(x)))
# $pbmc3k
# 
#   1   2   3   4   5   6   7   8   9  10
# 487 154 603 514  31 150 179 333 147  11
# 
# $pbmc4k
# 
#    1    2    3    4    5    6    7    8    9   10   11   12   13
#  497  185  569  786  373  232   44 1023   77  218   88   54   36
# 
# $pbmc8k
# 
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18
# 1004  759 1073 1543  367  150  201 2067   59  154  244   67   76  285   20   15   64    9

# 绘制t-SNE图
pbmc_tsne_all <- list()

for (n in names(sce_pbmc_all)) {
    pbmc_tsne_all[[n]] <- plotTSNE(
        sce_pbmc_all[[n]],
        colour_by = "label"
    ) + ggtitle(n)
}
do.call(gridExtra::grid.arrange, c(pbmc_tsne_all, list(ncol = 3)))


# 数据整合 ------------------------------------------------------------------------
# 查看三个数据集的gene数量
lapply(sce_pbmc_all, function(x) length(rownames(x)))
# $pbmc3k
# [1] 32738
# 
# $pbmc4k
# [1] 33694
# 
# $pbmc8k
# [1] 33694

# 找共有gene
universe <- Reduce(intersect, lapply(sce_pbmc_all, rownames))
length(universe)
# [1] 31232

# 对每个数据集取子集
sce_pbmc_all_2 <- lapply(sce_pbmc_all, "[", i = universe,)

lapply(sce_pbmc_all_2, dim)
# $pbmc3k
# [1] 31232  2609
# 
# $pbmc4k
# [1] 31232  4182
# 
# $pbmc8k
# [1] 31232  8157

# 对找HVGs的结果取子集
# 末尾的逗号不加会报错
pbmc_dec_all_2 <- lapply(pbmc_dec_all, "[", i = universe,)

# 把三个数据集当作三个batch，重新进行归一化
library(batchelor)
sce_pbmc_all_2_renorm <- do.call(multiBatchNorm, sce_pbmc_all_2)

# 根据重新归一化的结果，再次寻找HVGs
library(scran)
combined_dec <- do.call(combineVar, pbmc_dec_all_2)
combined_hvgs <- getTopHVGs(combined_dec, n = 5000)
length(combined_hvgs)
# [1] 5000

# 使用MNN进行批次效应矫正，同时还进行PCA降维
set.seed(1000101)
pbmc_merged <- do.call(
    fastMNN,
    c(
        sce_pbmc_all_2_renorm,
        list(
            subset.row = combined_hvgs,
            BSPARAM = BiocSingular::RandomParam()
        )
    )
)
pbmc_merged
# class: SingleCellExperiment
# dim: 5000 14948
# metadata(2): merge.info pca.info
# assays(1): reconstructed
# rownames(5000): ENSG00000090382 ENSG00000163220 ... ENSG00000122068
#   ENSG00000011132
# rowData names(1): rotation
# colnames: NULL
# colData names(1): batch
# reducedDimNames(1): corrected
# altExpNames(0):

table(pbmc_merged$batch)
# pbmc3k pbmc4k pbmc8k
#   2609   4182   8157

# 使用lost.var检查结果
metadata(pbmc_merged)$merge.info$lost.var
#            pbmc3k       pbmc4k      pbmc8k
# [1,] 7.003371e-03 3.125816e-03 0.000000000
# [2,] 7.137205e-05 5.125053e-05 0.003002682

# 重新进行聚类
library(scran)
snn_gr_merged <- buildSNNGraph(pbmc_merged, use.dimred = "corrected")
cluster_merged <- igraph::cluster_walktrap(snn_gr_merged)$membership
colLabels(pbmc_merged) <- factor(cluster_merged)
tab_merged <- table(
    Cluster = cluster_merged,
    Batch = pbmc_merged$batch
)
tab_merged
#        Batch
# Cluster pbmc3k pbmc4k pbmc8k
#      1     473    774   1576
#      2     596    523   1071
#      3     327    588   1123
#      4     300    541   1024
#      5      21     66    130
#      6      11     53    157
#      7     138    160    272
#      8      16     50     83
#      9     445   1010   2013
#      10      6      8     14
#      11    137     84    146
#      12      8     13     77
#      13     92    214    308
#      14      4      9     15
#      15      2      5     19
#      16      6      9     13
#      17      1     10      9
#      18      3     36     64
#      19     12     26     34
#      20     11      3      9

# 绘制t-SNE图
set.seed(10101010)
pbmc_merged <- runTSNE(pbmc_merged, dimred = "corrected")
# 查看3个数据混合后的聚类效果，以及是否存在批次效应
gridExtra::grid.arrange(
    plotTSNE(
        pbmc_merged,
        colour_by = "label",
        text_by = "label",
        text_colour = "red"
    ),
    plotTSNE(pbmc_merged, colour_by = "batch"),
    ncol = 2
)
