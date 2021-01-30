########################################
# practice 3, 10x, unfiltered PBMC
# date: 2021.01.30 - 01.30
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-3
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
####################10X Genomics PBMC 4K
# 数据下载
library(BiocFileCache)
# 在当前工作目录新建raw_data目录
bfc <- BiocFileCache("raw_data", ask = FALSE)
# 提供下载链接
# http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz
raw_path <- bfcrpath(
    bfc,
    file.path(
        "http://cf.10xgenomics.com/samples",
        "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
    )
)

# 解压数据至当前工作目录，并新建pbmc4k目录
untar(raw_path, exdir = file.path(getwd(), "pbmc4k"))

library(DropletUtils)
fname <- file.path(getwd(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce_pbmc_4k<- read10xCounts(fname, col.names = TRUE)
dim(sce_pbmc_4k)
# [1]  33694 737280

sce_pbmc_4k
# class: SingleCellExperiment
# dim: 33694 737280
# metadata(1): Samples
# assays(1): counts
# rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#   ENSG00000268674
# rowData names(2): ID Symbol
# colnames(737280): AAACCTGAGAAACCAT-1 AAACCTGAGAAACCGC-1 ... TTTGTCATCTTTAGTC-1
#   TTTGTCATCTTTCCTC-1
# colData names(2): Sample Barcode
# reducedDimNames(0):
# altExpNames(0):

# gene annotation -------------------------------------------------------------
# ID整合
library(scater)
rownames(sce_pbmc_4k) <- uniquifyFeatureNames(
    rowData(sce_pbmc_4k)$ID,
    rowData(sce_pbmc_4k)$Symbol
)

# 添加位置信息
# BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
location <- mapIds(
    EnsDb.Hsapiens.v86,
    keys = rowData(sce_pbmc_4k)$ID,
    column = "SEQNAME",
    keytype = "GENEID"
)


# detecting dropout -----------------------------------------------------------
set.seed(100)
e_out <- emptyDrops(counts(sce_pbmc_4k))
e_out
# DataFrame with 737280 rows and 5 columns
#                        Total   LogProb    PValue   Limited       FDR
#                    <integer> <numeric> <numeric> <logical> <numeric>
# AAACCTGAGAAACCAT-1         1        NA        NA        NA        NA
# AAACCTGAGAAACCGC-1         0        NA        NA        NA        NA
# AAACCTGAGAAACCTA-1         1        NA        NA        NA        NA
# AAACCTGAGAAACGAG-1         0        NA        NA        NA        NA
# AAACCTGAGAAACGCC-1         1        NA        NA        NA        NA
# ...                      ...       ...       ...       ...       ...
# TTTGTCATCTTTACAC-1         2        NA        NA        NA        NA
# TTTGTCATCTTTACGT-1        33        NA        NA        NA        NA
# TTTGTCATCTTTAGGG-1         0        NA        NA        NA        NA
# TTTGTCATCTTTAGTC-1         0        NA        NA        NA        NA
# TTTGTCATCTTTCCTC-1         1        NA        NA        NA        NA

sce_pbmc_4k_filtered <- sce_pbmc_4k[, which(e_out$FDR <= 0.001)]
dim(sce_pbmc_4k_filtered)
# [1] 33694  4300


# qc, especially for mitochondrial --------------------------------------------
stats <- perCellQCMetrics(
    sce_pbmc_4k_filtered,
    subsets = list(Mito = which(location == "MT"))
)
high_mito <- isOutlier(stats$subsets_Mito_percent, type = "higher")
summary(high_mito)
#    Mode   FALSE    TRUE
# logical    3985     315

sce_pbmc_4k_final <- sce_pbmc_4k_filtered[, !high_mito]
dim(sce_pbmc_4k_final)
# [1] 33694  3985

##### 使用qc标准对原数据作图
colData(sce_pbmc_4k_filtered) <- cbind(colData(sce_pbmc_4k_filtered), stats)
sce_pbmc_4k_filtered$discard <- high_mito

# 使用qc标准对原数据作图
gridExtra::grid.arrange(
    plotColData(sce_pbmc_4k_filtered, y = "sum", colour_by = "discard") +
    scale_y_log10() + ggtitle("Total count"),
    plotColData(sce_pbmc_4k_filtered, y = "detected", colour_by = "discard") +
    scale_y_log10() + ggtitle("Detected features"),
    plotColData(sce_pbmc_4k_filtered, y = "subsets_Mito_percent",
        colour_by = "discard") + ggtitle("Mito percent"),
    ncol = 3
)

# 查看线粒体含量与文库大小的关系
plotColData(sce_pbmc_4k_filtered, x = "sum", y = "subsets_Mito_percent",
    colour_by = "discard") + scale_x_log10() +
ggtitle("Mito percent to total count")


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(1000)
cluster_pbmc_4k <- quickCluster(sce_pbmc_4k_final)
sce_pbmc_4k_final <- computeSumFactors(
    sce_pbmc_4k_final,
    cluster = cluster_pbmc_4k
)
# logNormCounts()
sce_pbmc_4k_final <- logNormCounts(sce_pbmc_4k_final)

summary(sizeFactors(sce_pbmc_4k_final))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.00749  0.71207  0.87490  1.00000  1.09900 12.25412

# 比较librarySizeFactors()（常规方法）、
# computeSumFactors()（去卷积方法）进行归一化的差异
plot(
    librarySizeFactors(sce_pbmc_4k_final),
    sizeFactors(sce_pbmc_4k_final),
    pch = 16,
    xlab = "Library size factors",
    ylab = "Deconvolution factors",
    log = "xy"
)
abline(a = 0, b = 1, col = "red")


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
set.seed(1001)
dec_pbmc_4k_pois <- modelGeneVarByPoisson(sce_pbmc_4k_final)
top_hvgs_pbmc_4k <- getTopHVGs(dec_pbmc_4k_pois, prop = 0.1)
length(top_hvgs_pbmc_4k)
# [1] 1599

# 查看方差大小
plot(
    dec_pbmc_4k_pois$mean,
    dec_pbmc_4k_pois$total,
    main = "10X Genomics, unfiltered PBMC", pch = 16, cex = 0.5,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
cur_fit <- metadata(dec_pbmc_4k_pois)
# 蓝线指的是所有gene都会存在的一种偏差
curve(cur_fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


# dimension reduce, using three methods ---------------------------------------
##### PCA
set.seed(10000)
sce_pbmc_4k_final <- denoisePCA(
    sce_pbmc_4k_final,
    subset.row = top_hvgs_pbmc_4k,
    technical = dec_pbmc_4k_pois
)
dim(reducedDim(sce_pbmc_4k_final, "PCA"))
# [1] 3985    9

##### t-SNE
set.seed(100000)
sce_pbmc_4k_final <- runTSNE(sce_pbmc_4k_final, dimred = "PCA")
dim(reducedDim(sce_pbmc_4k_final, "TSNE"))
# [1] 3985    2

##### UMAP
set.seed(1000000)
sce_pbmc_4k_final <- runUMAP(sce_pbmc_4k_final, dimred = "PCA")
dim(reducedDim(sce_pbmc_4k_final, "UMAP"))
# [1] 3985    2


# clustering ------------------------------------------------------------------
##### graph-based clustering
# 使用PCA的前几个PCs构建SNNG
library(scran)
snn_gr_pbmc_4k_k10 <- buildSNNGraph(
    sce_pbmc_4k_final,
    k = 10,
    use.dimred = "PCA"
)

# 鉴定cluster
cluster_pbmc_4k_k10 <- igraph::cluster_walktrap(snn_gr_pbmc_4k_k10)$membership
table(cluster_pbmc_4k_k10)
# clust_pbmc_4k_k10
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
# 205 508 541  56 374 125  46 432 302 867  47 155 166  61  84  16

# 把cluster信息存成SingleCellExperiment对象的一个因子
library(scater)
colLabels(sce_pbmc_4k_final) <- factor(cluster_pbmc_4k_k10)

sce_pbmc_4k_final
# class: SingleCellExperiment
# dim: 33694 3985
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(33694): RP11-34P13.3 FAM138A ... AC213203.1 FAM231B
# rowData names(2): ID Symbol
# colnames(3985): AAACCTGAGAAGGCCT-1 AAACCTGAGACAGACC-1 ... TTTGTCAGTTAAGACA-1
#   TTTGTCATCCCAAGAT-1
# colData names(4): Sample Barcode sizeFactor label
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):

# 绘制t-SNE图
plotTSNE(sce_pbmc_4k_final, colour_by = "label")


# detecting markers -----------------------------------------------------------
markers_pbmc_4k_final <- findMarkers(
    sce_pbmc_4k_final,
    groups = colLabels(sce_pbmc_4k_final),
    pval.type = "some",
    direction = "up"
)

# use cluster 8 as an explaination
chosen_cluster <- "8"
markers_cluster_8 <- markers_pbmc_4k_final[[chosen_cluster]]

index <- grep(
    paste(c("CD14", "CD68", "MNDA"), collapse = "|"),
    rownames(markers_cluster_8)
)
index
# [1]  2 13 25
markers_cluster_8[index, 1:4]

plotExpression(
    sce_pbmc_4k_final,
    features = c("CD14", "CD68", "MNDA", "FCGR3A"),
    x = "label",
    colour_by = "label"
)
