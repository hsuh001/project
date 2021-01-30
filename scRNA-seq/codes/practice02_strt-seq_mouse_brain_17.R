########################################
# practice 2, STRT-seq, mouse brain
# date: 2021.01.29 - 01.30
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-2
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
library(scRNAseq)
load("sce_zeisel.RData")
# sce_zeisel <- ZeiselBrainData()
sce_zeisel <- sce.zeisel
sce_zeisel
# class: SingleCellExperiment
# dim: 20006 3005
# metadata(0):
# assays(1): counts
# rownames(20006): Tspan12 Tshz1 ... mt-Rnr1 mt-Nd4l
# rowData names(1): featureType
# colnames(3005): 1772071015_C02 1772071017_G12 ... 1772066098_A12 1772058148_F03
# colData names(10): tissue group # ... level1class level2class
# reducedDimNames(0):
# altExpNames(2): ERCC repeat

# ERCC数量为57
dim(altExp(sce_zeisel, "ERCC"))
# [1]   57 3005

# 已经标注好线粒体gene
table(rowData(sce_zeisel))
# endogenous       mito
#      19972         34


# preprocessing and gene annotation -------------------------------------------
library(scater)
# 将每个细胞的所有gene表达量加起来，得到每个细胞的文库大小
# 同时替换一些奇怪的gene名
sce_zeisel <- aggregateAcrossFeatures(
    sce_zeisel,
    id = sub("_loc[0-9]+$", "", rownames(sce_zeisel))
)
dim(sce_zeisel)
# [1] 19839  3005

# 有些gene会有多个loc
head(rownames(sce_zeisel)[grep("_loc[0-9]+$", rownames(sce_zeisel))])
# [1] "Syne1_loc2"              "Hist1h2ap_loc1"          "Inadl_loc1"
# [4] "OTTMUSG00000016609_loc4" "OTTMUSG00000016609_loc3" "Gm5643_loc2"

# 以Syne1为例，查看aggregateAcrossFeatures()的效果
length(grep("Syne1", rownames(test)))
# [1] 1
class(counts(sce_zeisel)[grep("Syne1", rownames(sce_zeisel)), ])
# [1] "matrix" "array"
counts(sce_zeisel)[grep("Syne1", rownames(sce_zeisel)), ][, 1:3]
#            1772071015_C02 1772071017_G12 1772071017_A05
# Syne1_loc2             11              2              4
# Syne1_loc1              0              0              4

test <- aggregateAcrossFeatures(
    sce_zeisel,
    id = sub("_loc[0-9]+$", "", rownames(sce_zeisel))
)
# 两个loc合为一个
length(grep("Syne1", rownames(test)))
# [1] 1
class(counts(test)[grep("Syne1", rownames(test)), ])
# [1] "numeric"
# 表达量也合二为一
counts(test)[grep("Syne1", rownames(test)), ][1:3]
# 1772071015_C02 1772071017_G12 1772071017_A05
#             11              2              8

# gene annotation, adding Ensembl ID
library(org.Mm.eg.db)
rowData(sce_zeisel)$Ensembl <- mapIds(
    org.Mm.eg.db,
    keys = rownames(sce_zeisel),
    keytype = "SYMBOL",
    column = "ENSEMBL"
)


# qc, 先perCellQCMetrics()，后quickPerCellQC() -----------------------------------
stats <- perCellQCMetrics(
    sce_zeisel,
    subsets = list(Mt = rowData(sce_zeisel)$featureType == "mito")
)
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent")
)
sce_zeisel_filtered <- sce_zeisel[, !qc$discard]

sum(qc$discard)
# [1] 189

dim(sce_zeisel_filtered)
# [1] 19839  2816

##### 使用qc标准对原数据作图
colData(sce_zeisel) <- cbind(colData(sce_zeisel), stats)
sce_zeisel$discard <- qc$discard

# 使用qc标准对原数据作图
gridExtra::grid.arrange(
    plotColData(sce_zeisel, y = "sum", colour_by = "discard") +
    scale_y_log10() + ggtitle("Total count"),
    plotColData(sce_zeisel, y = "detected", colour_by = "discard") +
    scale_y_log10() + ggtitle("Detected features"),
    plotColData(sce_zeisel, y = "subsets_Mt_percent",
        colour_by = "discard") + ggtitle("Mito percent"),
    plotColData(sce_zeisel, y = "altexps_ERCC_percent",
        colour_by = "discard") + ggtitle("ERCC percent"),
    ncol = 2
)

# 查看线粒体含量与文库大小、ERCC的关系
gridExtra::grid.arrange(
    plotColData(sce_zeisel, x = "sum", y = "subsets_Mt_percent",
        colour_by = "discard") + scale_x_log10() +
    ggtitle("Mito percent to total count"),
    plotColData(sce_zeisel, x = "altexps_ERCC_percent",
        y = "subsets_Mt_percent", colour_by = "discard") +
    ggtitle("Mito percent to ERCC percent"),
    ncol = 2
)

# 检查被过滤的gene
colSums(as.matrix(qc))
#              low_lib_size            low_n_features high_altexps_ERCC_percent
#                         0                         3                        65
#   high_subsets_Mt_percent                   discard
#                       128                       189


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(1000)
cluster_zeisel <- quickCluster(sce_zeisel_filtered)
sce_zeisel_filtered <- computeSumFactors(
    sce_zeisel_filtered,
    cluster = cluster_zeisel
)
# logNormCounts()
sce_zeisel_filtered <- logNormCounts(sce_zeisel_filtered)

summary(sizeFactors(sce_zeisel_filtered))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.1186  0.4860  0.8314  1.0000  1.3209  4.5090

# 比较librarySizeFactors()（常规方法）、
# computeSumFactors()（去卷积方法）进行归一化的差异
plot(
    librarySizeFactors(sce_zeisel_filtered),
    sizeFactors(sce_zeisel_filtered),
    pch = 16,
    xlab = "Library size factors",
    ylab = "Deconvolution factors",
    log = "xy"
)


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
library(scran)
set.seed(1001)
# with spike-in
dec_zeisel_spike <- modelGeneVarWithSpikes(
    sce_zeisel_filtered,
    "ERCC",
)
top_hvgs_zeisel <- getTopHVGs(dec_zeisel_spike, prop = 0.1)
length(top_hvgs_zeisel)
# [1] 1816

# 查看方差大小
plot(
    dec_zeisel_spike$mean,
    dec_zeisel_spike$total,
    main = "STRT-seq", pch = 16, cex = 0.5,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
cur_fit <- metadata(dec_zeisel_spike)
points(cur_fit$mean, cur_fit$var, col = "red", pch = 16)
curve(cur_fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


# dimension reduce, using three methods ---------------------------------------
##### PCA
library(BiocSingular)
set.seed(101011001)
sce_zeisel_filtered <- denoisePCA(
    sce_zeisel_filtered,
    technical = dec_zeisel_spike,
    subset.row = top_hvgs_zeisel
)
dim(reducedDim(sce_zeisel_filtered, "PCA"))
# [1] 2816   50

##### t-SNE
sce_zeisel_filtered <- runTSNE(sce_zeisel_filtered, dimred = "PCA")
dim(reducedDim(sce_zeisel_filtered, "TSNE"))
# [1] 2816   2
sce_zeisel_filtered
# class: SingleCellExperiment
# dim: 19839 2816
# metadata(0):
# assays(2): counts logcounts
# rownames(19839): 0610005C13Rik 0610007N19Rik ... Zzef1 Zzz3
# rowData names(2): featureType Ensembl
# colnames(2816): 1772071015_C02 1772071017_G12 ... 1772063068_D01 1772066098_A12
# colData names(11): tissue group # ... level2class sizeFactor
# reducedDimNames(2): PCA TSNE
# altExpNames(2): ERCC repeat


# clustering, graph-based -----------------------------------------------------
snn_gr_zeisel <- buildSNNGraph(sce_zeisel_filtered, use.dimred = "PCA")
# 鉴定cluster
cluster_zeisel <- igraph::cluster_walktrap(snn_gr_zeisel)$membership
table(cluster_zeisel)
# cluster_zeisel
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14
# 283 451 114 143 599 167 191 128 350  70 199  58  39  24

colLabels(sce_zeisel_filtered) <- factor(cluster_zeisel)

# 绘制t-SNE图
plotTSNE(sce_zeisel_filtered, colour_by = "label")


# detecting markers -----------------------------------------------------------
markers_zeisel <- findMarkers(
    sce_zeisel_filtered,
    groups = colLabels(sce_zeisel_filtered),
    direction = "up"
)

# use cluster 1 as an explaination
chosen_cluster <- "1"
markers_cluster_1 <- markers_zeisel[[chosen_cluster]]

colnames(markers_cluster_1)
# [1] "Top"           "p.value"       "FDR"           "summary.logFC" "logFC.2"      
# [6] "logFC.3"       "logFC.4"       "logFC.5"       "logFC.6"       "logFC.7"      
# [11] "logFC.8"       "logFC.9"       "logFC.10"      "logFC.11"      "logFC.12"     
# [16] "logFC.13"      "logFC.14"

head(markers_cluster_1[, 1:4], 5)
# DataFrame with 5 rows and 4 columns
#              Top      p.value          FDR summary.logFC
#        <integer>    <numeric>    <numeric>     <numeric>
# Atp1a3         1 1.45982e-282 7.24035e-279       3.45669
# Celf4          1 2.27030e-246 4.50404e-243       3.10465
# Gad1           1 7.44925e-232 1.34351e-228       4.57719
# Gad2           1 2.88086e-207 3.57208e-204       4.25393
# Mllt11         1 1.72982e-249 3.81309e-246       2.88363

# 使用cluster 1的top 10 gene绘制热图
markers_cluster_1_top_10 <- rownames(
    markers_cluster_1
)[markers_cluster_1$Top <= 10]
length(markers_cluster_1_top_10)
# [1] 58

plotHeatmap(
    sce_zeisel_filtered,
    features = markers_cluster_1_top_10,
    order_columns_by = "label"
)

# 选取markers_cluster_1的前50个gene，基于logFC绘制热图
library(pheatmap)
log_fc <- getMarkerEffects(markers_cluster_1[1:50, ])

log_fc[1:4, 1:4]
#                 2          3         4        5
# Atp1a3 0.03985679 0.08939429 1.2413877 3.456688
# Celf4  0.38867160 0.61450233 0.8693342 3.104649
# Gad1   4.53927507 4.30032797 4.0503051 4.472360
# Gad2   4.23224869 3.88846542 3.7695559 4.169020

pheatmap(
    log_fc,
    breaks = seq(-5, 5, length.out = 101)
)
