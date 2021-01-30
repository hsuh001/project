########################################
# practice 1, smart-seq2, mouse bone marrow
# date: 2021.01.29 - 01.29
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-1
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
##### 关于scRNAseq数据集
# 早期版本仅内置3个数据集Fluidig、Th2、Allen
library(scRNAseq)
fluidigm <- ReprocessedFluidigmData()       # provides 65 cells
th2 <- ReprocessedTh2Data()                 # provides 96 T helper cells
allen <- ReprocessedAllenData()             # provides 379 cells from the mouse visual cortex

##### 使用sce_416b数据
sce_416b <- LunSpikeInData(which = "416b")
dim(sce_416b)
# [1] 46604   192
sce_416b
# class: SingleCellExperiment
# dim: 46604 192
# metadata(0):
# assays(1): counts
# rownames(46604): ENSMUSG00000102693 ENSMUSG00000064842 ... ENSMUSG00000095742
#   CBFB-MYH11-mcherry
# rowData names(1): Length
# colnames(192): SLX-9555.N701_S502.C89V9ANXX.s_1.r_1
#   SLX-9555.N701_S503.C89V9ANXX.s_1.r_1 ... SLX-11312.N712_S508.H5H5YBBXX.s_8.r_1
#   SLX-11312.N712_S517.H5H5YBBXX.s_8.r_1
# colData names(9): Source Name cell line ... spike-in addition block
# reducedDimNames(0):
# altExpNames(2): ERCC SIRV

# 设置分组信息（来自2个96孔板），为后面的批处理做准备
sce_416b$block <- factor(sce_416b$block)
table(sce_416b$block)
# 20160113 20160325
#       96       96

# 当前关于行（即gene）的信息只有gene length
rowData(sce_416b)
# DataFrame with 46604 rows and 1 column
#                       Length
#                    <integer>
# ENSMUSG00000102693      1070
# ENSMUSG00000064842       110
# ENSMUSG00000051951      6094
# ENSMUSG00000102851       480
# ENSMUSG00000103377      2819
# ...                      ...
# ENSMUSG00000094621       121
# ENSMUSG00000098647        99
# ENSMUSG00000096730      3077
# ENSMUSG00000095742       243
# CBFB-MYH11-mcherry      2998


# gene annotation -------------------------------------------------------------
#BiocManager::install("AnnotationHub")
library(AnnotationHub)
ens_mm_v97 <- AnnotationHub()[["AH73905"]]
rowData(sce_416b)$ENSEMBL <- rownames(sce_416b)
rowData(sce_416b)$SYMBOL <- mapIds(
    ens_mm_v97,
    keys = rownames(sce_416b),
    keytype = "GENEID",
    column = "SYMBOL"
)

# 便于识别线粒体gene
rowData(sce_416b)$SEQNAME <- mapIds(
    ens_mm_v97,
    keys = rownames(sce_416b),
    keytype = "GENEID",
    column = "SEQNAME"
)
table(rowData(sce_416b)$SEQNAME == "MT")
# FALSE  TRUE
# 46004    37

library(scater)
rownames(sce_416b) <- uniquifyFeatureNames(
    rowData(sce_416b)$ENSEMBL,
    rowData(sce_416b)$SYMBOL
)


# qc, mitochondrial, ERCC, and batch ------------------------------------------
is_mito <- which(rowData(sce_416b)$SEQNAME == "MT")
# 计算各个细胞的qc指标
stats <- perCellQCMetrics(
    sce_416b,
    subsets = list(Mt = is_mito)
)
# 根据qc指标判断要去除的cells
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent"),
    batch = sce_416b$block
)
sce_416b_filtered <- sce_416b[, !qc$discard]
dim(sce_416b_filtered)
# [1] 46604   185

##### 使用qc标准对原数据作图
colData(sce_416b) <- cbind(colData(sce_416b), stats)
sce_416b$discard <- qc$discard

# 使用qc标准对原数据作图
gridExtra::grid.arrange(
    plotColData(sce_416b, x = "block", y = "sum", colour_by = "discard") +
    scale_y_log10() + ggtitle("Total count"),
    plotColData(sce_416b, x = "block", y = "detected", colour_by = "discard") +
    scale_y_log10() + ggtitle("Detected features"),
    plotColData(sce_416b, x = "block", y = "subsets_Mt_percent",
        colour_by = "discard") + ggtitle("Mito percent"),
    plotColData(sce_416b, x = "block", y = "altexps_ERCC_percent",
        colour_by = "discard") + ggtitle("ERCC percent"),
    nrow = 2, ncol = 2
)

# 查看线粒体含量与文库大小、ERCC的关系
gridExtra::grid.arrange(
    plotColData(sce_416b, x = "sum", y = "subsets_Mt_percent",
        colour_by = "discard") + scale_x_log10() +
    ggtitle("Mito percent to total count"),
    plotColData(sce_416b, x = "altexps_ERCC_percent", y = "subsets_Mt_percent",
        colour_by = "discard") +
    ggtitle("Mito percent to ERCC percent"),
    ncol = 2
)

# 检查被过滤的gene
colSums(as.matrix(qc))
#              low_lib_size            low_n_features high_altexps_ERCC_percent
#                         5                         0                         2
#   high_subsets_Mt_percent                   discard
#                         2                         7


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
# 未使用quickCluster()
sce_416b_filtered <- computeSumFactors(sce_416b_filtered)
# logNormCounts()
sce_416b_filtered <- logNormCounts(sce_416b_filtered)
summary(sizeFactors(sce_416b_filtered))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.3468  0.7111  0.9207  1.0000  1.1520  3.6037

# 比较librarySizeFactors()（常规方法）、
# computeSumFactors()（去卷积方法）进行归一化的差异
plot(
    librarySizeFactors(sce_416b_filtered),
    sizeFactors(sce_416b_filtered),
    pch = 16,
    xlab = "Library size factors",
    ylab = "Deconvolution factors",
    # 直接将不同类型赋予不同颜色
    # + 1将False(0)、True(1)变成1、2，就可以取颜色
    # 因此black是uninduced、red是induced
    col = c("black", "red")[grepl("induced", sce_416b_filtered$phenotype) + 1],
    log = "xy"
)


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
library(scran)
set.seed(1001)
# with spike-in
dec_416b_spike <- modelGeneVarWithSpikes(
    sce_416b_filtered,
    "ERCC",
    block = sce_416b_filtered$block
)
top_hvgs_416b <- getTopHVGs(dec_416b_spike, prop = 0.1)
length(top_hvgs_416b)
# [1] 1067

# 对两个batch分别作图
par(mfrow = c(1, 2))
blocked_state <- dec_416b_spike$per.block
for (i in colnames(blocked_state)) {
    current <- blocked_state[[i]]
    plot(current$mean, current$total,
        main = i, pch = 16, cex = 0.5,
        xlab = "Mean of log-expression",
        ylab = "Variance of log-expression"
    )
    cur_fit <- metadata(current)
    points(cur_fit$mean, cur_fit$var, col = "red", pch = 16)
    curve(cur_fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
}


# correcting batch effects by linear regression -------------------------------
library(limma)
assay(sce_416b_filtered, "corrected") <- removeBatchEffect(
    logcounts(sce_416b_filtered),
    design = model.matrix(~ sce_416b_filtered$phenotype),
    batch = sce_416b_filtered$block
)


# dimension reduce, using three methods ---------------------------------------
##### PCA
library(scater)
library(BiocSingular)
sce_416b_filtered <- runPCA(
    sce_416b_filtered,
    ncomponents = 10,
    subset_row = top_hvgs_416b,
    exprs_values = "corrected",
    BSPARAM = ExactParam()
)
dim(reducedDim(sce_416b_filtered, "PCA"))
# [1] 185    10

##### t-SNE
set.seed(1010)
sce_416b_filtered <- runTSNE(sce_416b_filtered, dimred = "PCA", perplexity = 10)
dim(reducedDim(sce_416b_filtered, "TSNE"))
# [1] 185    2
sce_416b_filtered
# class: SingleCellExperiment
# dim: 46604 185
# metadata(0):
# assays(3): counts logcounts corrected
# rownames(46604): 4933401J01Rik Gm26206 ... CAAA01147332.1 CBFB-MYH11-mcherry
# rowData names(4): Length ENSEMBL SYMBOL SEQNAME
# colnames(185): SLX-9555.N701_S502.C89V9ANXX.s_1.r_1
#   SLX-9555.N701_S503.C89V9ANXX.s_1.r_1 ... SLX-11312.N712_S507.H5H5YBBXX.s_8.r_1
#   SLX-11312.N712_S517.H5H5YBBXX.s_8.r_1
# colData names(10): Source Name cell line ... block sizeFactor
# reducedDimNames(2): PCA TSNE
# altExpNames(2): ERCC SIRV


# hierarchical clustering -----------------------------------------------------
# 利用前几个PCs计算细胞间的距离
dist_416b <- dist(reducedDim(sce_416b_filtered, "PCA"))

# 使用Ward's方法进行聚类
tree_416b <- hclust(dist_416b, method = "ward.D2")

## 对树进行修剪
# install.packages("dynamicTreeCut")
library(dynamicTreeCut)
cluster_416b <- cutreeDynamic(
    tree_416b,
    distM = as.matrix(dist_416b),
    minClusterSize = 10,
    verbose = 0
)
table(cluster_416b)
# cluster_416b
#  1  2  3  4
# 78 69 24 14

colLabels(sce_416b_filtered) <- factor(cluster_416b)

# 查看分群和批次之间的关系
table(
    Cluster = colLabels(sce_416b_filtered),
    Plate = sce_416b_filtered$block
)
#        Plate
# Cluster 20160113 20160325
#       1       40       38
#       2       37       32
#       3       10       14
#       4        6        8

# 查看分群和表型之间的关系
table(
    Cluster = colLabels(sce_416b_filtered),
    Phenotype = sce_416b_filtered$phenotype
)
#        Phenotype
# Cluster induced CBFB-MYH11 oncogene expression wild type phenotype
#       1                                     78                   0
#       2                                      0                  69
#       3                                      1                  23
#       4                                     14                   0

# 绘制t-SNE图
plotTSNE(sce_416b_filtered, colour_by = "label")

# evaluation of hierarchical clustering, using silhouette width
library(cluster)
sil <- silhouette(cluster_416b, dist = dist_416b)
colnames(sil)
# [1] "cluster"   "neighbor"  "sil_width"

# 设置每个cluster的颜色
cluster_col <- scater:::.get_palette("tableau10medium")
# 如果sil_width>0，就属于和自己接近，否则属于和其他亚群接近
sil_cols <- cluster_col[ifelse(sil[, 3] > 0, sil[, 1], sil[, 2])]
sil_cols <- sil_cols[order(-sil[, 1], sil[, 3])]

plot(
    sil,
    main = paste(length(unique(cluster_416b)), "clusters"),
    border = sil_cols,
    col = sil_cols,
    do.col.sort = FALSE
)


# detecting markers -----------------------------------------------------------
markers_416b <- findMarkers(
    sce_416b_filtered,
    groups = colLabels(sce_416b_filtered),
    block = sce_416b_filtered$block
)
# use cluster 1 as an explaination
chosen_cluster <- "1"
markers_cluster_1 <- markers_416b[[chosen_cluster]]
# cluster 1的top 10
head(markers_cluster_1, 10)
# DataFrame with 10 rows and 7 columns
#             Top     p.value         FDR summary.logFC   logFC.2   logFC.3    logFC.4
#       <integer>   <numeric>   <numeric>     <numeric> <numeric> <numeric>  <numeric>
# Ccna2         1 9.85422e-67 4.59246e-62      -7.13310  -7.13310  -2.20632 -7.3451052
# Cdca8         1 1.01449e-41 1.52514e-38      -7.26175  -6.00378  -2.03841 -7.2617478
# Pirb          1 4.16555e-33 1.95516e-30       5.87820   5.28149   5.87820  0.0352849
# Cks1b         2 2.98233e-40 3.23229e-37      -6.43381  -6.43381  -4.15385 -6.4385323
# Aurkb         2 2.41436e-64 5.62593e-60      -6.94063  -6.94063  -1.65534 -6.4162126
# Myh11         2 1.28865e-46 3.75353e-43       4.38182   4.38182   4.29290  0.9410499
# Mcm6          3 1.15877e-28 3.69887e-26      -5.44558  -5.44558  -5.82130 -3.5804973
# Cdca3         3 5.02047e-45 1.23144e-41      -6.22179  -6.22179  -2.10502 -7.0539510
# Top2a         3 7.25965e-61 1.12776e-56      -7.07811  -7.07811  -2.39123 -6.8297343
# Mcm2          4 1.50854e-33 7.98908e-31      -5.54197  -5.54197  -6.09178 -3.8238103

# 使用cluster 1的top 10 marker gene绘制热图
markers_cluster_1_top_10 <- rownames(
    markers_cluster_1
)[markers_cluster_1$Top <= 10]
length(markers_cluster_1_top_10)
# [1] 29

plotHeatmap(
    sce_416b_filtered,
    features = markers_cluster_1_top_10,
    order_columns_by = "label",
    colour_columns_by = c("label", "block", "phenotype"),
    center = TRUE,
    symmetric = TRUE,
    zlim = c(-5, 5)
)
