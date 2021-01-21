########################################
# infering cell cycle
# date: 2021.01.21 - 01.21
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-11
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
library(scRNAseq)
sce_416b <- LunSpikeInData(which = "416b")
sce_416b$block <- factor(sce_416b$block)
dim(sce_416b)
# [1] 46604   192


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
rowData(sce_416b)$SEQNAME <- mapIds(
    ens_mm_v97,
    keys = rownames(sce_416b),
    keytype = "GENEID",
    column = "SEQNAME"
)

library(scater)
rownames(sce_416b) <- uniquifyFeatureNames(
    rowData(sce_416b)$ENSEMBL,
    rowData(sce_416b)$SYMBOL
)


# qc, mitochondrial, ERCC, and batch ------------------------------------------
is_mito <- which(rowData(sce_416b)$SEQNAME == "MT")
stats <- perCellQCMetrics(
    sce_416b,
    subsets = list(Mt = is_mito)
)
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent"),
    batch = sce_416b$block
)
sce_416b_filterd <- sce_416b[, !qc$discard]
dim(sce_416b_filterd)
# [1] 46604   185


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
# 未使用quichCluster()
sce_416b_filterd <- computeSumFactors(sce_416b_filterd)
# logNormCounts()
sce_416b_filterd <- logNormCounts(sce_416b_filterd)


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
library(scran)
set.seed(1001)
# with spike-in
dec_416b_spike <- modelGeneVarWithSpikes(sce_416b_filterd, "ERCC")
top_hvgs_416b <- getTopHVGs(dec_416b_spike, prop = 0.1)
length(top_hvgs_416b)
# [1] 1109


# dimension reduce, using three methods ---------------------------------------
##### PCA
set.seed(101010011)
sce_416b_filterd <- denoisePCA(
    sce_416b_filterd,
    subset.row = top_hvgs_416b,
    technical = dec_416b_spike
)
dim(reducedDim(sce_416b_filterd, "PCA"))
# [1] 185    42

##### t-SNE
set.seed(1011101001)
sce_416b_filterd <- runTSNE(sce_416b_filterd, dimred = "PCA")
dim(reducedDim(sce_416b_filterd, "TSNE"))
# [1] 185    2


# hierarchical clustering -----------------------------------------------------
# 利用前几个PCs计算细胞间的距离
dist_416b <- dist(reducedDim(sce_416b_filterd, "PCA"))

# 使用Ward's方法进行聚类
tree_416b <- hclust(dist_416b, "ward.D2")

tree_416b$labels <- seq_along(tree_416b$labels)

# 对树进行修剪
# install.packages("dynamicTreeCut")
library(dynamicTreeCut)
clust_416b <- cutreeDynamic(
    tree_416b,
    distM = as.matrix(dist_416b),
    minClusterSize = 10,
    deepSplit = 1
)

table(clust_416b)
# clust_416b
#  1  2  3  4
# 79 65 28 13

colLabels(sce_416b_filterd) <- factor(clust_416b)

sce_416b_filterd
# class: SingleCellExperiment
# dim: 46604 185
# metadata(0):
# assays(2): counts logcounts
# rownames(46604): 4933401J01Rik Gm26206 ... CAAA01147332.1 CBFB-MYH11-mcherry
# rowData names(4): Length ENSEMBL SYMBOL SEQNAME
# colnames(185): SLX-9555.N701_S502.C89V9ANXX.s_1.r_1
#   SLX-9555.N701_S503.C89V9ANXX.s_1.r_1 ... SLX-11312.N712_S507.H5H5YBBXX.s_8.r_1
#   SLX-11312.N712_S517.H5H5YBBXX.s_8.r_1
# colData names(11): Source Name cell line ... sizeFactor label
# reducedDimNames(2): PCA TSNE
# altExpNames(2): ERCC SIRV


# infering cell cycle ---------------------------------------------------------
##### 使用细胞周期蛋白
## 先找细胞周期蛋白gene，再对照数据集
cyclin_genes_index <- grep("^Ccn[abde][0-9]$", rowData(sce_416b_filterd)$SYMBOL)
cyclin_genes <- rownames(sce_416b_filterd)[cyclin_genes_index]
cyclin_genes
# [1] "Ccnb3" "Ccna2" "Ccna1" "Ccne2" "Ccnd2" "Ccne1" "Ccnd1" "Ccnb2" "Ccnb1" "Ccnd3"

# 绘制热图
library(scater)
plotHeatmap(
    sce_416b_filterd,
    order_columns_by = "label",
    cluster_rows = FALSE,
    features = sort(cyclin_genes)
)

## 先找marker gene，再对应细胞周期蛋白
# 寻找marker gene
library(scran)
markers_416b <- findMarkers(
    sce_416b_filterd,
    subset.row = cyclin_genes,
    test.type = "wilcox",
    direction = "up"
)
# 4个cluster都有相应的细胞周期蛋白gene排名
markers_416b
# List of length 4
# names(4): 1 2 3 4

# 对应细胞周期蛋白
# 查看cluster 4为例
markers_416b[[4]]
# DataFrame with 10 rows and 7 columns
#             Top     p.value         FDR summary.AUC     AUC.1     AUC.2     AUC.3
#       <integer>   <numeric>   <numeric>   <numeric> <numeric> <numeric> <numeric>
# Ccna2         1 2.42516e-08 2.42516e-07    0.987342  0.987342  0.597633  0.923077
# Ccnd1         1 7.28308e-05 1.82077e-04    0.859172  0.395326  0.859172  0.777473
# Ccnb2         2 1.34077e-06 4.46923e-06    0.926972  0.926972  0.772781  0.870879
# Ccnb1         2 5.36875e-07 2.68437e-06    0.939630  0.939630  0.492308  0.923077
# Ccna1         4 2.24034e-02 4.48067e-02    0.538462  0.538462  0.495858  0.538462
# Ccne2         5 3.38530e-02 5.64217e-02    0.665531  0.665531  0.476923  0.446429
# Ccne1         6 3.97636e-01 5.68052e-01    0.585686  0.585686  0.379882  0.483516
# Ccnd2         8 9.99994e-01 1.00000e+00    0.280220  0.120740  0.286391  0.280220
# Ccnd3         8 9.96849e-01 1.00000e+00    0.389484  0.389484  0.278107  0.233516
# Ccnb3        10 1.00000e+00 1.00000e+00    0.500000  0.500000  0.500000  0.500000

##### 辅助参考数据
library(scRNAseq)
sce_ref <- BuettnerESCData()
sce_ref
# class: SingleCellExperiment
# dim: 38293 288
# metadata(0):
# assays(1): counts
# rownames(38293): ENSMUSG00000000001 ENSMUSG00000000003 ... ENSMUSG00000097934
#   ENSMUSG00000097935
# rowData names(3): EnsemblTranscriptID AssociatedGeneName GeneLength
# colnames(288): G1_cell1_count G1_cell2_count ... G2M_cell95_count G2M_cell96_count
# colData names(1): phase
# reducedDimNames(0):
# altExpNames(1): ERCC

table(sce_ref$phase)
# G1 G2M   S
# 96  96  96

# 再加入GO中与细胞周期相关的gene，然后二者与自身数据的gene取交集
library(org.Mm.eg.db)
cycle_anno <- select(
    org.Mm.eg.db,
    keytype = "GOALL",
    keys = "GO:0007049",
    columns = "ENSEMBL"
)[, "ENSEMBL"]
# 取交集
candidates <- Reduce(
    intersect,
    list(rownames(sce_ref), rowData(sce_416b_filterd)$ENSEMBL, cycle_anno)
)
str(candidates)
# chr [1:1605] "ENSMUSG00000000001" "ENSMUSG00000000028" ...

# 用三者交集的gene，带入到注释好的参考数据集中，继续寻找marker gene
sce_ref <- logNormCounts(sce_ref)
phase_stats <- pairwiseWilcox(
    logcounts(sce_ref),
    sce_ref$phase,
    direction = "up",
    subset.row = candidates
)
cycle_markers <- getTopMarkers(phase_stats[[1]], phase_stats[[2]])

# 将得到的marker gene提供给sce_416b_filtered数据
test_data <- logcounts(sce_416b_filterd)
# 由于sce_ref的gene ID是Ensembl，所以这里也使用Ensembl
rownames(test_data) <- rowData(sce_416b_filterd)$ENSEMBL

library(SingleR)
assignments <- SingleR(
    test_data,
    ref = sce_ref,
    label = sce_ref$phase,
    genes = cycle_markers
)
tab <- table(assignments$labels, colLabels(sce_416b_filterd))
tab
#        1  2  3  4
#   G1  71  7 19  1
#   G2M  2 60  1 13
#   S    5  2  4  0

# 再加一步检验，查看两个cluster之间细胞分布的差异大小
chisq.test(tab[, 1:2])
#  Pearson's Chi-squared test
# 
# data:  tab[, 1:2]
# X-squared = 107.91, df = 2, p-value < 2.2e-16

##### 用已发表的分类器
## cyclone分类器
# 加载内置数据
set.seed(100)
library(scran)
mm_pairs <- readRDS(
    system.file(
        "exdata",
        "mouse_cycle_markers.rds",
        package = "scran"
    )
)
class(mm_pairs)
# [1] "list"
names(mm_pairs)
# [1] "G1"  "S"   "G2M"
head(mm_pairs$G1)
#                first             second
# 1 ENSMUSG00000000001 ENSMUSG00000001785
# 2 ENSMUSG00000000001 ENSMUSG00000005470
# 3 ENSMUSG00000000001 ENSMUSG00000012443
# 4 ENSMUSG00000000001 ENSMUSG00000015120
# 5 ENSMUSG00000000001 ENSMUSG00000022033
# 6 ENSMUSG00000000001 ENSMUSG00000023015

# 进行比较
assignments <- cyclone(
    sce_416b_filterd,
    mm_pairs,
    gene.names = rowData(sce_416b_filterd)$ENSEMBL
)
class(assignments)
# [1] "list"
names(assignments)
# [1] "phases"            "scores"            "normalized.scores"
dim(assignments$scores)
# [1] 185   3
head(assignments$scores)
#      G1     S   G2M
# 1 0.585 0.998 0.000
# 2 0.976 0.455 0.000
# 3 0.055 0.518 0.149
# 4 0.000 0.762 1.000
# 5 0.919 0.431 0.002
# 6 0.901 0.948 0.000

# 查看G1期和G2 / M期
plot(
    assignments$scores$G1,
    assignments$scores$G2M,
    xlab = "G1 score",
    ylab = "G2 / M score",
    pch = 16
)

# 进行判断
table(assignments$phases, colLabels(sce_416b_filterd))
#       1  2  3  4
#  G1  74  5 23  0
#  G2M  2 48  0 12
#  S    3 12  5  1

##### 移除细胞周期导致的差异
# method_1, regressBatches()
library(batchelor)
sce_416b_filterd_nocycle_1 <- regressBatches(
    sce_416b_filterd,
    batch = assignments$phases
)

# method_2, 支持block参数的函数
sce_416b_filterd_nocycle_2 <- modelGeneVarWithSpikes(
    sce_416b_filterd,
    "ERCC",
    block = assignments$phases
)
markers_416b_nocycle <- findMarkers(
    sce_416b_filterd,
    block = assignments$phases
)
