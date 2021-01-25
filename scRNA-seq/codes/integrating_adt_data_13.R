########################################
# intergrating antibody-derived tag (ADT) data to analyse protein abundance
# date: 2021.01.23 - 01.25
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
# PBMC 10K
# 数据下载
library(BiocFileCache)
# 在当前工作目录新建raw_data目录
bfc <- BiocFileCache("raw_data", ask = FALSE)
# 提供下载链接
# http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz
raw_path <- bfcrpath(
    bfc,
    file.path(
        "http://cf.10xgenomics.com",
        "samples/cell-exp/3.0.0/pbmc_10k_protein_v3",
        "pbmc_10k_protein_v3_filtered_feature_bc_matrix.tar.gz"
    )
)

# 解压数据至当前工作目录，并新建pbmc10k目录
untar(raw_path, exdir = file.path(getwd(), "pbmc10k"))

library(DropletUtils)
fname <- file.path(getwd(), "pbmc10k/filtered_feature_bc_matrix")
sce_pbmc_10k <- read10xCounts(fname)
dim(sce_pbmc_10k)
# [1] 33555  7865
sce_pbmc_10k
# class: SingleCellExperiment
# dim: 33555 7865
# metadata(1): Samples
# assays(1): counts
# rownames(33555): ENSG00000243485 ENSG00000237613 ... IgG1 IgG2b
# rowData names(3): ID Symbol Type
# colnames: NULL
# colData names(2): Sample Barcode
# reducedDimNames(0):
# altExpNames(0):

names(rowData(sce_pbmc_10k))
# [1] "ID"     "Symbol" "Type"
table(rowData(sce_pbmc_10k)$Type)
# Antibody Capture  Gene Expression
#               17            33538


# 分离ADT数据 ---------------------------------------------------------------------
sce_pbmc_10k <- splitAltExps(
    sce_pbmc_10k,
    rowData(sce_pbmc_10k)$Type
)
altExpNames(sce_pbmc_10k)
# [1] "Antibody Capture"

# 一共17个ADT的信息
altExp(sce_pbmc_10k)
# class: SingleCellExperiment
# dim: 17 7865
# metadata(1): Samples
# assays(1): counts
# rownames(17): CD3 CD4 ... IgG1 IgG2b
# rowData names(3): ID Symbol Type
# colnames: NULL
# colData names(0):
# reducedDimNames(0):
# altExpNames(0):

# ADT表达量是一个稀疏矩阵
counts(altExp(sce_pbmc_10k))[1:3, 1:3]
# 3 x 3 sparse Matrix of class "dgCMatrix"
#                 
# CD3   18  30  18
# CD4  138 119 207
# CD8a  13  19  10

# 将稀疏矩阵转换成常规矩阵
counts(altExp(sce_pbmc_10k)) <- as.matrix(counts(altExp(sce_pbmc_10k)))


# 质控 --------------------------------------------------------------------------
# 根据线粒体水平去掉低质量细胞
library(scater)
is_mito <- grep("^MT-", rowData(sce_pbmc_10k)$Symbol)
stats <- perCellQCMetrics(
    sce_pbmc_10k,
    subsets = list(Mito = is_mito)
)
high_mito <- isOutlier(stats$subsets_Mito_percent, type = "higher")
summary(high_mito)
#    Mode   FALSE    TRUE
# logical    7569     296

# 由于ADT表达量的值一般很大，需要log2转换
# min_diff, minimum difference from the median to consider as an outlier
low_ab <- isOutlier(
    stats$`altexps_Antibody Capture_detected`,
    log = TRUE,
    type = "lower",
    min_diff = 1
)
summary(low_ab)
#    Mode   FALSE    TRUE
# logical    7864       1

# 整合过滤信息
discard <- high_mito | low_ab
table(discard)
discard
# FALSE  TRUE
#  7568   297

sce_pbmc_10k_filterd <- sce_pbmc_10k[, !discard]
dim(sce_pbmc_10k_filterd)
# [1] 33538  7568


# 归一化 -------------------------------------------------------------------------
##### 计算size factor
# method_1, 基于ADTs的文库大小的计算
sf_lib <- librarySizeFactors(altExp(sce_pbmc_10k_filterd))
summary(sf_lib)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.02668  0.52418  0.90485  1.00000  1.26625 22.81740

# method_2, 基于DESeq的计算
# pseudosample
ambient <- rowMeans(counts(altExp(sce_pbmc_10k_filterd)))
# DESeq-like size factors
sf_amb <- medianSizeFactors(altExp(sce_pbmc_10k_filterd), reference = ambient)
summary(sf_amb)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0030  0.5927  0.8282  1.0000  1.1439 41.7972

# method_3, 基于对照的计算
controls <- grep("^Ig", rownames(altExp(sce_pbmc_10k_filterd)))
controls
# [1] 15 16 17
sf_control <- librarySizeFactors(
    altExp(sce_pbmc_10k_filterd),
    subset_row = controls
)
summary(sf_control)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.6854  0.8757  1.0000  1.1423 44.0155

##### 推荐使用第二种基于DESeq的方法
sizeFactors(altExp(sce_pbmc_10k_filterd)) <- sf_amb
sce_pbmc_10k_filterd <- logNormCounts(
    sce_pbmc_10k_filterd,
    use_altexps = TRUE
)

assayNames(sce_pbmc_10k_filterd)
# [1] "counts"    "logcounts"
assayNames(altExp(sce_pbmc_10k_filterd))
# [1] "counts"    "logcounts"


# 聚类 -------------------------------------------------------------------------
# 设置d = NA可以跳过PCA，直接分群
g_adt <- buildSNNGraph(altExp(sce_pbmc_10k_filterd), d = NA)
clusters_adt <- igraph::cluster_walktrap(g_adt)$membership

# 可视化
set.seed(1010010)
altExp(sce_pbmc_10k_filterd) <- runTSNE(altExp(sce_pbmc_10k_filterd))
colLabels(altExp(sce_pbmc_10k_filterd)) <- factor(clusters_adt)
plotTSNE(
    altExp(sce_pbmc_10k_filterd),
    colour_by = "label",
    text_by = "label",
    text_col = "red"
)

# 绘制热图检查cluster
se_averaged <- sumCountsAcrossCells(
    altExp(sce_pbmc_10k_filterd),
    clusters_adt,
    exprs_values = "logcounts",
    average = TRUE
)
dim(se_averaged)
# [1] 17 28

assay(se_averaged)[1:3, 1:3]
#             1        2        3
# CD3  4.621559 9.556522 5.673089
# CD4  4.348560 9.795961 9.078986
# CD8a 4.797436 5.262174 5.306986

# 绘制热图
library(pheatmap)
averaged <- assay(se_averaged)
pheatmap(
    averaged - rowMeans(averaged),
    breaks = seq(-3, 3, length.out = 101)
)


# 利用转录组的分群找蛋白相关的marker gene ---------------------------------------------------
##### 继续分亚群
## 对ADT得到的cluster进行循环分亚群
set.seed(101010)
all_sce_pbmc_10k <- quickSubCluster(
    sce_pbmc_10k_filterd,
    clusters_adt,
    prepFUN = function(x) {
        dec <- modelGeneVar(x)
        top <- getTopHVGs(dec, prop = 0.1)
        x <- runPCA(x, subset_row = top, ncomponents = 25)
    },
    clusterFUN = function(x) {
        g_trans <- buildSNNGraph(x, use.dimred = "PCA")
        igraph::cluster_walktrap(g_trans)$membership
    }
)
# You're computing too large a percentage of total singular values, use a standard svd instead.
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :

all_sce_pbmc_10k
# List of length 28
# names(28): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28

# 查看cluster 3
all_sce_pbmc_10k[[3]]
# class: SingleCellExperiment
# dim: 33538 1828
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#   ENSG00000268674
# rowData names(3): ID Symbol Type
# colnames: NULL
# colData names(4): Sample Barcode sizeFactor subcluster
# reducedDimNames(1): PCA
# altExpNames(1): Antibody Capture

table(all_sce_pbmc_10k[[3]]$subcluster)
# 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8
#  90 772  25 178 695  44  19   5

## 查看每个ADT分群结果各自对应的亚群
# 每个ADT群的细胞数量
ncells <- sapply(all_sce_pbmc_10k, ncol, simplify = TRUE)
ncells
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18
#  888   92 1828 1186  512   55   71   85  999  571   70  148   73   39  409   21   39   69
#   19   20   21   22   23   24   25   26   27   28 
#   74  126   32   50   13   17   29   38   22   12

# 每个ADT群的亚群数量
nsubclusters <- sapply(
    all_sce_pbmc_10k,
    function(x) {
        length(unique(x$subcluster))
    },
    simplify = TRUE
)
nsubclusters
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
# 21  3  8  7  3  3  2  3  9  8  3  2  3  1 12  1  1  2  3  3  1  2  1  1  1  1  1  1

# 绘图查看
plot(
    ncells,
    nsubclusters,
    xlab = "Number of cells",
    type = "n",
    ylab = "Number of subclusters",
    log = "xy"
)
text(ncells, nsubclusters, names(all_sce_pbmc_10k))

# 以cluster 12为例
of_interest <- "12"
# 左边是GZMH，右边是GZHK
plotExpression(
    all_sce_pbmc_10k[[of_interest]],
    x = "subcluster",
    features = c("ENSG00000100450", "ENSG00000113088")
)

# 补充检查
of_interest <- "12"
sce_cd8 <- all_sce_pbmc_10k[[of_interest]]
plotExpression(
    altExp(sce_cd8),
    x = I(sce_cd8$subcluster),
    features = c("CD3", "CD8a")
)

##### 利用转录组的分群找蛋白相关的marker gene
## 对转录组数据进行快速分群
# sce_pbmc_10k_filterd <- logNormCounts(sce_pbmc_10k_filterd)
dec_pbmc_10k <- modelGeneVar(sce_pbmc_10k_filterd)
top_hvg_pbmc_10k <- getTopHVGs(dec_pbmc_10k, prop = 0.1)
length(top_hvg_pbmc_10k)
# [1] 1387

set.seed(1001010)
sce_pbmc_10k_filterd <- runPCA(
    sce_pbmc_10k_filterd,
    subset_row = top_hvg_pbmc_10k,
    ncomponents = 25
)

g_pbmc_10k <- buildSNNGraph(sce_pbmc_10k_filterd, use.dimred = "PCA")
clusters_pbmc_10k <- igraph::cluster_walktrap(g_pbmc_10k)$membership
colLabels(sce_pbmc_10k_filterd) <- factor(clusters_pbmc_10k)

set.seed(1000010)
sce_pbmc_10k_filterd <- runTSNE(sce_pbmc_10k_filterd, dimred = "PCA")
plotTSNE(
    sce_pbmc_10k_filterd,
    colour_by = "label",
    text_by = "label"
)

## ADT数据 + 分群结果，找marker gene
markers_pbmc_10k <- findMarkers(
    altExp(sce_pbmc_10k_filterd),
    colLabels(sce_pbmc_10k_filterd)
)
markers_pbmc_10k
# List of length 18
# names(18): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

of_interest_markers <- markers_pbmc_10k[[16]]
pheatmap(
    getMarkerEffects(of_interest_markers),
    breaks = seq(-3, 3, length.out = 101)
)
