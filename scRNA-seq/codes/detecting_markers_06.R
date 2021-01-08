########################################
# detect marker genes
# date: 2021.01.06 - 01.08
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-6
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
#################### PBMC
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
sce_pbmc <- read10xCounts(fname, col.names = TRUE)
dim(sce_pbmc)
# [1]  33694 737280

# gene annotation -------------------------------------------------------------
# ID整合
library(scater)
rownames(sce_pbmc) <- uniquifyFeatureNames(
    rowData(sce_pbmc)$ID,
    rowData(sce_pbmc)$Symbol
)

# 添加位置信息
# BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
location <- mapIds(
    EnsDb.Hsapiens.v86,
    keys = rowData(sce_pbmc)$ID,
    column = "SEQNAME",
    keytype = "GENEID"
)


# detect dropout -----------------------------------------------------------
set.seed(100)
e_out <- emptyDrops(counts(sce_pbmc))
sce_pbmc_filtered <- sce_pbmc[, which(e_out$FDR <= 0.001)]
dim(sce_pbmc_filtered)
# [1] 33694  4300


# qc, especially for mitochondrial --------------------------------------------
stats <- perCellQCMetrics(
    sce_pbmc_filtered,
    subsets = list(Mito = which(location == "MT"))
)
high_mito <- isOutlier(stats$subsets_Mito_percent, type = "higher")
sce_pbmc_final <- sce_pbmc_filtered[, !high_mito]
dim(sce_pbmc_final)
# [1] 33694  3985


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(1000)
clust_pbmc <- quickCluster(sce_pbmc_final)
sce_pbmc_final <- computeSumFactors(
    sce_pbmc_final,
    cluster = clust_pbmc
)
# logNormCounts()
sce_pbmc_final <- logNormCounts(sce_pbmc_final)


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
set.seed(1001)
dec_pbmc_pois <- modelGeneVarByPoisson(sce_pbmc_final)
top_hvgs_pbmc <- getTopHVGs(dec_pbmc_pois, prop = 0.1)
length(top_hvgs_pbmc)
# [1] 1599


# dimension reduce, using three methods ---------------------------------------
##### PCA
set.seed(10000)
sce_pbmc_final <- denoisePCA(
    sce_pbmc_final,
    subset.row = top_hvgs_pbmc,
    technical = dec_pbmc_pois
)
dim(reducedDim(sce_pbmc_final, "PCA"))
# [1] 3985    9

##### t-SNE
set.seed(100000)
sce_pbmc_final <- runTSNE(sce_pbmc_final, dimred = "PCA")
dim(reducedDim(sce_pbmc_final, "TSNE"))
# [1] 3985    2

##### UMAP
set.seed(1000000)
sce_pbmc_final <- runUMAP(sce_pbmc_final, dimred = "PCA")
dim(reducedDim(sce_pbmc_final, "UMAP"))
# [1] 3985    2


# clustering ------------------------------------------------------------------
##### graph-based clustering
# 使用PCA的前几个PCs构建SNNG
library(scran)
g_pbmc_k10 <- buildSNNGraph(sce_pbmc_final, k = 10, use.dimred = "PCA")

# 鉴定cluster
clust_pbmc_k10 <- igraph::cluster_walktrap(g_pbmc_k10)$membership
table(clust_pbmc_k10)
# clust_pbmc_k10
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
# 205 508 541  56 374 125  46 432 302 867  47 155 166  61  84  16

# 把cluster信息存成SingleCellExperiment对象的一个因子
library(scater)
colLabels(sce_pbmc_final) <- factor(clust_pbmc_k10)

sce_pbmc_final
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


# detecting marker gene -------------------------------------------------------
##### method_1, t-test
library(scran)
markers_pbmc <- findMarkers(
    sce_pbmc_final,
    groups = colLabels(sce_pbmc_final),
    test.type = "t"
)
markers_pbmc
# List of length 16
# names(16): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16

# use cluster 9 as an explaination
chosen_cluster <- "9"
interesting_markers <- markers_pbmc[[chosen_cluster]]
colnames(interesting_markers)
# [1] "Top"           "p.value"       "FDR"           "summary.logFC" "logFC.1"
# [6] "logFC.2"       "logFC.3"       "logFC.4"       "logFC.5"       "logFC.6"
# [11] "logFC.7"       "logFC.8"       "logFC.10"      "logFC.11"      "logFC.12"
# [16] "logFC.13"      "logFC.14"      "logFC.15"      "logFC.16"

interesting_markers[1:10, 1:6]
# DataFrame with 10 rows and 6 columns
#                Top      p.value          FDR summary.logFC   logFC.1   logFC.2
#          <integer>    <numeric>    <numeric>     <numeric> <numeric> <numeric>
# MNDA             1 9.77406e-142 3.57964e-139       2.13463  2.134626  2.094336
# CD74             1 1.04576e-296 2.34906e-293       2.88749  2.673481  2.596252
# HLA-DQA1         1  2.53408e-97  4.49386e-95       1.76140  1.761396  1.678937
# FCN1             1 1.53584e-199 1.29371e-196       2.70777  2.707774  2.693735
# MALAT1           1 1.46693e-155 6.67928e-153      -1.28299 -0.560604 -0.817246
# LYZ              1  0.00000e+00  0.00000e+00       4.72716  4.978850  4.742220
# TYROBP           1  0.00000e+00  0.00000e+00       3.75792  1.473897  3.805600
# FTL              1  0.00000e+00  0.00000e+00       3.12210  3.652643  3.284101
# LGALS2           1 4.11545e-137 1.37293e-134       2.21572  2.215722  2.187566
# RPS27            2 1.96725e-240 2.76185e-237      -1.60059  0.355781 -1.119805

# extract all Top6 genes
top6_markers <- interesting_markers[interesting_markers$Top <= 6, ]
dim(top6_markers)
# [1] 46 19
# 提取这些gene在与各个cluster比较时的logFC
log_fcs <- getMarkerEffects(top6_markers)
dim(log_fcs)
# [1] 46 15
log_fcs[1:4, 1:4]
#                 1        2          3         4
# MNDA     2.134626 2.094336  2.0907383 0.2249112
# CD74     2.673481 2.596252 -1.5832557 0.8555704
# HLA-DQA1 1.761396 1.678937 -0.2155989 0.9291779
# FCN1     2.707774 2.693735  2.7332584 0.8423725

library(pheatmap)
pheatmap(log_fcs, breaks = seq(-5, 5, length.out = 101))

## 只关注上调gene
markers_pbmc_up <- findMarkers(
    sce_pbmc_final,
    groups = colLabels(sce_pbmc_final),
    test.type = "t",
    direction = "up"
)
chosen_cluster <- "9"
interesting_markers_up <- markers_pbmc_up[[chosen_cluster]]
interesting_markers_up[1:10, 1:4]
# DataFrame with 10 rows and 4 columns
#                Top      p.value          FDR summary.logFC
#          <integer>    <numeric>    <numeric>     <numeric>
# MNDA             1 4.88703e-142 4.33325e-139       2.13463
# CD74             1 5.22880e-297 1.95755e-293       2.88749
# HLA-DQA1         1  1.26704e-97  5.08234e-95       1.76140
# FCN1             1 7.67919e-200 1.17610e-196       2.70777
# LYZ              1  0.00000e+00  0.00000e+00       4.72716
# CST3             1  0.00000e+00  0.00000e+00       3.75562
# TYROBP           1  0.00000e+00  0.00000e+00       3.75792
# FTL              1  0.00000e+00  0.00000e+00       3.12210
# LGALS2           1 2.05773e-137 1.69105e-134       2.21572
# PRELID1          2 5.52807e-108 3.15700e-105       1.61240

# 参数lfc设置logFC阈值进行过滤
markers_pbmc_up_2 <- findMarkers(
    sce_pbmc_final,
    groups = colLabels(sce_pbmc_final),
    test.type = "t",
    direction = "up",
    lfc = 1
)
chosen_cluster <- "9"
interesting_markers_up_2 <- markers_pbmc_up_2[[chosen_cluster]]
interesting_markers_up_2[1:10, 1:4]
# DataFrame with 10 rows and 4 columns
#                Top      p.value          FDR summary.logFC
#          <integer>    <numeric>    <numeric>     <numeric>
# HLA-DRB1         1 3.27481e-176 1.10342e-172       3.07301
# HLA-DPB1         1  2.41337e-93  3.12754e-90       2.47239
# FCN1             1 1.31080e-129 2.32452e-126       2.81479
# LYZ              1 1.11607e-280 9.40121e-277       4.72716
# TYROBP           1  0.00000e+00  0.00000e+00       3.75792
# FTL              1 7.54749e-310 1.27153e-305       3.12210
# CTSS             2 9.70475e-265 6.53984e-261       3.28621
# MNDA             2  1.46900e-71  1.49989e-68       2.17427
# PRELID1          2  1.15869e-29  7.09837e-27       1.61240
# AIF1             2 2.68067e-148 6.94787e-145       2.88474

# 绘制根据direction = "up" + lfc = 1过滤后的cluster9的marker gene热图
top5_markers <- interesting_markers_up_2[interesting_markers_up_2$Top <= 5, ]
dim(top5_markers)
# [1] 29 19
# 提取这些gene在与各个cluster比较时的logFC
log_fcs_2 <- getMarkerEffects(top5_markers)
dim(log_fcs_2)
# [1] 29 15
log_fcs_2[1:4, 1:4]
#                 1        2          3         4
# HLA-DRB1 3.103054 2.932936 -0.1971619 1.0687736
# HLA-DPB1 2.334938 2.350922 -0.6860820 0.9979761
# FCN1     2.707774 2.693735  2.7332584 0.8423725
# LYZ      4.978850 4.742220  4.7271646 0.5008614

library(pheatmap)
pheatmap(log_fcs_2, breaks = seq(-5, 5, length.out = 101))

## 寻找cluster特异的marker gene
# pval.type = "all" + direction = "up"
markers_pbmc_up_3 <- findMarkers(
    sce_pbmc_final,
    groups = colLabels(sce_pbmc_final),
    test.type = "t",
    pval.type = "all",
    direction = "up",
)
chosen_cluster <- "9"
interesting_markers_up_3 <- markers_pbmc_up_3[[chosen_cluster]]
interesting_markers_up_3[1:10, 1:4]
# DataFrame with 10 rows and 4 columns
#             p.value       FDR summary.logFC   logFC.1
#           <numeric> <numeric>     <numeric> <numeric>
# LGALS2  5.86255e-05         1     0.2751608 2.2157218
# PID1    3.83746e-04         1     0.1637791 0.4152300
# FPR3    4.21811e-03         1     0.0248670 0.0272235
# ATP6V1F 1.06877e-02         1     0.1321922 0.5028556
# GSTP1   1.51436e-02         1     0.1403132 1.4078303
# FOLR2   2.18023e-02         1     0.0456450 0.1098009
# LGALS3  2.25770e-02         1     0.1806593 1.1678046
# TNFAIP2 2.56018e-02         1     0.0926097 0.5456399
# ODF3B   2.86300e-02         1     0.0872068 0.4883289
# DENND6B 4.07048e-02         1     0.0770453 0.2351185

##### method_2, 使用Wilcoxon检验寻找marker gene
markers_pbmc_wilcox <- findMarkers(
    sce_pbmc_final,
    groups = colLabels(sce_pbmc_final),
    test.type = "wilcox",
    direction = "up",
)
markers_pbmc_wilcox
# List of length 16
# names(16): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16

chosen_cluster <- "9"
interesting_markers_wilcox <- markers_pbmc_wilcox[[chosen_cluster]]
interesting_markers_wilcox[1:10, 1:4]
# DataFrame with 10 rows and 4 columns
#              Top      p.value          FDR summary.AUC
#        <integer>    <numeric>    <numeric>   <numeric>
# CSTA           1 5.37447e-216 1.81087e-211    0.982569
# CD74           1 1.04715e-146 7.67012e-144    0.998785
# RPL10          1  1.95068e-47  1.96197e-45    0.881925
# FCN1           1 2.34234e-201 1.97307e-197    0.994065
# LYZ            1 3.81471e-154 3.78038e-151    0.999568
# CST3           1 3.14477e-161 4.23839e-158    0.998709
# TYROBP         1 3.27636e-176 6.13298e-173    0.999183
# FTL            1 1.01900e-146 7.62981e-144    0.999263
# CTSS           2 5.88112e-161 7.62147e-158    0.999141
# S100A4         2 6.32516e-146 4.53447e-143    0.996051

##### method_3, 使用binomial检验寻找marker gene
markers_pbmc_binom <- findMarkers(
    sce_pbmc_final,
    groups = colLabels(sce_pbmc_final),
    test.type = "binom",
    direction = "up"
)
markers_pbmc_binom
# List of length 16
# names(16): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16

chosen_cluster <- "9"
interesting_markers_binom <- markers_pbmc_binom[[chosen_cluster]]
interesting_markers_binom[1:10, 1:4]
# DataFrame with 10 rows and 4 columns
#                Top      p.value          FDR summary.logFC
#          <integer>    <numeric>    <numeric>     <numeric>
# RBP7             1  1.01977e-78  6.36296e-76       5.69911
# MNDA             1 6.83896e-117 3.29188e-113       4.00500
# CSTA             1 8.17807e-126 1.37776e-121       4.37953
# HLA-DQA1         1  6.70291e-82  4.80527e-79       3.07999
# FGL2             1 2.91867e-117 1.63903e-113       4.41722
# FCN1             1  6.75826e-94  9.48803e-91       3.13098
# CLEC7A           1 1.51091e-120 1.01817e-116       5.26471
# IFI30            1 3.92909e-128 1.32387e-123       6.06525
# TYROBP           1  2.24733e-52  5.40869e-50       1.96480
# LILRB4           1  1.31740e-34  1.68138e-32       4.53379

library(scater)
top_markers_binom <- head(rownames(interesting_markers_binom))
plotExpression(sce_pbmc_final, x = "label", features = top_markers_binom)

##### 整合t检验、Wilcoxon检验和binomial检验
combined_markers <- multiMarkerStats(
    t = findMarkers(
        sce_pbmc_final,
        groups = colLabels(sce_pbmc_final),
        test.type = "t",
        direction = "up"
    ),
    wilcox = findMarkers(
        sce_pbmc_final,
        groups = colLabels(sce_pbmc_final),
        test.type = "wilcox",
        direction = "up"
    ),
    binom = findMarkers(
        sce_pbmc_final,
        groups = colLabels(sce_pbmc_final),
        test.type = "binom",
        direction = "up"
    )
)
combined_markers
# List of length 16
# names(16): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16

combined_markers[[9]][1:4, 1:15]
# DataFrame with 4 rows and 15 columns
#              Top      p.value          FDR     t.Top wilcox.Top binom.Top    t.p.value
#        <integer>    <numeric>    <numeric> <integer>  <integer> <integer>    <numeric>
# FCN1           1  6.75826e-94  2.07012e-90         1          1         1 7.67919e-200
# TYROBP         1  2.24733e-52  1.08174e-49         1          1         1  0.00000e+00
# MNDA           2 6.83896e-117 7.68106e-113         1          2         1 4.88703e-142
# LGALS2         3 4.22567e-113 3.55949e-109         1          3         1 2.05773e-137
#        wilcox.p.value binom.p.value        t.FDR   wilcox.FDR    binom.FDR t.summary.logFC
#             <numeric>     <numeric>    <numeric>    <numeric>    <numeric>       <numeric>
# FCN1     2.34234e-201   6.75826e-94 1.17610e-196 1.97307e-197  9.48803e-91         2.70777
# TYROBP   3.27636e-176   2.24733e-52  0.00000e+00 6.13298e-173  5.40869e-50         3.75792
# MNDA     1.29794e-207  6.83896e-117 4.33325e-139 2.18663e-203 3.29188e-113         2.13463
# LGALS2   2.01724e-201  4.22567e-113 1.69105e-134 1.97307e-197 1.42380e-109         2.21572
#        wilcox.summary.AUC binom.summary.logFC
#                 <numeric>           <numeric>
# FCN1             0.994065             3.13098
# TYROBP           0.999183             1.96480
# MNDA             0.978838             4.00500
# LGALS2           0.969257             3.94797

##### block一些不重要因素
markers_416b_without_batch <- findMarkers(
    sce_416b_filterd,
    groups = colLabels(sce_416b_filterd),
    test.type = "t",
    direction = "up",
    block = sce_416b_filterd$block
)
