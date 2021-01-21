########################################
# differential expression and differential abundance analysis
# date: 2021.01.20 - 01.21
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-10
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
library(scRNAseq)
sce_mam <- BachMammaryData(samples = "G_1")


# gene annotation -------------------------------------------------------------
library(scater)
rownames(sce_mam) <- uniquifyFeatureNames(
    rowData(sce_mam)$Ensembl,
    rowData(sce_mam)$Symbol
)

library(AnnotationHub)
ens_mm_v97 <- AnnotationHub()[["AH73905"]]
rowData(sce_mam)$SEQNAME <- mapIds(
    ens_mm_v97,
    keys = rowData(sce_mam)$Ensembl,
    keytype = "GENEID",
    column = "SEQNAME"
)


# quality control -------------------------------------------------------------
is_mito <- rowData(sce_mam)$SEQNAME == "MT"
stats <- perCellQCMetrics(
    sce_mam,
    subsets = list(Mito = which(is_mito))
)
qc <- quickPerCellQC(stats, percent_subsets = "subsets_Mito_percent")
sce_mam_filtered <- sce_mam[, !qc$discard]


# normalization by deconvolution ----------------------------------------------
library(scran)
set.seed(101000110)
cluster_mam <- quickCluster(sce_mam_filtered)
sce_mam_final <- computeSumFactors(sce_mam_filtered, clusters = cluster_mam)
sce_mam_final <- logNormCounts(sce_mam_final)


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
set.seed(00010101)
dec_mam_pois <- modelGeneVarByPoisson(sce_mam_final)
top_hvgs_mam <- getTopHVGs(dec_mam_pois, prop = 0.1)
length(top_hvgs_mam)
# [1] 1366


# dimension reduce ------------------------------------------------------------
library(BiocSingular)
set.seed(101010011)
##### PCA
sce_mam_final <- denoisePCA(
    sce_mam_final,
    technical = dec_mam_pois,
    subset.row = top_hvgs_mam
)

##### t-SNE
sce_mam_final <- runTSNE(sce_mam_final, dimred = "PCA")


# clustering ------------------------------------------------------------------
##### graph-based clustering
# 使用PCA的前几个PCs构建SNNG
library(scran)
snn_gr_mam <- buildSNNGraph(sce_mam_final, k = 25, use.dimred = "PCA")
colLabels(sce_mam_final) <- factor(
    igraph::cluster_walktrap(snn_gr_mam)$membership
)
sce_mam_final
# class: SingleCellExperiment
# dim: 27998 2772
# metadata(0):
# assays(2): counts logcounts
# rownames(27998): Xkr4 Gm1992 ... Vmn2r122 CAAA01147332.1
# rowData names(3): Ensembl Symbol SEQNAME
# colnames: NULL
# colData names(5): Barcode Sample Condition sizeFactor label
# reducedDimNames(2): PCA TSNE
# altExpNames(0):


# detecting doublets -----------------------------------------------------------
##### method_1，基于分群结果的检测
library(scran)
table(sce_mam_final$label)
#   1   2   3   4   5   6   7   8   9  10
# 550 799 716 452  24  84  52  39  32  24

# sce_mam_final一共分了10个cluster，因此，下面的结果也是10行，每一群都做了一次检测
doublet_out <- doubletCluster(sce_mam_final)
dim(doublet_out)
# [1] 10  9

doublet_out
# DataFrame with 10 rows and 9 columns
#        source1     source2         N        best     p.value lib.size1 lib.size2       prop
#    <character> <character> <integer> <character>   <numeric> <numeric> <numeric>  <numeric>
# 6            2           1        13       Pcbp2 1.28336e-03  0.811531  0.516399 0.03030303
# 2           10           3       109        Pigr 4.34790e-21  0.619865  1.411579 0.28823954
# 4            6           5       111       Cotl1 1.09709e-08  1.540751  0.688651 0.16305916
# 5           10           7       139        Gde1 9.30195e-12  1.125474  1.167854 0.00865801
# 10           8           5       191       Krt18 5.54539e-20  0.888432  0.888514 0.00865801
# 7            8           5       270    AF251705 3.29661e-24  0.856192  0.856271 0.01875902
# 9            8           5       295       Fabp4 2.21523e-32  0.655624  0.655685 0.01154401
# 8           10           9       388      Col1a1 6.82664e-32  1.125578  1.525264 0.01406926
# 1            8           6       513       Acta2 1.07294e-24  0.865449  1.936489 0.19841270
# 3            6           5       530      Sapcd2 6.08574e-16  0.872951  0.390173 0.25829726
#                                    all.pairs
#                              <DataFrameList>
# 6       2:1:13:...,4:1:14:...,3:1:28:...,...
# 2   10:3:109:...,5:3:194:...,8:3:205:...,...
# 4   6:5:111:...,10:6:146:...,6:2:209:...,...
# 5  10:7:139:...,9:7:153:...,10:9:168:...,...
# 10   8:5:191:...,9:5:224:...,6:5:229:...,...
# 7   8:5:270:...,9:5:337:...,10:5:362:...,...
# 9   8:5:295:...,10:8:323:...,5:1:338:...,...
# 8  10:9:388:...,9:7:400:...,10:7:436:...,...
# 1    8:6:513:...,9:6:609:...,6:5:854:...,...
# 3    6:5:530:...,5:2:593:...,5:4:700:...,...

# 主要看N的数量。因此，可以按N对结果排序
library(scater)
is_doublet <- isOutlier(doublet_out$N, type = "lower", log = TRUE)
chosen_doublet <- rownames(doublet_out)[is_doublet]
chosen_doublet
# [1] "6"

# 寻找cluster 6的marker gene
markers_mam <- findMarkers(sce_mam_final, direction = "up")
doublet_markers <- markers_mam[[chosen_doublet]]
dim(doublet_markers)
# [1] 27998    13

# 提取Top10的gene
library(scater)
top10_doublet_markers <- rownames(doublet_markers)[doublet_markers$Top <= 10]
length(top10_doublet_markers)
# [1] 43

# 绘制Top10的gene的表达谱热图
plotHeatmap(
    sce_mam_final,
    order_columns_by = "label",
    features = top10_doublet_markers,
    center = TRUE,
    symmetric = TRUE,
    zlim = c(-5, 5)
)

# 绘制Acta2、Csn2的在各个cluster中的表达量情况
plotExpression(
    sce_mam_final,
    features = c("Acta2", "Csn2"),
    x = "label",
    colour_by = "label"
)

##### method_2，基于模拟推断的检测
library(BiocSingular)
set.seed(100)

doublet_density <- doubletCells(
    sce_mam_final,
    subset.row = top_hvgs_mam,
    d = ncol(reducedDim(sce_mam_final))
)
summary(doublet_density)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#     0.00     7.63    21.04   395.42    49.30 39572.52

# 绘制doublet score
sce_mam_final$DoubletScore <- log10(doublet_density + 1)
plotTSNE(sce_mam_final, colour_by = "DoubletScore")

# 绘制cluster的t-SNE结果
plotTSNE(sce_mam_final, colour_by = "label")
