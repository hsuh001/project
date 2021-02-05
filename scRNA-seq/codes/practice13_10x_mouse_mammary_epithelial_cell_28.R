########################################
# practice 13, 10X Genomics, mouse mammary epithelial cells
# date: 2021.02.04 - 02.05
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-13
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
sce_mam
# class: SingleCellExperiment
# dim: 27998 2915
# metadata(0):
# assays(1): counts
# rownames(27998): ENSMUSG00000051951 ENSMUSG00000089699 ... ENSMUSG00000096730
#   ENSMUSG00000095742
# rowData names(2): Ensembl Symbol
# colnames: NULL
# colData names(3): Barcode Sample Condition
# reducedDimNames(0):
# altExpNames(0):

# 样本信息
sapply(
    names(colData(sce_mam)), function(x) head(colData(sce_mam)[, x])
)
#      Barcode              Sample Condition
# [1,] "AAACCTGAGGATGCGT-1" "G_1"  "Gestation"
# [2,] "AAACCTGGTAGTAGTA-1" "G_1"  "Gestation"
# [3,] "AAACCTGTCAGCATGT-1" "G_1"  "Gestation"
# [4,] "AAACCTGTCGTCCGTT-1" "G_1"  "Gestation"
# [5,] "AAACGGGCACGAAATA-1" "G_1"  "Gestation"
# [6,] "AAACGGGCAGACGCTC-1" "G_1"  "Gestation"


# gene annotation -------------------------------------------------------------
library(scater)
rownames(sce_mam) <- uniquifyFeatureNames(
    rowData(sce_mam)$Ensembl,
    rowData(sce_mam)$Symbol
)

library(AnnotationHub)
ens_mm_v97 <- AnnotationHub(localHub = TRUE)[["AH73905"]]
rowData(sce_mam)$SEQNAME <- mapIds(
    ens_mm_v97,
    keys = rowData(sce_mam)$Ensembl,
    keytype = "GENEID",
    column = "SEQNAME"
)

# 线粒体gene个数
sum(grepl("MT", rowData(sce_mam)$SEQNAME))
# [1] 13


# quality control -------------------------------------------------------------
is_mito <- rowData(sce_mam)$SEQNAME == "MT"
stats <- perCellQCMetrics(
    sce_mam,
    subsets = list(Mito = which(is_mito))
)
qc <- quickPerCellQC(stats, percent_subsets = "subsets_Mito_percent")

colSums(as.matrix(qc), na.rm = TRUE)
#             low_lib_size            low_n_features high_subsets_Mito_percent
#                        0                         0                       143
#                  discard
#                      143

sce_mam_filtered <- sce_mam[, !qc$discard]

##### 使用qc标准对原数据作图
colData(sce_mam) <- cbind(colData(sce_mam), stats)
sce_mam$discard <- qc$discard

# 使用qc标准对原数据作图
gridExtra::grid.arrange(
    plotColData(sce_mam, y = "sum", colour_by = "discard") +
    scale_y_log10() + ggtitle("Total count"),
    plotColData(sce_mam, y = "detected", colour_by = "discard") +
    scale_y_log10() + ggtitle("Detected features"),
    plotColData(sce_mam, y = "subsets_Mito_percent",
        colour_by = "discard") + ggtitle("Mito percent"),
    ncol = 3
)

# 查看线粒体含量与文库大小的关系
plotColData(sce_mam, x = "sum", y = "subsets_Mito_percent",
    colour_by = "discard") +
scale_y_log10() + ggtitle("Mito percent to total count")


# normalization by deconvolution ----------------------------------------------
library(scran)
set.seed(101000110)
cluster_mam <- quickCluster(sce_mam_filtered)
sce_mam_final <- computeSumFactors(sce_mam_filtered, clusters = cluster_mam)
sce_mam_final <- logNormCounts(sce_mam_final)
sce_mam_final
# class: SingleCellExperiment
# dim: 27998 2772
# metadata(0):
# assays(2): counts logcounts
# rownames(27998): Xkr4 Gm1992 ... Vmn2r122 CAAA01147332.1
# rowData names(3): Ensembl Symbol SEQNAME
# colnames: NULL
# colData names(4): Barcode Sample Condition sizeFactor
# reducedDimNames(0):
# altExpNames(0):

dim(sce_mam)
# [1] 27998  2915

dim(sce_mam_filtered)
# [1] 27998  2772

dim(sce_mam_final)
# [1] 27998  2772


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
set.seed(00010101)
dec_mam_pois <- modelGeneVarByPoisson(sce_mam_final)
top_hvgs_mam <- getTopHVGs(dec_mam_pois, prop = 0.1)
length(top_hvgs_mam)
# [1] 1366

# 查看方差大小
plot(
    dec_mam_pois$mean,
    dec_mam_pois$total,
    main = "10X_mouse_mammary epithelial cells",
    pch = 16, cex = 0.5,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
cur_fit <- metadata(dec_mam_pois)
curve(cur_fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


# dimension reduce ------------------------------------------------------------
library(BiocSingular)
set.seed(101010011)
##### PCA
sce_mam_final <- denoisePCA(
    sce_mam_final,
    technical = dec_mam_pois,
    subset.row = top_hvgs_mam
)

ncol(reducedDim(sce_mam_final, "PCA"))
# [1] 15

##### t-SNE
sce_mam_final <- runTSNE(sce_mam_final, dimred = "PCA")


# clustering ------------------------------------------------------------------
##### graph-based clustering
# 使用PCA的前几个PCs构建SNNG
library(scran)
snn_gr_mam <- buildSNNGraph(sce_mam_final, k = 25, use.dimred = "PCA")
cluster_mam <- igraph::cluster_walktrap(snn_gr_mam)$membership
colLabels(sce_mam_final) <- factor(cluster_mam)

table(cluster_mam)
# cluster_mam
#   1   2   3   4   5   6   7   8   9  10
# 550 799 716 452  24  84  52  39  32  24

plotTSNE(sce_mam_final, colour_by = "label")
