########################################
# practice 10, CEL-seq, mouse haematopoietic stem cell
# date: 2021.02.03 - 02.04
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-10
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
library(scRNAseq)
sce_grun_hsc <- GrunHSCData(ensembl = TRUE)
sce_grun_hsc
# class: SingleCellExperiment
# dim: 21817 1915
# metadata(0):
# assays(1): counts
# rownames(21817): ENSMUSG00000109644 ENSMUSG00000007777 ... ENSMUSG00000055670
#   ENSMUSG00000039068
# rowData names(3): symbol chr originalName
# colnames(1915): JC4_349_HSC_FE_S13_ JC4_350_HSC_FE_S13_ ... JC48P6_1203_HSC_FE_S8_
#   JC48P6_1204_HSC_FE_S8_
# colData names(2): sample protocol
# reducedDimNames(0):
# altExpNames(0):


# gene annotation -------------------------------------------------------------
library(AnnotationHub)
ens_mm_v97 <- AnnotationHub(localHub = TRUE)[["AH73905"]]
anno <- select(
    ens_mm_v97,
    keys = rownames(sce_grun_hsc),
    keytype = "GENEID",
    columns = c("SYMBOL", "SEQNAME")
)

# 全部对应
sum(is.na(anno$SYMBOL))
# [1] 0

sum(is.na(anno$SEQNAME))
# [1] 0

# 接下来只需要匹配顺序即可
rowData(sce_grun_hsc) <- anno[match(rownames(sce_grun_hsc), anno$GENEID),]


# qc --------------------------------------------------------------------------
# 没有线粒体
grep("MT", rowData(sce_grun_hsc)$SEQNAME)
# integer(0)

# protocol表示细胞提取方法
table(sce_grun_hsc$protocol)
#          micro-dissected cells sorted hematopoietic stem cells
#                           1546                             369

# 查看sample与protocol的关系
table(
    sce_grun_hsc$protocol,
    sce_grun_hsc$sample
)
#                                   JC20 JC21 JC26 JC27 JC28 JC30 JC32 JC35 JC36 JC37 JC39 JC4
#   micro-dissected cells             87   96   85   91   80   96   93   96   80   87   93   0
#   sorted hematopoietic stem cells    0    0    0    0    0    0    0    0    0    0    0  84
#                                  
#                                   JC40 JC41 JC43 JC44 JC45 JC46 JC48P4 JC48P6 JC48P7
#   micro-dissected cells             96   94   92   94   90   96      0      0      0
#   sorted hematopoietic stem cells    0    0    0    0    0    0     95     96     94

# 不将micro-dissected纳入考虑条件
library(scater)
stats <- perCellQCMetrics(sce_grun_hsc)
qc <- quickPerCellQC(
    stats,
    batch = sce_grun_hsc$protocol,
    subset = grepl("sorted", sce_grun_hsc$protocol)
)

colSums(as.matrix(qc), na.rm = TRUE)
#   low_lib_size low_n_features        discard
#            465            482            488

sce_grun_hsc_filtered <- sce_grun_hsc[, !qc$discard]
dim(sce_grun_hsc_filtered)
# [1] 21817  1427

##### 使用qc标准对原数据作图
colData(sce_grun_hsc) <- cbind(colData(sce_grun_hsc), stats)
sce_grun_hsc$discard <- qc$discard

# 使用qc标准对原数据作图
gridExtra::grid.arrange(
    plotColData(sce_grun_hsc, x = "sample", y = "sum",
        colour_by = "discard", other_fields = "protocol") +
    scale_y_log10() + ggtitle("Total count") + facet_wrap(~ protocol),
    plotColData(sce_grun_hsc, x = "sample", y = "detected",
    colour_by = "discard", other_fields = "protocol") +
    scale_y_log10() + ggtitle("Detected features") + facet_wrap(~ protocol),
    ncol = 1
)


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(101000110)
cluster_grun_hsc <- quickCluster(sce_grun_hsc_filtered)
sce_grun_hsc_filtered <- computeSumFactors(
    sce_grun_hsc_filtered,
    cluster = cluster_grun_hsc
)
# logNormCounts()
sce_grun_hsc_filtered <- logNormCounts(sce_grun_hsc_filtered)
summary(sizeFactors(sce_grun_hsc_filtered))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.02746  0.28972  0.60315  1.00000  1.20111 16.43327


# measure the degree of change ------------------------------------------------
# and HVGs selection by proportion
# 没有指定批次信息
dec_grun_hsc <- modelGeneVarByPoisson(sce_grun_hsc_filtered)
top_hvgs_grun_hsc <- getTopHVGs(dec_grun_hsc, prop = 0.1)
length(top_hvgs_grun_hsc)
# [1] 881

# 绘图查看方差
plot(
    dec_grun_hsc$mean,
    dec_grun_hsc$total,
    pch = 16,
    cex = 0.5,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
cur_fit <- metadata(dec_grun_hsc)
curve(cur_fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


# dimension reduce ------------------------------------------------------------
set.seed(101010011)
sce_grun_hsc_filtered <- denoisePCA(
    sce_grun_hsc_filtered,
    subset.row = top_hvgs_grun_hsc,
    technical = dec_grun_hsc
)
# 检查PCs的数量
ncol(reducedDim(sce_grun_hsc_filtered, "PCA"))
# [1] 9

sce_grun_hsc_filtered <- runTSNE(sce_grun_hsc_filtered, dimred = "PCA")


# clustering, graph-based -----------------------------------------------------
snn_gr_grun_hsc <- buildSNNGraph(sce_grun_hsc_filtered, use.dimred = "PCA")
# 鉴定cluster
cluster_grun_hsc <- igraph::cluster_walktrap(snn_gr_grun_hsc)$membership
colLabels(sce_grun_hsc_filtered) <- factor(cluster_grun_hsc)

table(cluster_grun_hsc)
# cluster_grun_hsc
#   1   2   3   4   5   6   7   8   9  10  11  12
# 257 148 219 151 102  48 105 148 103  63  65  18

# 绘制t-SNE图查看两种技术的差异
techs <- ifelse(
    grepl("micro", sce_grun_hsc_filtered$protocol),
    "micro",
    "sorted"
)
gridExtra::grid.arrange(
    plotTSNE(sce_grun_hsc_filtered, colour_by = "label"),
    plotTSNE(sce_grun_hsc_filtered, colour_by = I(techs)),
    ncol = 2
)


# detecting markers -----------------------------------------------------------
markers_grun_hsc <- findMarkers(
    sce_grun_hsc_filtered,
    groups = colLabels(sce_grun_hsc_filtered),
    test.type = "wilcox",
    direction = "up",
    row.data = rowData(sce_grun_hsc_filtered)[, "SYMBOL", drop = FALSE]
)

# use cluster 2 as an explaination
chosen_cluster <- "2"
markers_cluster_2 <- markers_grun_hsc[[chosen_cluster]]
# cluster 1的top 10
interest_markers <- markers_cluster_2[markers_cluster_2$Top <= 10, ]
length(interest_markers)
# [1] 16

# 提取cluster 2与其他clusters对比的AUC结果
aucs <- getMarkerEffects(interest_markers, prefix = "AUC")
rownames(aucs) <- interest_markers$SYMBOL

library(pheatmap)
pheatmap(aucs, color = viridis::plasma(100))
