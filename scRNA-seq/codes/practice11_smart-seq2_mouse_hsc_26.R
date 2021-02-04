########################################
# practice 11, Smart-seq2, mouse haematopoietic stem cell
# date: 2021.02.04 - 02.04
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-11
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
library(scRNAseq)
sce_nest_hsc <- NestorowaHSCData()
sce_nest_hsc
# class: SingleCellExperiment
# dim: 46078 1920
# metadata(0):
# assays(1): counts
# rownames(46078): ENSMUSG00000000001 ENSMUSG00000000003 ... ENSMUSG00000107391
#   ENSMUSG00000107392
# rowData names(0):
# colnames(1920): HSPC_007 HSPC_013 ... Prog_852 Prog_810
# colData names(2): cell.type FACS
# reducedDimNames(1): diffusion
# altExpNames(1): ERCC


# gene annotation -------------------------------------------------------------
library(AnnotationHub)
ens_mm_v97 <- AnnotationHub(localHub = TRUE)[["AH73905"]]
anno <- select(
    ens_mm_v97,
    keys = rownames(sce_nest_hsc),
    keytype = "GENEID",
    columns = c("SYMBOL", "SEQNAME")
)


# 全部对应
sum(is.na(anno$SYMBOL))
# [1] 0

sum(is.na(anno$SEQNAME))
# [1] 0

# 接下来只需要匹配顺序即可
rowData(sce_nest_hsc) <- anno[match(rownames(sce_nest_hsc), anno$GENEID),]
sce_nest_hsc
# class: SingleCellExperiment
# dim: 46078 1920
# metadata(0):
# assays(1): counts
# rownames(46078): ENSMUSG00000000001 ENSMUSG00000000003 ... ENSMUSG00000107391
#   ENSMUSG00000107392
# rowData names(3): GENEID SYMBOL SEQNAME
# colnames(1920): HSPC_007 HSPC_013 ... Prog_852 Prog_810
# colData names(2): cell.type FACS
# reducedDimNames(1): diffusion
# altExpNames(1): ERCC


# qc --------------------------------------------------------------------------
# 其实有线粒体gene，也有其他类型的如CHR_MG4151_PATCH
grep("MT", rowData(sce_nest_hsc)$SEQNAME)
#  [1] 17989 17990 17991 17992 17993 17994 17995 17996 17997 17998 17999 18000 18001 18002
# [15] 18003 18004 18005 18006 18007 18008 18009 18010 18011 18012 18013 18014 18015 18016
# [29] 18017 18018 18019 18020 18021 18022 18023 18024 18974

library(scater)
stats <- perCellQCMetrics(sce_nest_hsc)
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent")
)

colSums(as.matrix(qc), na.rm = TRUE)
#             low_lib_size            low_n_features high_altexps_ERCC_percent
#                      146                        28                       241
#                  discard
#                      264

sce_nest_hsc_filtered <- sce_nest_hsc[, !qc$discard]
dim(sce_nest_hsc_filtered)
# [1] 46078  1656

##### 使用qc标准对原数据作图
colData(sce_nest_hsc) <- cbind(colData(sce_nest_hsc), stats)
sce_nest_hsc$discard <- qc$discard

# 使用qc标准对原数据作图
gridExtra::grid.arrange(
    plotColData(sce_nest_hsc, y = "sum", colour_by = "discard") +
    scale_y_log10() + ggtitle("Total count"),
    plotColData(sce_nest_hsc, y = "detected", colour_by = "discard") +
    scale_y_log10() + ggtitle("Detected features"),
    plotColData(sce_nest_hsc, y = "altexps_ERCC_percent",
        colour_by = "discard") + ggtitle("ERCC percent"),
    ncol = 3
)


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(101000110)
cluster_nest_hsc <- quickCluster(sce_nest_hsc_filtered)
sce_nest_hsc_filtered <- computeSumFactors(
    sce_nest_hsc_filtered,
    cluster = cluster_nest_hsc
)
# logNormCounts()
sce_nest_hsc_filtered <- logNormCounts(sce_nest_hsc_filtered)
summary(sizeFactors(sce_nest_hsc_filtered))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.04368  0.42180  0.74844  1.00000  1.24926 15.92737


# measure the degree of change ------------------------------------------------
# and HVGs selection by proportion
dec_nest_hsc_spike <- modelGeneVarWithSpikes(
    sce_nest_hsc_filtered,
    spikes = "ERCC"
)
top_hvgs_nest_hsc <- getTopHVGs(dec_nest_hsc_spike, prop = 0.1)
length(top_hvgs_nest_hsc)
# [1] 384

# 查看方差大小
plot(
    dec_nest_hsc_spike$mean,
    dec_nest_hsc_spike$total,
    main = "Smart-seq2_mouse_hsc", pch = 16, cex = 0.5,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
cur_fit <- metadata(dec_nest_hsc_spike)
points(cur_fit$mean, cur_fit$var, col = "red", pch = 16)
curve(cur_fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

names(cur_fit)
# [1] "mean"    "var"     "trend"   "std.dev"

# 一共92个ERCC spike-in
length(unique(names(cur_fit$mean)))
# [1] 92


# dimension reduce ------------------------------------------------------------
set.seed(101010011)
sce_nest_hsc_filtered <- denoisePCA(
    sce_nest_hsc_filtered,
    subset.row = top_hvgs_nest_hsc,
    technical = dec_nest_hsc_spike
)
# 检查PCs的数量
ncol(reducedDim(sce_nest_hsc_filtered, "PCA"))
# [1] 9

sce_nest_hsc_filtered <- runTSNE(sce_nest_hsc_filtered, dimred = "PCA")


# clustering, graph-based -----------------------------------------------------
snn_gr_nest_hsc <- buildSNNGraph(sce_nest_hsc_filtered, use.dimred = "PCA")
# 鉴定cluster
cluster_nest_hsc <- igraph::cluster_walktrap(snn_gr_nest_hsc)$membership
colLabels(sce_nest_hsc_filtered) <- factor(cluster_nest_hsc)

table(cluster_nest_hsc)
# cluster_nest_hsc
#   1   2   3   4   5   6   7   8   9
# 203 472 258 175 142 229  20  83  74

# 绘制t-SNE图查看分群
plotTSNE(sce_nest_hsc_filtered, colour_by = "label")


# detecting markers -----------------------------------------------------------
markers_nest_hsc <- findMarkers(
    sce_nest_hsc_filtered,
    groups = colLabels(sce_nest_hsc_filtered),
    test.type = "wilcox",
    direction = "up",
    lfc = 0.5,
    row.data = rowData(sce_nest_hsc_filtered)[, "SYMBOL", drop = FALSE]
)

# use cluster 8 as an explaination
chosen_cluster <- "8"
markers_cluster_8 <- markers_nest_hsc[[chosen_cluster]]
# cluster 8的top 10
interest_markers <- markers_cluster_8[markers_cluster_8$Top <= 10, ]
length(interest_markers)
# [1] 13

# 提取cluster 8与其他clusters对比的AUC结果
aucs <- getMarkerEffects(interest_markers, prefix = "AUC")
rownames(aucs) <- interest_markers$SYMBOL

library(pheatmap)
pheatmap(aucs, color = viridis::plasma(100))


# annotating cell type --------------------------------------------------------
library(SingleR)
mm_ref <- MouseRNAseqData()
mm_ref
# class: SummarizedExperiment
# dim: 21214 358
# metadata(0):
# assays(1): logcounts
# rownames(21214): Xkr4 Rp1 ... LOC100039574 LOC100039753
# rowData names(0):
# colnames(358): ERR525589Aligned ERR525592Aligned ... SRR1044043Aligned
#   SRR1044044Aligned
# colData names(3): label.main label.fine label.ont

# 进行转换
renamed <- sce_nest_hsc_filtered
# mm_ref使用的事symbol name，需要进行转换
rownames(renamed) <- uniquifyFeatureNames(
    rownames(renamed),
    rowData(sce_nest_hsc_filtered)$SYMBOL
)
# 在参考数据集中找cell对应的细胞类型
predict <- SingleR(
    test = renamed,
    ref = mm_ref,
    labels = mm_ref$label.fine
)

table(predict$labels)
#          B cells Endothelial cells      Erythrocytes      Granulocytes       Macrophages
#               61                 1              1005                 1                 2
#        Monocytes          NK cells           T cells
#              500                 1                85

tab <- table(
    Pred = predict$labels,
    Cluster = sce_nest_hsc_filtered$label
)

pheatmap::pheatmap(
    log10(tab + 10),
    color = viridis::viridis(100)
)
