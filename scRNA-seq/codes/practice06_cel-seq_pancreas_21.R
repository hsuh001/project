########################################
# practice 6, CEL-seq, human pancreatic cells
# date: 2021.02.01 - 02.01
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-6
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
library(scRNAseq)
sce_muraro <- MuraroPancreasData()

sce_muraro
# class: SingleCellExperiment
# dim: 19059 3072
# metadata(0):
# assays(1): counts
# rownames(19059): A1BG-AS1__chr19 A1BG__chr19 ... ZZEF1__chr17 ZZZ3__chr1
# rowData names(2): symbol chr
# colnames(3072): D28-1_1 D28-1_2 ... D30-8_95 D30-8_96
# colData names(3): label donor plate
# reducedDimNames(0):
# altExpNames(1): ERCC

# 有4个donor
table(sce_muraro$donor)
# D28 D29 D30 D31
# 768 768 768 768


# gene annotation -------------------------------------------------------------
# 将没有匹配的NA去掉，并且去掉重复的行
head(rownames(sce_muraro))
# [1] "A1BG-AS1__chr19" "A1BG__chr19"     "A1CF__chr10"
# [4] "A2M-AS1__chr12"  "A2ML1__chr12"    "A2M__chr12"

library(AnnotationHub)
edb <- AnnotationHub(localHub = TRUE)[["AH73881"]]
gene_symbols <- sub("__chr.*", "", rownames(sce_muraro))
gene_ids <- mapIds(
    edb,
    keys = gene_symbols,
    keytype = "SYMBOL",
    column = "GENEID"
)

keep_row <- !is.na(gene_ids) & !duplicated(gene_ids)
sum(is.na(gene_ids))
# [1] 2110

sum(duplicated(gene_ids))
# [1] 2118

sce_muraro_filtered <- sce_muraro[keep_row, ]
rownames(sce_muraro_filtered) <- gene_ids[keep_row]
dim(sce_muraro_filtered)
# [1] 16940  3072


# qc --------------------------------------------------------------------------
# 备份数据用于质控探索
unfiltered <- sce_muraro_filtered

# 不含有线粒体gene
table(rowData(sce_muraro_filtered)$chr)
#  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22
#  1710   659   981   894   333   540   530   716   996   241  1186  1116   460   179   392
#  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX  chrY
#   988   656   788   922   815   564   637   613    24

stats <- perCellQCMetrics(sce_muraro_filtered)
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent"),
    batch = sce_muraro_filtered$donor
)

colSums(as.matrix(qc), na.rm = TRUE)
#              low_lib_size            low_n_features high_altexps_ERCC_percent
#                       250                       285                       307
#                   discard
#                       338

# 使用qc标准对原数据作图
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
    plotColData(unfiltered, x = "donor", y = "sum", colour_by = "discard") +
    scale_y_log10() + ggtitle("Total count"),
    plotColData(unfiltered, x = "donor", y = "detected",
        colour_by = "discard") + scale_y_log10() +
    ggtitle("Detected features"),
    plotColData(unfiltered, x = "donor", y = "altexps_ERCC_percent",
        colour_by = "discard") + ggtitle("ERCC percent"),
    ncol = 3
)

# 重新计算过滤条件
# 注意，这里用的subset，而非subsets
qc2 <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent"),
    batch = sce_muraro_filtered$donor,
    subset = sce_muraro_filtered$donor != "D28"
)
colSums(as.matrix(qc2), na.rm = TRUE)
#             low_lib_size            low_n_features high_altexps_ERCC_percent
#                      663                       700                       738
#                  discard
#                      773

# 重新作图
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc2$discard

gridExtra::grid.arrange(
    plotColData(unfiltered, x = "donor", y = "sum", colour_by = "discard") +
    scale_y_log10() + ggtitle("Total count"),
    plotColData(unfiltered, x = "donor", y = "detected",
        colour_by = "discard") + scale_y_log10() +
    ggtitle("Detected features"),
    plotColData(unfiltered, x = "donor", y = "altexps_ERCC_percent",
        colour_by = "discard") + ggtitle("ERCC percent"),
    ncol = 3
)

# 进行过滤
sce_muraro_final <- sce_muraro_filtered[, !qc2$discard]
dim(sce_muraro)
# [1] 19059  3072

dim(sce_muraro_filtered)
# [1] 16940  3072

dim(sce_muraro_final)
# [1] 16940  2299


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(1000)
cluster_muraro <- quickCluster(sce_muraro_final)
sce_muraro_final <- computeSumFactors(
    sce_muraro_final,
    cluster = cluster_muraro
)
# logNormCounts()
sce_muraro_final <- logNormCounts(sce_muraro_final)
summary(sizeFactors(sce_muraro_final))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.08782  0.54109  0.82081  1.00000  1.21079 13.98692


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
library(scran)
# 结合plate、donor，一起作为batch
table(sce_muraro_final$plate)
#   1   2   3   4   5   6   7   8
# 281 292 292 295 282 285 283 289

table(sce_muraro_final$donor)
# D28 D29 D30 D31
# 333 601 676 689

muraro_block <- paste0(sce_muraro_final$plate, "_", sce_muraro_final$donor)

# with spike-in
dec_muraro_spike <- modelGeneVarWithSpikes(
    sce_muraro_final,
    spikes = "ERCC",
    block = muraro_block
)
top_hvgs_muraro <- getTopHVGs(dec_muraro_spike, prop = 0.1)
length(top_hvgs_muraro)
# [1] 1557


# correcting batch effect -----------------------------------------------------
# 使用MNN进行批次效应矫正，同时还进行PCA降维
library(batchelor)
set.seed(1001010)
sce_muraro_megered <- fastMNN(
    sce_muraro_final,
    subset.row = top_hvgs_muraro,
    batch = sce_muraro_final$donor
)
sce_muraro_megered
# class: SingleCellExperiment
# dim: 1557 2299
# metadata(2): merge.info pca.info
# assays(1): reconstructed
# rownames(1557): ENSG00000089199 ENSG00000115263 ... ENSG00000163346
#   ENSG00000134278
# rowData names(1): rotation
# colnames(2299): D28-1_1 D28-1_2 ... D30-8_93 D30-8_94
# colData names(1): batch
# reducedDimNames(1): corrected
# altExpNames(0):

# 使用lost.var检查结果
metadata(sce_muraro_megered)$merge.info$lost.var
#              D28         D29         D30       D31
# [1,] 0.060846866 0.024121138 0.000000000 0.0000000
# [2,] 0.002645897 0.003017710 0.062420804 0.0000000
# [3,] 0.003448618 0.002640773 0.002598054 0.0816244


# dimension reduce ------------------------------------------------------------
set.seed(100111)
sce_muraro_megered <- runTSNE(sce_muraro_megered, dimred = "corrected")


# clustering, graph-based -----------------------------------------------------
snn_gr_muraro <- buildSNNGraph(sce_muraro_megered, use.dimred = "corrected")
# 鉴定cluster
cluster_muraro <- igraph::cluster_walktrap(snn_gr_muraro)$membership
colLabels(sce_muraro_megered) <- factor(cluster_muraro)

table(cluster_muraro)
# cluster_muraro
#   1   2   3   4   5   6   7   8   9  10
# 279 254 194 423 839 108  50  19 114  19

# 查看分群和细胞类型之间的关系
tab_merged <- table(
    Cluster = cluster_muraro,
    CellType = sce_muraro_final$label
)
tab_merged
#        CellType
# Cluster acinar alpha beta delta duct endothelial epsilon mesenchymal  pp unclear
#      1     216     1    3     0    7           0       0           0   2       0
#      2       0     5    4     0  212           0       0           1   0       4
#      3       0     0    0   184    0           0       0           0   0       0
#      4       0    17  383     3    0           0       0           0   0       0
#      5       0   770    2     1    0           0       3           0   0       0
#      6       1     0    2     0    6           1       0          79   0       0
#      7       0     1   41     1    3           0       0           0   0       0
#      8       0     0    1     0   11           0       0           0   0       0
#      9       0     1    6     0    0           0       0           0  94       0
#      10      0     0    0     0    0          17       0           0   0       0

library(pheatmap)
pheatmap(
    log10(tab_merged + 10),
    color = viridis::viridis(100)
)

# 绘制t-SNE图，查看批次效应的矫正效果
gridExtra::grid.arrange(
    plotTSNE(sce_muraro_megered, colour_by = "label"),
    plotTSNE(sce_muraro_megered, colour_by = "batch"),
    ncol = 2
)
