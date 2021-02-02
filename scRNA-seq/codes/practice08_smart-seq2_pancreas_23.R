########################################
# practice 8, Smart-seq2, human pancreatic cells
# date: 2021.02.02 - 02.02
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-8
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
library(scRNAseq)
sce_seger <- SegerstolpePancreasData()

sce_seger
# class: SingleCellExperiment
# dim: 26179 3514
# metadata(0):
# assays(1): counts
# rownames(26179): SGIP1 AZIN2 ... BIVM-ERCC5 eGFP
# rowData names(2): symbol refseq
# colnames(3514): HP1502401_N13 HP1502401_D14 ... HP1526901T2D_O11 HP1526901T2D_A8
# colData names(8): Source Name individual ... age body mass index
# reducedDimNames(0):
# altExpNames(1): ERCC


# gene annotation -------------------------------------------------------------
library(AnnotationHub)
edb <- AnnotationHub(localHub = TRUE)[["AH73881"]]
gene_symbols <- rowData(sce_seger)$symbol
gene_ids <- mapIds(
    edb,
    keys = gene_symbols,
    keytype = "SYMBOL",
    column = "GENEID"
)

# 将NA替换成symbol，去除重复行
gene_ids <- ifelse(is.na(gene_ids), gene_symbols, gene_ids)
keep_row <- !duplicated(gene_ids)

sce_seger_filtered <- sce_seger[keep_row, ]
rownames(sce_seger_filtered) <- gene_ids[keep_row]
dim(sce_seger_filtered)
# [1] 25454  3514

# 编辑样本信息，只保留3列
keep_coldata <- c("cell type", "individual", "single cell well quality")
emtab_meta <- colData(sce_seger_filtered)[, keep_coldata]
colnames(emtab_meta) <- c("cell_type", "donor", "quality")
colData(sce_seger_filtered) <- emtab_meta
sce_seger_filtered
# class: SingleCellExperiment
# dim: 25454 3514
# metadata(0):
# assays(1): counts
# rownames(25454): ENSG00000118473 ENSG00000142920 ... ENSG00000278306 eGFP
# rowData names(2): symbol refseq
# colnames(3514): HP1502401_N13 HP1502401_D14 ... HP1526901T2D_O11 HP1526901T2D_A8
# colData names(3): cell_type donor quality
# reducedDimNames(0):
# altExpNames(1): ERCC

# 去掉cell_type中的cell字符，并将首字母大写
sce_seger_filtered$cell_type <- gsub(" cell", "", sce_seger_filtered$cell_type)
sce_seger_filtered$cell_type <- paste0(
    toupper(substr(sce_seger_filtered$cell_type, 1, 1)),
    substring(sce_seger_filtered$cell_type, 2)
)


# qc --------------------------------------------------------------------------
# 备份数据用于质控探索
unfiltered <- sce_seger_filtered

# 作者标注的细胞质量信息
table(sce_seger_filtered$quality)
# control, 2-cell well  control, empty well     low quality cell                   OK
#                   32                   96                 1177                 2209

library(scater)
stats <- perCellQCMetrics(sce_seger_filtered)
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent"),
    batch = sce_seger_filtered$donor
)

colSums(as.matrix(qc), na.rm = TRUE)
#             low_lib_size            low_n_features high_altexps_ERCC_percent
#                      594                       698                       934
#                  discard
#                     1047

##### 使用qc标准对原数据作图
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

# 使用qc标准对原数据作图
gridExtra::grid.arrange(
    plotColData(unfiltered, x = "donor", y = "sum",
        colour_by = "discard") +
    scale_y_log10() +
    ggtitle("Total count") +
    theme(axis.text.x = element_text(angle = 90)),
    plotColData(unfiltered, x = "donor", y = "detected",
    colour_by = "discard") +
    scale_y_log10() +
    ggtitle("Detected features") +
    theme(axis.text.x = element_text(angle = 90)),
    plotColData(unfiltered, x = "donor", y = "altexps_ERCC_percent",
        colour_by = "discard") +
    ggtitle("ERCC percent") +
    theme(axis.text.x = element_text(angle = 90)),
    ncol = 3
)

# 重新计算过滤条件
library(scater)
qc2 <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent"),
    batch = sce_seger_filtered$donor,
    subset = !sce_seger_filtered$donor %in% c("HP1504901", "HP1509101")
)

##### 重新作图
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc2$discard

# 使用qc标准对原数据作图
gridExtra::grid.arrange(
    plotColData(unfiltered, x = "donor", y = "sum",
        colour_by = "discard") +
    scale_y_log10() +
    ggtitle("Total count") +
    theme(axis.text.x = element_text(angle = 90)),
    plotColData(unfiltered, x = "donor", y = "detected",
    colour_by = "discard") +
    scale_y_log10() +
    ggtitle("Detected features") +
    theme(axis.text.x = element_text(angle = 90)),
    plotColData(unfiltered, x = "donor", y = "altexps_ERCC_percent",
        colour_by = "discard") +
    ggtitle("ERCC percent") +
    theme(axis.text.x = element_text(angle = 90)),
    ncol = 3
)

colSums(as.matrix(qc2))
#             low_lib_size            low_n_features high_altexps_ERCC_percent
#                      788                      1056                      1031
#                  discard
#                     1246

# 将qc + 作者标注的低质量，一起过滤
low_equal <- sce_seger_filtered$quality == "low quality cell"
length(low_equal)
# [1] 3514

sce_seger_final <- sce_seger_filtered[, !(qc2$discard | low_equal)]

dim(sce_seger)
# [1] 26179  3514

dim(sce_seger_filtered)
# [1] 25454  3514

dim(sce_seger_final)
# [1] 25454  2090


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
sce_seger_test_1 <- sce_seger_final
sce_seger_test_1 <- computeSpikeFactors(sce_seger_test_1, "ERCC")
sce_seger_test_1 <- logNormCounts(sce_seger_test_1)
# Error in .local(x, ...) : size factors should be positive

# method_1，去掉这几个细胞
table(colSums(counts(altExp(sce_seger_final))) == 0)
# FALSE  TRUE
#  2087     3

sce_seger_test_2 <- sce_seger_final[, !colSums(counts(altExp(sce_seger_final))) == 0]
sce_seger_test_2 <- computeSpikeFactors(sce_seger_test_2, "ERCC")
sce_seger_test_2 <- logNormCounts(sce_seger_test_2)

# method_2，使用computeSumFactors()
cluster_seger <- quickCluster(sce_seger_final)
sce_seger_final <- computeSumFactors(
    sce_seger_final,
    cluster = cluster_seger
)
# logNormCounts()
sce_seger_final <- logNormCounts(sce_seger_final)
summary(sizeFactors(sce_seger_final))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.01406  0.39040  0.70823  1.00000  1.33182 11.18150


# measure the degree of change ------------------------------------------------
# and HVGs selection by number

# 批次AZ相对其他batch的细胞数量少很多，去除
table(sce_seger_final$donor)
#          AZ    HP1502401 HP1504101T2D    HP1504901    HP1506401    HP1507101 HP1508501T2D
#          65          249          247          139          230          203          276
#   HP1509101 HP1525301T2D HP1526901T2D
#         109          274          298

# 同时还要去除没有ERCC的细胞
keep_col <- librarySizeFactors(altExp(sce_seger_final)) > 0 & sce_seger_final$donor != "AZ"
length(keep_col)
# [1] 2090

sce_seger_final_new <- sce_seger_final[, keep_col]

dec_seger_spike <- modelGeneVarWithSpikes(
    sce_seger_final_new,
    spikes = "ERCC",
    block = sce_seger_final_new$donor
)
top_hvgs_seger <- getTopHVGs(dec_seger_spike, n = 2000)
length(top_hvgs_seger)
# [1] 2000

# 查看方差大小
par(mfrow = c(3, 3))
blcok_stats <- dec_seger_spike$per.block
for (i in colnames(blcok_stats)) {
    current <- blcok_stats[[i]]
    plot(
        current$mean,
        current$total,
        pch = 16,
        cex = 0.5,
        main = i,
        xlab = "Mean of log-expression",
        ylab = "Variance of log-expression"
        )
        cur_fit <- metadata(current)
        points(cur_fit$mean, cur_fit$var, col = "red", pch = 16)
        # 蓝线指的是所有gene都会存在的一种偏差
        curve(cur_fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
}


# dimension reduce ------------------------------------------------------------
library(BiocSingular)
set.seed(101011001)
sce_seger_final_new <- runPCA(
    sce_seger_final_new,
    subset_row = top_hvgs_seger,
    ncomponents = 25
)
sce_seger_final_new <- runTSNE(sce_seger_final_new, dimred = "PCA")


# clustering, graph-based -----------------------------------------------------
snn_gr_seger <- buildSNNGraph(sce_seger_final_new, use.dimred = "PCA")
# 鉴定cluster
cluster_seger <- igraph::cluster_walktrap(snn_gr_seger)$membership
colLabels(sce_seger_final_new) <- factor(cluster_seger)

table(cluster_seger)
# cluster_seger
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22
#  99 137 150 135 275  68 177 147  72 100  47  42  34  83  26 185  15  65  16  76  44  32

# 查看分群和批次之间的关系
cluster_batch_tab <- table(
    Cluster = cluster_seger,
    Donor = sce_seger_final_new$donor
)

library(pheatmap)
pheatmap(
    log10(cluster_batch_tab + 10),
    color = viridis::viridis(100)
)

# 绘制t-SNE，查看批次效应
gridExtra::grid.arrange(
    plotTSNE(sce_seger_final_new, colour_by = "label"),
    plotTSNE(sce_seger_final_new, colour_by = "donor"),
    ncol = 2
)


# correcting batch effect -----------------------------------------------------
# 进行批次效应矫正
library(batchelor)
set.seed(1001010)
sce_seger_megered <- fastMNN(
    sce_seger_final_new,
    subset.row = top_hvgs_seger,
    batch = sce_seger_final_new$donor
)
sce_seger_megered
# class: SingleCellExperiment
# dim: 2000 2025
# metadata(2): merge.info pca.info
# assays(1): reconstructed
# rownames(2000): ENSG00000115263 ENSG00000118271 ... ENSG00000145495
#   ENSG00000070540
# rowData names(1): rotation
# colnames(2025): HP1502401_H13 HP1502401_J14 ... HP1526901T2D_N8 HP1526901T2D_A8
# colData names(1): batch
# reducedDimNames(1): corrected
# altExpNames(0):

# 使用lost.var检查结果
metadata(sce_seger_megered)$merge.info$lost.var
#        HP1502401 HP1504101T2D   HP1504901    HP1506401   HP1507101 HP1508501T2D   HP1509101
# [1,] 0.013849337  0.021962945 0.000000000 0.0000000000 0.000000000 0.0000000000 0.000000000
# [2,] 0.012936219  0.008721311 0.023408664 0.0000000000 0.000000000 0.0000000000 0.000000000
# [3,] 0.011344055  0.008272761 0.006292893 0.0336172395 0.000000000 0.0000000000 0.000000000
# [4,] 0.002286099  0.001960364 0.001705354 0.0054620488 0.029222844 0.0000000000 0.000000000
# [5,] 0.004847777  0.003995802 0.004238750 0.0055952986 0.005600396 0.0437603456 0.000000000
# [6,] 0.001961449  0.001615226 0.001384978 0.0016600023 0.002045468 0.0013627051 0.040699573
# [7,] 0.002352900  0.002165374 0.001516092 0.0008315320 0.001692425 0.0009056075 0.002149675
# [8,] 0.001033656  0.001349923 0.001785929 0.0007330182 0.002067223 0.0012553188 0.001506425
#      HP1525301T2D HP1526901T2D
# [1,]  0.000000000   0.00000000
# [2,]  0.000000000   0.00000000
# [3,]  0.000000000   0.00000000
# [4,]  0.000000000   0.00000000
# [5,]  0.000000000   0.00000000
# [6,]  0.000000000   0.00000000
# [7,]  0.039095755   0.00000000
# [8,]  0.001235738   0.04090019


# again, dimension reduce -----------------------------------------------------
library(BiocSingular)
set.seed(101011001)
sce_seger_megered <- runTSNE(sce_seger_megered, dimred = "corrected")


# again, clustering, graph-based ----------------------------------------------
snn_gr_seger <- buildSNNGraph(sce_seger_megered, use.dimred = "corrected")
# 鉴定cluster
cluster_seger <- igraph::cluster_walktrap(snn_gr_seger)$membership
colLabels(sce_seger_megered) <- factor(cluster_seger)

table(cluster_seger)
# cluster_seger
#   1   2   3   4   5   6   7   8   9
# 384 158 276 842 169  26  54 100  16

# 查看分群和批次之间的关系
cluster_batch_tab <- table(
    Cluster = cluster_seger,
    Donor = sce_seger_megered$batch
)

library(pheatmap)
pheatmap(
    log10(cluster_batch_tab + 10),
    color = viridis::viridis(100)
)

# 绘制t-SNE，查看批次效应
gridExtra::grid.arrange(
    plotTSNE(sce_seger_megered, colour_by = "label"),
    plotTSNE(sce_seger_megered, colour_by = "batch"),
    ncol = 2
)
