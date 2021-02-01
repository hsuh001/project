########################################
# practice 5, CEL-seq2, human pancreatic cells
# date: 2021.01.31 - 02.01
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-5
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
library(scRNAseq)
sce_grun <- GrunPancreasData()

sce_grun
# class: SingleCellExperiment
# dim: 20064 1728
# metadata(0):
# assays(1): counts
# rownames(20064): A1BG-AS1__chr19 A1BG__chr19 ... ZZEF1__chr17 ZZZ3__chr1
# rowData names(2): symbol chr
# colnames(1728): D2ex_1 D2ex_2 ... D17TGFB_95 D17TGFB_96
# colData names(2): donor sample
# reducedDimNames(0):
# altExpNames(1): ERCC

rowData(sce_grun)
# DataFrame with 20064 rows and 2 columns
#                      symbol         chr
#                 <character> <character>
# A1BG-AS1__chr19    A1BG-AS1       chr19
# A1BG__chr19            A1BG       chr19
# A1CF__chr10            A1CF       chr10
# A2M-AS1__chr12      A2M-AS1       chr12
# A2ML1__chr12          A2ML1       chr12
# ...                     ...         ...
# ZYG11A__chr1         ZYG11A        chr1
# ZYG11B__chr1         ZYG11B        chr1
# ZYX__chr7               ZYX        chr7
# ZZEF1__chr17          ZZEF1       chr17
# ZZZ3__chr1             ZZZ3        chr1

# batch information
table(sce_grun$donor)
# D10 D17  D2  D3  D7
# 288 480  96 480 384


# gene annotation -------------------------------------------------------------
# symbol to Ensembl ID
library(org.Hs.eg.db)
gene_ids <- mapIds(
    org.Hs.eg.db,
    keys = rowData(sce_grun)$symbol,
    keytype = "SYMBOL",
    column = "ENSEMBL"
)

# method_1，将Ensembl ID添加到原来的数据中
library(scater)
rowData(sce_grun)$ensembl <- uniquifyFeatureNames(
    rowData(sce_grun)$symbol,
    gene_ids
)
rowData(sce_grun)
# DataFrame with 20064 rows and 3 columns
#                      symbol         chr         ensembl
#                 <character> <character>     <character>
# A1BG-AS1__chr19    A1BG-AS1       chr19 ENSG00000268895
# A1BG__chr19            A1BG       chr19 ENSG00000121410
# A1CF__chr10            A1CF       chr10 ENSG00000148584
# A2M-AS1__chr12      A2M-AS1       chr12 ENSG00000245105
# A2ML1__chr12          A2ML1       chr12 ENSG00000166535
# ...                     ...         ...             ...
# ZYG11A__chr1         ZYG11A        chr1 ENSG00000203995
# ZYG11B__chr1         ZYG11B        chr1 ENSG00000162378
# ZYX__chr7               ZYX        chr7 ENSG00000159840
# ZZEF1__chr17          ZZEF1       chr17 ENSG00000074755
# ZZZ3__chr1             ZZZ3        chr1 ENSG00000036549

# method_2， 将没有匹配的NA去掉，并且去掉重复的行
keep_row <- !is.na(gene_ids) & !duplicated(gene_ids)
sum(is.na(gene_ids))
# [1] 2570

sum(duplicated(gene_ids))
# [1] 2589

sce_grun_filtered <- sce_grun[keep_row, ]
rownames(sce_grun_filtered) <- gene_ids[keep_row]
dim(sce_grun_filtered)
# [1] 17474  1728


# qc --------------------------------------------------------------------------
# 备份数据用于质控探索
unfiltered <- sce_grun_filtered

# 不含有线粒体gene
table(rowData(sce_grun_filtered)$chr)
# chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22
# 1755   669  1064   941   327   550   533   738  1049   257  1240  1139   479   188   399
# chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX  chrY
# 1007   689   791   897   823   603   660   656    20

stats <- perCellQCMetrics(sce_grun_filtered)
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent"),
    batch = sce_grun_filtered$donor
)

colSums(as.matrix(qc), na.rm = TRUE)
#              low_lib_size            low_n_features high_altexps_ERCC_percent
#                        85                        93                        88
#                   discard
#                       129

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
    batch = sce_grun_filtered$donor,
    subset = sce_grun_filtered$donor %in% c("D17", "D7", "D2")
)
colSums(as.matrix(qc2), na.rm = TRUE)
#             low_lib_size            low_n_features high_altexps_ERCC_percent
#                      450                       512                       606
#                  discard
#                      665

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
sce_grun_final <- sce_grun_filtered[, !qc2$discard]
dim(sce_grun)
# [1] 20064  1728

dim(sce_grun_filtered)
# [1] 17474  1728

dim(sce_grun_final)
# [1] 17474  1063


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(1000)
cluster_grun <- quickCluster(sce_grun_final)
sce_grun_final <- computeSumFactors(
    sce_grun_final,
    cluster = cluster_grun
)
# logNormCounts()
sce_grun_final <- logNormCounts(sce_grun_final)
summary(sizeFactors(sce_grun_final))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.09863 0.51067 0.79625 1.00000 1.23058 8.83814

# 比较librarySizeFactors()（常规方法）、
# computeSumFactors()（去卷积方法）进行归一化的差异
plot(
    librarySizeFactors(sce_grun_final),
    sizeFactors(sce_grun_final),
    pch = 16,
    xlab = "Library size factors",
    ylab = "Deconvolution factors",
    col = as.integer(factor(sce_grun_final$donor)),
    log = "xy"
)
abline(a = 0, b = 1, col = "red")


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
library(scran)
grun_block <- paste0(sce_grun_final$sample, "_", sce_grun_final$donor)

# with spike-in
dec_grun_spike <- modelGeneVarWithSpikes(
    sce_grun_final,
    spikes = "ERCC",
    block = grun_block
)
top_hvgs_grun <- getTopHVGs(dec_grun_spike, prop = 0.1)
length(top_hvgs_grun)
# [1] 1461


# correcting batch effect -----------------------------------------------------
# 使用MNN进行批次效应矫正，同时还进行PCA降维
library(batchelor)
set.seed(1001010)
sce_grun_megered <- fastMNN(
    sce_grun_final,
    subset.row = top_hvgs_grun,
    batch = sce_grun_final$donor
)
sce_grun_megered
# class: SingleCellExperiment
# dim: 1461 1063
# metadata(2): merge.info pca.info
# assays(1): reconstructed
# rownames(1461): REG1B__chr2 REG3A__chr2 ... ANXA3__chr4 XPR1__chr1
# rowData names(1): rotation
# colnames(1063): D2ex_1 D2ex_2 ... D17TGFB_94 D17TGFB_95
# colData names(1): batch
# reducedDimNames(1): corrected
# altExpNames(0):

# 使用lost.var检查结果
metadata(sce_grun_megered)$merge.info$lost.var
#              D10         D17          D2         D3         D7
# [1,] 0.030626488 0.032122881 0.000000000 0.00000000 0.00000000
# [2,] 0.007150824 0.011372000 0.036091011 0.00000000 0.00000000
# [3,] 0.003904948 0.005134747 0.007729219 0.05238854 0.00000000
# [4,] 0.011861752 0.014642939 0.013594152 0.01235050 0.05386529


# dimension reduce ------------------------------------------------------------
set.seed(100111)
sce_grun_megered <- runTSNE(sce_grun_megered, dimred = "corrected")


# clustering, graph-based -----------------------------------------------------
snn_gr_grun <- buildSNNGraph(sce_grun_megered, use.dimred = "corrected")
# 鉴定cluster
cluster_grun <- igraph::cluster_walktrap(snn_gr_grun)$membership
colLabels(sce_grun_megered) <- factor(cluster_grun)

table(cluster_grun)
# cluster_grun
#   1   2   3   4   5   6   7   8   9  10  11  12
# 241 120 187  17 185  22  53  17 119  28  18  56

tab_merged <- table(
    Cluster = cluster_grun,
    Donor = sce_grun_megered$batch
)
tab_merged
#        Donor
# Cluster D10 D17  D2  D3  D7
#      1   32  70  31  80  28
#      2   14  34   3   2  67
#      3   12  71  31   2  71
#      4    5   4   2   4   2
#      5   11 119   0   0  55
#      6    2   8   3   3   6
#      7    3  40   0   0  10
#      8    1   9   0   0   7
#      9   15  36  12  11  45
#      10   5  13   0   0  10
#      11   4  13   0   0   1
#      12   5  17   0   1  33

# 绘制t-SNE图，查看批次效应的矫正效果
gridExtra::grid.arrange(
    plotTSNE(sce_grun_megered, colour_by = "label"),
    plotTSNE(sce_grun_megered, colour_by = "batch"),
    ncol = 2
)
