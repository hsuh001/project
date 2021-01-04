########################################
# normalization
# date: 2020.12.28 - 12.29
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-2
########################################


# rm all objects ---------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load packages ---------------------------------------------------------------
library(scRNAseq)
sce_zeisel <- ZeiselBrainData()

library(scater)
# 将每个细胞的所有gene表达量加起来，得到每个细胞的文库大小
# 同时替换一些奇怪的gene名
sce_zeisel <- aggregateAcrossFeatures(
    sce_zeisel,
    id = sub("_loc[0-9]+$", "", rownames(sce_zeisel))
)

# gene annotation
library(org.Mm.eg.db)
rowData(sce_zeisel)$Ensembl <- mapIds(
    org.Mm.eg.db,
    keys = rownames(sce_zeisel),
    keytype = "SYMBOL",
    column = "ENSEMBL"
)

# qc, 先perCellQCMetrics()，后quickPerCellQC()
stats <- perCellQCMetrics(
    sce_zeisel,
    subsets = list(Mt = rowData(sce_zeisel)$featureType == "mito")
)

qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent")
)

sce_zeisel_filtered <- sce_zeisel[, !qc$discard]

sce_zeisel
# class: SingleCellExperiment
# dim: 19839 3005
# metadata(0):
# assays(1): counts
# rownames(19839): 0610005C13Rik 0610007N19Rik ... Zzef1 Zzz3
# rowData names(2): featureType Ensembl
# colnames(3005): 1772071015_C02 1772071017_G12 ... 1772066098_A12
#   1772058148_F03
# colData names(10): tissue group # ... level1class level2class
# reducedDimNames(0):
# spikeNames(0):
# altExpNames(2): ERCC repeat

sce_zeisel_filtered
# class: SingleCellExperiment
# dim: 19839 2816
# metadata(0):
# assays(1): counts
# rownames(19839): 0610005C13Rik 0610007N19Rik ... Zzef1 Zzz3
# rowData names(2): featureType Ensembl
# colnames(2816): 1772071015_C02 1772071017_G12 ... 1772063068_D01
#   1772066098_A12
# colData names(10): tissue group # ... level1class level2class
# reducedDimNames(0):
# spikeNames(0):
# altExpNames(2): ERCC repeat


# library size factor ---------------------------------------------------------
lib_sf_zeisel <- librarySizeFactors(sce_zeisel_filtered)
length(lib_sf_zeisel)
# [1] 2816

ncol(sce_zeisel_filtered)
# [1] 2816

summary(lib_sf_zeisel)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.1757  0.5680  0.8680  1.0000  1.2783  4.0839

# 直方图展示library size factor
hist(log10(lib_sf_zeisel), xlab = "Log10(size factor)", col = "grey80")


# normalization by deconvolution ----------------------------------------------
library(scran)
set.seed(100)
clust_zeisel <- quickCluster(sce_zeisel_filtered)
table(clust_zeisel)
# clust_zeisel
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14
# 170 254 441 178 393 148 219 240 189 123 112 103 135 111

# 去卷积
deconv_sf_zeisel <- calculateSumFactors(
    sce_zeisel_filtered,
    cluster = clust_zeisel
)
summary(deconv_sf_zeisel)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.1186  0.4860  0.8314  1.0000  1.3209  4.5090

# 数据包含的细胞类型
table(sce_zeisel_filtered$level1class)
# astrocytes_ependymal    endothelial-mural         interneurons
#                  179                  160                  290
#            microglia     oligodendrocytes        pyramidal CA1
#                   78                  774                  938
#         pyramidal SS 
#                  397

# 比较去卷积与常规方法得到的size factor
plot(lib_sf_zeisel, deconv_sf_zeisel,
    xlab = "Library size factor",
    ylab = "Deconvolution size factor",
    log = "xy",
    pch = 16,
    col = as.integer(factor(sce_zeisel_filtered$level1class))
) +
abline(a = 0, b = 1, col = "red")   # 截距为0，斜率为1


# calculate spike-in size factor ----------------------------------------------
library(scRNAseq)
sce_richard <- RichardTCellData()
sce_richard_filtered <- sce_richard[, sce_richard$`single cell quality` == "OK"]
sce_richard_filtered
# class: SingleCellExperiment
# dim: 46603 528
# metadata(0):
# assays(1): counts
# rownames(46603): ENSMUSG00000102693 ENSMUSG00000064842 ...
#   ENSMUSG00000096730 ENSMUSG00000095742
# rowData names(0):
# colnames(528): SLX-12611.N701_S502. SLX-12611.N702_S502. ...
#   SLX-12612.i712_i522. SLX-12612.i714_i522.
# colData names(13): age individual ... stimulus time
# reducedDimNames(0):
# spikeNames(0):
# altExpNames(1): ERCC

# 计算所有细胞的spike-in size factor
sce_richard_filtered <- computeSpikeFactors(sce_richard_filtered, "ERCC")
summary(sizeFactors(sce_richard_filtered))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.1247  0.4282  0.6274  1.0000  1.0699 23.3161


# plot the results of spike-in and deconvolution ------------------------------
to_plot <- data.frame(
    deconv_sf = calculateSumFactors(sce_richard_filtered),
    spike_sf = sizeFactors(sce_richard_filtered),
    stimulus = sce_richard_filtered$stimulus,
    time = sce_richard_filtered$time
)

ggplot(
    to_plot,
    aes(x = deconv_sf, y = spike_sf, color = time)
) +
geom_point() +
facet_wrap(~stimulus) +
scale_x_log10() +
scale_y_log10() +
geom_abline(intercept = 0, slope = 1, color = "red")


# log-transformation ----------------------------------------------------------
set.seed(100)
# normalization by deconvolution
clust_zeisel <- quickCluster(sce_zeisel_filtered)
sce_zeisel_filtered <- computeSumFactors(
    sce_zeisel_filtered,
    cluster = clust_zeisel,
    min.mean = 0.1
)
# logNormCounts()
sce_zeisel_filtered <- logNormCounts(sce_zeisel_filtered)
assayNames(sce_zeisel_filtered)
# [1] "counts"    "logcounts"
