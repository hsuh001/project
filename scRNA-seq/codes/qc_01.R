########################################
# quality control for scRNA-seq
# date: 2020.12.26 - 12.28
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-1
########################################


# rm all objects ---------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load packages ---------------------------------------------------------------
# BiocManager::install("scRNAseq")
library(scRNAseq)


# download test data ----------------------------------------------------------
sce_416b <- LunSpikeInData(which = "416b")

table(sce_416b$block)
# 两个96孔板的细胞
# 20160113 20160325
#       96       96

sce_416b$block <- factor(sce_416b$block)
# class: SingleCellExperiment
# dim: 46604 192
# metadata(0):
# assays(1): counts
# rownames(46604): ENSMUSG00000102693 ENSMUSG00000064842 ...
#   ENSMUSG00000095742 CBFB-MYH11-mcherry
# rowData names(1): Length
# colnames(192): SLX-9555.N701_S502.C89V9ANXX.s_1.r_1
#   SLX-9555.N701_S503.C89V9ANXX.s_1.r_1 ...
#   SLX-11312.N712_S508.H5H5YBBXX.s_8.r_1
#   SLX-11312.N712_S517.H5H5YBBXX.s_8.r_1
# colData names(9): Source Name cell line ... spike-in addition block
# reducedDimNames(0):
# spikeNames(0):
# altExpNames(2): ERCC SIRV


# measurements of qc ----------------------------------------------------------
# library size
# the number of expressed features in each cell
# the proportion of reads mapped to spike-in transcripts
# the proportion of reads mapped to genes in the mitochondrial


# state the information of mitochondrial --------------------------------------
# method_1, AnnotationHub()
library(AnnotationHub)
ah <- AnnotationHub()
ens_mm_v97 <- ah[["AH73905"]]
location <- mapIds(ens_mm_v97, keys = rownames(sce.416b),
    keytype = "GENEID", column = "SEQNAME")
is_mito <- which(location == "MT")

# method_2,
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene,
    keys = rownames(sce_416b),
    keytype = "GENEID", column = "CDSCHROM")
is_mito <- which(location == "MT")
# named integer(0)
summary(location == "chrM")
#    Mode   FALSE    TRUE    NA's
# logical   22428      13   24163

# method_3, 调用rowRanges()获取基因组坐标，再结合seqnames()找到"MT"
location <- rowRanges(sce_416b)
is_mito <- any(seqnames(location) == "MT")
table(is_mito)
# is_mito
# FALSE
# 46604


# calculate metrics of quality control ----------------------------------------
# method_1, perCellQCMetrics()计算
library(scater)
df <- perCellQCMetrics(sce_416b, subsets = list(Mito = is_mito))
df
# DataFrame with 192 rows and 16 columns
#           sum  detected   percent_top_50  percent_top_100  percent_top_200
#     <integer> <integer>        <numeric>        <numeric>        <numeric>
# 1      865936      7618 26.7218362557972 32.2773276546997 39.7208338722492
# 2     1076277      7521 29.4043262097025 35.0354044544295  42.258080401235
# 3     1180138      8306 27.3453613052033 32.4769645583822 39.3295529844815
# ...       ...       ...              ...              ...              ...
# 191     46731      6649 32.2997581904945 37.9148744944469 44.5999443624147
# 192   1866692     10964 26.6632095707273 31.2583972074665 37.5607759608977
#      percent_top_500 subsets_Mito_sum subsets_Mito_detected
#            <numeric>        <integer>             <integer>
# 1   52.9037942757894                0                     0
# 2   55.7454075484285                0                     0
# 3   51.9336721637639                0                     0
# ...              ...              ...                   ...
# 191 56.5235068798014                0                     0
# 192 48.9489428357758                0                     0
#     subsets_Mito_percent altexps_ERCC_sum altexps_ERCC_detected
#                <numeric>        <integer>             <integer>
# 1                      0            65278                    39
# 2                      0            74748                    40
# 3                      0            60878                    42
# ...                  ...              ...                   ...
# 191                    0             7580                    44
# 192                    0            48664                    39
#     altexps_ERCC_percent altexps_SIRV_sum altexps_SIRV_detected
#                <numeric>        <integer>             <integer>
# 1       6.80658407035354            27828                     7
# 2       6.28029958040595            39173                     7
# 3       4.78949297995239            30058                     7
# ...                  ...              ...                   ...
# 191      13.488984589102             1883                     7
# 192     2.51930349520745            16289                     7
#     altexps_SIRV_percent     total
#                <numeric> <integer>
# 1        2.9016456005055    959042
# 2       3.29130111124368   1190198
# 3       2.36477183861837   1271074
# ...                  ...       ...
# 191     3.35089155425846     56194
# 192    0.843270890872805   1931645

# method_2, addPerCellQC()计算
library(scater)
sce_416b <- addCellQC(sce_416b, subsets = list(Mito = is_mito))
colnames(colData(sce_416b))
#  [1] "Source Name"              "cell line"
#  [3] "cell type"                "single cell well quality"
#  [5] "genotype"                 "phenotype"
#  [7] "strain"                   "spike-in addition"
#  [9] "block"                    "sum"
# [13] "percent_top_100"          "percent_top_200"
# [15] "percent_top_500"          "subsets_Mito_sum"
# [17] "subsets_Mito_detected"    "subsets_Mito_percent"
# [19] "altexps_ERCC_sum"         "altexps_ERCC_detected"
# [21] "altexps_ERCC_percent"     "altexps_SIRV_sum"
# [23] "altexps_SIRV_detected"    "altexps_SIRV_percent"
# [25] "total"


# filter out low quality cells ------------------------------------------------
# method_1, 使用固定阈值
# 基于之前perCellQCMetrics()计算的QC结果df
qc_lib <- df$sum < 1e5
qc_nexprs <- df$detected < 5e3
qc_spike <- df$altexps_ERCC_percent > 10
qc_mito <- df$subsets_Mito_percent > 10
# 都设定好以后传给discard
discard <- qc_lib | qc_nexprs | qc_spike | qc_mito

# 统计每种指标过滤的细胞数量
filter_out_cells <- DataFrame(
    LibSize = sum(qc_lib),
    NExprs = sum(qc_nexprs),
    SpikeProp = sum(qc_spike),
    MitoProp = sum(qc_mito),
    Total = sum(discard))
filter_out_cells
# DataFrame with 1 row and 5 columns
#     LibSize    NExprs SpikeProp  MitoProp     Total
#   <integer> <integer> <integer> <integer> <integer>
# 1         3         0        19         0        21
# 这里看到线粒体部分按照阈值为10%，是不能对细胞进行过滤的，
# 因为根据TxDb.Mmusculus.UCSC.mm10.ensGene一共才找到13个线粒体基因；
# 而作者的教程中可能找到的线粒体基因比较多，因此他最终过滤掉了14个细胞

# 另外看到Total这里并不是简单地将四项指标相加，这是因为可能一个细胞同时符合两项或多项过滤条件，最终只计算一次
# 比如下面就是查看：既符合qc_lib过滤又符合qc_spike过滤条件的细胞，第191个细胞
intersect(which(qc_lib == 1), which(qc_spike == 1))
# [1] 191

# method_2, 使用相对阈值，离群点
# isOutlier()鉴定离群点
qc_lib_2 <- isOutlier(df$sum, log = TRUE, type = "lower")
qc_nexprs_2 <- isOutlier(df$detected, log = TRUE, type = "lower")
qc_spike_2 <- isOutlier(df$altexps_ERCC_percent, type = "higher")
qc_mito_2 <- isOutlier(df$subsets_Mito_percent, type = "higher")

# 查看结果
attr(qc_lib_2, "thresholds")
#    lower   higher
# 434082.9      Inf

attr(qc_nexprs_2, "thresholds")
#    lower   higher
# 5231.468      Inf

attr(qc_spike_2, "thresholds")
#    lower   higher
#     -Inf 14.15371

attr(qc_mito_2, "thresholds")
#  lower higher
#   -Inf      0

# 使用相对阈值对过滤的细胞数量进行统计
discard_2 <- qc_lib_2 | qc_nexprs_2 | qc_spike_2 | qc_mito_2
filter_out_cells_2 <- DataFrame(
    LibSize = sum(qc_lib_2),
    NExprs = sum(qc_nexprs_2),
    SpikeProp = sum(qc_spike_2),
    MitoProp = sum(qc_mito_2),
    Total = sum(discard_2))
filter_out_cells_2
# DataFrame with 1 row and 5 columns
#     LibSize    NExprs SpikeProp  MitoProp     Total
#   <integer> <integer> <integer> <integer> <integer>
# 1         4         0         1         0         4

# quickPerCellQC()一步整合上面的操作
reasons <- quickPerCellQC(
    df,
    percent_subsets = c("subsets_Mito_percent", "altexps_ERCC_percent")
)
colSums(as.matrix(reasons))
#              low_lib_size            low_n_features
#                         4                         0
# high_subsets_Mito_percent high_altexps_ERCC_percent
#                         0                         1
#                   discard
#                         4

# considering batch effect during detecting outliers
# 一个是实验批次，两个细胞板
table(sce_416b$block)
# 20160113 20160325
#       96       96

# 一个是细胞表型，两种生物状态
table(sce_416b$phenotype)
# induced CBFB-MYH11 oncogene expression
#                                     96
#                    wild type phenotype
#                                     96

# 设定batch信息
batch <- paste0(sce_416b$phenotype, "-", sce_416b$block)
table(batch)
# batch
# induced CBFB-MYH11 oncogene expression-20160113
#                                              48
# induced CBFB-MYH11 oncogene expression-20160325
#                                              48
#                    wild type phenotype-20160113
#                                              48
#                    wild type phenotype-20160325
#                                              48
batch_reasons <- quickPerCellQC(
    df,
    percent_subsets = c("subsets_Mito_percent", "altexps_ERCC_percent"),
    batch = batch
)
colSums(as.matrix(batch_reasons))
#              low_lib_size            low_n_features
#                         5                         4
# high_subsets_Mito_percent high_altexps_ERCC_percent
#                         0                         6
#                   discard
#                         9

# refering to otherr batches to filter out problematical batches
library(scRNAseq)
sce_grun <- GrunPancreasData()
# batch information
table(sce_grun$donor)
# D10 D17  D2  D3  D7
# 288 480  96 480 384

sce_grun <- addPerCellQC(sce_grun)
# 分批次进行检测
discard_ercc <- isOutlier(
    sce_grun$altexps_ERCC_percent,
    type = "higher",
    batch = sce_grun$donor
)
with_blocking <- plotColData(
    sce_grun,
    x = "donor",
    y = "altexps_ERCC_percent",
    colour_by = I(discard_ercc)
)

# 先算其他几个批次的离群点
discard_ercc_2 <- isOutlier(
    sce_grun$altexps_ERCC_percent,
    type = "higher",
    batch = sce_grun$donor,
    subset = sce_grun$donor %in% c("D17", "D2", "D7")
)
without_blocking <- plotColData(
    sce_grun,
    x = "donor",
    y = "altexps_ERCC_percent",
    colour_by = I(discard_ercc_2)
)

# 拼图
gridExtra::grid.arrange(with_blocking, without_blocking, ncol = 2)


# identify problematical batches ----------------------------------------------
ercc_thresholds <- attr(discard_ercc, "thresholds")["higher", ]
ercc_thresholds
#        D10        D17         D2         D3         D7
#  73.610696   7.599947   6.010975 113.105828  15.216956

names(ercc_thresholds)[isOutlier(ercc_thresholds, type = "higher")]
# [1] "D10" "D3"


# use robustbase package ------------------------------------------------------
stats <- cbind(
    log10(df$sum),
    log10(df$detected),
    df$subsets_Mito_percent,
    df$altexps_ERCC_percent
)

library(robustbase)
outlying <- adjOutlyingness(fullRank(stats), only.outlyingness = TRUE)
multi_outlier <- isOutlier(outlying, type = "higher")

attr(multi_outlier, "thresholds")
#    lower   higher
#     -Inf 2.081908

summary(multi_outlier)
#    Mode   FALSE    TRUE
# logical     183       9


# make a plot to check --------------------------------------------------------
# The 1st figure, metrics of qc in different batches
library(scRNAseq)
sce_416b <- LunSpikeInData(which = "416b")
dim(colData(sce_416b))
# [1] 192   9

library(scater)
df <- perCellQCMetrics(sce_416b, subsets = list(Mito = is_mito))
# [1] 192  16

# 整合qc指标信息和初始细胞信息
colData(sce_416b) <- cbind(colData(sce_416b), df)
names(colData(sce_416b))
#  [1] "Source Name"              "cell line"
#  [3] "cell type"                "single cell well quality"
#  [5] "genotype"                 "phenotype"
#  [7] "strain"                   "spike-in addition"
#  [9] "block"                    "sum"
# [11] "detected"                 "percent_top_50"
# [13] "percent_top_100"          "percent_top_200"
# [15] "percent_top_500"          "subsets_Mito_sum"
# [17] "subsets_Mito_detected"    "subsets_Mito_percent"
# [19] "altexps_ERCC_sum"         "altexps_ERCC_detected"
# [21] "altexps_ERCC_percent"     "altexps_SIRV_sum"
# [23] "altexps_SIRV_detected"    "altexps_SIRV_percent"
# [25] "total"

# 修改整合后的信息
sce_416b$block <- factor(sce_416b$block)
sce_416b$phenotype <- ifelse(
    grepl("induced", sce_416b$phenotype), "induced", "wild type"
)

# 添加过滤的结果
reasons <- quickPerCellQC(
    df,
    percent_subsets = c("subsets_Mito_percent", "altexps_ERCC_percent")
)
sce_416b$discard <- reasons$discard

# 作图
gridExtra::grid.arrange(
    plotColData(sce_416b, x = "block", y = "sum",
        colour_by = "discard", other_fields = "phenotype") +
        facet_wrap(~phenotype) +
        scale_y_log10() +
        ggtitle("Total count"),
    plotColData(sce_416b, x = "block", y = "detected",
        colour_by = "discard", other_fields = "phenotype") +
        facet_wrap(~phenotype) +
        scale_y_log10() +
        ggtitle("Detected features"),
    plotColData(sce_416b, x = "block", y = "subsets_Mito_percent",
        colour_by = "discard", other_fields = "phenotype") +
        facet_wrap(~phenotype) +
        ggtitle("Mito percent"),
    plotColData(sce_416b, x = "block", y = "altexps_ERCC_percent",
        colour_by = "discard", other_fields = "phenotype") +
        facet_wrap(~phenotype) +
        ggtitle("ERCC percent"),
    ncol = 2
)

# The 2nd figure, mito_percent and library size
plotColData(
    sce_416b, x = "sum", y = "subsets_Mito_percent",
    colour_by = "discard", other_fields = c("block", "phenotype")
) +
facet_grid(block ~ phenotype) +
theme(panel.border = element_rect(color = "grey"))

# The 3rd figure, ERCC and mito_percent
plotColData(
    sce_416b, x = "altexps_ERCC_percent", y = "subsets_Mito_percent",
    colour_by = "discard", other_fields = c("block", "phenotype")
) +
facet_grid(block ~ phenotype) +
theme(panel.border = element_rect(color = "grey"))


# qc for droplet-based data, such as 10X Genomics -----------------------------
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

BiocManager::install("DropletUtils")
library(DropletUtils)
library(Matrix)
fname <- file.path(getwd(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce_pbmc <- read10xCounts(fname, col.names = TRUE)
sce_pbmc
# class: SingleCellExperiment 
# dim: 33694 737280 
# metadata(1): Samples
# assays(1): counts
# rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#   ENSG00000268674
# rowData names(2): ID Symbol
# colnames(737280): AAACCTGAGAAACCAT-1 AAACCTGAGAAACCGC-1 ...
#   TTTGTCATCTTTAGTC-1 TTTGTCATCTTTCCTC-1
# colData names(2): Sample Barcode
# reducedDimNames(0):
# spikeNames(0):
# altExpNames(0):

# dispaly distribution of all barcodes' library size
bcrank <- barcodeRanks(counts(sce_pbmc))
# 为了作图速度，每个rank只展示一个点（去重复）
uniq <- !duplicated(bcrank$rank)
barcodes_lib_distribution_plot <- plot(
    bcrank$rank[uniq],
    bcrank$total[uniq],
    log = "xy",
    xlab = "Rank",
    ylab = "Total UMI count",
    cex.lab = 1.2
) +
abline(h = metadata(bcrank)$inflection, col = "darkgreen", lty = 2) +
abline(h = metadata(bcrank)$knee, col = "dodgerblue", lty = 2) +
legend(
    "bottomleft", legend = c("Inflection", "Knee"),
    col = c("darkgreen", "dodgerblue"),
    lty = 2, cex = 1.2
)

# detect dropouts
set.seed(100)
e_out <- emptyDrops(counts(sce_pbmc))
summary(e_out$FDR <= 0.001)
# Mode   FALSE    TRUE    NA's
# logical    1056    4233  731991

# 保存结果
sce_pbmc_final <- sce_pbmc[, which(e_out$FDR <= 0.001)]

ncol(counts(sce_pbmc))
# [1] 737280
ncol(counts(sce_pbmc_final))
# [1] 4233

# combining emptyDrops() and qc metrics
is_mito <- grep("^MT-", rowData(sce_pbmc_final)$Symbol)
# [1] 13
pbmc_qc_final <- perCellQCMetrics(sce_pbmc_final, subsets = list(MT = is_mito))
discard_mito <- isOutlier(pbmc_qc_final$subsets_MT_percent, type = "higher")
summary(discard_mito)
#    Mode   FALSE    TRUE
# logical    3922     311

plot(pbmc_qc_final$sum, pbmc_qc_final$subsets_MT_percent,
    log = "x",
    xlab = "Total count",
    ylab = "Mitochondrial %"
) +
abline(h = attr(discard_mito, "thresholds")["higher"], col = "red")


# filter out cells with low quality -------------------------------------------
# 直接去除
sce_416b_filterd <- sce_416b[, !reasons$discard]

# 先做标记
# example 1
# 计算舍弃组和保留组的平均表达量及logFC
kept <- calculateAverage(counts(sce_416b)[, !discard])
lost <- calculateAverage(counts(sce_416b)[, discard])

library(edgeR)
# 计算log(cpm + 2)
logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)

abundance <- rowMeans(logged)
log_fc <- logged[, 1] - logged[, 2]
plot(abundance, log_fc,
    xlab = "Average count", ylab = "log_lost_kept", pch = 16)
points(abundance[is_mito], log_fc[is_mito], col = "dodgerblue", pch = 16)

# example 2
alt_discard <- colSums(counts(sce_pbmc_final)) < 500
kept <- calculateAverage(counts(sce_pbmc_final)[, !alt_discard])
lost <- calculateAverage(counts(sce_pbmc_final)[, alt_discard])
logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
abundance <- rowMeans(logged)
log_fc <- logged[, 1] - logged[, 2]
plot(abundance, log_fc,
    xlab = "Average count", ylab = "log_lost_kept", pch = 16)

platelet <- c("PF4", "PPBP", "CAVIN2")
library(org.Hs.eg.db)
ids <- mapIds(
    org.Hs.eg.db,
    keys = platelet, column = "ENSEMBL", keytype = "SYMBOL"
)
points(abundance[ids], log_fc[ids], col = "orange", pch = 16)

# 只做标记
marked <- sce_416b
marked$discard <- batch_reasons$discard
