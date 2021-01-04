########################################
# dimension reduce
# date: 2021.01.01 - 01.04
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-4
########################################


# rm all objects ---------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
#################### sce_zeisel
library(scRNAseq)
load("sce_zeisel.RData")
# sce_zeisel <- ZeiselBrainData()
sce_zeisel <- sce.zeisel

library(scater)
# 将每个细胞的所有gene表达量加起来，得到每个细胞的文库大小
# 同时替换一些奇怪的gene名
sce_zeisel <- aggregateAcrossFeatures(
    sce_zeisel,
    id = sub("_loc[0-9]+$", "", rownames(sce_zeisel))
)
dim(sce_zeisel)
# [1] 19839  3005

# gene annotation
library(org.Mm.eg.db)
rowData(sce_zeisel)$Ensembl <- mapIds(
    org.Mm.eg.db,
    keys = rownames(sce_zeisel),
    keytype = "SYMBOL",
    column = "ENSEMBL"
)

# qc, 先perCellQCMetrics()，后quickPerCellQC() -----------------------------------
stats <- perCellQCMetrics(
    sce_zeisel,
    subsets = list(Mt = rowData(sce_zeisel)$featureType == "mito")
)
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent")
)
sce_zeisel_filtered <- sce_zeisel[, !qc$discard]
dim(sce_zeisel_filtered)
# [1] 19839  2816


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(1000)
clust_zeisel <- quickCluster(sce_zeisel_filtered)
sce_zeisel_filtered <- computeSumFactors(
    sce_zeisel_filtered,
    cluster = clust_zeisel
)
# logNormCounts()
sce_zeisel_filtered <- logNormCounts(sce_zeisel_filtered)


# measure the degree of change using log-counts, with spike-in ----------------
# and HVGs selection by proportion, top 10%
dec_zeisel <- modelGeneVarWithSpikes(sce_zeisel_filtered, "ERCC")
top_hvgs_zeisel <- getTopHVGs(dec_zeisel, prop = 0.1)
length(top_hvgs_zeisel)
# [1] 1816

sce_zeisel_filtered
# class: SingleCellExperiment
# dim: 19839 2816
# metadata(0):
# assays(2): counts logcounts
# rownames(19839): 0610005C13Rik 0610007N19Rik ... Zzef1 Zzz3
# rowData names(2): featureType Ensembl
# colnames(2816): 1772071015_C02 1772071017_G12 ... 1772063068_D01 1772066098_A12
# colData names(11): tissue group # ... level2class sizeFactor
# reducedDimNames(0):
# altExpNames(2): ERCC repeat





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

sce_pbmc_final
# class: SingleCellExperiment
# dim: 33694 3985
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(33694): RP11-34P13.3 FAM138A ... AC213203.1 FAM231B
# rowData names(2): ID Symbol
# colnames(3985): AAACCTGAGAAGGCCT-1 AAACCTGAGACAGACC-1 ... TTTGTCAGTTAAGACA-1
#   TTTGTCATCCCAAGAT-1
# colData names(3): Sample Barcode sizeFactor
# reducedDimNames(0):
# altExpNames(0):



# dimension reduce ------------------------------------------------------------
##### PCA
library(scran)
top_hvgs_zeisel_new <- getTopHVGs(dec_zeisel, n = 2000)


library(scater)
set.seed(100)   # PCA是随机的
sce_zeisel_filtered <- runPCA(
    sce_zeisel_filtered,
    subset_row = top_hvgs_zeisel_new
)
reducedDimNames(sce_zeisel_filtered)
# [1] "PCA"
dim(reducedDim(sce_zeisel_filtered, "PCA")) # 2816个细胞，50个PCs
# [1] 2816   50

## SVD, 使用奇异值分解
library(BiocSingular)
set.seed(1000)
sce_zeisel_filtered <- runPCA(
    sce_zeisel_filtered,
    subset_row = top_hvgs_zeisel_new,
    BSPARAM = RandomParam(),
    name = "IRLBA"
)
reducedDimNames(sce_zeisel_filtered)
# [1] "PCA"   "IRLBA"
dim(reducedDim(sce_zeisel_filtered, "IRLBA"))
# [1] 2816   50


# PCs selection ---------------------------------------------------------------
### method_1, using elbow point
# 首先提取各个PCs对整体差异的贡献比例
percent_var <- attr(reducedDim(sce_zeisel_filtered), "percentVar")

# 辅助选择
BiocManager::install("PCAtools")
chosen_elbow <- PCAtools::findElbowPoint(percent_var)
chosen_elbow
# [1] 7
# 画图查看
plot(percent_var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen_elbow, col = "red")

### method_2, using technical noise
library(scran)
set.seed(111001001)
denoised_pbmc_pois <- denoisePCA(
    sce_pbmc_final,
    technical = dec_pbmc_pois,
    subset.row = top_hvgs_pbmc
)
ncol(reducedDim(denoised_pbmc_pois))
# [1] 9

set.seed(001001001)
denoised_seizel <- denoisePCA(
    sce_zeisel_filtered,
    technical = dec_zeisel,
    subset.row = top_hvgs_zeisel_new
)
ncol(reducedDim(denoised_seizel))
# [1] 50

# modelGeneVar() and denoisePCA()
dec_pbmc_var <- modelGeneVar(sce_pbmc_final)
denoised_pbmc_var <- denoisePCA(
    sce_pbmc_final,
    technical = dec_pbmc_var,
    subset.row = top_hvgs_pbmc
)
ncol(reducedDim(denoised_pbmc_var))
# [1] 5

### method_3, based on cell clustering
pcs <- reducedDim(sce_zeisel_filtered)
choices <- getClusteredPCs(pcs)
metadata(choices)$chosen
# [1] 17

plot(
    choices$n.pcs, choices$n.clusters,
    xlab = "The number of PCs",
    ylab = "The number of clusters"
)
abline(a = 1, b = 1, col = "red")
abline(v = metadata(choices)$chosen, col = "grey80", lty = 2)

### method_4, self-determined PCs’ number
# 直接在PCA结果上修改
reducedDim(sce_zeisel_filtered, "PCA") <- reducedDim(sce_zeisel_filtered,
    "PCA")[, 1:20]
ncol(reducedDim(sce_zeisel_filtered, "PCA"))
# [1] 20

# 赋值新的变量名，新增一个PCA结果
reducedDim(sce_zeisel_filtered, "PCA_20") <- reducedDim(sce_zeisel_filtered,
    "PCA")[, 1:20]
reducedDimNames(sce_zeisel_filtered)
# [1] "PCA"    "IRLBA"  "PCA_20"


# NMF, 非负矩阵分解 -----------------------------------------------------------------
# 注意，需要先安装NMF包
install.packages("NMF")
library(NMF)

library(scater)
set.seed(101001)
nmf_zeisel <- runNMF(
    sce_zeisel_filtered,
    ncomponents = 10,
    subset_row = top_hvgs_zeisel_new
)

# 提取结果
nmf_out <- reducedDim(nmf_zeisel, "NMF")
nmf_basis <- attr(nmf_out, "basis")
colnames(nmf_out) <- colnames(nmf_basis) <- 1:10

# 画每个细胞与因子的关系
per_cell <- pheatmap::pheatmap(
    nmf_out,
    silent = TRUE,
    main = "By cell",
    show_rownames = FALSE,
    color = rev(viridis::magma(100)),
    cluster_cols = FALSE
)

# 画每个gene与因子的关系
per_gene <- pheatmap::pheatmap(
    nmf_basis,
    silent = TRUE,
    main = "By gene",
    show_rownames = FALSE,
    color = rev(viridis::magma(100)),
    cluster_cols = FALSE
)

# 组图
gridExtra::grid.arrange(per_cell[[4]], per_gene[[4]], ncol = 2)


# plot for dim reduce ---------------------------------------------------------
##### visualization of PCA
plotReducedDim(
    sce_zeisel_filtered,
    dimred = "PCA",
    colour_by = "level1class"
)

# 绘制多个PCs
plotReducedDim(
    sce_zeisel_filtered,
    dimred = "PCA",
    ncomponents = 4,
    colour_by = "level1class"
)

##### visualization of t-SNE
set.seed(00101001101)
install.packages("Rtsne")
library(Rtsne)
# runTSNE()的结果也存储在reducedDims()中
sce_zeisel_filtered <- runTSNE(
    sce_zeisel_filtered,
    dimred = "PCA"
)
plotReducedDim(
    sce_zeisel_filtered,
    dimred = "TSNE",
    colour_by = "level1class"
)

### different values of perplexity argument
# perplexity = 5
set.seed(100)
sce_zeisel_filtered_5 <- runTSNE(
    sce_zeisel_filtered, dimred = "PCA", perplexity = 5)
out_5 <- plotReducedDim(
    sce_zeisel_filtered_5,
    dimred = "TSNE", colour_by = "level1class"
) + ggtitle("perplexity = 5")

# perplexity = 20
set.seed(100)
sce_zeisel_filtered_20 <- runTSNE(
    sce_zeisel_filtered, dimred = "PCA", perplexity = 20)
out_20 <- plotReducedDim(
    sce_zeisel_filtered_20,
    dimred = "TSNE", colour_by = "level1class"
) + ggtitle("perplexity = 20")

# perplexity = 80
set.seed(100)
sce_zeisel_filtered_80 <- runTSNE(
    sce_zeisel_filtered, dimred = "PCA", perplexity = 80)
out_80 <- plotReducedDim(
    sce_zeisel_filtered_80,
    dimred = "TSNE", colour_by = "level1class"
) + ggtitle("perplexity = 80")

gridExtra::grid.arrange(out_5, out_20, out_80, ncol = 3)

##### visualization of UMAP
install.packages("uwot")
library(uwot)
set.seed(1100101001)
sce_zeisel_filtered <- runUMAP(sce_zeisel_filtered, dimred = "PCA")
plotReducedDim(
    sce_zeisel_filtered,
    dimred = "UMAP", colour_by = "level1class"
)
