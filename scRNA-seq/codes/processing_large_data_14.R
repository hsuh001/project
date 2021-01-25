########################################
# processing large-scale data
# date: 2021.01.25 - 01.25
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-14
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
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


# dimension reduce, using three methods ---------------------------------------
##### PCA
set.seed(10000)
sce_pbmc_final <- denoisePCA(
    sce_pbmc_final,
    subset.row = top_hvgs_pbmc,
    technical = dec_pbmc_pois
)
dim(reducedDim(sce_pbmc_final, "PCA"))
# [1] 3985    9

##### t-SNE
set.seed(100000)
sce_pbmc_final <- runTSNE(sce_pbmc_final, dimred = "PCA")
dim(reducedDim(sce_pbmc_final, "TSNE"))
# [1] 3985    2

##### UMAP
set.seed(1000000)
sce_pbmc_final <- runUMAP(sce_pbmc_final, dimred = "PCA")
dim(reducedDim(sce_pbmc_final, "UMAP"))
# [1] 3985    2


# clustering ------------------------------------------------------------------
##### graph-based clustering
# 使用PCA的前几个PCs构建SNNG
library(scran)
g_pbmc_k10 <- buildSNNGraph(sce_pbmc_final, k = 10, use.dimred = "PCA")

# 鉴定cluster
clust_pbmc_k10 <- igraph::cluster_walktrap(g_pbmc_k10)$membership
table(clust_pbmc_k10)
# clust_pbmc_k10
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
# 205 508 541  56 374 125  46 432 302 867  47 155 166  61  84  16

# 把cluster信息存成SingleCellExperiment对象的一个因子
library(scater)
colLabels(sce_pbmc_final) <- factor(clust_pbmc_k10)

sce_pbmc_final
# class: SingleCellExperiment
# dim: 33694 3985
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(33694): RP11-34P13.3 FAM138A ... AC213203.1 FAM231B
# rowData names(2): ID Symbol
# colnames(3985): AAACCTGAGAAGGCCT-1 AAACCTGAGACAGACC-1 ... TTTGTCAGTTAAGACA-1
#   TTTGTCATCCCAAGAT-1
# colData names(4): Sample Barcode sizeFactor label
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):


##### 快速估算
## 使用近邻搜素
library(scran)
BiocManager::install("BiocNeighbors")
library(BiocNeighbors)
snn_gr <- buildSNNGraph(
    sce_pbmc_final,
    BNPARAM = AnnoyParam(),
    use.dimred = "PCA"
)

# 鉴定cluster
clusters <- igraph::cluster_walktrap(snn_gr)$membership
table(Exact = colLabels(sce_pbmc_final), Approx = clusters)
#      Approx
# Exact   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
#    1  205   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
#    2    0 479   0   0   2   0   0   0   0  27   0   0   0   0   0   0
#    3    0   0 540   0   1   0   0   0   0   0   0   0   0   0   0   0
#    4    0   0   0  55   0   0   0   0   0   1   0   0   0   0   0   0
#    5    0  25   0   0 349   0   0   0   0   0   0   0   0   0   0   0
#    6    0   0   0   0   0 125   0   0   0   0   0   0   0   0   0   0
#    7    0   0   0   0   0   0  46   0   0   0   0   0   0   0   0   0
#    8    0   0   0   0   0   0   0 432   0   0   0   0   0   0   0   0
#    9    0   0   0   1   0   0   0  10 291   0   0   0   0   0   0   0
#    10   0  28   0   0   0   0   0   0   0 839   0   0   0   0   0   0
#    11   0   0   0   0   0   0   0   0   0   0  47   0   0   0   0   0
#    12   0   0   0   0   0   0   0   0   0   0   0 155   0   0   0   0
#    13   0   0   0   0   0   0   0   0   0   0   0   0 166   0   0   0
#    14   0   0   0   0   0   0   0   0   0   0   0   0   0  61   0   0
#    15   0   0   0   0   0   0   0   0   0   0   0   0   0   0  84   0
#    16   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  16

## 奇异值分解
library(scater)
library(BiocSingular)

# method_1, randomized SVD (RSVD)
set.seed(101000)
r_out <- runPCA(
    sce_pbmc_final,
    ncomponents = 20,
    BSPARAM = RandomParam()
)
str(reducedDim(r_out))
#  num [1:3985, 1:20] 15.05 13.43 -8.67 -7.74 6.45 ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : chr [1:3985] "AAACCTGAGAAGGCCT-1" "AAACCTGAGACAGACC-1" "AAACCTGAGGCATGGT-1" "AAACCTGCAAGGTTCT-1" ...
#   ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...
#  - attr(*, "varExplained")= num [1:20] 85.36 40.43 23.22 8.99 6.66 ...
#  - attr(*, "percentVar")= num [1:20] 19.85 9.4 5.4 2.09 1.55 ...
#  - attr(*, "rotation")= num [1:500, 1:20] 0.203 0.1834 0.1779 0.1063 0.0647 ...
#   ..- attr(*, "dimnames")=List of 2
#   .. ..$ : chr [1:500] "LYZ" "S100A9" "S100A8" "HLA-DRA" ...
#   .. ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...

# method_2, IRLBA
set.seed(101001)
i_out <- runPCA(
    sce_pbmc_final,
    ncomponents = 20,
    BSPARAM = IrlbaParam()
)
str(reducedDim(i_out))
#  num [1:3985, 1:20] 15.05 13.43 -8.67 -7.74 6.45 ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : chr [1:3985] "AAACCTGAGAAGGCCT-1" "AAACCTGAGACAGACC-1" "AAACCTGAGGCATGGT-1" "AAACCTGCAAGGTTCT-1" ...
#   ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...
#  - attr(*, "varExplained")= num [1:20] 85.36 40.43 23.22 8.99 6.66 ...
#  - attr(*, "percentVar")= num [1:20] 19.85 9.4 5.4 2.09 1.55 ...
#  - attr(*, "rotation")= num [1:500, 1:20] 0.203 0.1834 0.1779 0.1063 0.0647 ...
#   ..- attr(*, "dimnames")=List of 2
#   .. ..$ : chr [1:500] "LYZ" "S100A9" "S100A8" "HLA-DRA" ...
#   .. ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...


# 并行计算 ------------------------------------------------------------------------
library(BiocParallel)
##### 查看支持加速的类型。最顶上的是默认的。
registered()
# $SnowParam
# class: SnowParam
#   bpisup: FALSE; bpnworkers: 6; bptasks: 0; bpjobname: BPJOB
#   bplog: FALSE; bpthreshold: INFO; bpstopOnError: TRUE
#   bpRNGseed: ; bptimeout: 2592000; bpprogressbar: FALSE
#   bpexportglobals: TRUE
#   bplogdir: NA
#   bpresultdir: NA
#   cluster type: SOCK
# 
# $SerialParam
# class: SerialParam
#   bpisup: FALSE; bpnworkers: 1; bptasks: 0; bpjobname: BPJOB
#   bplog: FALSE; bpthreshold: INFO; bpstopOnError: TRUE
#   bpRNGseed: ; bptimeout: 2592000; bpprogressbar: FALSE
#   bpexportglobals: TRUE
#   bplogdir: NA
#   bpresultdir: NA

##### 更改加速的类型
default <- registered()

# 改为BatchtoolsParam
register(BatchtoolsParam(workers = 10), default = TRUE)
names(registered())

# 恢复原来的设置
for(param in rev(default)) {
    register(param)
}

##### 使用方法
# 不同方法
dec_pbmc_mc <- modelGeneVar(sce_pbmc_final, BPPARAM = MulticoreParam(2))
dec_pbmc_snow <- modelGeneVar(sce_pbmc_final, BPPARAM = SnowParam(5))

##### 在SLURM HPC中使用BatchtoolsParam()
# 设置每个任务2h、8G内存、1CPU、总共10个任务
bpp <- BatchtoolsParam(
    10,
    cluster = "slurm",
    resources = list(
        walltime = 7200,
        memory = 8000,
        ncpus = 1
    )
)


# 可能会遇到内存不足 -------------------------------------------------------------------
# 从130万个脑细胞数据选取2万个
BiocManager::install("TENxBrainData")
library(TENxBrainData)
sce_brain <- TENxBrainData20k()
sce_brain
# class: SingleCellExperiment
# dim: 27998 20000
# metadata(0):
# assays(1): counts
# rownames: NULL
# rowData names(2): Ensembl Symbol
# colnames: NULL
# colData names(4): Barcode Sequence Library Mouse
# reducedDimNames(0):
# altExpNames(0):

# 查看表达矩阵
counts(sce_brain)
# <27998 x 20000> matrix of class HDF5Matrix and type "integer":
#              [,1]     [,2]     [,3]     [,4] ... [,19997] [,19998] [,19999] [,20000]
#     [1,]        0        0        0        0   .        0        0        0        0
#     [2,]        0        0        0        0   .        0        0        0        0
#     [3,]        0        0        0        0   .        0        0        0        0
#     [4,]        0        0        0        0   .        0        0        0        0
#     [5,]        0        0        0        0   .        0        0        0        0
#      ...        .        .        .        .   .        .        .        .        .
# [27994,]        0        0        0        0   .        0        0        0        0
# [27995,]        0        0        0        1   .        0        2        0        0
# [27996,]        0        0        0        0   .        0        1        0        0
# [27997,]        0        0        0        0   .        0        0        0        0
# [27998,]        0        0        0        0   .        0        0        0        0

# 查看大小
object.size(counts(sce_brain))
# 2560 bytes

# 实际大小
file.info(path(counts(sce_brain)))$size
# [1] 76264332
