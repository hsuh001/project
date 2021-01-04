########################################
# understanding the SingleCellExperiment object
# date: 2020.10.07
# author: Jing Xiao
########################################


# set work directory ----------------------------------------------------------
work_dir <- "/home1/jxiao/project/scRNA-seq"
setwd(work_dir)


# rm all objects ---------------------------------------------------------------
rm(list = ls())


# load packages ---------------------------------------------------------------
library(SingleCellExperiment)


###############################################################################
# hsuh，写于2020.10.07，更新于2020.10.07
###############################################################################
######################################### 1.1 数据结构：SingleCellExperiment对象

# build sce object ------------------------------------------------------------
# 1、构建表达矩阵。行是gene、列是样本。一定要是matrix
counts_df <- data.frame(
    cell_1 = rpois(10, 10),
    cell_2 = rpois(10, 10),
    cell_3 = rpois(10, 10)
)
rownames(counts_df) <- paste0("gene_", 1:10)
counts_matrix <- as.matrix(counts_df)

# 2、构建assays，assays是一个list，这个list可以包括任意个矩阵。再构建sce对象
sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

# 查看sce的内容
# > sce
## class: SingleCellExperiment
## dim: 10 3
## metadata(0):
## assays(1): counts
## rownames(10): gene_1 gene_2 ... gene_9 gene_10
## rowData names(0):
## colnames(3): cell_1 cell_2 cell_3
## colData names(0):
## reducedDimNames(0):
## spikeNames(0):
## altExpNames(0):

# 3、提取这个sce对象中的矩阵，有两种方法：
# ① 使用assay(sce, "counts")。最通用的办法，这里的"counts"可以换成其他名称，只要是出现在之前的list中都可以
assay(sce, "counts")

# ② 使用counts(sce)。只能提取名称为"counts"的矩阵
counts(sce)

##         cell_1 cell_2 cell_3
## gene_1      11      8      6
## gene_2       9     12     10
## gene_3       9      9     11
## gene_4       8     10      8
## gene_5       7      9      8
## gene_6      12      9      6
## gene_7       9     14     11
## gene_8      13      9     13
## gene_9       7     14     12
## gene_10      8      5     11

# 4、初始assays只有原始表达矩阵，可以多向扩展，也有两种方式：
# ① 使用一些R包的标准函数扩展，如扩展到归一化矩阵
sce <- scran::computeSumFactors(sce)
sce <- scater::normalize(sce)

# > sce
## class: SingleCellExperiment
## dim: 10 3
## metadata(1): log.exprs.offset
## assays(2): counts logcounts
## rownames(10): gene_1 gene_2 ... gene_9 gene_10
## rowData names(0):
## colnames(3): cell_1 cell_2 cell_3
## colData names(0):
## reducedDimNames(0):
## spikeNames(0):
## altExpNames(0):

# 使用assay(sce, "logcounts")或logcounts(sce)提取名称为"logcounts"的矩阵
assay(sce, "logcounts")
logcounts(sce)

# 使用assays(sce)查看结果
assays(sce)
## List of length 2
## names(2): counts logcounts

# ② 使用assay()自定义扩展assays，由个人需求根据原始矩阵得到新的矩阵
counts_100 <- assay(sce, "counts") + 100
assay(sce, "counts_100") <- counts_100

# 使用assays(sce)查看结果
assays(sce)
## List of length 3
## names(3): counts logcounts counts_100


# build cols annotation: colData ----------------------------------------------
# 1、手动构建细胞批次注释信息
cell_metadata <- data.frame(batch = c(1, 1, 2))
rownames(cell_metadata) <- paste0("cell_", 1:3)

# 2、将注释信息加入sce对象，两种方式：
# ① 直接构建
sce <- SingleCellExperiment(
    assays = list(counts = counts_matrix),
    colData = cell_metadata
)

# ② 后续添加
sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
colData(sce) <- DataFrame(cell_metadata)

# > sce
## class: SingleCellExperiment
## dim: 10 3
## metadata(0):
## assays(1): counts
## rownames(10): gene_1 gene_2 ... gene_9 gene_10
## rowData names(0):
## colnames(3): cell_1 cell_2 cell_3
## colData names(1): batch
## reducedDimNames(0):
## spikeNames(0):
## altExpNames(0):

# 3、获取加入sce对象以后的batch信息，两种方式：
# ① 使用colData()
colData(sce)
## DataFrame with 3 rows and 1 column
##            batch
##        <numeric>
## cell_1         1
## cell_2         1
## cell_3         2

# ② sce$batch直接查看结果
sce$batch
## [1] 1 1 2

# 4、使用一些R包构建注释信息。例如scater包的calculateQCMetrics()会计算几十项细胞的质量信息，
# 结果依然使用colData调用注释信息
sce <- SingleCellExperiment(
    assays = list(counts = counts_matrix),
    colData = cell_metadata
)
sce <- scater::addPerCellQC(sce)

# > sce
## class: SingleCellExperiment
## dim: 10 3
## metadata(0):
## assays(1): counts
## rownames(10): gene_1 gene_2 ... gene_9 gene_10
## rowData names(0):
## colnames(3): cell_1 cell_2 cell_3
## colData names(8): batch sum ... percent_top_500 total
## reducedDimNames(0):
## spikeNames(0):
## altExpNames(0):

colData(sce)[, 1:5]
## DataFrame with 3 rows and 5 columns
##            batch       sum  detected percent_top_50 percent_top_100
##        <numeric> <integer> <integer>      <numeric>       <numeric>
## cell_1         1        93        10            100             100
## cell_2         1        99        10            100             100
## cell_3         2        96        10            100             100

# 也可以任意添加
sce$more_stuff <- runif(ncol(sce))  # ncol(sce)值为3
sce$more_stuff
## [1] 0.3879090 0.2464490 0.1110965

colnames(colData(sce))
## [1] "batch"  "sum"   "detected"  "percent_top_50"    "percent_top_100"
## [6] "percent_top_200"    "percent_top_500"   "total" "more_stuff"

# 5、选取colData中的一部分
# 提取前3列注释信息
colData(sce)[, 1:3]
## DataFrame with 3 rows and 3 columns
##            batch       sum  detected
##        <numeric> <integer> <integer>
## cell_1         1        93        10
## cell_2         1        99        10
## cell_3         2        96        10

# 提取batch信息，数据库取子集的思路，两种方式：
colData(sce)$batch

sce$batch

# 提取batch中为1的信息，数据框的思路取子集
sce_1 <- sce[, sce$batch == 1]

# > sce_1
## class: SingleCellExperiment
## dim: 10 2
## metadata(0):
## assays(1): counts
## rownames(10): gene_1 gene_2 ... gene_9 gene_10
## rowData names(0):
## colnames(2): cell_1 cell_2
## colData names(9): batch sum ... total more_stuff
## reducedDimNames(0):
## spikeNames(0):
## altExpNames(0):


# build rows annotation: rowData / rowRanges-----------------------------------
# 1、rowData，是一个数据框
sce <- SingleCellExperiment(
    assays = list(counts = counts_matrix),
    colData = cell_metadata
)
rowData(sce)    # 一开始是空的
## DataFrame with 10 rows and 0 columns

# 使用R包构建行注释信息
sce <- scater::addPerFeatureQC(sce)
rowData(sce)
## DataFrame with 10 rows and 2 columns
##                     mean  detected
##                <numeric> <numeric>
## gene_1  8.33333333333333       100
## gene_2  10.3333333333333       100
## gene_3  9.66666666666667       100
## gene_4  8.66666666666667       100
## gene_5                 8       100
## gene_6                 9       100
## gene_7  11.3333333333333       100
## gene_8  11.6666666666667       100
## gene_9                11       100
## gene_10                8       100

# 2、rowRanges，是一个GRange对象，存储基因坐标信息：chr、start、end等
rowRanges(sce)  # 一开始是空的

## GRangesList object of length 10:
## $gene_1
## GRanges object with 0 ranges and 0 metadata columns:
##    seqnames    ranges strand
##       <Rle> <IRanges>  <Rle>
##   -------
##   seqinfo: no sequences
## 
## $gene_2
## GRanges object with 0 ranges and 0 metadata columns:
##    seqnames    ranges strand
##       <Rle> <IRanges>  <Rle>
##   -------
##   seqinfo: no sequences
## 
## $gene_3
## GRanges object with 0 ranges and 0 metadata columns:
##    seqnames    ranges strand
##       <Rle> <IRanges>  <Rle>
##   -------
##   seqinfo: no sequences
## 
## ...
## <7 more elements>

# 3、按行取子集，可以按位置、或者按名称
sce_2 <- sce[c("gene_1", "gene_4"), ]

sce_2 <- sce[c(1, 4), ]

# > sce_2
## class: SingleCellExperiment
## dim: 2 3
## metadata(0):
## assays(1): counts
## rownames(2): gene_1 gene_4
## rowData names(2): mean detected
## colnames(3): cell_1 cell_2 cell_3
## colData names(1): batch
## reducedDimNames(0):
## spikeNames(0):
## altExpNames(0):


# 对样本进行归一化，sizeFactors() ------------------------------------------------------
sce <- SingleCellExperiment(
    assays = list(counts = counts_matrix),
    colData = cell_metadata
)
sce <- scran::computeSumFactors(sce)
sizeFactors(sce)
## [1] 0.96875 1.03125 1.00000

# 或者手动添加
sizeFactors(sce) <- scater::librarySizeFactors(sce)
sizeFactors(sce)
## cell_1  cell_2  cell_3
## 0.96875 1.03125 1.00000


# 降维结果，reducedDims ------------------------------------------------------------
# 1、使用scater计算PCA
sce <- scater::logNormCounts(sce)
sce <- scater::runPCA(sce)
# 这个算法利用了sce对象的归一化结果logcounts(sce)
reducedDim(sce, "PCA")

##               PC1        PC2
## cell_1 -0.8880790 -0.1107840
## cell_2  0.5497601 -0.6425679
## cell_3  0.3383189  0.7533519
## attr(,"percentVar")
## [1] 54.83772 45.16228

# 2、使用scater计算tSNE
sce <- scater::runTSNE(sce, perplexity = 0.1)
reducedDim(sce, "TSNE")

##              [,1]      [,2]
## cell_1  5309.0777 -2051.699
## cell_2 -4417.5601 -3569.787
## cell_3  -891.5176  5621.486

# 3、使用reducedDims(sce)查看全部结果。这里是复数形式，与assays()查看结果一样
reducedDims(sce)

## List of length 2
## names(2): PCA TSNE

# 4、手动添加自己计算的降维结果，而不用scater包
# 进行UMAP降维，可以用scater::runUMAP()，但依然可以自己处理
# 使用uwot包，不过这个包只能计算，需要手动将结果添加到sce对象
u <- uwot::umap(t(logcounts(sce)), n_neighbors = 2)
reducedDim(sce, "UMAP_uwot") <- u

# >u
##              [,1]       [,2]
## cell_1 -0.0645930  0.6433452
## cell_2  0.3503579 -0.5297200
## cell_3 -0.2857649 -0.1136252
## attr(,"scaled:center")
## [1] -2.457387 10.183555

# 使用reducedDims(sce)查看全部结果
reducedDims(sce)

## List of length 3
## names(3): PCA TSNE UMAP_uwot


# metadata接口 ------------------------------------------------------------------
# 许多DIY的接口，不能直接导入，但仍需要这些信息。metadata接口是一个list，
# 可以存储任意类型的数据，只要给它一个名字。
# 例如，有几个感兴趣的基因（比如是高变化基因），要把它保存在`sce`对象中，便于以后使用：
my_genes <- c("gene_4", "gene_8")
metadata(sce) <- list(interest_genes = my_genes)

metadata(sce)

## $interest_genes
## [1] "gene_4" "gene_8"

# metadata接口是一个list，就意味着支持扩展
your_genes <- c("gene_1", "gene_5")
metadata(sce)$other_genes <- your_genes

metadata(sce)

## $interest_genes
## [1] "gene_4" "gene_8"
## 
## $other_genes
## [1] "gene_1" "gene_5"


# spike-in --------------------------------------------------------------------
# 还是可以继续使用isSpike来添加，虽然会显示'isSpike' is deprecated
isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))

# 结果就会在sce中添加：
## spikeNames(1): ERCC

spikeNames(sce)
## [1] "ERCC"

# 查看spike-in个数
table(isSpike(sce, "ERCC"))
## FALSE 
##    10

# 现在推出了alternative Experiments
spike_counts <- cbind(
    cell_1 = rpois(5, 10),
    cell_2 = rpois(5, 10),
    cell_3 = rpois(5, 10)
)
rownames(spike_counts) <- paste0("spike_", 1:5)
spike_se <- SummarizedExperiment(list(counts = spike_counts))

# > spike_se
## class: SummarizedExperiment
## dim: 5 3
## metadata(0):
## assays(1): counts
## rownames(5): spike_1 spike_2 spike_3 spike_4 spike_5
## rowData names(0):
## colnames(3): cell_1 cell_2 cell_3
## colData names(0):

# 然后通过altExp()加入进来，操作和assay()、reducedDim()类似
altExp(sce, "spike") <- spike_se
altExps(sce)

## List of length 1
## names(1): spike

# 提取数据也类似
sub <- sce[, 1:2]
altExp(sub, "spike")

## class: SummarizedExperiment
## dim: 5 2
## metadata(0):
## assays(1): counts
## rownames(5): spike_1 spike_2 spike_3 spike_4 spike_5
## rowData names(0):
## colnames(2): cell_1 cell_2
## colData names(0):


# 对列添加label -------------------------------------------------------------------
# 使用colLabels()，尤其在非监督聚类过程中对细胞添加label，进行分组
# 例如，scran::findMarkers()就是通过colLabels()提取细胞信息
# 注意，colLabels()在v1.8.0并不存在，在v1.9.3才存在
colLabels(sce) <- LETTERS[1:3]
colLabels(sce)
## [1] "A" "B" "C"



###############################################################################
# # hsuh，写于2020.10.07，更新于2020.10.08
###############################################################################
######################################################## 1.2 总览 | 从实验到分析

################################################## Example for scRNAseq package
# set work directory ----------------------------------------------------------
work_dir <- "/home1/jxiao/project/scRNA-seq_practice"
setwd(work_dir)


# rm all objects ---------------------------------------------------------------
rm(list = ls())


# load packages ---------------------------------------------------------------
library(scRNAseq)
library(scater)
library(scran)


# 质控 --------------------------------------------------------------------------
sce <- MacoskoRetinaData()
# 筛选以"MT-"开头的gene
is_mito <- grepl("^MT-", rownames(sce))
qc_states <- perCellQCMetrics(sce, subsets = list(mito = is_mito))
filterd <- quickPerCellQC(qc_states, percent_subsets = "subsets_mito_percent")
sce_filterd <- sce[, !filterd$discard]


# 归一化 -------------------------------------------------------------------------
sce_filterd <- logNormCounts(sce_filterd)


# 挑选高变化基因 ------------------------------------------------------------------
dec <- modelGeneVar(sce_filterd)
hvg <- getTopHVGs(dec, prop = 0.1)


# PCA、UMAP降维 ------------------------------------------------------------------
set.seed(1234)
sce_filterd <- runPCA(sce_filterd, ncomponents = 25, subset_row = hvg)
sce_filterd <- runUMAP(sce_filterd, dimred = "PCA", external_neighbors = TRUE)


# 聚类 --------------------------------------------------------------------------
g <- buildSNNGraph(sce_filterd, use.dimred = "PCA")
sce_filterd$clusters <- factor(igraph::cluster_louvain(g)$membership)


# 可视化 -------------------------------------------------------------------------
umap_plot <- plotUMAP(sce_filterd, colour_by = "clusters")
pdf("figs/1-2_umap_plot.pdf")
print(umap_plot)
dev.off()
png("figs/1-2_umap_plot.png")
print(umap_plot)
dev.off()



###############################################################################
# hsuh，写于2020.10.08，更新于2020.10.08
###############################################################################
############################################# 2.1 综述 | 单细胞转录组分析最佳思路
