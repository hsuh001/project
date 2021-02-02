########################################
# practice 7, SMARTer, human pancreatic cells
# date: 2021.02.01 - 02.02
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-7
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
library(scRNAseq)
sce_lawlor <- LawlorPancreasData()

sce_lawlor
# class: SingleCellExperiment
# dim: 26616 638
# metadata(0):
# assays(1): counts
# rownames(26616): ENSG00000229483 ENSG00000232849 ... ENSG00000251576
#   ENSG00000082898
# rowData names(0):
# colnames(638): 10th_C10_S104 10th_C11_S96 ... 9th-C96_S81 9th-C9_S13
# colData names(8): title age ... race Sex
# reducedDimNames(0):
# altExpNames(0):


# gene annotation -------------------------------------------------------------
library(AnnotationHub)
edb <- AnnotationHub(localHub = TRUE)[["AH73881"]]
anno <- select(
    edb,
    keys = rownames(sce_lawlor),
    keytype = "GENEID",
    column = c("SYMBOL", "SEQNAME")
)
rowData(sce_lawlor) <- anno[match(rownames(sce_lawlor), anno[, 1]), -1]
rowData(sce_lawlor)
# DataFrame with 26616 rows and 2 columns
#                      SYMBOL     SEQNAME
#                 <character> <character>
# ENSG00000229483   LINC00362          13
# ENSG00000232849   LINC00363          13
# ENSG00000229558    SACS-AS1          13
# ENSG00000232977   LINC00327          13
# ENSG00000227893   LINC00352          13
# ...                     ...         ...
# ENSG00000232746   LINC02022           3
# ENSG00000150867     PIP4K2A          10
# ENSG00000255021  AC093496.1           3
# ENSG00000251576   LINC01267           3
# ENSG00000082898        XPO1           2


# qc --------------------------------------------------------------------------
# 检查是否含有线粒体gene和批次信息
table(rowData(sce_lawlor)$SEQNAME == "MT")
# FALSE  TRUE
# 25269    13

table(sce_lawlor$`islet unos id`)
#  ACCG268 ACCR015A ACEK420A  ACEL337  ACHY057  ACIB065  ACIW009  ACJV399
#      136       57       45      103       39       57       93      108

# 进行质控
stats <- perCellQCMetrics(
    sce_lawlor,
    subsets = list(Mito = which(rowData(sce_lawlor)$SEQNAME == "MT"))
)
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("subsets_Mito_percent"),
    batch = sce_lawlor$`islet unos id`
)

table(qc$discard)
# FALSE  TRUE
#   604    34

colSums(as.matrix(qc), na.rm = TRUE)
#             low_lib_size            low_n_features high_subsets_Mito_percent
#                        9                         5                        25
#                  discard
#                       34

sce_lawlor_filtered <- sce_lawlor[, !qc$discard]
dim(sce_lawlor_filtered)
# [1] 26616   604

##### 使用qc标准对原数据作图
colData(sce_lawlor) <- cbind(colData(sce_lawlor), stats)
sce_lawlor$discard <- qc$discard

# 使用qc标准对原数据作图
gridExtra::grid.arrange(
    plotColData(sce_lawlor, x = "islet unos id", y = "sum",
        colour_by = "discard") +
    scale_y_log10() +
    ggtitle("Total count") +
    theme(axis.text.x = element_text(angle = 90)),
    plotColData(sce_lawlor, x = "islet unos id", y = "detected",
    colour_by = "discard") +
    scale_y_log10() +
    ggtitle("Detected features") +
    theme(axis.text.x = element_text(angle = 90)),
    plotColData(sce_lawlor, x = "islet unos id", y = "subsets_Mito_percent",
        colour_by = "discard") +
    ggtitle("Mito percent") +
    theme(axis.text.x = element_text(angle = 90)),
    ncol = 3
)

# 查看线粒体含量与文库大小的关系
plotColData(sce_lawlor, x = "sum", y = "subsets_Mito_percent",
    colour_by = "discard") + scale_x_log10() +
ggtitle("Mito percent to total count")


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(1000)
cluster_lawlor <- quickCluster(sce_lawlor_filtered)
sce_lawlor_filtered <- computeSumFactors(
    sce_lawlor_filtered,
    cluster = cluster_lawlor
)
# logNormCounts()
sce_lawlor_filtered <- logNormCounts(sce_lawlor_filtered)
summary(sizeFactors(sce_lawlor_filtered))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.2955  0.7807  0.9633  1.0000  1.1820  2.6287


# measure the degree of change ------------------------------------------------
# and HVGs selection by number
# 没有ERCC、没有UMI，因此使用最基础的方法`modelGeneVar()`，指定批次信息
dec_lawlor <- modelGeneVar(
    sce_lawlor_filtered,
    block = sce_lawlor_filtered$`islet unos id`
)
top_hvgs_lawlor <- getTopHVGs(dec_lawlor, n = 2000)
length(top_hvgs_lawlor)
# [1] 2000


# correcting batch effect -----------------------------------------------------
table(sce_lawlor_filtered$`islet unos id`)
# ACCG268 ACCR015A ACEK420A  ACEL337  ACHY057  ACIB065  ACIW009  ACJV399
#     132       54       37       96       36       55       88      106

# 尝试进行批次效应矫正
library(batchelor)
set.seed(1001010)
sce_lawlor_megered <- fastMNN(
    sce_lawlor_filtered,
    subset.row = top_hvgs_lawlor,
    batch = sce_lawlor_filtered$`islet unos id`
)
sce_lawlor_megered
# class: SingleCellExperiment
# dim: 2000 604
# metadata(2): merge.info pca.info
# assays(1): reconstructed
# rownames(2000): ENSG00000115263 ENSG00000254647 ... ENSG00000164211
#   ENSG00000143384
# rowData names(1): rotation
# colnames(604): 10th_C11_S96 10th_C13_S61 ... 9th-C96_S81 9th-C9_S13
# colData names(1): batch
# reducedDimNames(1): corrected
# altExpNames(0):

# 使用lost.var检查结果
metadata(sce_lawlor_megered)$merge.info$lost.var
#          ACCG268    ACCR015A    ACEK420A     ACEL337     ACHY057     ACIB065    ACIW009
# [1,] 0.019621021 0.026312727 0.000000000 0.000000000 0.000000000 0.000000000 0.00000000
# [2,] 0.033516619 0.038240736 0.066943490 0.000000000 0.000000000 0.000000000 0.00000000
# [3,] 0.021763444 0.015673182 0.012069398 0.086400534 0.000000000 0.000000000 0.00000000
# [4,] 0.017190903 0.014664258 0.008104988 0.022556177 0.053556280 0.000000000 0.00000000
# [5,] 0.025003743 0.024054493 0.011861787 0.027338154 0.006991160 0.109056341 0.00000000
# [6,] 0.008414106 0.005546284 0.005505436 0.005774706 0.003433287 0.006163555 0.17838597
# [7,] 0.015579037 0.016210843 0.008480433 0.017508841 0.005196159 0.014656537 0.01550881
#       ACJV399
# [1,] 0.000000
# [2,] 0.000000
# [3,] 0.000000
# [4,] 0.000000
# [5,] 0.000000
# [6,] 0.000000
# [7,] 0.182398


# dimension reduce ------------------------------------------------------------
library(BiocSingular)
set.seed(101011001)
sce_lawlor_filtered <- runPCA(
    sce_lawlor_filtered,
    subset_row = top_hvgs_lawlor,
    ncomponents = 25
)
sce_lawlor_filtered <- runTSNE(sce_lawlor_filtered, dimred = "PCA")


# clustering, graph-based -----------------------------------------------------
snn_gr_lawlor <- buildSNNGraph(sce_lawlor_filtered, use.dimred = "PCA")
# 鉴定cluster
cluster_lawlor <- igraph::cluster_walktrap(snn_gr_lawlor)$membership
colLabels(sce_lawlor_filtered) <- factor(cluster_lawlor)

table(cluster_lawlor)
# cluster_lawlor
#   1   2   3   4   5   6   7   8
#  34  78 165  26 181  22  75  23

# 查看分群和细胞类型之间的关系
cluster_celltype_tab <- table(
    Cluster = cluster_lawlor,
    CellType = sce_lawlor_filtered$`cell type`
)
cluster_celltype_tab
#        CellType
# Cluster Acinar Alpha Beta Delta Ductal Gamma/PP None/Other Stellate
#       1      1     0    0    13      2       16          2        0
#       2      0     1   76     1      0        0          0        0
#       3      0   161    1     0      0        1          2        0
#       4      0     1    0     1      0        0          5       19
#       5      0     0  175     4      1        0          1        0
#       6     22     0    0     0      0        0          0        0
#       7      0    75    0     0      0        0          0        0
#       8      0     0    0     1     20        0          2        0

library(pheatmap)
pheatmap(
    log10(cluster_celltype_tab + 10),
    color = viridis::viridis(100)
)

# 查看分群和批次之间的关系
cluster_batch_tab <- table(
    Cluster = cluster_lawlor,
    Batch = sce_lawlor_filtered$`islet unos id`
)
cluster_batch_tab
#        Batch
# Cluster ACCG268 ACCR015A ACEK420A ACEL337 ACHY057 ACIB065 ACIW009 ACJV399
#       1       8        2        2       4       4       4       9       1
#       2      14        3        2      33       3       2       4      17
#       3      36       23       14      13      14      14      21      30
#       4       7        1        0       1       0       4       9       4
#       5      34       10        4      39       7      23      24      40
#       6       0        2       13       0       0       0       5       2
#       7      32       12        0       5       6       7       4       9
#       8       1        1        2       1       2       1      12       3

library(pheatmap)
pheatmap(
    log10(cluster_batch_tab + 10),
    color = viridis::viridis(100)
)

# 绘制t-SNE图，查看批次效应的矫正效果
gridExtra::grid.arrange(
    plotTSNE(sce_lawlor_filtered, colour_by = "label"),
    plotTSNE(sce_lawlor_filtered, colour_by = "islet unos id"),
    ncol = 2
)
