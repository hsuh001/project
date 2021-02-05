########################################
# practice 14, 10X Genomics, HCA计划的38万骨髓细胞
# date: 2021.02.05 - 02.05
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-14
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
BiocManager::install("HCAData")
library(HCAData)
sce_bone <- HCAData("ica_bone_marrow")
sce_bone
# class: SingleCellExperiment
# dim: 33694 378000
# metadata(0):
# assays(1): counts
# rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#   ENSG00000268674
# rowData names(2): ID Symbol
# colnames(378000): MantonBM1_HiSeq_1-AAACCTGAGCAGGTCA-1
#   MantonBM1_HiSeq_1-AAACCTGCACACTGCG-1 ... MantonBM8_HiSeq_8-TTTGTCATCTGCCAGG-1
#   MantonBM8_HiSeq_8-TTTGTCATCTTGAGAC-1
# colData names(1): Barcode
# reducedDimNames(0):
# altExpNames(0):

# 数据大小为41M
object.size(counts(sce_bone))
# 41741432 bytes

# 实际文件大小为700+M
file.info(path(counts(sce_bone)))$size
# [1] 769046295

# 包含供体信息
head(sce_bone$Barcode)
# [1] "MantonBM1_HiSeq_1-AAACCTGAGCAGGTCA-1" "MantonBM1_HiSeq_1-AAACCTGCACACTGCG-1"
# [3] "MantonBM1_HiSeq_1-AAACCTGCACCGGAAA-1" "MantonBM1_HiSeq_1-AAACCTGCATAGACTC-1"
# [5] "MantonBM1_HiSeq_1-AAACCTGCATCGATGT-1" "MantonBM1_HiSeq_1-AAACCTGCATTCCTGC-1"

# 提取供体信息
sce_bone$donor <- sub("_.*", "", sce_bone$Barcode)
table(sce_bone$donor)
# MantonBM1 MantonBM2 MantonBM3 MantonBM4 MantonBM5 MantonBM6 MantonBM7 MantonBM8
#     48000     48000     48000     48000     48000     42000     48000     48000


# gene annotation -------------------------------------------------------------
library(EnsDb.Hsapiens.v86)
# 获取染色体信息
rowData(sce_bone)$chr <- mapIds(
    EnsDb.Hsapiens.v86,
    keys = rownames(sce_bone),
    keytype = "GENEID",
    column = "SEQNAME"
)

# 线粒体gene个数
table()(grepl("MT", rowData(sce_bone)$chr))
# FALSE  TRUE
# 33681    13

rownames(sce_bone) <- uniquifyFeatureNames(
    rowData(sce_bone)$ID,
    names = rowData(sce_bone)$Symbol
)


# quality control -------------------------------------------------------------
library(BiocParallel)
start <- Sys.time()

is_mito <- rowData(sce_bone)$chr == "MT"
stats <- perCellQCMetrics(
    sce_bone,
    BPPARAM = SnowParam(4),     # 调用4线程
    subsets = list(Mito = which(is_mito))
)
qc <- quickPerCellQC(stats, percent_subsets = "subsets_Mito_percent")

end <- Sys.time()
end - start
# Time difference of 4.841791 mins

colSums(as.matrix(qc), na.rm = TRUE)
#             low_lib_size            low_n_features high_subsets_Mito_percent
#                    35419                     42961                     44937
#                  discard
#                    62149

sce_bone_filtered <- sce_bone[, !qc$discard]
dim(sce_bone_filtered)
# [1]  33694 315851

##### 使用qc标准对原数据作图
colData(sce_bone) <- cbind(colData(sce_bone), stats)
sce_bone$discard <- qc$discard

# 使用qc标准对原数据作图
gridExtra::grid.arrange(
    plotColData(sce_bone, x = "donor", y = "sum",
        colour_by = "discard") +
    scale_y_log10() + ggtitle("Total count"),
    plotColData(sce_bone, x = "donor", y = "detected",
        colour_by = "discard") +
    scale_y_log10() + ggtitle("Detected features"),
    plotColData(sce_bone, x = "donor", y = "subsets_Mito_percent",
        colour_by = "discard") + ggtitle("Mito percent"),
    ncol = 3
)

# 查看线粒体含量与文库大小的关系
plotColData(sce_bone, x = "sum", y = "subsets_Mito_percent",
    colour_by = "discard") +
scale_x_log10() + ggtitle("Mito percent to total count")


# normalization by deconvolution ----------------------------------------------
sce_bone_filtered <- logNormCounts(
    sce_bone_filtered,
    size_factors = sce_bone_filtered$sum
)
summary(sizeFactors(sce_bone_filtered))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.09402  0.47036  0.64737  1.00000  0.88834 42.29495


# measure the degree of change ---------------------------
# and HVGs selection by number
library(scran)
start <- Sys.time()

dec_bone <- modelGeneVar(
    sce_bone_filtered,
    block = sce_bone_filtered$donor,
    BPPARAM = SnowParam(4)     # 调用4线程
)
end <- Sys.time()
end - start
# Time difference of 10.8131 mins

top_hvgs_bone <- getTopHVGs(dec_bone, n = 5000)


# correcting batch effects ----------------------------------------------------
library(batchelor)
library(BiocNeighbors)

set.seed(1010001)
sce_bone_merged <- fastMNN(
    sce_bone_filtered,
    batch = sce_bone_filtered$donor,
    subset.row = top_hvgs_bone,
    BSPARAM = BiocSingular::RandomParam(deferred = TRUE),
    BNPARAM = AnnoyParam(),
    BPPARAM = SnowParam(4)
)

reducedDim(sce_bone_filtered, "MNN") <- reducedDim(
    sce_bone_merged, "corrected"
)

# 查看lost.var值
metadata(sce_bone_merged)$merge.info$lost.var
#        MantonBM1   MantonBM2   MantonBM3   MantonBM4   MantonBM5   MantonBM6   MantonBM7
# [1,] 0.011673815 0.010618772 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
# [2,] 0.008204251 0.008318387 0.025851198 0.000000000 0.000000000 0.000000000 0.000000000
# [3,] 0.005832440 0.004156914 0.005306508 0.026059097 0.000000000 0.000000000 0.000000000
# [4,] 0.003626484 0.003026335 0.002719067 0.003106281 0.029682575 0.000000000 0.000000000
# [5,] 0.005410018 0.003743349 0.003345411 0.005458235 0.006387179 0.035140708 0.000000000
# [6,] 0.005507858 0.005863778 0.005171567 0.004774945 0.005838170 0.006308313 0.038122413
# [7,] 0.002885792 0.002483701 0.003354281 0.003203810 0.002203076 0.003087428 0.002640433
#       MantonBM8
# [1,] 0.00000000
# [2,] 0.00000000
# [3,] 0.00000000
# [4,] 0.00000000
# [5,] 0.00000000
# [6,] 0.00000000
# [7,] 0.04405415


# dimension reduce ------------------------------------------------------------
set.seed(01010100)
##### UMAP
sce_bone_filtered <- runUMAP(
    sce_bone_filtered,
    dimred = "MNN",
    BNPARAM = AnnoyParam(),
    BPPARAM = SnowParam(4),
    n_threads = bpnworkers(SnowParam(4))
)


# clustering ------------------------------------------------------------------
##### graph-based clustering
set.seed(1000)
cluster_bone <- clusterSNNGraph(
    sce_bone_filtered,
    use.dimred = "MNN",
    use.kmeans = TRUE,
    kmeans.centers = 1000
)
colLabels(sce_bone_filtered) <- factor(cluster_bone)

table(cluster_bone)
# cluster_bone
#     1     2     3     4     5     6     7     8     9    10
# 59828 27480 41724 36645 74230 15673  9753 19555 13966 16997

tab <- table(
    Cluster = cluster_bone,
    Donor = sce_bone_filtered$donor
)
pheatmap::pheatmap(
    log10(tab + 10),
    color = viridis::viridis(100)
)


# dectecting markers ----------------------------------------------------------
markers_bone <- findMarkers(
    sce_bone_filtered,
    block = sce_bone_filtered$donor,
    direction = "up",
    lfc = 1,
    BPPARAM = SnowParam(4)
)

# use cluster 1 as an explaination
chosen_cluster <- "1"
markers_cluster_1 <- markers_bone[[chosen_cluster]]
# cluster 1的top 10
interest_markers <- markers_cluster_1[markers_cluster_1$Top <= 10, ]
length(interest_markers)
# [1] 13

# 提取cluster 1与其他clusters对比的logFC结果
lfcs <- getMarkerEffects(interest_markers)

library(pheatmap)
pheatmap(lfcs, breaks = seq(-5, 5, length.out = 101))


# annotating cell type --------------------------------------------------------
sce_bone_aggregated <- sumCountsAcrossCells(
    sce_bone_filtered,
    id = colLabels(sce_bone_filtered)
)
sce_bone_aggregated
# class: SummarizedExperiment
# dim: 33694 10
# metadata(0):
# assays(1): sum
# rownames(33694): RP11-34P13.3 FAM138A ... AC213203.1 FAM231B
# rowData names(0):
# colnames(10): 1 2 ... 9 10
# colData names(2): ids ncells

library(SingleR)
hpc_ref <- HumanPrimaryCellAtlasData()
hpc_ref
# class: SummarizedExperiment
# dim: 19363 713
# metadata(0):
# assays(1): logcounts
# rownames(19363): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
# rowData names(0):
# colnames(713): GSM112490 GSM112491 ... GSM92233 GSM92234
# colData names(3): label.main label.fine label.ont

# 在参考数据集中找cell对应的细胞类型
predict <- SingleR(
    test = sce_bone_aggregated,
    ref = hpc_ref,
    labels = hpc_ref$label.main,
    assay.type.test = "sum"
)

predict
# DataFrame with 10 rows and 5 columns
#                            scores     first.labels       tuning.scores           labels
#                          <matrix>      <character>         <DataFrame>      <character>
# 1  0.310398:0.650347:0.759637:... Pre-B_cell_CD34- 0.491239: 0.1326720         Monocyte
# 2  0.339986:0.765679:0.671786:...           B_cell 0.712581:-0.2524437           B_cell
# 3  0.325049:0.688650:0.654647:...          NK_cell 0.574427: 0.4547246          T_cells
# 4  0.346716:0.627450:0.600886:...          T_cells 0.641036:-0.0218996          T_cells
# 5  0.326947:0.677199:0.650322:...          T_cells 0.690266: 0.2010081          T_cells
# 6  0.385569:0.642941:0.733746:...              CMP 0.569197: 0.3772705              CMP
# 7  0.379836:0.709464:0.701503:... Pro-B_cell_CD34+ 0.815505: 0.7349840 Pro-B_cell_CD34+
# 8  0.320753:0.685264:0.665057:...          NK_cell 0.826341: 0.7544158          NK_cell
# 9  0.404920:0.605420:0.697091:...              MEP 0.412257: 0.3132766       BM & Prog.
# 10 0.310629:0.770565:0.667568:...           B_cell 0.721718:-0.2705010           B_cell
#       pruned.labels
#         <character>
# 1          Monocyte
# 2            B_cell
# 3                NA
# 4           T_cells
# 5           T_cells
# 6               CMP
# 7  Pro-B_cell_CD34+
# 8           NK_cell
# 9        BM & Prog.
# 10           B_cell
