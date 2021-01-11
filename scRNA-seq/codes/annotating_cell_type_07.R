########################################
# annotating cell type
# date: 2021.01.08 - 01.11
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-7
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

#################### seger
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
edb <- AnnotationHub()[["AH73881"]]
symbols <- rowData(sce_seger)$symbol
ens_id <- mapIds(
    edb,
    keys = symbols,
    keytype = "SYMBOL",
    column = "GENEID"
)
ens_id <- ifelse(is.na(ens_id), symbols, ens_id)
# 去重id重复的行
ens_id_keep <- !duplicated(ens_id)
sce_seger_keep <- sce_seger[ens_id_keep, ]
sce_seger_keep
# class: SingleCellExperiment
# dim: 25454 3514
# metadata(0):
# assays(1): counts
# rownames(25454): SGIP1 AZIN2 ... KIR2DS3 eGFP
# rowData names(2): symbol refseq
# colnames(3514): HP1502401_N13 HP1502401_D14 ... HP1526901T2D_O11 HP1526901T2D_A8
# colData names(8): Source Name individual ... age body mass index
# reducedDimNames(0):
# altExpNames(1): ERCC


# sample annotation -----------------------------------------------------------
emtab_meta <- colData(sce_seger_keep)[, c("cell type", "individual", "single cell well quality")]
colnames(emtab_meta) <- c("cell_type", "donor", "quality")
colData(sce_seger_keep) <- emtab_meta

sce_seger_keep$cell_type <- gsub(" cell", "", sce_seger_keep$cell_type)
sce_seger_keep$cell_type <- paste0(
    toupper(substr(sce_seger_keep$cell_type, 1, 1)),
    substring(sce_seger_keep$cell_type, 2)
)


# qc --------------------------------------------------------------------------
low_quality <- sce_seger_keep$quality == "low quality cell"

library(scater)
stats <- perCellQCMetrics(sce_seger_keep)
qc <- quickPerCellQC(
    stats,
    percent_subsets = "altexps_ERCC_percent",
    batch = sce_seger_keep$donor,
    subset = !sce_seger_keep$donor %in% c("HP1504901", "HP1509101")
)
sce_seger_final <- sce_seger_keep[, !(qc$discard | low_quality)]
sce_seger_final
# class: SingleCellExperiment
# dim: 25454 2090
# metadata(0):
# assays(1): counts
# rownames(25454): SGIP1 AZIN2 ... KIR2DS3 eGFP
# rowData names(2): symbol refseq
# colnames(2090): HP1502401_H13 HP1502401_J14 ... HP1526901T2D_N8 HP1526901T2D_A8
# colData names(3): cell_type donor quality
# reducedDimNames(0):
# altExpNames(1): ERCC


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
clust_seger <- quickCluster(sce_seger_final)
sce_seger_final <- computeSumFactors(
    sce_seger_final,
    cluster = clust_seger
)
sce_seger_final <- logNormCounts(sce_seger_final)
sce_seger_final
# class: SingleCellExperiment
# dim: 25454 2090
# metadata(0):
# assays(2): counts logcounts
# rownames(25454): SGIP1 AZIN2 ... KIR2DS3 eGFP
# rowData names(2): symbol refseq
# colnames(2090): HP1502401_H13 HP1502401_J14 ... HP1526901T2D_N8 HP1526901T2D_A8
# colData names(4): cell_type donor quality sizeFactor
# reducedDimNames(0):
# altExpNames(1): ERCC


#################### 参考数据集
##### 使用内置的参考注释数据
BiocManager::install("SingleR")
BiocManager::install("celldex")
library(SingleR)
library(celldex)
ref_data <- celldex::BlueprintEncodeData()
ref_data
# class: SummarizedExperiment
# dim: 19859 259
# metadata(0):
# assays(1): logcounts
# rownames(19859): TSPAN6 TNMD ... LINC00550 GIMAP1-GIMAP5
# rowData names(0):
# colnames(259): mature.neutrophil CD14.positive..CD16.negative.classical.monocyte
#   ... epithelial.cell.of.umbilical.artery.1
#   dermis.lymphatic.vessel.endothelial.cell.1
# colData names(3): label.main label.fine label.ont

# 将自测的细胞在参考数据集中找对应的细胞类型
# 返回的predict是一个数据框，每行是自测数据的一个细胞
predict <- SingleR(
    test = sce_pbmc_final,
    ref = ref_data,
    labels = ref_data$label.main
)
table(predict$labels)
# B-cells CD4+ T-cells CD8+ T-cells           DC  Eosinophils Erythrocytes          HSC
#     549          773         1274            1            1            5           14
#    Monocytes     NK cells
#         1117          251

plotScoreHeatmap(predict)

sum(is.na(predict$pruned.labels))
# [1] 76

# 以下函数绘制三种情况：match（黄色assigned）、mismatch（灰色other）、unclear（橙色pruned）
plotScoreDistribution(predict)

# 查看cluster与label的对应关系
tab <- table(
    assigned = predict$pruned.labels,
    cluster = colLabels(sce_pbmc_final)
)
pheatmap::pheatmap(
    log2(tab + 10),
    color = colorRampPalette(c("white", "blue"))(101)
)

##### 使用自定义的注释数据
sce_muraro
# class: SingleCellExperiment
# dim: 16940 2299
# metadata(0):
# assays(2): counts logcounts
# rownames(16940): ENSG00000268895 ENSG00000121410 ... ENSG00000159840
#   ENSG00000074755
# rowData names(2): symbol chr
# colnames(2299): D28-1_1 D28-1_2 ... D30-8_93 D30-8_94
# colData names(4): label donor plate sizeFactor
# reducedDimNames(0):
# altExpNames(1): ERCC

# 去掉不明确的细胞类型
sce_muraro_filtered <- sce_muraro[, !is.na(sce_muraro$label) & sce_muraro$label != "unclear"]
table(sce_muraro_filtered$label)
#      acinar       alpha        beta       delta        duct endothelial
#         217         795         442         189         239          18
#     epsilon mesenchymal          pp
#           3          80          96

# 以sce_muraro_filtered`作为参考数据集，sce_seger_final作为自测数据进行注释
predict_seger <- SingleR(
    test = sce_seger_final,
    ref = sce_muraro_filtered,
    labels = sce_muraro_filtered$label,
    de.method = "wilcox"
)
table(predict_seger$labels)
#      acinar       alpha        beta       delta        duct endothelial
#         188         889         279         105         385          17
#     epsilon mesenchymal          pp
#           5          53         169

# 将predict_seger$pruned.labels与sce_seger_final$cell_type进行比较
tab <- table(predict_seger$pruned.labels, sce_seger_final$cell_type)
pheatmap::pheatmap(
    log2(tab + 10),
    color = colorRampPalette(c("white", "blue"))(101)
)

##### 使用marker基因集
## 参考数据集, sce_zeisel_filtered
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

table(sce_zeisel_filtered$level1class)
# astrocytes_ependymal    endothelial-mural         interneurons            microglia
#                  179                  160                  290                   78
#     oligodendrocytes        pyramidal CA1         pyramidal SS
#                  774                  938                  397

library(scran)
wilcox_zeisel <- pairwiseWilcox(
    sce_zeisel_filtered,
    sce_zeisel_filtered$level1class,
    lfc = 1,
    direction = "up"
)
markers_zeisel <- getTopMarkers(
    wilcox_zeisel$statistics,
    wilcox_zeisel$pairs,
    pairwise = FALSE,
    n = 50
)
lengths(markers_zeisel)
# astrocytes_ependymal    endothelial-mural         interneurons            microglia
#                   79                   83                  118                   69
#     oligodendrocytes        pyramidal CA1         pyramidal SS
#                   81                  125                  149

## 自测数据集, sce_tasic
library(scRNAseq)
sce_tasic <- TasicBrainData()
sce_tasic
# class: SingleCellExperiment
# dim: 24058 1809
# metadata(0):
# assays(1): counts
# rownames(24058): 0610005C13Rik 0610007C21Rik ... mt_X57780 tdTomato
# rowData names(0):
# colnames(1809): Calb2_tdTpositive_cell_1 Calb2_tdTpositive_cell_2 ...
#   Rbp4_CTX_250ng_2 Trib2_CTX_250ng_1
# colData names(13): sample_title mouse_line ... secondary_type aibs_vignette_id
# reducedDimNames(0):
# altExpNames(1): ERCC

# 使用参考数据集的markers_zeisel，构建一个GeneSetCollection对象
names(markers_zeisel)
# [1] "astrocytes_ependymal" "endothelial-mural"    "interneurons"
# [4] "microglia"            "oligodendrocytes"     "pyramidal CA1"
# [7] "pyramidal SS"

library(GSEABase)
all_sets <- lapply(
    names(markers_zeisel),
    function(x) {
        GeneSet(markers_zeisel[[x]], setName = x)
    }
)
all_sets <- GeneSetCollection(all_sets)
all_sets
# GeneSetCollection
#   names: astrocytes_ependymal, endothelial-mural, ..., pyramidal SS (7 total)
#   unique identifiers: Apoe, Clu, ..., Gabrb2 (555 total)
#   types in collection:
#     geneIdType: NullIdentifier (1 total)
#     collectionType: NullCollection (1 total)

# 对自测数据集sce_tasic的gene表达量进行排序，然后统计得到AUC指标
BiocManager::install("AUCell")
library(AUCell)
# gene表达量排序
gene_rank <- AUCell_buildRankings(
    counts(sce_tasic),
    plotStats = FALSE,
    verbose = FALSE
)
# 计算AUC，然后转置
cell_auc <- AUCell_calcAUC(all_sets, gene_rank)
# Genes in the gene sets NOT available in the dataset:
# 	endothelial-mural: 	7 (8% of 83)
# 	interneurons: 	1 (1% of 118)
# 	oligodendrocytes: 	2 (2% of 81)
# 	pyramidal CA1: 	4 (3% of 125)
# 	pyramidal SS: 	4 (3% of 149)
cell_auc
# AUC for 7 gene sets (rows) and 1809 cells (columns).
#
# Top-left corner of the AUC matrix:
#                       cells
# gene sets              Calb2_tdTpositive_cell_1 Calb2_tdTpositive_cell_2
#   astrocytes_ependymal               0.13865276               0.13662832
#   endothelial-mural                  0.04264310               0.04884635
#   interneurons                       0.53060935               0.45383569
#   microglia                          0.04845394               0.02682648
#   oligodendrocytes                   0.13179577               0.12112934
#   pyramidal CA1                      0.23176680               0.20625697
#                       cells
# gene sets              Calb2_tdTpositive_cell_3 Calb2_tdTpositive_cell_4
#   astrocytes_ependymal               0.10873233               0.13216583
#   endothelial-mural                  0.07269892               0.04993108
#   interneurons                       0.34589983               0.51126651
#   microglia                          0.03583482               0.05387632
#   oligodendrocytes                   0.15672040               0.14808929
#   pyramidal CA1                      0.32187260               0.25467138
#                       cells
# gene sets              Calb2_tdTpositive_cell_5
#   astrocytes_ependymal               0.15128922
#   endothelial-mural                  0.07161420
#   interneurons                       0.49297711
#   microglia                          0.06655747
#   oligodendrocytes                   0.13857658
#   pyramidal CA1                      0.20884775

results <- t(assay(cell_auc))

# 与参考数据中的cluster类型比较
new_label <- colnames(results)[max.col(results)]
table(new_label, sce_tasic$broad_type)
# new_label              Astrocyte Endothelial Cell GABA-ergic Neuron Glutamatergic Neuron
#   astrocytes_ependymal        43                2                 0                    0
#   endothelial-mural            0               27                 0                    0
#   interneurons                 0                0               759                    2
#   microglia                    0                0                 0                    0
#   oligodendrocytes             0                0                 1                    0
#   pyramidal SS                 0                0                 1                  810
#                       
# new_label              Microglia Oligodendrocyte Oligodendrocyte Precursor Cell Unclassified
#   astrocytes_ependymal         0               0                             20            4
#   endothelial-mural            0               0                              0            2
#   interneurons                 0               0                              0           15
#   microglia                   22               0                              0            1
#   oligodendrocytes             0              38                              2            0
#   pyramidal SS                 0               0                              0           60

# 计算AUC最大值减去中位数
library(scater)
library(DelayedArray)
deltas <- rowMaxs(results) - rowMedians(results)
discard <- isOutlier(deltas, type = "lower", batch = new_label)
table(new_label[discard])
# astrocytes_ependymal    endothelial-mural         interneurons     oligodendrocytes
#                   24                    1                    7                   10
#         pyramidal SS
#                   16

# 查看被过滤掉的细胞在整体分布中的位置
par(mar = c(10, 4, 1, 1))
boxplot(split(deltas, new_label), las = 2)
points(attr(discard, "thresholds")[1, ], col = "red", pch = 4, cex = 2)

## 从MsigDb下载gene列表
library(BiocFileCache)
bfc <- BiocFileCache(ask = FALSE)
scsig_path <- bfcrpath(
    bfc,
    file.path("http://software.broadinstitute.org",
        "gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt")
)
scsigs <- getGmt(scsig_path)
scsigs
# GeneSetCollection
#   names: Zheng_Cord_Blood_C1_Putative_Megakaryocyte_Progenitor, Zheng_Cord_Blood_C2_Putative_Basophil_Eosinophil_Mast_cell_Progenitor, ..., Muraro_Pancreas_Endothelial_Cell (256 total)
#   unique identifiers: ABCC3, ABCC4, ..., SEC31A (18656 total)
#   types in collection:
#     geneIdType: NullIdentifier (1 total)
#     collectionType: NullCollection (1 total)

## 使用不同的gene列表作为marker基因集
# 对自测数据集的gene进行排序
muraro_mat <- counts(sce_muraro_filtered)
rownames(muraro_mat) <- rowData(sce_muraro_filtered)$symbol
muraro_rank <- AUCell_buildRankings(
    muraro_mat,
    plotStats = FALSE,
    verbose = FALSE
)

# 将MsigDb的scsigs作为marker基因集用于自测数据
scsig_auc <- AUCell_calcAUC(scsigs, muraro_rank)
scsig_results <- t(assay(scsig_auc))
# 绘图
full_label <- colnames(scsig_results)[max.col(scsig_results)]
full_tab <- table(full_label, sce_muraro_filtered$label)
full_heat <- pheatmap::pheatmap(
    log10(full_tab + 10),
    color = viridis::viridis(100),
    silent = TRUE
)

# 将MsigDb的scsigs的一个子集（pancreas相关）作为marker基因集用于自测数据
scsigs_sub <- scsigs[grep("Pancreas", names(scsigs))]
scsig_sub_auc <- AUCell_calcAUC(scsigs_sub, muraro_rank)
scsig_sub_results <- t(assay(scsig_sub_auc))
# 绘图
sub_label <- colnames(scsig_sub_results)[max.col(scsig_sub_results)]
sub_tab <- table(sub_label, sce_muraro_filtered$label)
sub_heat <- pheatmap::pheatmap(
    log10(sub_tab + 10),
    color = viridis::viridis(100),
    silent = TRUE
)

gridExtra::grid.arrange(full_heat[[4]], sub_heat[[4]])


#################### 基于marker gene的富集分析
# load data -------------------------------------------------------------------
library(scRNAseq)
sce_mam <- BachMammaryData(samples = "G_1")


# gene annotation -------------------------------------------------------------
library(scater)
rownames(sce_mam) <- uniquifyFeatureNames(
    rowData(sce_mam)$Ensembl,
    rowData(sce_mam)$Symbol
)

library(AnnotationHub)
ens_mm_v97 <- AnnotationHub()[["AH73905"]]
rowData(sce_mam)$SEQNAME <- mapIds(
    ens_mm_v97,
    keys = rowData(sce_mam)$Ensembl,
    keytype = "GENEID",
    column = "SEQNAME"
)


# quality control -------------------------------------------------------------
is_mito <- rowData(sce_mam)$SEQNAME == "MT"
stats <- perCellQCMetrics(
    sce_mam,
    subsets = list(Mito = which(is_mito))
)
qc <- quickPerCellQC(stats, percent_subsets = "subsets_Mito_percent")
sce_mam_filtered <- sce_mam[, !qc$discard]


# normalization by deconvolution ----------------------------------------------
library(scran)
set.seed(101000110)
cluster_mam <- quickCluster(sce_mam_filtered)
sce_mam_final <- computeSumFactors(sce_mam_filtered, clusters = cluster_mam)
sce_mam_final <- logNormCounts(sce_mam_final)
sce_mam_final
# class: SingleCellExperiment
# dim: 27998 2772
# metadata(0):
# assays(2): counts logcounts
# rownames(27998): Xkr4 Gm1992 ... Vmn2r122 CAAA01147332.1
# rowData names(3): Ensembl Symbol SEQNAME
# colnames: NULL
# colData names(4): Barcode Sample Condition sizeFactor
# reducedDimNames(0):
# altExpNames(0):


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
set.seed(00010101)
dec_mam_pois <- modelGeneVarByPoisson(sce_mam_final)
top_hvgs_mam <- getTopHVGs(dec_mam_pois, prop = 0.1)
length(top_hvgs_mam)
# [1] 1366


# dimension reduce ------------------------------------------------------------
library(BiocSingular)
set.seed(101010011)
##### PCA
sce_mam_final <- denoisePCA(
    sce_mam_final,
    technical = dec_mam_pois,
    subset.row = top_hvgs_mam
)

##### t-SNE
sce_mam_final <- runTSNE(sce_mam_final, dimred = "PCA")


# clustering ------------------------------------------------------------------
##### graph-based clustering
# 使用PCA的前几个PCs构建SNNG
library(scran)
snn_gr_mam <- buildSNNGraph(sce_mam_final, k = 25, use.dimred = "PCA")
colLabels(sce_mam_final) <- factor(
    igraph::cluster_walktrap(snn_gr_mam)$membership
)


# detecting markers -----------------------------------------------------------
markers_mam <- findMarkers(sce_mam_final, direction = "up", lfc = 1)
markers_mam
# List of length 10
# names(10): 1 2 3 4 5 6 7 8 9 10


# GO analysis for cluster 2 ---------------------------------------------------
chosen_cluster <- "2"
chosen_markers <- markers_mam[[chosen_cluster]]
is_de <- chosen_markers$FDR <= 0.05
summary(is_de)
#    Mode   FALSE    TRUE
# logical   27819     179

# GO analysis
library(org.Mm.eg.db)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
gn <- unique(rownames(chosen_markers)[is_de])
trans <- bitr(
    gn,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
)
go_out <- enrichGO(
    gene = trans$ENTREZID,
    keyType = "ENTREZID",
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)
View(head(summary(go_out), 20))

# 检查Csn2和Csn3的表达情况
plotExpression(
    sce_mam_final,
    features = c("Csn2", "Csn3"),
    x = "label",
    colour_by = "label"
)

#看cluster 2中某个感兴趣通路中的gene表现情况
library(stringr)
adhesion <- go_out["GO:0022408", ]
adhesion_gn <- str_split(adhesion$geneID, "/", simplify = TRUE)
head(chosen_markers[rownames(chosen_markers) %in% adhesion_gn, 1:4], 10)
# DataFrame with 10 rows and 4 columns
#               Top     p.value         FDR summary.logFC
#         <integer>   <numeric>   <numeric>     <numeric>
# Spint2         11 3.28234e-34 1.37163e-31       2.39280
# Epcam          17 8.86978e-94 7.09531e-91       2.32968
# Cebpb          21 6.76957e-16 2.03800e-13       1.80192
# Cd24a          21 3.24195e-33 1.29669e-30       1.72318
# Btn1a1         24 2.16574e-13 6.12488e-11       1.26343
# Cd9            51 1.41373e-11 3.56592e-09       2.73785
# Ceacam1        52 1.66948e-38 7.79034e-36       1.56912
# Sdc4           59 9.15001e-07 1.75467e-04       1.84014
# Anxa1          68 2.58840e-06 4.76777e-04       1.29724
# Cdh1           69 1.73658e-07 3.54897e-05       1.31265
