########################################
# practice 9, different techs, human pancreatic cells
# date: 2021.02.02 - 02.03
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/04/04-9
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


#################### 整合sce_muraro、sce_grun
# load data -------------------------------------------------------------------
########## CEL-seq, sce_muraro
# load data -------------------------------------------------------------------
library(scRNAseq)
sce_muraro <- MuraroPancreasData()

sce_muraro
# class: SingleCellExperiment
# dim: 19059 3072
# metadata(0):
# assays(1): counts
# rownames(19059): A1BG-AS1__chr19 A1BG__chr19 ... ZZEF1__chr17 ZZZ3__chr1
# rowData names(2): symbol chr
# colnames(3072): D28-1_1 D28-1_2 ... D30-8_95 D30-8_96
# colData names(3): label donor plate
# reducedDimNames(0):
# altExpNames(1): ERCC

# 有4个donor
table(sce_muraro$donor)
# D28 D29 D30 D31
# 768 768 768 768


# gene annotation -------------------------------------------------------------
# 将没有匹配的NA去掉，并且去掉重复的行
head(rownames(sce_muraro))
# [1] "A1BG-AS1__chr19" "A1BG__chr19"     "A1CF__chr10"
# [4] "A2M-AS1__chr12"  "A2ML1__chr12"    "A2M__chr12"

library(AnnotationHub)
edb <- AnnotationHub(localHub = TRUE)[["AH73881"]]
gene_symbols <- sub("__chr.*", "", rownames(sce_muraro))
gene_ids <- mapIds(
    edb,
    keys = gene_symbols,
    keytype = "SYMBOL",
    column = "GENEID"
)

keep_row <- !is.na(gene_ids) & !duplicated(gene_ids)
sum(is.na(gene_ids))
# [1] 2110

sum(duplicated(gene_ids))
# [1] 2118

sce_muraro_filtered <- sce_muraro[keep_row, ]
rownames(sce_muraro_filtered) <- gene_ids[keep_row]
dim(sce_muraro_filtered)
# [1] 16940  3072


# qc --------------------------------------------------------------------------
# 备份数据用于质控探索
unfiltered <- sce_muraro_filtered

# 不含有线粒体gene
table(rowData(sce_muraro_filtered)$chr)
#  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22
#  1710   659   981   894   333   540   530   716   996   241  1186  1116   460   179   392
#  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX  chrY
#   988   656   788   922   815   564   637   613    24

stats <- perCellQCMetrics(sce_muraro_filtered)
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent"),
    batch = sce_muraro_filtered$donor
)

colSums(as.matrix(qc), na.rm = TRUE)
#              low_lib_size            low_n_features high_altexps_ERCC_percent
#                       250                       285                       307
#                   discard
#                       338

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
    batch = sce_muraro_filtered$donor,
    subset = sce_muraro_filtered$donor != "D28"
)
colSums(as.matrix(qc2), na.rm = TRUE)
#             low_lib_size            low_n_features high_altexps_ERCC_percent
#                      663                       700                       738
#                  discard
#                      773

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
sce_muraro_final <- sce_muraro_filtered[, !qc2$discard]
dim(sce_muraro)
# [1] 19059  3072

dim(sce_muraro_filtered)
# [1] 16940  3072

dim(sce_muraro_final)
# [1] 16940  2299


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(1000)
cluster_muraro <- quickCluster(sce_muraro_final)
sce_muraro_final <- computeSumFactors(
    sce_muraro_final,
    cluster = cluster_muraro
)
# logNormCounts()
sce_muraro_final <- logNormCounts(sce_muraro_final)
summary(sizeFactors(sce_muraro_final))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.08782  0.54109  0.82081  1.00000  1.21079 13.98692

sce_muraro_final
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


########## CEL-seq2, sce_grun
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

sce_grun_final
# class: SingleCellExperiment
# dim: 17474 1063
# metadata(0):
# assays(2): counts logcounts
# rownames(17474): ENSG00000268895 ENSG00000121410 ... ENSG00000074755
#   ENSG00000036549
# rowData names(3): symbol chr ensembl
# colnames(1063): D2ex_1 D2ex_2 ... D17TGFB_94 D17TGFB_95
# colData names(3): donor sample sizeFactor
# reducedDimNames(0):
# altExpNames(1): ERCC


# 取两个数据的交集子集 ------------------------------------------------------------------
# 首先获取交集gene
universe <- intersect(
    rownames(sce_muraro_final),
    rownames(sce_grun_final)
)
length(universe)
# [1] 15914

nrow(sce_muraro_final)
# [1] 16940

nrow(sce_grun_final)
# [1] 17474

# 对数据集取子集
sce_muraro_final_sub <- sce_muraro_final[universe, ]
sce_grun_final_sub <- sce_grun_final[universe, ]


# 整合后的归一化normalization -------------------------------------------
# 查看各自原本的文库大小
summary(colSums(logcounts(sce_muraro_final_sub)))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    4425    8018    8663    8538    9172   10566

summary(colSums(logcounts(sce_grun_final_sub)))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    2319    4156    4574    4526    4938    7268

# 使用batchelor::multiBatchNorm()进行处理
library(batchelor)
normed_muraro_grun <- multiBatchNorm(
    sce_muraro_final_sub,
    sce_grun_final_sub
)
sce_muraro_final_sub_2 <- normed_muraro_grun[[1]]
sce_grun_final_sub_2 <- normed_muraro_grun[[2]]
summary(colSums(logcounts(sce_muraro_final_sub_2)))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    3063    4716    4975    4936    5185    5845

summary(colSums(logcounts(sce_grun_final_sub_2)))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    2319    4156    4574    4526    4938    7268


# 整合后的挑选HVGs ------------------------------------------------------------------
# 首先对表达量变化模型取子集
muraro_block <- paste0(sce_muraro_final$plate, "_", sce_muraro_final$donor)
# with spike-in
dec_muraro_spike <- modelGeneVarWithSpikes(
    sce_muraro_final,
    spikes = "ERCC",
    block = muraro_block
)
dec_muraro_spike_sub <- dec_muraro_spike[universe, ]


library(scran)
grun_block <- paste0(sce_grun_final$sample, "_", sce_grun_final$donor)
# with spike-in
dec_grun_spike <- modelGeneVarWithSpikes(
    sce_grun_final,
    spikes = "ERCC",
    block = grun_block
)
dec_grun_spike_sub <- dec_grun_spike[universe, ]

# 使用combineVar()组合两组的结果
library(scran)
dec_combined_pan <- combineVar(dec_muraro_spike_sub, dec_grun_spike_sub)
# 挑选更有可能代表生物差异的gene，用于下游的PCA和聚类
chosen_genes <- rownames(dec_combined_pan)[dec_combined_pan$bio > 0]
length(chosen_genes)
# [1] 14999


# 矫正批次效应 ----------------------------------------------------------------------
# 使用batchelor::rescaleBatches()进行线性回归矫正
library(scater)
sce_rescaled_pan <- rescaleBatches(
    sce_muraro_final_sub_2,
    sce_grun_final_sub_2
)

set.seed(100101)
sce_rescaled_pan <- runPCA(
    sce_rescaled_pan,
    subset_row = chosen_genes,
    exprs_values = "corrected"
)

sce_rescaled_pan <- runTSNE(sce_rescaled_pan, dimred = "PCA")

plotTSNE(sce_rescaled_pan, colour_by = "batch")

# 使用fastMNN()进行矫正
set.seed(1011011)
sce_mnn_pan <- fastMNN(
    sce_muraro_final_sub_2,
    sce_grun_final_sub_2,
    subset.row = chosen_genes
)

snn_gr <- buildSNNGraph(sce_mnn_pan, use.dimred = "corrected")
clusters <- igraph::cluster_walktrap(snn_gr)$membership
tab <- table(
    Cluster = clusters,
    Batch = sce_mnn_pan$batch
)
tab
#        Batch
# Cluster   1   2
#      1  281 244
#      2  419 114
#      3  592 136
#      4  252 306
#      5  195  56
#      6    2  37
#      7  243  63
#      8  108  23
#      9  117  19
#      10  54  60
#      11  17   0
#      12  19   5

sce_mnn_pan <- runTSNE(sce_mnn_pan, dimred = "corrected")

plotTSNE(sce_mnn_pan, colour_by = "batch")


#################### 整合sce_muraro、sce_grun、sce_lawlor、sce_seger
# load data -------------------------------------------------------------------
########## CEL-seq, sce_muraro_final
########## CEL-seq2, sce_grun_final

########## SMARTer, sce_lawlor_final
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

sce_lawlor_filtered
# class: SingleCellExperiment
# dim: 26616 604
# metadata(0):
# assays(2): counts logcounts
# rownames(26616): ENSG00000229483 ENSG00000232849 ... ENSG00000251576
#   ENSG00000082898
# rowData names(2): SYMBOL SEQNAME
# colnames(604): 10th_C11_S96 10th_C13_S61 ... 9th-C96_S81 9th-C9_S13
# colData names(9): title age ... Sex sizeFactor
# reducedDimNames(0):
# altExpNames(0):


########## Smart-seq2, sce_seger
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

sce_seger_final
# class: SingleCellExperiment
# dim: 25454 2090
# metadata(0):
# assays(2): counts logcounts
# rownames(25454): ENSG00000118473 ENSG00000142920 ... ENSG00000278306 eGFP
# rowData names(2): symbol refseq
# colnames(2090): HP1502401_H13 HP1502401_J14 ... HP1526901T2D_N8 HP1526901T2D_A8
# colData names(4): cell_type donor quality sizeFactor
# reducedDimNames(0):
# altExpNames(1): ERCC


# 取四个数据的交集子集 ------------------------------------------------------------------
# 首先获取交集gene
sce_all <- list(
    muraro = sce_muraro_final,
    grun = sce_grun_final,
    lawlor = sce_lawlor_filtered,
    seger = sce_seger_final
)

universe <- Reduce(intersect, lapply(sce_all, rownames))
length(universe)
# [1] 15165

nrow(sce_muraro_final)
# [1] 16940

nrow(sce_grun_final)
# [1] 17474

nrow(sce_lawlor_filtered)
# [1] 26616

nrow(sce_seger_final)
# [1] 25454

# 分别取子集
sce_all_sub <- lapply(sce_all, "[", i = universe, )


# 整合后的归一化normalization -------------------------------------------
# 使用batchelor::multiBatchNorm()进行处理
library(batchelor)
normed_4_pancreas <- do.call(multiBatchNorm, sce_all_sub)


# 整合后的挑选HVGs ------------------------------------------------------------------
# 首先对表达量变化模型取子集
# sce_muraro_final
muraro_block <- paste0(sce_muraro_final$plate, "_", sce_muraro_final$donor)
# with spike-in
dec_muraro_spike <- modelGeneVarWithSpikes(
    sce_muraro_final,
    spikes = "ERCC",
    block = muraro_block
)

# sce_grun_final
library(scran)
grun_block <- paste0(sce_grun_final$sample, "_", sce_grun_final$donor)
# with spike-in
dec_grun_spike <- modelGeneVarWithSpikes(
    sce_grun_final,
    spikes = "ERCC",
    block = grun_block
)

# sce_lawlor_filtered, 没有ERCC、没有UMI，使用modelGeneVar()，指定批次信息
dec_lawlor <- modelGeneVar(
    sce_lawlor_filtered,
    block = sce_lawlor_filtered$`islet unos id`
)

# sce_seger_final, 去除批次AZ + 去除没有ERCC的细胞
keep_col <- librarySizeFactors(altExp(sce_seger_final)) > 0 & sce_seger_final$donor != "AZ"
length(keep_col)
# [1] 2090

sce_seger_final_new <- sce_seger_final[, keep_col]

dec_seger_spike <- modelGeneVarWithSpikes(
    sce_seger_final_new,
    spikes = "ERCC",
    block = sce_seger_final_new$donor
)

dec_all <- list(
    muraro = dec_muraro_spike,
    grun = dec_grun_spike,
    lawlor = dec_lawlor,
    seger = dec_seger_spike
)

# 取表达量变化模型的子集
dec_all_sub <- lapply(dec_all, "[", i = universe, )

# 使用combineVar()组合四组的结果
library(scran)
dec_combined_pan <- do.call(combineVar, dec_all_sub)
# 挑选更有可能代表生物差异的gene，用于下游的PCA和聚类
chosen_genes <- rownames(dec_combined_pan)[dec_combined_pan$bio > 0]
length(chosen_genes)
# [1] 9222


# 矫正批次效应 ----------------------------------------------------------------------
# 使用fastMNN()进行矫正
set.seed(1011011)
sce_mnn_4pan <- fastMNN(normed_4_pancreas)


# clustering, 聚类 --------------------------------------------------------------
snn_gr <- buildSNNGraph(sce_mnn_4pan, use.dimred = "corrected", k = 20)
clusters <- factor(igraph::cluster_walktrap(snn_gr)$membership)
tab <- table(
    Cluster = clusters,
    Batch = sce_mnn_4pan$batch
)
tab
#        Batch
# Cluster grun lawlor muraro seger
#      1   173    255    463   209
#      2   248     22    282   187
#      3     0      0      0   202
#      4   200    237    846   337
#      5   306     27    253   230
#      6     0      0      0   156
#      7     0      0      0   154
#      8    56     17    194   107
#      9    19     18    117   171
#      10   32      0      1     0
#      11   24     19    108    55
#      12    0      0      0    45
#      13    0      0      0   185
#      14    5      8     19    17
#      15    0      1     16     4
#      16    0      0      0    31

# 绘制t-SNE图
sce_mnn_4pan <- runTSNE(sce_mnn_4pan, dimred = "corrected")

gridExtra::grid.arrange(
    plotTSNE(sce_mnn_4pan, colour_by = "batch", text_by = I(clusters)),
    plotTSNE(sce_mnn_4pan, colour_by = I(clusters), text_by = I(clusters)),
    ncol = 2
)


# 对来自供体种类的batch effect的检查 -----------------------------------------------------
# 查看各个数据的供体信息
table(normed_4_pancreas$muraro$donor)
# D28 D29 D30 D31
# 333 601 676 689

table(normed_4_pancreas$grun$donor)
# D10 D17  D2  D3  D7
# 109 434  82 103 335

table(normed_4_pancreas$lawlor$`islet unos id`)
#  ACCG268 ACCR015A ACEK420A  ACEL337  ACHY057  ACIB065  ACIW009  ACJV399
#      132       54       37       96       36       55       88      106

table(normed_4_pancreas$seger$donor)
#           AZ    HP1502401 HP1504101T2D    HP1504901    HP1506401    HP1507101 HP1508501T2D
#           65          249          247          139          230          203          276
#    HP1509101 HP1525301T2D HP1526901T2D
#          109          274          298

# 作图验证
donors_all <- c(
    normed_4_pancreas$muraro$donor,
    normed_4_pancreas$grun$donor,
    normed_4_pancreas$lawlor$`islet unos id`,
    normed_4_pancreas$seger$donor
)

muraro_donors <- donors_all
muraro_donors[sce_mnn_4pan$batch != "muraro"] <- "NA"

grun_donors <- donors_all
grun_donors[sce_mnn_4pan$batch != "grun"] <- "NA"

lawlor_donors <- donors_all
lawlor_donors[sce_mnn_4pan$batch != "lawlor"] <- "NA"

seger_donors <- donors_all
seger_donors[sce_mnn_4pan$batch != "seger"] <- "NA"

gridExtra::grid.arrange(
    plotTSNE(sce_mnn_4pan, colour_by = I(muraro_donors)) +
    ggtitle("muraro_donors"),
    plotTSNE(sce_mnn_4pan, colour_by = I(grun_donors)) +
    ggtitle("grun_donors"),
    plotTSNE(sce_mnn_4pan, colour_by = I(lawlor_donors)) +
    ggtitle("lawlor_donors"),
    plotTSNE(sce_mnn_4pan, colour_by = I(seger_donors)) +
    ggtitle("seger_donors"),
    ncol = 2
)


# correcting batch effect -----------------------------------------------------
# 使用noCorrrect()简单地聚合4个分开的数据集
sce_combined_all <- noCorrect(normed_4_pancreas)
assayNames(sce_combined_all)
# [1] "merged"
# merged重命名为logcounts
assayNames(sce_combined_all) <- "logcounts"
sce_combined_all$donor <- donors_all
sce_combined_all
# class: SingleCellExperiment
# dim: 15165 6056
# metadata(0):
# assays(1): logcounts
# rownames(15165): ENSG00000121410 ENSG00000148584 ... ENSG00000159840
#   ENSG00000074755
# rowData names(0):
# colnames(6056): D28-1_1 D28-1_2 ... HP1526901T2D_N8 HP1526901T2D_A8
# colData names(2): batch donor
# reducedDimNames(0):
# altExpNames(0):

# 划分数据批次、供体批次
donors_per_batch <- split(
    sce_combined_all$donor,
    sce_combined_all$batch
)
donors_per_batch <- lapply(donors_per_batch, unique)
donors_per_batch
# $grun
# [1] "D2"  "D3"  "D7"  "D10" "D17"
# 
# $lawlor
# [1] "ACIW009"  "ACJV399"  "ACCG268"  "ACCR015A" "ACEK420A" "ACEL337"  "ACHY057"  "ACIB065"
# 
# $muraro
# [1] "D28" "D29" "D31" "D30"
# 
# $seger
#  [1] "HP1502401"    "HP1504101T2D" "AZ"           "HP1508501T2D" "HP1506401"
#  [6] "HP1507101"    "HP1509101"    "HP1504901"    "HP1525301T2D" "HP1526901T2D"

# 使用fastMNN()进行矫正
set.seed(1010100)
# batch使用全部的供体
# weights：可以理解为哪些供体属于哪个数据集
sce_multiout <- fastMNN(
    sce_combined_all,
    batch = sce_combined_all$donor,
    subset.row = chosen_genes,
    weights = donors_per_batch
)

# 重新记录数据集、供体
sce_multiout$dataset <- sce_combined_all$batch
sce_multiout$donor <- sce_multiout$batch


# again, dimension reduce -----------------------------------------------------
sce_multiout <- runTSNE(sce_multiout, dimred = "corrected")


# again, clustering, graph-based ----------------------------------------------
snn_gr_multiout <- buildSNNGraph(sce_multiout, use.dimred = 1, k = 20)
# 鉴定cluster
cluster_multiout <- igraph::cluster_walktrap(snn_gr_multiout)$membership

# 查看分群和数据集之间的关系
table(
    Cluster = cluster_multiout,
    Dataset = sce_multiout$dataset
)
#        Dataset
# Cluster grun lawlor muraro seger
#       1  247     20    278   186
#       2  338     28    259   388
#       3  173    256    474   292
#       4  200    239    835   869
#       5   57     17    193   108
#       6   24     17    108    55
#       7    5      9     19    17
#       8    0      1     17     4
#       9   19     17    116   171

# 绘制t-SNE，查看批次效应
gridExtra::grid.arrange(
    plotTSNE(sce_multiout,
        colour_by = "dataset",
        text_by = I(cluster_multiout)
    ),
    plotTSNE(sce_multiout, colour_by = I(seger_donors)),
    ncol = 2
)

# 获取已注释的细胞类型信息
proposed <- c(
    sce_muraro_final$label,
    rep(NA, ncol(sce_grun_final)),
    sce_lawlor_filtered$`cell type`,
    sce_seger_final$cell_type
)

proposed <- tolower(proposed)

# 根据原文修改细胞类型名称
proposed[proposed == "gamma/pp"] <- "gamma"
proposed[proposed == "pp"] <- "gamma"
proposed[proposed == "duct"] <- "ductal"
proposed[proposed == "psc"] <- "stellate"

table(proposed)
proposed
proposed
#                acinar                  alpha                   beta          co-expression
#                   425                   1887                    940                     39
#                 delta                 ductal            endothelial                epsilon
#                   316                    640                     34                      8
#                 gamma            mesenchymal           mhc class ii                   nana
#                   282                     80                      4                     26
#            none/other               stellate           unclassified unclassified endocrine
#                    12                     72                      2                      6
#               unclear
#                     4

table(proposed, cluster_multiout)
#                         cluster_multiout
# proposed                    1    2    3    4    5    6    7    8    9
#   acinar                  421    2    0    1    0    0    1    0    0
#   alpha                     1    8   25 1851    1    0    1    0    0
#   beta                      3    4  919    3    0    2    1    1    7
#   co-expression             0    0   38    1    0    0    0    0    0
#   delta                     0    2    5    2  306    0    1    0    0
#   ductal                    6  613    4    0    0    6    0   10    1
#   endothelial               0    0    0    0    0    1   33    0    0
#   epsilon                   0    0    1    3    0    0    0    0    4
#   gamma                     2    0    0    1    0    0    0    0  279
#   mesenchymal               0    1    0    0    0   79    0    0    0
#   mhc class ii              0    0    0    0    0    0    0    4    0
#   nana                      2    8    2   12    1    0    1    0    0
#   none/other                0    3    1    3    0    0    4    1    0
#   stellate                  0    1    0    0    0   70    1    0    0
#   unclassified              0    0    0    0    0    2    0    0    0
#   unclassified endocrine    0    0    3    3    0    0    0    0    0
#   unclear                   0    4    0    0    0    0    0    0    0
