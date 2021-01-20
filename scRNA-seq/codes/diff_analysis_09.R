########################################
# differential expression and differential abundance analysis
# date: 2021.01.13 - 01.20
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-8
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


# load data -------------------------------------------------------------------
# 使用小鼠E8.5时期的嵌合体胚胎数据
BiocManager::install("MouseGastrulationData")
library(MouseGastrulationData)
sce_chimera <- WTChimeraData(samples = 5:10)
sce_chimera
# class: SingleCellExperiment
# dim: 29453 20935
# metadata(0):
# assays(1): counts
# rownames(29453): ENSMUSG00000051951 ENSMUSG00000089699 ... ENSMUSG00000095742
#   tomato-td
# rowData names(2): ENSEMBL SYMBOL
# colnames(20935): cell_9769 cell_9770 ... cell_30702 cell_30703
# colData names(11): cell barcode ... doub.density sizeFactor
# reducedDimNames(2): pca.corrected.E7.5 pca.corrected.E8.5
# altExpNames(0):


# gene annotation -------------------------------------------------------------
library(scater)
rownames(sce_chimera) <- uniquifyFeatureNames(
    rowData(sce_chimera)$ENSEMBL,
    rowData(sce_chimera)$SYMBOL
)


# quality control -------------------------------------------------------------
drop <- sce_chimera$celltype.mapped %in% c("stripped", "Doublet")
sce_chimera_filtered <- sce_chimera[, !drop]
sce_chimera_filtered
# class: SingleCellExperiment
# dim: 29453 19426
# metadata(0):
# assays(1): counts
# rownames(29453): Xkr4 Gm1992 ... CAAA01147332.1 tomato-td
# rowData names(2): ENSEMBL SYMBOL
# colnames(19426): cell_9769 cell_9770 ... cell_30701 cell_30702
# colData names(11): cell barcode ... doub.density sizeFactor
# reducedDimNames(2): pca.corrected.E7.5 pca.corrected.E8.5
# altExpNames(0):


# normalization -------------------------------------------------------------
sce_chimera_filtered <- logNormCounts(sce_chimera_filtered)


# variance modeling -----------------------------------------------------------
library(scran)
dec_chimera <- modelGeneVar(
    sce_chimera_filtered,
    block = sce_chimera_filtered$sample
)
chosen_hvgs_chimera <- dec_chimera$bio > 0


# merging ---------------------------------------------------------------------
library(batchelor)
set.seed(01001001)
sce_chimera_merged <- correctExperiments(
    sce_chimera_filtered,
    batch = sce_chimera_filtered$sample,
    subset.row = chosen_hvgs_chimera,
    PARAM = FastMnnParam(
        merge.order = list(
            list(1, 3, 5),  # WT (3 replicates)
            list(2, 4, 6)   # td-Tomato (3 replicates)
        )
    )
)
sce_chimera_merged
# class: SingleCellExperiment
# dim: 14699 19426
# metadata(2): merge.info pca.info
# assays(3): reconstructed counts logcounts
# rownames(14699): Xkr4 Rp1 ... Vmn2r122 CAAA01147332.1
# rowData names(3): rotation ENSEMBL SYMBOL
# colnames(19426): cell_9769 cell_9770 ... cell_30701 cell_30702
# colData names(12): batch cell ... doub.density sizeFactor
# reducedDimNames(3): corrected pca.corrected.E7.5 pca.corrected.E8.5
# altExpNames(0):


# clustering ------------------------------------------------------------------
snn_gr_chimera <- buildSNNGraph(sce_chimera_merged, use.dimred = "corrected")
clusters_chimera <- igraph::cluster_louvain(snn_gr_chimera)
colLabels(sce_chimera_merged) <- factor(clusters_chimera$membership)


# dimensionality reduce -------------------------------------------------------
sce_chimera_merged <- runTSNE(
    sce_chimera_merged,
    dimred = "corrected",
    external_neighbors = TRUE
)
sce_chimera_merged <- runUMAP(
    sce_chimera_merged,
    dimred = "corrected",
    external_neighbors = TRUE
)


# 初步认识数据 ----------------------------------------------------------------------
sce_chimera_merged
# class: SingleCellExperiment
# dim: 14699 19426
# metadata(2): merge.info pca.info
# assays(3): reconstructed counts logcounts
# rownames(14699): Xkr4 Rp1 ... Vmn2r122 CAAA01147332.1
# rowData names(3): rotation ENSEMBL SYMBOL
# colnames(19426): cell_9769 cell_9770 ... cell_30701 cell_30702
# colData names(13): batch cell ... sizeFactor label
# reducedDimNames(5): corrected pca.corrected.E7.5 pca.corrected.E8.5 TSNE UMAP
# altExpNames(0):

# 包含的样本信息
names(colData(sce_chimera_merged))
#  [1] "batch"           "cell"            "barcode"         "sample"
#  [5] "stage"           "tomato"          "pool"            "stage.mapped"   
#  [9] "celltype.mapped" "closest.cell"    "doub.density"    "sizeFactor"     
# [13] "label"

# 查看三个批次
table(sce_chimera_merged$pool)
#     3     4     5
#  3324  5644 10458

# 查看两种处理
table(sce_chimera_merged$tomato)
# FALSE  TRUE 
# 10331  9095

# 查看分群结果
table(colLabels(sce_chimera_merged))
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18
#  947  112  868  680  606  507 2208  424 1259  252  104  727   58  423 1044  872  432 1264
#   19   20   21   22   23   24   25   26
#  454 1022  211  160  156 1640  860 2136

# 将三个批次和分群结果放在一起看
table(colLabels(sce_chimera_merged), sce_chimera_merged$pool)

# 将两个处理和分群结果放在一起看
table(colLabels(sce_chimera_merged), sce_chimera_merged$tomato)

# 绘制t-SNE图查看数据的批次效应处理效果
gridExtra::grid.arrange(
    plotTSNE(sce_chimera_merged, colour_by = "tomato", text_by = "label"),
    plotTSNE(sce_chimera_merged, colour_by = data.frame(
        pool = factor(sce_chimera_merged$pool))
    ),
    ncol = 2
)

# 绘制热图查看使用已有研究的细胞类型注释数据在自己数据的效果
by_label <- table(
    colLabels(sce_chimera_merged),
    sce_chimera_merged$celltype.mapped
)
pheatmap::pheatmap(
    log2(by_label + 1),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    color = viridis::viridis(101)
)


# 不同处理间差异基因表达分析 DE analysis ---------------------------------------------------
##### 构造一个“拟bulk转录组”的样本
colData(sce_chimera_merged)[, c("celltype.mapped", "sample")]
# DataFrame with 19426 rows and 2 columns
#                   celltype.mapped    sample
#                       <character> <integer>
# cell_9769              Mesenchyme         5
# cell_9770             Endothelium         5
# cell_9771               Allantois         5
# cell_9772              Erythroid3         5
# cell_9773              Erythroid1         5
# ...                           ...       ...
# cell_30698             Erythroid2        10
# cell_30699             Erythroid3        10
# cell_30700       Surface ectoderm        10
# cell_30701 Forebrain/Midbrain/H..        10
# cell_30702             Erythroid3        10

# 根据细胞类型注释信息和分群信息，将各自的组内细胞表达量分别加起来
sce_chimera_summed <- aggregateAcrossCells(
    sce_chimera_merged,
    id = colData(sce_chimera_merged)[, c("celltype.mapped", "sample")]
)

dim(sce_chimera_merged)
# [1] 14699 19426

dim(sce_chimera_summed)
# [1] 14699   186

##### 进行差异分析
### bulk转录组的edgeR差异分析
options(stringsAsFactors = FALSE)

load("airway.exprSet.Rdata")
library(edgeR)

# step1：创建edgeR对象，需要表达矩阵和分组信息
e <- DGEList(counts = expr, group = factor(grp))

# step2：预处理，过滤 + 校正
keep <- rowSums(cpm(e) > 1) >= 2
e_keep <- e[keep, , keep.lib.sizes = FALSE]

e_keep$samples$lib.size <- colSums(e_keep$counts)
e_keep <- calcNormFactors(e_keep)
e_keep$samples

# step3：创建模型分析，得到分组矩阵（实验设计矩阵）
DEG <- e_keep
design_mat <- model.matrix(~ 0 + factor(grp))
rownames(design_mat) <- colnames(DEG)
colnames(design_mat) <- levels(factor(grp))

DEG <- estimateGLMCommonDisp(DEG, design_mat)
DEG <- estimateGLMTrendedDisp(DEG, design_mat)
DEG <- estimateGLMTagwiseDisp(DEG, design_mat)

fit <- glmFit(DEG, design_mat)
# edgeR User guide (Page. 30 => "GLM Approach")
lrt <- glmLRT(fit, contrast = c(-1, 1)) # accoding to design to modify
# or we can use another way (glmQLFTest):
# CasevsCtrl <- makeContrasts(Case-Ctrl=trt-untrt, levels=design)
# lrt <- glmQLFTest(fit,contrast=CasevsCtrl)

nrDEG <- topTags(lrt, n = nrow(DEG))
nrDEG <- as.data.frame(nrDEG)

### scRNA-seq数据的edgeR差异分析
# 指定的细胞类型是间质细胞
label <- "Mesenchyme"
sce_chimera_current <- sce_chimera_summed[, label == sce_chimera_summed$celltype.mapped]
dim(sce_chimera_current)
# [1] 14699     6
# 从原来的186组样本被再次过滤为6组样本

## step1：创建edgeR对象
library(edgeR)
deg_chimera_current <- DGEList(
    counts(sce_chimera_current),
    samples = colData(sce_chimera_current)
)
# deg_chimera_current
# An object of class "DGEList"
# $counts
#        Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
# Xkr4         2       0       0       0       3       0
# Rp1          0       0       1       0       0       0
# Sox17        7       0       3       0      14       9
# Mrpl15    1420     271    1009     379    1578     749
# Rgs20        3       0       1       1       0       0
# 14694 more rows ...
# 
# $samples
#         group lib.size norm.factors batch cell barcode sample stage tomato pool stage.mapped
# Sample1     1  4607053            1     5 <NA>    <NA>      5  E8.5   TRUE    3         <NA>
# Sample2     1  1064970            1     6 <NA>    <NA>      6  E8.5  FALSE    3         <NA>
# Sample3     1  2494010            1     7 <NA>    <NA>      7  E8.5   TRUE    4         <NA>
# Sample4     1  1028668            1     8 <NA>    <NA>      8  E8.5  FALSE    4         <NA>
# Sample5     1  4290221            1     9 <NA>    <NA>      9  E8.5   TRUE    5         <NA>
# Sample6     1  1950840            1    10 <NA>    <NA>     10  E8.5  FALSE    5         <NA>
#         celltype.mapped closest.cell doub.density sizeFactor label celltype.mapped.1
# Sample1      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
# Sample2      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
# Sample3      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
# Sample4      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
# Sample5      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
# Sample6      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
#         sample.1 ncells
# Sample1        5    286
# Sample2        6     55
# Sample3        7    243
# Sample4        8    134
# Sample5        9    478
# Sample6       10    299

## step2：预处理，过滤 + 矫正
# 去除低质量样本
discarded <- sce_chimera_current$ncells < 20
summary(discarded)
#    Mode   FALSE
# logical       6

deg_chimera_current_filtered <- deg_chimera_current[, !discarded]

# 去除低质量gene
keep_gene <- filterByExpr(
    deg_chimera_current_filtered,
    group = sce_chimera_current$tomato
)
summary(keep_gene)
#    Mode   FALSE    TRUE
# logical    9011    5688

deg_chimera_current_final <- deg_chimera_current_filtered[keep_gene, ]
dim(deg_chimera_current_final)
# [1] 5688    6

# 归一化矫正，使用TMM（trimmed mean of M-values）方法
deg_chimera_current_final <- calcNormFactors(deg_chimera_current_final)
deg_chimera_current_final$samples
#         group lib.size norm.factors batch cell barcode sample stage tomato pool stage.mapped
# Sample1     1  4607053    1.0683392     5 <NA>    <NA>      5  E8.5   TRUE    3         <NA>
# Sample2     1  1064970    1.0487418     6 <NA>    <NA>      6  E8.5  FALSE    3         <NA>
# Sample3     1  2494010    0.9582296     7 <NA>    <NA>      7  E8.5   TRUE    4         <NA>
# Sample4     1  1028668    0.9774156     8 <NA>    <NA>      8  E8.5  FALSE    4         <NA>
# Sample5     1  4290221    0.9707300     9 <NA>    <NA>      9  E8.5   TRUE    5         <NA>
# Sample6     1  1950840    0.9816914    10 <NA>    <NA>     10  E8.5  FALSE    5         <NA>
#         celltype.mapped closest.cell doub.density sizeFactor label celltype.mapped.1
# Sample1      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
# Sample2      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
# Sample3      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
# Sample4      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
# Sample5      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
# Sample6      Mesenchyme         <NA>           NA         NA  <NA>        Mesenchyme
#         sample.1 ncells
# Sample1        5    286
# Sample2        6     55
# Sample3        7    243
# Sample4        8    134
# Sample5        9    478
# Sample6       10    299

## step3：构建模型分析，得到分组矩阵（实验设计矩阵）
design_mat <- model.matrix(
    ~ factor(pool) + factor(tomato),
    deg_chimera_current_final$samples
)
design_mat
#         (Intercept) factor(pool)4 factor(pool)5 factor(tomato)TRUE
# Sample1           1             0             0                  1
# Sample2           1             0             0                  0
# Sample3           1             1             0                  1
# Sample4           1             1             0                  0
# Sample5           1             0             1                  1
# Sample6           1             0             1                  0
# attr(,"assign")
# [1] 0 1 1 2
# attr(,"contrasts")
# attr(,"contrasts")$`factor(pool)`
# [1] "contr.treatment"
# 
# attr(,"contrasts")$`factor(tomato)`
# [1] "contr.treatment"

# estimateDisp()估计负二项分布negative binomial（NB）
deg_chimera_current_final <- estimateDisp(
    deg_chimera_current_final,
    design_mat
)

# 使用glmQLFit()估计quasi-likelihood分布
fit <- glmQLFit(deg_chimera_current_final, design_mat, robust = TRUE)

# 使用glmQLFTest()查看表达量差异
res <- glmQLFTest(fit, coef = ncol(design_mat))
summary(decideTests(res))
#        factor(tomato)TRUE
# Down                    8
# NotSig               5672
# Up                      8

# 最后得到一个整合的结果
deg_results <- topTags(res)
deg_results <- as.data.frame(deg_results)
deg_results
#               logFC   logCPM          F       PValue          FDR
# Phlda2   -4.3873546 9.934130 1638.59469 1.812293e-16 1.030832e-12
# Erdr1     2.0690698 8.832662  356.36590 1.060836e-11 3.017019e-08
# Mid1      1.5190728 6.931325  120.14656 1.844351e-08 3.496889e-05
# H13      -1.0596020 7.540121   80.79795 2.373001e-07 2.526790e-04
# Kcnq1ot1  1.3762700 7.241651   83.30701 2.392052e-07 2.526790e-04
# Akr1e1   -1.7205826 5.127802   79.31386 2.665390e-07 2.526790e-04
# Zdbf2     1.8008336 6.797367   83.66324 6.808994e-07 5.532794e-04
# Asb4     -0.9234911 7.340648   53.44578 2.918297e-06 2.074909e-03
# Impact    0.8516300 7.353208   50.31429 4.145416e-06 2.619903e-03
# Lum      -0.6031413 9.274529   41.67104 1.204523e-05 6.851324e-03

##### 循环操作
## step1：过滤细胞，基于"拟bulk转录组"的数据
sce_chimera_summed_filtered <- sce_chimera_summed[, sce_chimera_summed$ncells >= 20]
dim(sce_chimera_summed_filtered)
# [1] 14699   126

## step2：过滤细胞，基于"拟bulk转录组"的数据
# 限制为6个处理，获得每个处理出现第一次的细胞
targets <- colData(sce_chimera_merged)[!duplicated(sce_chimera_merged$sample), ]
dim(targets)
# [1]  6 13

# 再根据这个分组信息构建矩阵
design_mat <- model.matrix(
    ~ factor(pool) + factor(tomato),
    data = targets
)
rownames(design_mat) <- targets$sample
design_mat
#    (Intercept) factor(pool)4 factor(pool)5 factor(tomato)TRUE
# 5            1             0             0                  1
# 6            1             0             0                  0
# 7            1             1             0                  1
# 8            1             1             0                  0
# 9            1             0             1                  1
# 10           1             0             1                  0
# attr(,"assign")
# [1] 0 1 1 2
# attr(,"contrasts")
# attr(,"contrasts")$`factor(pool)`
# [1] "contr.treatment"

## step3：使用pseudoBulkDGE()
library(scran)
de_results <- pseudoBulkDGE(
    sce_chimera_summed_filtered,
    sample = sce_chimera_summed_filtered$sample,
    label = sce_chimera_summed_filtered$celltype.mapped,
    design = design_mat,
    coef = ncol(design_mat),
    condition = targets$tomato
)
# Error in .compute_offsets_by_lfc(design = curdesign, coef = coef, contrast = contrast,  :
#   argument "null.lfc" is missing, with no default
# In addition: Warning message:
# matrix arguments for 'design=' are deprecated.
# Use functions or formulas instead.

# 查看结果
dim(de_results)

# 查看5列的内容
de_results$Allantois[1:3, 1:5]
# DataFrame with 3 rows and 5 columns
#           logFC    logCPM         F    PValue       FDR
#       <numeric> <numeric> <numeric> <numeric> <numeric>
# Xkr4         NA        NA        NA        NA        NA
# Rp1          NA        NA        NA        NA        NA
# Sox17 -0.139844     4.909  0.335911  0.562224  0.891031

##### 查看循环操作的结果
## 查看每个细胞类型中FDR <= 0.05的差异gene数量
is_de <- decideTestsPerLabel(de_results, threshold = 0.05)
dim(is_de)
# [1] 14699    26

# 记录每个gene在每个细胞类型中是否差异表达
is_de[1:3, 1:3]
#       Allantois Blood progenitors 2 Cardiomyocytes
# Xkr4         NA                  NA             NA
# Rp1          NA                  NA             NA
# Sox17         0                  NA             NA

summarizeTestsPerLabel(is_de)
#                                -1    0  1    NA
# Allantois                      69 5048 66  9516
# Blood progenitors 2             1 2472  2 12224
# Cardiomyocytes                  6 4361  5 10327
# Caudal epiblast                 0    0  0 14699
# Caudal Mesoderm                 0    0  0 14699
# Def. endoderm                   0    0  0 14699
# Endothelium                     3 3222  6 11468
# Erythroid1                     12 3035 25 11627
# Erythroid2                      5 3389  8 11297
# Erythroid3                     13 5048 16  9622
# ExE ectoderm                    0    0  0 14699
# ExE mesoderm                    2 5097 10  9590
# Forebrain/Midbrain/Hindbrain    8 6226 11  8454
# Gut                             5 4482  6 10206
# Haematoendothelial progenitors  7 4347 17 10328
# Intermediate mesoderm           6 3256  8 11429
# Mesenchyme                      8 5672  8  9011
# Neural crest                    6 3311  8 11374
# NMP                             6 4107 10 10576
# Paraxial mesoderm               4 4756  5  9934
# Parietal endoderm               0    0  0 14699
# Pharyngeal mesoderm             2 5082  9  9606
# Rostral neurectoderm            0    0  0 14699
# Somitic mesoderm                7 2948 13 11731
# Spinal cord                     7 4591  7 10094
# Surface ectoderm                9 5556  8  9126

## 查看哪些上、下调gene在各种细胞类型中比较多
# 上调
up_de <- is_de > 0 & !is.na(is_de)
head(sort(rowMeans(up_de), decreasing = TRUE), 10)
#   Mid1    Erdr1   Impact    Mcts2     Nnat Kcnq1ot1  Slc38a4    Zdbf2
# 0.7692   0.6538   0.5385   0.5000   0.5000   0.5000   0.3846   0.3462
#   Hopx     Peg3
# 0.3462   0.2308

# 下调
down_de <- is_de < 0 & !is.na(is_de)
head(sort(rowMeans(down_de), decreasing = TRUE), 10)
#  Akr1e1          Xist        Cdkn1c        Phlda2           H13
# 0.61538       0.57692       0.57692       0.57692       0.46154
#   Wfdc2         Hbb-y         Grb10 B930036N10Rik         Pink1
# 0.19231       0.11538       0.11538       0.07692       0.07692

## 查看细胞类型特异的差异gene
# 复制一份差异分析的结果
remotely_de <- decideTestsPerLabel(de_results, threshold = 0.05)
# 挑选非差异的gene
not_de <- remotely_de == 0 | is.na(remotely_de)

# 如果只关注Allantois这个细胞类型，就把剩余细胞类型设为"其他"
other_labels <- setdiff(colnames(not_de), "Allantois")
# 从差异gene中找Allantois的，is.de[, "Allantois"] != 0；
# 同时还需要排除那些也存在于其他细胞类型的gene，即选择非差异gene中其他细胞类型均为TRUE的
# 即，rowMeans(not.de[, other.labels]) == 1
unique_degs <- is.de[, "Allantois"] != 0 & rowMeans(not_de[, other_labels]) == 1
# 提取gene名
unique_degs_name <- names(which(unique_degs))

# 只看Allantois中存在的是14699个，特有的差异gene是44个
length(is_de[, "Allantois"] != 0)
# [1] 14699
length(unique_degs_name)
# [1] 44

## 对细胞类型特异的差异gene排序，例如下面按PValue排序
de_allantois <- de_results$Allantois
de_allantois <- de_allantois[order(de_allantois$PValue), ]
de_allantois <- de_allantois[rownames(de_allantois) %in% unique_degs_name, ]

head(de_allantois)
# DataFrame with 6 rows and 5 columns
#              logFC    logCPM         F      PValue         FDR
#          <numeric> <numeric> <numeric>   <numeric>   <numeric>
# Slc22a18 -4.629906   3.87065  118.8863 2.19531e-27 1.42228e-24
# Eif2s3y   1.757650   5.84561   74.0124 1.01564e-17 4.04928e-15
# Rbp4      1.978988   4.38295   32.7443 1.11016e-08 1.74363e-06
# Cfc1     -0.884383   5.66713   23.7376 1.13690e-06 1.33922e-04
# Pdgfa    -0.423246   8.70475   23.4935 1.28991e-06 1.48569e-04
# H3f3b     0.269312  12.07841   22.8530 1.79717e-06 2.02494e-04

# 提取由于没有重复或对照而不能进行差异分析的细胞类型
metadata(de_results)$failed
# [1] "Caudal epiblast"      "Caudal Mesoderm"      "Def. endoderm"
# [4] "ExE ectoderm"         "Parietal endoderm"    "Rostral neurectoderm"


# 差异细胞丰度分析 DA analysis --------------------------------------------------------
# 将细胞按细胞类型 + 处理两个维度归类
abundances <- table(sce_chimera_merged$celltype.mapped, sce_chimera_merged$sample)

class(abundances)
# [1] "table"
dim(abundances)
# [1] 34  6
head(abundances)
#                         5   6   7   8   9  10
#   Allantois            97  15 139 127 318 259
#   Blood progenitors 1   6   3  16   6   8  17
#   Blood progenitors 2  31   8  28  21  43 114
#   Cardiomyocytes       85  21  79  31 174 211
#   Caudal epiblast       2   2   0   0  22  45
#   Caudal Mesoderm      10  10   9   3  10  29

##### 基于edgeR进行DA analysis
## 为列（即样本） 添加一些信息
extra_info <- colData(sce_chimera_merged)[match(colnames(abundances), sce_chimera_merged$sample), ]
dim(extra_info)
# [1]  6 13

# 创建edgeR对象
deg_chimera_ab <- DGEList(abundances, samples = extra_info)
# deg_chimera_ab
# An object of class "DGEList"
# $counts
# 
#                         5   6   7   8   9  10
#   Allantois            97  15 139 127 318 259
#   Blood progenitors 1   6   3  16   6   8  17
#   Blood progenitors 2  31   8  28  21  43 114
#   Cardiomyocytes       85  21  79  31 174 211
#   Caudal epiblast       2   2   0   0  22  45
# 29 more rows ...
# 
# $samples
#    group lib.size norm.factors batch       cell          barcode sample stage tomato pool
# 5      1     2298            1     5  cell_9769 AAACCTGAGACTGTAA      5  E8.5   TRUE    3
# 6      1     1026            1     6 cell_12180 AAACCTGCAGATGGCA      6  E8.5  FALSE    3
# 7      1     2740            1     7 cell_13227 AAACCTGAGACAAGCC      7  E8.5   TRUE    4
# 8      1     2904            1     8 cell_16234 AAACCTGCAAACCCAT      8  E8.5  FALSE    4
# 9      1     4057            1     9 cell_19332 AAACCTGCAACGATCT      9  E8.5   TRUE    5
# 10     1     6401            1    10 cell_23875 AAACCTGAGGCATGTG     10  E8.5  FALSE    5
#    stage.mapped              celltype.mapped closest.cell doub.density sizeFactor label
# 5         E8.25                   Mesenchyme   cell_24159   0.02985045  1.6348759    19
# 6         E8.25             Somitic mesoderm   cell_63247   0.29191572  2.5980769     6
# 7          E8.5             Somitic mesoderm   cell_25454   0.60173979  1.5939009    17
# 8         E8.25                 ExE mesoderm  cell_139075   0.00473259  0.8707367     9
# 9          E8.0                 ExE mesoderm  cell_116116   0.07941464  0.8932525    15
# 10         E8.5 Forebrain/Midbrain/Hindbrain   cell_39343   0.04074704  0.3947355     1

# 按行根据细胞数量进行过滤
keep_row <- filterByExpr(
    deg_chimera_ab,
    group = deg_chimera_ab$samples$tomato
)
summary(keep_row)
#    Mode   FALSE    TRUE
# logical      10      24

deg_chimera_ab_final <- deg_chimera_ab[keep_row, ]
# Error in subsetListOfArrays(object, i, j, IJ = IJ, IX = IX, I = I, JX = JX) : 
#   argument "j" is missing, with no default
dim(deg_chimera_ab_final)
# [1] 24    6

# 构建模型分析，得到分组矩阵（实验设计矩阵）
design_mat <- model.matrix(
    ~ factor(pool) + factor(tomato),
    deg_chimera_ab_final$samples
)
design_mat

# estimateDisp()估计负二项分布negative binomial（NB）
# 由于细胞数量远少于gene数量，点太少就不需要绘制趋势线
deg_chimera_ab_final <- estimateDisp(
    deg_chimera_ab_final,
    design_mat,
    trend = "none"
)

# 使用glmQLFit()估计quasi-likelihood分布
fit_ab <- glmQLFit(
    deg_chimera_ab_final,
    design_mat,
    robust = TRUE,
    abundance.trend = FALSE
)

# 使用glmQLFTest()查看表达量差异
res_ab <- glmQLFTest(fit_ab, coef = ncol(design_mat))
summary(decideTests(res_ab))
#        factor(tomato)TRUE
# Down                    1
# NotSig                 22
# Up                      1

# 最后得到一个整合的结果
res_ab_results <- topTags(res_ab)
res_ab_results <- as.data.frame(res_ab_results)
res_ab_results
# Coefficient:  factor(tomato)TRUE 
#                                  logFC logCPM      F    PValue       FDR
# ExE ectoderm                   -6.5663  13.02 66.267 1.352e-10 3.245e-09
# Mesenchyme                      1.1652  16.29 11.291 1.535e-03 1.841e-02
# Allantois                       0.8345  15.51  5.312 2.555e-02 1.621e-01
# Cardiomyocytes                  0.8484  14.86  5.204 2.701e-02 1.621e-01
# Neural crest                   -0.7706  14.76  4.106 4.830e-02 2.149e-01
# Endothelium                     0.7519  14.29  3.912 5.371e-02 2.149e-01
# Erythroid3                     -0.6431  17.28  3.604 6.367e-02 2.183e-01
# Haematoendothelial progenitors  0.6581  14.72  3.124 8.351e-02 2.505e-01
# ExE mesoderm                    0.3805  15.68  1.181 2.827e-01 6.258e-01
# Pharyngeal mesoderm             0.3793  15.72  1.169 2.850e-01 6.258e-01

# 关于实验设计矩阵
# pool是哪两个样本属于重复
deg_chimera_ab_final$samples$pool
# [1] 3 3 4 4 5 5

# tomato指处理还是对照
> deg_chimera_ab_final$samples$tomato
# [1]  TRUE FALSE  TRUE FALSE  TRUE FALSE

##### 处理细胞组成的影响
## method_1, 假设大部分细胞类型在各个处理中都存在
# 使用calcNormFactors()计算
deg_chimera_ab_final_1 <- calcNormFactors(deg_chimera_ab_final)
deg_chimera_ab_final_1$samples$norm.factors
# estimateDisp()估计负二项分布negative binomial（NB）
deg_chimera_ab_final_1 <- estimateDisp(
    deg_chimera_ab_final_1,
    design_mat,
    trend = "none"    
)
# 使用glmQLFit()估计quasi-likelihood分布
fit_ab_1 <- glmQLFit(
    deg_chimera_ab_final_1,
    design_mat,
    robust = TRUE,
    abundance.trend = FALSE
)

# 使用glmQLFTest()查看表达量差异
res_ab_1 <- glmQLFTest(fit_ab_1, coef = ncol(design_mat))

# 最后得到一个整合的结果
res_ab_results_1 <- topTags(res_ab_1, n = 10)
res_ab_results_1 <- as.data.frame(res_ab_results_1)
res_ab_results_1
# Coefficient:  factor(tomato)TRUE 
#                                  logFC logCPM      F    PValue       FDR
# ExE ectoderm                   -6.9215  13.17 70.364 5.738e-11 1.377e-09
# Mesenchyme                      0.9513  16.27  6.787 1.219e-02 1.143e-01
# Neural crest                   -1.0032  14.78  6.464 1.429e-02 1.143e-01
# Erythroid3                     -0.8504  17.35  5.517 2.299e-02 1.380e-01
# Cardiomyocytes                  0.6400  14.84  2.735 1.047e-01 4.809e-01
# Allantois                       0.6054  15.51  2.503 1.202e-01 4.809e-01
# Forebrain/Midbrain/Hindbrain   -0.4943  16.55  1.928 1.713e-01 5.178e-01
# Endothelium                     0.5482  14.27  1.917 1.726e-01 5.178e-01
# Erythroid2                     -0.4818  16.00  1.677 2.015e-01 5.373e-01
# Haematoendothelial progenitors  0.4262  14.73  1.185 2.818e-01 6.240e-01

## method_2, 移除细胞数量很多的细胞类型
offenders <- "ExE ectoderm"
deg_chimera_ab_final_2 <- deg_chimera_ab_final[setdiff(rownames(deg_chimera_ab_final), offenders), , keep.lib.sizes = FALSE]
deg_chimera_ab_final_2 <- estimateDisp(
    deg_chimera_ab_final_2,
    design_mat,
    trend = "none"    
)
# 使用glmQLFit()估计quasi-likelihood分布
fit_ab_2 <- glmQLFit(
    deg_chimera_ab_final_2,
    design_mat,
    robust = TRUE,
    abundance.trend = FALSE
)

# 使用glmQLFTest()查看表达量差异
res_ab_2 <- glmQLFTest(fit_ab_2, coef = ncol(design_mat))

# 最后得到一个整合的结果
res_ab_results_2 <- topTags(res_ab_2, n = 10)
res_ab_results_2 <- as.data.frame(res_ab_results_2)
res_ab_results_2
# Coefficient:  factor(tomato)TRUE 
#                                  logFC logCPM      F   PValue     FDR
# Mesenchyme                      1.1274  16.32 11.501 0.001438 0.03308
# Allantois                       0.7950  15.54  5.231 0.026836 0.18284
# Cardiomyocytes                  0.8104  14.90  5.152 0.027956 0.18284
# Neural crest                   -0.8085  14.80  4.903 0.031798 0.18284
# Erythroid3                     -0.6808  17.32  4.387 0.041743 0.19202
# Endothelium                     0.7151  14.32  3.830 0.056443 0.21636
# Haematoendothelial progenitors  0.6189  14.76  2.993 0.090338 0.29683
# Def. endoderm                   0.4911  12.43  1.084 0.303347 0.67818
# ExE mesoderm                    0.3419  15.71  1.036 0.314058 0.67818
# Pharyngeal mesoderm             0.3407  15.76  1.025 0.316623 0.67818
