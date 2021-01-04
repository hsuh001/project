########################################
# selection of highly variable genes
# date: 2020.12.29 - 12.31
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-3
########################################


# rm all objects ---------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)

##### PBMC data
# load packages and PBMC data -------------------------------------------------
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
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
location <- mapIds(
    EnsDb.Hsapiens.v86,
    keys = rowData(sce_pbmc)$ID,
    column = "SEQNAME",
    keytype = "GENEID"
)
head(location)
# ENSG00000243485 ENSG00000237613 ENSG00000186092 ENSG00000238009 
#             "1"             "1"             "1"             "1" 
# ENSG00000239945 ENSG00000239906 
#             "1"             "1"
table(location)
# location
#          1         10         11         12         13         14         15 
#       3125       1243       1902       1774        668       1391       1145 
#         16         17         18         19          2         20         21 
#       1535       1901        661       1997       2284        891        506 
#         22          3          4          5          6          7          8 
#        836       1721       1360       1658       1634       1535       1336 
#          9 GL000009.2 GL000194.1 GL000195.1 GL000205.2 GL000213.1 GL000218.1 
#       1227          1          2          2          1          1          1 
# GL000219.1 KI270711.1 KI270713.1 KI270721.1 KI270726.1 KI270727.1 KI270728.1 
#          1          1          2          1          2          4          6 
# KI270731.1 KI270734.1         MT          X          Y
#          1          3         13       1067        111


# detect dropout -----------------------------------------------------------
set.seed(100)
e_out <- emptyDrops(counts(sce_pbmc))
sce_pbmc_filtered <- sce_pbmc[, which(e_out$FDR <= 0.001)]
dim(sce_pbmc_filtered)
# [1] 33694  4233

# qc, especially for mitochondrial --------------------------------------------
stats <- perCellQCMetrics(
    sce_pbmc_filtered,
    subsets = list(Mito = which(location == "MT"))
)

high_mito <- isOutlier(stats$subsets_Mito_percent, type = "higher")
sce_pbmc_final <- sce_pbmc_filtered[, !high_mito]
dim(sce_pbmc_final)
# [1] 33694  3922


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
set.seed(100)
clust_pbmc <- quickCluster(sce_pbmc_final)
sce_pbmc_final <- computeSumFactors(
    sce_pbmc_final,
    cluster = clust_pbmc
)
# logNormCounts()
sce_pbmc_final <- logNormCounts(sce_pbmc_final)
sce_pbmc_final
# class: SingleCellExperiment
# dim: 33694 3922
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(33694): RP11-34P13.3 FAM138A ... AC213203.1 FAM231B
# rowData names(2): ID Symbol
# colnames(3922): AAACCTGAGAAGGCCT-1 AAACCTGAGACAGACC-1 ...
#   TTTGTCACAGGTCCAC-1 TTTGTCATCCCAAGAT-1
# colData names(2): Sample Barcode
# reducedDimNames(0):
# spikeNames(0):
# altExpNames(0):


# measure the degree of change ------------------------------------------------
# using log-counts, based on mean-variance relationship model
library(scran)
dec_pbmc <- modelGeneVar(sce_pbmc_final)

# 蓝线指所有的gene都存在的一种偏差
fit_pbmc <- metadata(dec_pbmc)
plot(fit_pbmc$mean, fit_pbmc$var,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
curve(fit_pbmc$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
dec_pbmc[order(dec_pbmc$bio, decreasing = TRUE), ]

# using CV
dec_pbmc_cv2 <- modelGeneCV2(sce_pbmc_final)
fit_pbmc_cv2 <- metadata(dec_pbmc_cv2)
plot(fit_pbmc_cv2$mean, fit_pbmc_cv2$cv2, log = "xy")
curve(fit_pbmc_cv2$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
dec_pbmc_cv2[order(dec_pbmc_cv2$ratio, decreasing = TRUE), ]


# without spike-in
set.seed(0010101)
dec_416b_pois <- modelGeneVarByPoisson(sce_pbmc_final)
dec_416b_pois <- dec_416b_pois[order(dec_416b_pois$bio, decreasing = TRUE), ]
head(dec_416b_pois)
# DataFrame with 6 rows and 6 columns
#                     mean            total              tech              bio
#                <numeric>        <numeric>         <numeric>        <numeric>
# LYZ     1.97769673707568  5.1159500611095 0.621547116723033 4.49440294438646
# S100A9  1.94950920785793 4.58859014766227 0.627306158817932 3.96128398884434
# S100A8  1.71827632030156 4.45723403089262 0.669427971460391 3.78780605943223
# HLA-DRA 2.09694376970247  3.7268983984442 0.596371862311291 3.13052653613291
# CD74    2.89839753243019 3.30912153093723 0.422623998345508 2.88649753259172
# CST3    1.49284939321847 2.97369478511548 0.695366743996103 2.27832804111937
#           p.value       FDR
#         <numeric> <numeric>
# LYZ             0         0
# S100A9          0         0
# S100A8          0         0
# HLA-DRA         0         0
# CD74            0         0
# CST3            0         0

# plot
plot(dec_416b_pois$mean, dec_416b_pois$total, pch = 16,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
curve(metadata(dec_416b_pois)$trend(x), col = "dodgerblue", add = TRUE)


# select the HVGs -------------------------------------------------------------
# based on the maximum variance
# for modelGeneVar() and modelGeneVarWithSpikes()
hvg_pbmc_var <- getTopHVGs(dec_pbmc, n = 1000)
str(hvg_pbmc_var)
#  chr [1:1000] "LYZ" "S100A9" "S100A8" "HLA-DRA" "CD74" "CST3" "TYROBP" ...

# for modelGeneCV2() and modelGeneCV2WithSpikes()
hvg_pbmc_cv2 <- getTopHVGs(dec_pbmc_cv2, var.field = "ratio", n = 1000)
str(hvg_pbmc_cv2)
# chr [1:1000] "HIST1H2AC" "GNG11" "PRTFDC1" "TNNC2" "PF4" "HGD" "PPBP" ...

# based on hypothesis testing
hvg_pbmc_var_2 <- getTopHVGs(dec_pbmc, fdr.threshold = 0.05)
length(hvg_pbmc_var_2)
# [1] 651

# keeping genes above the fitted curve
# for modelGeneVar()
hvg_pbmc_var_3 <- getTopHVGs(dec_pbmc, var.threshold = 0)
length(hvg_pbmc_var_3)
# [1] 12791

# for modelGeneCV2()
hvg_pbmc_cv2_3 <- getTopHVGs(dec_pbmc_cv2,
    var.field = "ratio", var.threshold = 1)
length(hvg_pbmc_cv2_3)
# [1] 9295


# deal with genes other than HVGs ---------------------------------------------
# 选取前10%的HVGs
dec_pbmc <- modelGeneVar(sce_pbmc_final)
chosen_hvgs <- getTopHVGs(dec_pbmc, prop = 0.1)
str(chosen_hvgs)
# chr [1:1279] "LYZ" "S100A9" "S100A8" "HLA-DRA" "CD74" "CST3" "TYROBP" ...

# method_1, discard the rest genes
sce_pbmc_hvg <- sce_pbmc_final[chosen_hvgs, ]
dim(sce_pbmc_hvg)
# [1] 1279 3922

# method_2, 原数据不改变，下游分析基于小部分HVGs
library(scater)
sce_pbmc_runpca <- runPCA(sce_pbmc_final, subset_row = chosen_hvgs)
reducedDimNames(sce_pbmc_runpca)
# [1] "PCA"

# method_3, 备份原数据，下游分析基于原数据，最后再复原
# 备份原数据
altExp(sce_pbmc_hvg, "original") <- sce_pbmc_final
altExpNames(sce_pbmc_hvg)

# 按方法一得到处理的数据，进行分析。此时无需指定参数subset_row
sce_pbmc_hvg <- runPCA(sce_pbmc_hvg)

# 复原
sce_pbmc_original <- altExp(sce_pbmc_hvg, "original", withColData = TRUE)


# using genes with prior knowledge --------------------------------------------
# method_1, 选取部分先验gene
# 首先选取gene集合
library(msigdbr)
c7_sets <- msigdbr(species = "Homo sapiens", category = "C7")
head(unique(c7_sets$gs_name))
# [1] "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN"
# [2] "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP"
# [3] "GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_DN"
# [4] "GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP"
# [5] "GOLDRATH_NAIVE_VS_MEMORY_CD8_TCELL_DN"
# [6] "GOLDRATH_NAIVE_VS_MEMORY_CD8_TCELL_UP"

# 利用Goldrath基因集区分CD8细胞亚型
cd8_sets <- c7_sets[grep("GOLDRATH", c7_sets$gs_name), ]
cd8_genes <- rowData(sce_pbmc_final)$Symbol %in% cd8_sets$human_gene_symbol
summary(cd8_genes)
#    Mode   FALSE    TRUE
# logical   32869     825

# 利用GSE11924基因集区分T helper细胞亚型
th_sets <- c7_sets[grep("GSE11924", c7_sets$gs_name), ]
th_genes <- rowData(sce_pbmc_final)$Symbol %in% th_sets$human_gene_symbol
summary(th_genes)
#    Mode   FALSE    TRUE
# logical   31785    1909

# 利用GSE11961基因集区分B细胞亚型
b_sets <- c7_sets[grep("GSE11961", c7_sets$gs_name), ]
b_genes <- rowData(sce_pbmc_final)$Symbol %in% b_sets$human_gene_symbol
summary(b_genes)
#    Mode   FALSE    TRUE
# logical   28192    5502

# method_2, 剔除某些先验gene
# 剔除核糖体蛋白gene或线粒体gene
ribo_discard <- grepl("^RP[SL]\\d+", rownames(sce_pbmc_final))
sum(ribo_discard)
# [1] 99

# 另一种更准确的方法
c2_sets <- msigdbr(species = "Homo sapiens", category = "C2")
ribo_set <- c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ]$human_gene_symbol
ribo_discard <- rownames(sce_pbmc_final) %in% ribo_set
sum(ribo_discard)
# [1] 87

# 对于免疫细胞数据集，可以剔除免疫球蛋白gene和T细胞受体gene
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
anno <- select(
    edb,
    keys = rowData(sce_pbmc_final)$ID,
    keytype = "GENEID",
    columns = "TXBIOTYPE"
)

# 去除免疫球蛋白gene
igv_set <- anno$GENEID[anno$TXBIOTYPE %in% c("IG_V_gene", "IG_V_pseudogene")]
igv_discard <- rowData(sce_pbmc_final)$ID %in% igv_set
sum(igv_discard)
## [1] 326

# 去除T细胞受体gene
tcr_set <- anno$GENEID[anno$TXBIOTYPE %in% c("TR_V_gene", "TR_V_pseudogene")]
tcr_discard <- rowData(sce_pbmc_final)$ID %in% tcr_set
sum(tcr_discard)
## [1] 138


##### 416B data
# load packages and 416B data -------------------------------------------------
# 数据下载
library(scRNAseq)
sce_416b <- LunSpikeInData(which = "416b")
sce_416b$block <- factor(sce_416b$block)
dim(sce_416b)
# [1] 46604   192


# gene annotation -------------------------------------------------------------
BiocManager::install("AnnotationHub")
library(AnnotationHub)
ens_mm_v97 <- AnnotationHub()[["AH73905"]]
rowData(sce_416b)$ENSEMBL <- rownames(sce_416b)
rowData(sce_416b)$SYMBOL <- mapIds(
    ens_mm_v97,
    keys = rownames(sce_416b),
    keytype = "GENEID",
    column = "SYMBOL"
)
rowData(sce_416b)$SEQNAME <- mapIds(
    ens_mm_v97,
    keys = rownames(sce_416b),
    keytype = "GENEID",
    column = "SEQNAME"
)

library(scater)
rownames(sce_416b) <- uniquifyFeatureNames(
    rowData(sce_416b)$ENSEMBL,
    rowData(sce_416b)$SYMBOL
)


# qc, mitochondrial, ERCC, and batch ------------------------------------------
is_mito <- which(rowData(sce_416b)$SEQNAME == "MT")
stats <- perCellQCMetrics(
    sce_416b,
    subsets = list(Mt = is_mito)
)
qc <- quickPerCellQC(
    stats,
    percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent"),
    batch = sce_416b$block
)
sce_416b_filterd <- sce_416b[, !qc$discard]
dim(sce_416b_filterd)
# [1] 46604   185


# 归一化normalization by deconvolution -------------------------------------------
library(scran)
# 未使用quichCluster()
sce_416b_filterd <- computeSumFactors(sce_416b_filterd)
# logNormCounts()
sce_416b_filterd <- logNormCounts(sce_416b_filterd)
sce_416b_filterd
# class: SingleCellExperiment
# dim: 46604 185
# metadata(0):
# assays(2): counts logcounts
# rownames(46604): 4933401J01Rik Gm26206 ... CAAA01147332.1
#   CBFB-MYH11-mcherry
# rowData names(4): Length ENSEMBL SYMBOL SEQNAME
# colnames(185): SLX-9555.N701_S502.C89V9ANXX.s_1.r_1
#   SLX-9555.N701_S503.C89V9ANXX.s_1.r_1 ...
#   SLX-11312.N712_S507.H5H5YBBXX.s_8.r_1
#   SLX-11312.N712_S517.H5H5YBBXX.s_8.r_1
# colData names(9): Source Name cell line ... spike-in addition block
# reducedDimNames(0):
# spikeNames(0):
# altExpNames(2): ERCC SIRV


# measure the degree of change ------------------------------------------------
# using log-counts, based on mean-variance relationship model
library(scran)
# with spike-in
dec_416b_spike <- modelGeneVarWithSpikes(sce_416b_filterd, "ERCC")
dec_416b_spike[order(dec_416b_spike$bio, decreasing = TRUE), ]

# plot
plot(dec_416b_spike$mean, dec_416b_spike$total,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression")
fit_416b_spike <- metadata(dec_416b_spike)
points(fit_416b_spike$mean, fit_416b_spike$var, col = "red", pch = 16)
curve(fit_416b_spike$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

# using log-counts, based on mean-variance relationship model
# adding to batch information
dec_416b_block <- modelGeneVarWithSpikes(
    sce_416b_filterd,
    "ERCC", block = sce_416b_filterd$block
)
head(dec_416b_block[order(dec_416b_block$bio, decreasing = TRUE), 1:6])
# DataFrame with 6 rows and 6 columns
#                     mean            total             tech              bio
#                <numeric>        <numeric>        <numeric>        <numeric>
# Lyz2    6.61235092956153 13.8618988024144 1.58416440831797 12.2777343940965
# Ccl9    6.67841214065115 13.2598761518269 1.44553397941952 11.8143421724073
# Top2a   5.81274666129111 14.0191605357462 2.74571164062005 11.2734488951261
# Cd200r3 4.83305175110888 15.5908569933105 4.31892122081348 11.2719357724971
# Ccnb2   5.97999269432625 13.0256084334992 2.46646680213646 10.5591416313627
# Hbb-bt  4.91682531222784 14.6538670496416 4.12156476433127 10.5323022853104
#                       p.value                   FDR
#                     <numeric>             <numeric>
# Lyz2                        0                     0
# Ccl9                        0                     0
# Top2a   3.89854999267269e-137 8.43398128869382e-135
# Cd200r3  1.17783157586757e-54  7.00721450273016e-53
# Ccnb2   1.20380077734747e-151 2.98404657276435e-149
# Hbb-bt     2.526385935163e-49  1.34197335042576e-47

par(mfrow = c(1, 2))
blocked_state <- dec_416b_block$per.block
for(i in colnames(blocked_state)) {
    current <- blocked_state[[i]]
    plot(current$mean, current$total,
        main = i, pch = 16, cex = 0.5,
        xlab = "Mean of log-expression",
        ylab = "Variance of log-expression"
    )
    cur_fit <- metadata(current)
    points(cur_fit$mean, cur_fit$var, col = "red", pch = 16)
    curve(cur_fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
}


# adding to batch design matrix
table(colData(sce_416b_filterd)$phenotype)
# induced CBFB-MYH11 oncogene expression
#                                     93
#                    wild type phenotype
#                                     92

design <- model.matrix(~factor(block) + phenotype, colData(sce_416b_filterd))
dec_416b_design <- modelGeneVarWithSpikes(
    sce_416b_filterd, "ERCC", design = design
)
dec_416b_design[order(dec_416b_design$bio, decreasing = TRUE), ]
