########################################
# infering trajectory
# date: 2021.01.21 - 01.23
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-12
########################################


# rm all objects --------------------------------------------------------------
rm(list = ls())


# set work directory ----------------------------------------------------------
# work_dir <- "/home1/jxiao/project/scRNA-seq/data/test_data"
work_dir <- "D:/JLab/project/scRNA-seq/data/test_data"
setwd(work_dir)


############### Slingshot
# 数据准备 ------------------------------------------------------------------------
means <- rbind(
    # non-DE genes
    matrix(
        rep(rep(c(0.1, 0.5, 1, 2, 3), each = 300), 100),
        ncol = 300,
        byrow = TRUE
    ),
    # early deactivation
    matrix(
        rep(exp(atan(((300:1) - 200) / 50)), 50),
        ncol = 300,
        byrow = TRUE
    ),
    # late deactivation
    matrix(
        rep(exp(atan(((300:1) - 100) / 50)), 50),
        ncol = 300,
        byrow = TRUE
    ),
    # early activation
    matrix(
        rep(exp(atan(((1:300) - 100) / 50)), 50),
        ncol = 300,
        byrow = TRUE
    ),
    # late activation
    matrix(
        rep(exp(atan(((1:300) - 200) / 50)), 50),
        ncol = 300,
        byrow = TRUE
    ),
    # transient
    matrix(
        rep(exp(atan(c((1:100) / 33, rep(3, 100), (100:1) / 33))), 50),
        ncol = 300,
        byrow = TRUE
    )
)
dim(means)
# [1] 750 300

counts <- apply(means, 2, function(cell_means) {
    total <- rnbinom(1, mu = 7500, size = 4)
    rmultinom(1, total, cell_means)
})
dim(counts)
# [1] 750 300

rownames(counts) <- paste0("G", 1:750)
colnames(counts) <- paste0("c", 1:300)

counts[1:4, 1:4]
#    c1 c2 c3 c4
# G1  0  0  0  4
# G2  4  3  5  5
# G3 17 13  9 12
# G4 26 21 18 26

# 构建一个sce对象
sce_sim <- SingleCellExperiment(assays = List(counts = counts))
sce_sim
# class: SingleCellExperiment
# dim: 750 300
# metadata(0):
# assays(1): counts
# rownames(750): G1 G2 ... G749 G750
# rowData names(0):
# colnames(300): c1 c2 ... c299 c300
# colData names(0):
# reducedDimNames(0):
# altExpNames(0):


# gene过滤 ----------------------------------------------------------------------
gene_filter <- apply(assays(sce_sim)$counts, 1, function(x) {
    sum(x >= 3) >= 10
})
sce_sim_filtered <- sce_sim[gene_filter, ]
# 过滤掉10个gene
dim(sce_sim_filtered)
# [1] 740 300

# 归一化 -------------------------------------------------------------------------
fq_norm <- function(counts) {
    rk <- apply(counts, 2, rank, ties.method = "min")
    counts_sort <- apply(counts, 2, sort)
    refdist <- apply(counts_sort, 1, median)
    norm <- apply(rk, 2, function(r) { refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
}
assays(sce_sim_filtered)$norm <- fq_norm(assays(sce_sim_filtered)$counts)

# 降维 --------------------------------------------------------------------------
## method_1, PCA
pca <- prcomp(t(log1p(assays(sce_sim_filtered)$norm)), scale. = FALSE)
rd1 <- pca$x[, 1:2]
plot(rd1, col = rgb(0, 0, 0, 0.5), pch = 16, asp = 1)

## method_2, diffusion maps
BiocManager::install("destiny")
library(destiny, quietly = TRUE)
dm <- DiffusionMap(t(log1p(assays(sce_sim_filtered)$norm)))
rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)
plot(rd2, col = rgb(0, 0, 0, 0.5), pch = 16, asp = 1)

## 保存两种结果
reducedDims(sce_sim_filtered) <- SimpleList(PCA = rd1, DiffMap = rd2)

# 聚类 --------------------------------------------------------------------------
## method_1, gaussian mixture modeling based on PCA results
install.packages("mclust")
library(mclust, quietly = TRUE)
cl_1 <- Mclust(rd1)$classification
colData(sce_sim_filtered)$GMM <- cl_1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9, "Set1")[cl_1], pch = 16, asp = 1)

## method_2, k-means based on PCA results
cl_2 <- kmeans(rd1, centers = 4)$cluster
colData(sce_sim_filtered)$kmeans <- cl_2

plot(rd1, col = brewer.pal(9, "Set1")[cl_2], pch = 16, asp = 1)

#使用slingshot ------------------------------------------------------------------
BiocManager::install("slingshot")
library(slingshot)
sce_sim_filtered <- slingshot(
    sce_sim_filtered,
    clusterLabels = "GMM",
    reducedDim = "PCA"
)
summary(sce_sim_filtered$slingPseudotime_1)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.000   8.331  20.893  21.130  33.956  43.449

# 对结果可视化
colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)
plot_color <- colors[cut(sce_sim_filtered$slingPseudotime_1, breaks = 100)]
plot(reducedDims(sce_sim_filtered)$PCA, col = plot_color, pch = 16, asp = 1)
lines(SlingshotDataSet(sce_sim_filtered), lwd = 2, col = "black")

############### TSCAN
# 数据准备 ------------------------------------------------------------------------
BiocManager::install("TSCAN")
library(TSCAN)
data(lpsdata)
procdata <- preprocess(lpsdata)

# 构建拟时序 -----------------------------------------------------------------------
lps_clust <- exprmclust(procdata)
# 查看结果
plotmclust(lps_clust)

# 获取排序
lps_order <- TSCANorder(lps_clust)
head(lps_order)
# [1] "Unstimulated_S60" "Unstimulated_S38" "Unstimulated_S28" "Unstimulated_S93"
# [5] "Unstimulated_S66" "Unstimulated_S61"

# 基于找到的排序检测差异gene -------------------------------------------------------------
# 使用difftest()函数
diffval <- difftest(procdata, lps_order)
class(diffval)
# [1] "data.frame"
dim(diffval)
# [1] 534   2
head(diffval)
#                   pval          qval
# STX6      3.679845e-02  3.800845e-02
# MRPL28    3.804289e-68  6.156032e-67
# CUTA      1.120286e-60  1.424363e-59
# AI413582 1.992715e-115 9.673727e-114
# SNRPC     8.639473e-08  1.065468e-07
# MTCH1     6.339235e-37  3.489847e-36

# 根据q值找差异gene
length(rownames(diffval)[diffval$qval < 0.05])
# [1] 517
head(rownames(diffval)[diffval$qval < 0.05])
# [1] "STX6"     "MRPL28"   "CUTA"     "AI413582" "SNRPC"    "MTCH1"

# 对其中某个差异gene作图，以MTCH1为例
mtch1_expr <- log2(lpsdata["MTCH1", ] + 1)
singlegeneplot(
    mtch1_expr,
    TSCANorder(lps_clust, flip = TRUE, orderonly = FALSE)
)

############### monocle2
BiocManager::install("monocle")
library(monocle)
# 数据准备 ------------------------------------------------------------------------
library(scRNAseq)
fluidigm <- ReprocessedFluidigmData()
fluidigm
# class: SingleCellExperiment
# dim: 26255 130
# metadata(3): sample_info clusters which_qc
# assays(4): tophat_counts cufflinks_fpkm rsem_counts rsem_tpm
# rownames(26255): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
# rowData names(0):
# colnames(130): SRR1275356 SRR1274090 ... SRR1275366 SRR1275261
# colData names(28): NREADS NALIGNED ... Cluster1 Cluster2
# reducedDimNames(0):
# altExpNames(0):

# 使用newCellDataSet()创建CellDataSet对象 -------------------------------------------
# ① RSEM表达矩阵（ct = count）
assay(fluidigm) <- assays(fluidigm)$rsem_counts
ct <- floor(assays(fluidigm)$rsem_counts)
ct[1:4, 1:4]
#          SRR1275356 SRR1274090 SRR1275251 SRR1275287
# A1BG              0          0          0          0
# A1BG-AS1          0          0          0          0
# A1CF              0          0          0          0
# A2M               0          0          0         33

# ② 临床信息
sample_ann <- as.data.frame(colData(fluidigm))
dim(sample_ann)
# [1] 130  28

# ③ gene注释信息
gene_ann <- data.frame(
    gene_short_name = rownames(ct),
    row.names = rownames(ct)
)
dim(gene_ann)
# [1] 26255     2

# 转换为AnnotatedDataFrame对象
pd <- new("AnnotatedDataFrame", data = sample_ann)
pd
# An object of class 'AnnotatedDataFrame'
#   rowNames: SRR1275356 SRR1274090 ... SRR1275261 (130 total)
#   varLabels: NREADS NALIGNED ... Cluster2 (28 total)
#   varMetadata: labelDescription

fd <- new("AnnotatedDataFrame", data = gene_ann)
fd
# An object of class 'AnnotatedDataFrame'
#   rowNames: 1 2 ... 26255 (26255 total)
#   varLabels: gene_short_name row_names
#   varMetadata: labelDescription

# 构建CDS对象
sc_cds <- newCellDataSet(
    ct,
    phenoData = pd,
    featureData = fd,
    expressionFamily = negbinomial.size(),
    lowerDetectionLimit = 1
)
sc_cds
# CellDataSet (storageMode: environment)
# assayData: 26255 features, 130 samples
#   element names: exprs
# protocolData: none
# phenoData
#   sampleNames: SRR1275356 SRR1274090 ... SRR1275261 (130 total)
#   varLabels: NREADS NALIGNED ... Size_Factor (29 total)
#   varMetadata: labelDescription
# featureData
#   featureNames: A1BG A1BG-AS1 ... ZZZ3 (26255 total)
#   fvarLabels: gene_short_name
#   fvarMetadata: labelDescription
# experimentData: use 'experimentData(object)'
# Annotation:

# 使用detectGenes()进行质控过滤 -------------------------------------------------------
cds <- sc_cds
# 设置gene表达量的过滤阈值，结果会在cds@featureData@data中新增一列num_cells_expressed
# 记录该gene在多少细胞中表达
cds <- detectGenes(cds, min_expr = 0.1)
head(cds@featureData@data)
#          gene_short_name num_cells_expressed
# A1BG                A1BG                  10
# A1BG-AS1        A1BG-AS1                   2
# A1CF                A1CF                   1
# A2M                  A2M                  21
# A2M-AS1          A2M-AS1                   3
# A2ML1              A2ML1                   9

# 使用subset()进行gene过滤，num_cells_expressed >= 5
expressed_genes <- row.names(
    subset(
        cds@featureData@data,
        num_cells_expressed >= 5
    )
)
length(expressed_genes)
cds <- cds[expressed_genes, ]
dim(cds)
# Features  Samples
#    13385      130

# 进行细胞层面的过滤
cds@phenoData@data[1:4, 1:4]
#              NREADS NALIGNED  RALIGN TOTAL_DUP
# SRR1275356 10554900  7555880 71.5862   58.4931
# SRR1274090   196162   182494 93.0323   14.5122
# SRR1275251  8524470  5858130 68.7213   65.0428
# SRR1275287  7229920  5891540 81.4884   49.7609

# 比如，查看细胞注释的第一个NREADS信息
tmp <- pData(cds)
fivenum(tmp[, 1])
# [1]    91616   232899   892209  8130850 14477100

# 利用subset()过滤细胞，但这里的数据不会对细胞过滤
valid_cells <- rownames(cds@phenoData@data)
cds <- cds[, valid_cells]
cds
# CellDataSet (storageMode: environment)
# assayData: 13385 features, 130 samples
#   element names: exprs
# protocolData: none
# phenoData
#   sampleNames: SRR1275356 SRR1274090 ... SRR1275261 (130 total)
#   varLabels: NREADS NALIGNED ... num_genes_expressed (30 total)
#   varMetadata: labelDescription
# featureData
#   featureNames: A1BG A2M ... ZZZ3 (13385 total)
#   fvarLabels: gene_short_name num_cells_expressed
#   fvarMetadata: labelDescription
# experimentData: use 'experimentData(object)'
# Annotation:

# 聚类 --------------------------------------------------------------------------
# 不使用marker gene聚类，使用clusterCells()
# step1: dispersionTable()
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
# 挑选有差异的
disp_table <- dispersionTable(cds)
dim(disp_table)
# [1] 13229     4
head(disp_table)
#   gene_id mean_expression dispersion_fit dispersion_empirical
# 1    A1BG      0.36606347       36.79574            16.377992
# 2     A2M     22.98234518       12.04768            42.439210
# 3   A2ML1      0.04807691      203.13183             0.000000
# 4    AAAS     22.38633764       12.05835             9.373876
# 5    AACS     14.94809736       12.26298            14.569791
# 6   AADAT      0.80661151       23.06028            15.831734


# 挑选表达量不太低的
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
dim(unsup_clustering_genes)
# [1] 12911     4
head(unsup_clustering_genes)
#   gene_id mean_expression dispersion_fit dispersion_empirical
# 1    A1BG       0.3660635       36.79574            16.377992
# 2     A2M      22.9823452       12.04768            42.439210
# 4    AAAS      22.3863376       12.05835             9.373876
# 5    AACS      14.9480974       12.26298            14.569791
# 6   AADAT       0.8066115       23.06028            15.831734
# 7   AAGAB       4.2098687       13.83388            12.738739

# 准备聚类gene list
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
# 图中黑点是被标记出来等会儿要进行聚类的gene
plot_ordering_genes(cds)

# step2：plot_pc_variance_explained(), 选主成分
# norm_method = "log"
plot_pc_variance_explained(cds, return_all = FALSE)

# step3：选合适的PCs数量
# 进行降维
cds <- reduceDimension(
    cds,
    max_components = 2,
    num_dim = 6,
    reduction_method = "tSNE",
    verbose = TRUE
)
# 进行聚类
cds <- clusterCells(cds, num_clusters = 4)
# Distance cutoff calculated to 0.627547
plot_cell_clusters(cds, 1, 2, color = "Biological_Condition")
table(cds@phenoData@data$Biological_Condition)
#   GW16   GW21 GW21+3    NPC
#     52     16     32     30

# 使用前16个PCs进行降维
cds <- reduceDimension(
    cds,
    max_components = 2,
    num_dim = 16,
    reduction_method = "tSNE",
    verbose = TRUE
)
# 进行聚类
cds <- clusterCells(cds, num_clusters = 4)
# Distance cutoff calculated to 1.137356
plot_cell_clusters(cds, 1, 2, color = "Biological_Condition")

# 降维（reduceDimension()）时也能排除批次效应等干扰因素
# 去除本来的生物学意义
if(F) {
    cds <- reduceDimension(
        cds,
        max_components = 2,
        num_dim = 6,
        reduction_method = "tSNE",
        residualModelFormulaStr = "~ Biological_Condition + num_genes_expressed",
        verbose = TRUE
    )
    cds <- clusterCells(cds, num_clusters = 4)
    # Distance cutoff calculated to 0.8343124
    plot_cell_clusters(cds, 1, 2, color = "Biological_Condition")
}

# 去除生物学意义以外的效应
cds <- reduceDimension(
    cds,
    max_components = 2,
    num_dim = 6,
    reduction_method = "tSNE",
    residualModelFormulaStr = "~ NREADS + num_genes_expressed",
    verbose = TRUE
)
cds <- clusterCells(cds, num_clusters = 4)
# Distance cutoff calculated to 0.57421
plot_cell_clusters(cds, 1, 2, color = "Biological_Condition")


# 差异分析 ------------------------------------------------------------------------
start <- Sys.time()
diff_test_res <- differentialGeneTest(
    cds,
    fullModelFormulaStr = "~ Biological_Condition"
)
end <- Sys.time()
end - start # 运行时间在几分钟至十几分钟不等
# Time difference of 2.302013 mins

# 获取差异gene
sig_genes <- subset(diff_test_res, qval < 0.1)
head(sig_genes[, c("gene_short_name", "pval", "qval")])
#       gene_short_name         pval         qval
# A1BG             A1BG 4.112065e-04 1.460722e-03
# A2M               A2M 4.251744e-08 4.266086e-07
# AACS             AACS 2.881832e-03 8.275761e-03
# AADAT           AADAT 1.069794e-02 2.621123e-02
# AAGAB           AAGAB 1.156771e-07 1.021331e-06
# AAMP             AAMP 7.626789e-05 3.243869e-04

# 绘图，需要将gene变成character
cg <- as.character(head(sig_genes$gene_short_name))
plot_genes_jitter(
    cds[cg, ],
    grouping = "Biological_Condition",
    color_by = "Biological_Condition",
    nrow = 3,
    ncol = NULL
)

# 自定义根据A1BG的表达量差异和分组信息进行作图
boxplot(
    log10(cds@assayData$exprs["A1BG", ] + 1) ~ cds@phenoData@data$Biological_Condition
)


# 推断发育轨迹 ----------------------------------------------------------------------
# 从差异分析结果中选合适gene
ordering_genes <- rownames(subset(diff_test_res, qval < 0.01))
length(ordering_genes)
# [1] 4773
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

# 降维，默认使用DDRTree方法
cds <- reduceDimension(
    cds,
    max_components = 2,
    method = "DDRTree"
)

# 细胞排序
cds <- orderCells(cds)

# 可视化
plot_cell_trajectory(
    cds,
    color_by = "Biological_Condition"
)

# 绘制gene在不同细胞中的表达量变化
plot_genes_in_pseudotime(
    cds[cg, ],
    color_by = "Biological_Condition"
)
