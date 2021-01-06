########################################
# clustering
# date: 2021.01.05 - 01.06
# author: Jing Xiao
# ref: https://jieandze1314.osca.top/03/03-5
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
dim(reducedDim(sce_pbmc_final, "TSNE"))
# [1] 3985    2

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


# clustering ------------------------------------------------------------------
##### graph-based clustering
# 使用PCA的前几个PCs构建SNNG
library(scran)
g_k10 <- buildSNNGraph(sce_pbmc_final, k = 10, use.dimred = "PCA")

# 鉴定cluster
clust_k10 <- igraph::cluster_walktrap(g_k10)$membership
table(clust_k10)
# clust
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
# 205 508 541  56 374 125  46 432 302 867  47 155 166  61  84  16

# 把cluster信息存成SingleCellExperiment对象的一个因子
library(scater)
colLabels(sce_pbmc_final) <- factor(clust_k10)
plotReducedDim(sce_pbmc_final, "TSNE", colour_by = "label")

### set different k values
# k = 50，得到的cluster数量比k = 10少，但每个cluster平均的细胞数更多
g_k50 <- buildSNNGraph(sce_pbmc_final, k = 50, use.dimred = "PCA")
clust_k50 <- igraph::cluster_walktrap(g_k50)$membership
table(clust_k50)
# clust_k50
#   1   2   3   4   5   6   7   8   9  10
# 869 514 194 478 539 944 138 175  89  45

# k = 5，得到的cluster数量比k = 10多
g_k5 <- buildSNNGraph(sce_pbmc_final, k = 5, use.dimred = "PCA")
clust_k5 <- igraph::cluster_walktrap(g_k5)$membership
table(clust_k5)
# clust_k5
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22
# 523 302 125  45 172 573 249 439 293  95 772 142  38  18  62  38  30  16  15   9  16  13

### plot using force-directed
set.seed(2000)
reducedDim(sce_pbmc_final, "force") <- igraph::layout_with_fr(g_k10)
plotReducedDim(sce_pbmc_final, dimred = "force", colour_by = "label")

### other arguments
## type参数设置计算权重的方法
g_num <- buildSNNGraph(sce_pbmc_final, use.dimred = "PCA", type = "number")
g_jaccard <- buildSNNGraph(sce_pbmc_final, use.dimred = "PCA", type = "jaccard")

## 检测community的方法，调用igraph包的各种方法
clust_walktrap <- igraph::cluster_walktrap(g_k10)$membership
clust_louvain <- igraph::cluster_louvain(g_k10)$membership
clust_infomap <- igraph::cluster_infomap(g_k10)$membership
clust_fast <- igraph::cluster_fast_greedy(g_k10)$membership
clust_labprop <- igraph::cluster_label_prop(g_k10)$membership
clust_eigen <- igraph::cluster_leading_eigen(g_k10)$membership

# 不同检测community方法的比较
library(pheatmap)
# 比较infomap和walktrap
# Using a large pseudo-count for a smoother color transition
# between 0 and 1 cell in each 'tab'.
tab_ivw <- table(
    paste("Infomap", clust_infomap),
    paste("Walktrap", clust_walktrap)
)
ivw_plot <- pheatmap(
    log10(tab_ivw + 10),
    main = "Infomap vs Walktrap",
    color = viridis::viridis(100),
    silent = TRUE
)

# Fast-greedy和Walktrap
tab_fvw <- table(
    paste("Fast", clust_fast),
    paste("Walktrap", clust_walktrap)
)
fvw_plot <- pheatmap(
    log10(tab_fvw + 10),
    main = "Fast-greedy vs Walktrap",
    color = viridis::viridis(100),
    silent = TRUE
)
gridExtra::grid.arrange(ivw_plot[[4]], fvw_plot[[4]], ncol = 2)

## 调用cut_at()指定cluster的数量
community_walktrap <- igraph::cluster_walktrap(g_k10)
# 指定5个cluster
table(igraph::cut_at(community_walktrap, n = 5))
# 1    2    3    4    5
# 3546  221  125   46   47

# 指定20个cluster
table(igraph::cut_at(community_walktrap, n = 20))
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
# 462 374 125 437  46 432 302 173 867  47 155 166 104  40  61  84  46  32  16  16


### evaluation of SNNG
library(scran)
g_k10 <- buildSNNGraph(sce_pbmc_final, k = 10, use.dimred = "PCA")
# 鉴定cluster
clust_k10 <- igraph::cluster_walktrap(g_k10)$membership
table(clust_k10)
# clust
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
# 205 508 541  56 374 125  46 432 302 867  47 155 166  61  84  16

ratio <- clusterModularity(g_k10, clust_k10, as.ratio = TRUE)
# Warning message:
# 'clusterModularity' is deprecated.
# Use 'bluster::pairwiseModularity' instead.
# See help("Deprecated")
dim(ratio)
# [1] 16 16
# 是一个矩阵，行 / 列均表示一个cluster，ratio值表示观察到的权重与预期权重之比
library(pheatmap)
pheatmap(
    log2(ratio + 1),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    color = colorRampPalette(c("white", "blue"))(100)
)

## cluster_gr, based on ratio
cluster_gr <- igraph::graph_from_adjacency_matrix(
    log2(ratio + 1),
    mode = "upper",
    weighted = TRUE,
    diag = FALSE
)
set.seed(11001010)
# 权重（weight）越大，线越明显
plot(
    cluster_gr,
    edge.width = igraph::E(cluster_gr)$weight * 5,
    layout = igraph::layout_with_lgl
)

##### k-means
## k = 10
set.seed(100)
clust_kmeans_10 <- kmeans(reducedDim(sce_pbmc_final, "PCA"), centers = 10)
table(clust_kmeans_10$cluster)
#  1   2   3   4   5   6   7   8   9  10
# 548  46 408 270 539 199 148 783 163 881

# 作图
colLabels(sce_pbmc_final) <- factor(clust_kmeans_10$cluster)
kmeans_10_tsne <- plotReducedDim(sce_pbmc_final, "TSNE", colour_by = "label")

## k = 20
set.seed(1001)
clust_kmeans_20 <- kmeans(reducedDim(sce_pbmc_final, "PCA"), centers = 20)
table(clust_kmeans_20$cluster)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
# 159 289 229 277 165 212 287 127 546 134 110 257 538 128  46 131  57  28 220  45

# 作图
colLabels(sce_pbmc_final) <- factor(clust_kmeans_20$cluster)
kmeans_20_tsne <- plotTSNE(
    sce_pbmc_final,
    colour_by = "label",
    text_by = "label"
)

gridExtra::grid.arrange(kmeans_10_tsne, kmeans_20_tsne, ncol = 2)

### evaluation of k-means, using WCSS or RMSD
ncells <- tabulate(clust_kmeans_20$cluster)
tab <- data.frame(
    wcss = clust_kmeans_20$withinss,
    ncells = ncells
)
tab$rmsd <- sqrt(tab$wcss / tab$ncells)
tab
#         wcss ncells      rmsd
# 1   2813.338    159  4.206418
# 2   2093.957    289  2.691751
# 3   3252.230    229  3.768538
# 4   2099.800    277  2.753271
# 5   3711.348    165  4.742681
# 6   3592.732    212  4.116655
# 7   1978.072    287  2.625307
# 8   2891.033    127  4.771167
# 9  13251.813    546  4.926532
# 10  2011.248    134  3.874185
# 11  3595.168    110  5.716935
# 12  6605.020    257  5.069563
# 13 12308.960    538  4.783211
# 14  3635.571    128  5.329437
# 15  6530.806     46 11.915286
# 16  2291.118    131  4.182039
# 17  1378.934     57  4.918519
# 18  2836.975     28 10.065809
# 19  3615.990    220  4.054173
# 20  1204.931     45  5.174576

# hierarchical clustering to check the relationship between clusters
cent_tree <- hclust(dist(clust_kmeans_20$centers), "ward.D2")
plot(cent_tree)

### k-means + graph-based clustering
set.seed(0101010)
# kmeans.centers is for k-means, and k is for graph-based clustering
kgraph_clusters <- clusterSNNGraph(
    sce_pbmc_final,
    use.dimred = "PCA",
    use.kmeans = TRUE,
    kmeans.centers = 1000,
    k = 5
)
table(kgraph_clusters)
# kgraph_clusters
#   1   2   3   4   5   6   7   8   9  10  11  12
# 191 854 506 541 541 892  46 120  29 132  47  86

plotTSNE(sce_pbmc_final, colour_by = I(kgraph_clusters))



##### hierarchical clustering, using 416B data
# load packages and 416B data -------------------------------------------------
# 数据下载
library(scRNAseq)
sce_416b <- LunSpikeInData(which = "416b")
sce_416b$block <- factor(sce_416b$block)
dim(sce_416b)
# [1] 46604   192


# gene annotation -------------------------------------------------------------
#BiocManager::install("AnnotationHub")
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


# measure the degree of change by data distribution ---------------------------
# and HVGs selection by proportion
library(scran)
set.seed(1001)
# with spike-in
dec_416b_spike <- modelGeneVarWithSpikes(sce_416b_filterd, "ERCC")
top_hvgs_416b <- getTopHVGs(dec_416b_spike, prop = 0.1)
length(top_hvgs_416b)
# [1] 1109


# dimension reduce, using three methods ---------------------------------------
##### PCA
set.seed(101010011)
sce_416b_filterd <- denoisePCA(
    sce_416b_filterd,
    subset.row = top_hvgs_416b,
    technical = dec_416b_spike
)
dim(reducedDim(sce_416b_filterd, "PCA"))
# [1] 185    42

##### t-SNE
set.seed(1011101001)
sce_416b_filterd <- runTSNE(sce_416b_filterd, dimred = "PCA")
dim(reducedDim(sce_416b_filterd, "TSNE"))
# [1] 185    2
sce_416b_filterd
# class: SingleCellExperiment
# dim: 46604 185
# metadata(0):
# assays(2): counts logcounts
# rownames(46604): 4933401J01Rik Gm26206 ... CAAA01147332.1 CBFB-MYH11-mcherry
# rowData names(4): Length ENSEMBL SYMBOL SEQNAME
# colnames(185): SLX-9555.N701_S502.C89V9ANXX.s_1.r_1
#   SLX-9555.N701_S503.C89V9ANXX.s_1.r_1 ... SLX-11312.N712_S507.H5H5YBBXX.s_8.r_1
#   SLX-11312.N712_S517.H5H5YBBXX.s_8.r_1
# colData names(10): Source Name cell line ... block sizeFactor
# reducedDimNames(2): PCA TSNE
# altExpNames(2): ERCC SIRV


# hierarchical clustering -----------------------------------------------------
# 利用前几个PCs计算细胞间的距离
dist_416b <- dist(reducedDim(sce_416b_filterd, "PCA"))

# 使用Ward's方法进行聚类
tree_416b <- hclust(dist_416b, "ward.D2")

## 绘制聚类树
install.packages("dendextend")
library(dendextend)
tree_416b$labels <- seq_along(tree_416b$labels)
dend <- as.dendrogram(tree_416b, hang = 0.1)
# 合并两种批次信息
combined_fac <- paste0(
    sce_416b_filterd$block, ".",
    sub(" .*", "", sce_416b_filterd$phenotype)
)
table(combined_fac)
# combined_fac
# 20160113.induced    20160113.wild 20160325.induced    20160325.wild 
#               46               47               47               45

# 按照不同处理设置蓝色、红色，随后按照批次设置淡蓝色、淡红色
labels_colors(dend) <- c(
    `20160113.wild` = "blue",
    `20160113.induced` = "red",
    `20160325.wild` = "dodgerblue",
    `20160325.induced` = "salmon"
)[combined_fac][order.dendrogram(dend)]

plot(dend)

## 对树进行修剪
install.packages("dynamicTreeCut")
library(dynamicTreeCut)
clust_416b <- cutreeDynamic(
    tree_416b,
    distM = as.matrix(dist_416b),
    minClusterSize = 10,
    deepSplit = 1
)
table(clust_416b)
# clust_416b
#  1  2  3  4
# 79 65 28 13

labels_colors <- clust_416b[order.dendrogram(dend)]
plot(dend)

colLabels(sce_416b_filterd) <- factor(clust_416b)
plotReducedDim(sce_416b_filterd, "TSNE", colour_by = "label")

### evaluation of hierarchical clustering, using silhouette width
library(cluster)
sil <- silhouette(clust_416b, dist = dist_416b)
colnames(sil)
# [1] "cluster"   "neighbor"  "sil_width"

# 设置每个cluster的颜色
clust_col <- scater:::.get_palette("tableau10medium")
# 如果sil_width>0，就属于和自己接近，否则属于和其他亚群接近
sil_cols <- clust_col[ifelse(sil[, 3] > 0, sil[, 1], sil[, 2])]
sil_cols <- sil_cols[order(-sil[, 1], sil[, 3])]

plot(
    sil,
    main = paste(length(unique(clust_416b)), "clusters"),
    border = sil_cols,
    col = sil_cols,
    do.col.sort = FALSE
)


##### 评价分群的稳定性
# 例如使用graph-based clustering
my_cluster_fun <- function(x) {
    g <- buildSNNGraph(x, use.dimred = "PCA", type = "jaccard")
    igraph::cluster_louvain(g)$membership
}

originals <- my_cluster_fun(sce_pbmc_final)

# 进行自助判断
set.seed(0010010100)
coassign <- bootstrapCluster(
    sce_pbmc_final,
    FUN = my_cluster_fun,
    clusters = originals)
dim(coassign)
# [1] 16 16

# 绘制热图
pheatmap(
    coassign,
    cluster_row = FALSE,
    cluster_col = FALSE,
    color = rev(viridis::magma(100))
)

##### subclustering
g_pbmc_full <- buildSNNGraph(sce_pbmc_final, use.dimred = "PCA")
clust_pbmc_full <- igraph::cluster_walktrap(g_pbmc_full)$membership

# 绘制几个marker gene的表达量，使用log-normalized过的表达量
plotExpression(
    sce_pbmc_final,
    features = c("CD3E", "CCR7", "CD69", "CD44"),
    x = I(factor(clust_pbmc_full)),
    colour_by = I(factor(clust_pbmc_full))
)

# 推测cluster10为记忆T细胞，获取CD4+、CD8+的亚亚群
memory_t <- 10
sce_pbmc_memory <- sce_pbmc_final[, clust_pbmc_full == memory_t]
# 重新挑一遍marker gene，进行PCA
dec_pbmc_memory <- modelGeneVarByPoisson(sce_pbmc_memory)
sce_pbmc_memory <- denoisePCA(
    sce_pbmc_memory,
    technical = dec_pbmc_memory,
    subset.row = getTopHVGs(dec_pbmc_memory, prop = 0.1)
)
# 重新进行聚类
g_pbmc_memory <- buildSNNGraph(sce_pbmc_memory, use.dimred = "PCA")
clust_pbmc_memory <- igraph::cluster_walktrap(g_pbmc_memory)$membership
plotExpression(
    sce_pbmc_memory,
    features = c("CD8A", "CD4"),
    x = I(factor(clust_pbmc_memory)),
    colour_by = I(factor(clust_pbmc_memory))
)
