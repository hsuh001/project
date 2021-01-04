########################################
# packges for scRNA-seq data analysis
# date: 2021.01.04 - 01.04
# author: Jing Xiao
########################################


install.packages("tidyverse")
library(tidyverse)

install.packages("robustbase")
library(robustbase)

install.packages("BiocManager")
library(BiocManager)

install.packages("pheatmap")
library(pheatmap)

install.packages("remotes")
library(remotes)

install.packages("NMF")
library(NMF)

install.packages("Rtsne")
library(Rtsne)

install.packages("uwot")
library(uwot)

bio_pkgs <- c("scran", "scater", "scRNAseq", "org.Mm.eg.db",
    "BiocFileCache", "DropletUtils", "EnsDb.Hsapiens.v86",
    "BiocSingular", "AnnotationHub", "TxDb.Mmusculus.UCSC.mm10.ensGene",
    "org.Hs.eg.db", "msigdbr"
)
BiocManager::install(bio_pkgs)
