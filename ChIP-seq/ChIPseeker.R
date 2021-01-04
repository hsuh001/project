######
# ChIPseeker
######

require(tidyverse)



source("https://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")
# biocLite("org.Mm.eg.db")	#小鼠
# biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
# biocLite("clusterProfiler")
# biocLite("ReactomePA")
# biocLite("DOSE")

library(ChIPseeker)
# library(org.Mm.eg.db)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# library(clusterProfiler)

# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# 直接读文件出图
files <- getSampleFiles()
covplot(files[[4]])

# 支持GRanges对象，同时可以多个文件或者GRangeslist
peak <- GenomicRanges::GRangesList(
	CBX6 = readPeakFile(files[[4]]),
	CBX7 = readPeakFile(files[[5]])
	)
covplot(peak, weightCol="V5") + facet_grid(chr ~ .id)

# 支持可视化某个窗口
covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4e7, 5e7)) + facet_grid(chr ~ .id)






##########
# peak annotation
##########

# ChIPseeker自带5个BED文件
# getSampleFiles()拿到文件的全路径，返回names list
# 以第四个BED文件为例
f <- getSampleFiles()[[4]]

# 获取人类全基因组的注释信息
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# annotationPeak()进行peak注释
x <- annotationPeak(f, tssRegion = c(-1000, 1000), TxDb = txdb)
# x存储ChIPseq的位点落在基因组上什么样的区域，分布情况如何

# as.GRanges()查看具体信息
as.GRanges(x) %>%
	head(3)

# 第三种注释。看peak上下游某个范围内（比如说-5k到5k的距离）都有什么基因
x <- annotatePeak(f[[4]], 
	tssRegion = c(-1000, 1000), 
	TxDb = txdb, 
	addFlankGeneInfo = TRUE, 
	flankDistance = 5000)

# 基因注释。传入annoDb参数即可加入gene id
# 如果TxDb的基因ID是Entrez，会转成ENSEMBL；反之亦然。
# TxDb的基因ID类型是哪种，都会给出SYMBOL、描述性的GENENAME
x <- annotationPeak(f, tssRegion = c(-1000, 1000), TxDb = txdb, annoDb = "org.Hs.eg.db")





### 通过UCSC在线制作TxDb对象
require(GenomicFeatures)
hg19.refseq.db <- makeTxDbFromUCSC(genome = "hg19", table = "refGene")

### 通过解析GFF文件制作TxDb
download.file("ftp://ftp.ebi.ac.uk/pub/databases/pombase/pombe/Chromosome_Dumps/gff3/schizosaccharomyces_pombe.chr.gff3", "schizosaccharomyces_pombe.chr.gff3")
require(GenomicFeatures)
spombe <- makeTxDbFromGFF("schizosaccharomyces_pombe.chr.gff3")

