1、example workflow
(1) http://jura.wi.mit.edu/bio/education/hot_topics/QC_HTP/QC_HTP.pdf
(2) http://jura.wi.mit.edu/bio/education/hot_topics/RNAseq/RNAseaDE_Dec2011.pdf
(3) https://f1000research.com/articles/4-1070/v1
(4) https://f1000research.com/articles/5-1438/v1
(5) https://www.bioconductor.org/help/workflows/rnaseaGene/
(6) 一个RNA-seq实战-超级简单-2小时搞定！（http://www.bio-info-trainee.com/2218.html）
(7) 一个植物转录组项目的实战（http://www.bio-info-trainee.com/2809.html）

2、outline
(1) 参考基因组及参考转录组
    gtf、genome.fa
(2) 质控
    fastqc、multiqc、trimmomatic、cutadapt、trim_galore
(3) 比对
    star、hisat2、tophat、bowtie2、bwa、subread
(4) 计数
    htseq、bedtools、deeptools、salmon
(5) 归一化、差异分析等
    DEseq2、edgeR、limma（voom）

3、准备工作
(1) 安装软件
(2) 下载人类参考基因组
## gtf存放参考基因组注释
$ mkdir -p ~/public/reference/gtf/hg38
$ cd ~/public/reference/gtf/hg38
# ensembl
# $ wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.chr.gtf.gz
# ucsc
$ wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.knownGene.gtf.gz


## genome存放参考基因组
$ mkdir -p ~/public/reference/genome/hg38
$ cd ~/public/reference/genome/hg38
# ensembl
# $ wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
# ucsc
$ wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz
$ gunzip -c hg38.fa.gz > hg38.fa

# 依次更改染色体{1..22}、X、Y、MT
$ wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz

(3) 构建index
## index存放索引
$ mkdir -p ~/public/reference/index/{bowtie2,bwa,hisat2,STAR,subread}
## 构建索引index

# bowtie2
$ cd ~/public/reference/index/bowtie2
$ nohup bowtie2-build ~/public/reference/genome/hg38/hg38.fa hg38 &
### root共享的index: /home/share/bowtie2_index/

# bwa
### root共享的index: /home/share/bwa_index/

# hisat2
### root共享的index: /home/share/hisat2_index/

# STAR
$ cd ~/public/reference/index/STAR
$ STAR --runMode genomeGenerate \
--genomeDir ~/reference/index/star \
--genomeFastaFiles ~/reference/genome/hg38/hg38.fa \
-sjdbGTFfile ~/reference/gtf/hg38/hg38.knownGene.gtf.gz \
--sjdbOverhang 149--runThreadN 4
### root共享的index: /home/share/bowtie2_STAR_hg38_indexindex/

4、分析示例
(1) 创建目录
$ mkdir -p ~/project/RNA/work1/{sra,fastq,fastqc_out,multiqc_out,clean}
$ cd ~/project/RNA/work1

(2) 下载sra数据
## 方法一，使用 wget
$ cd sra
$ for((id=677;id<=680;id++)); do ( nohup wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR957$id/SRR957$id.1 -O SRR957$id.sra & ); done
## 方法二，使用 prefetch
# 将sra数据的编号保存至id.txt
# 下载的sra数据目录：~/ncbi/public/sra
prefetch="/home1/jxiao/biosoft/sratoolkit/sratoolkit.2.10.5-centos_linux64/bin/prefetch"
$ cat id.txt | while read id; do ( nohup $prefetch $id & ); done
$ mv ~/ncbi/public/sra/* ./

(3) sra2fastq，得到的.fastq文件存在 ../fastq目录
$ ls *.sra | while read id; do ( nohup fastq-dump --gzip --split-3 -O ../fastq/ $id & ); done

(4) 使用fastqc进行qc
$ cd ../fastq
$ ls *.gz | while read id; do ( nohup fastqc $id -o ../fastqc_out/ & ); done 
#或 $ ls ./fastq/*.gz | xargs $fastqc -o ./fastqc_out/

(5) 使用multiqc将fastqc生成的多个报告整合成一个报告
$ cd ~/project/RNA/work1
$ multiqc ./fastqc_out -o ./multiqc_out

(6) 过滤低质量reads、去除adaptor等
$ cd ~/project/RNA/work1/clean
$ ls ../fastq/*_1.fastq.gz > 1
$ ls ../fastq/*_2.fastq.gz > 2
$ paste 1 2 > config
$ trim_galore=
$ cat config | while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
nohup $trim_galore -q 25 --phread33 --length 36 -e 0.1 --stringency 3 --paired -o ./ $fq1 $fq2 &
done

(7) 将fastq比对至参考基因组
## 比对
$ cd ~/project/RNA/work1/clean
## 取前1000行为例
# basename $id ".gz"可以去除文件名的.gz后缀。防止软件因后缀错误识别文件格式
#$ ls ./clean/*gz | while read id; do (zcat $id | head-1000 > $(basename $id ".gz")); done

$ ls *gz | cut -d "_" -f 1 | sort -u | while read id;
do
# 使用hisat2
nohup hisat2 -p 10 -x /public/reference/index/hisat/hg38 -1 ${id}_1_val_1.fq.gz -2 ${id}_2_val_2.fq.gz -S ${id}.hisat.sam &

# 使用subjunc
nohup subjunc -T 5 -i /public/reference/index/subread/hg38 -r ${id}_1_val_1.fq.gz -R ${id}_2_val_2.fq.gz -o ${id}.subjunc.sam &

# 使用bowtie2
nohup bowtie2 -p 10 -x /public/reference/index/bowtie/hg38 -1 ${id}_1_val_1.fq.gz -2 ${id}_2_val_2.fq.gz -S ${id}.bowtie2.sam &
# nohup bowtie2 -p 2 -x /home/share/bowtie2_index/hg38 -U ${id}.fastq.gz -S ${id}.bowtie2.sam &

# 使用bwa
nohup bwa mem -t 5 -M /public/reference/index/bwa/hg38 ${id}_1_val_1.fq.gz ${id}_2_val_2.fq.gz > ${id}.bwa.sam &
done


(8) sam2bam
$ cd ~/project/RNA/work1/clean
ls *.sam | while read id; do (nohup samtools sort -O bam -@ 2 -o $(basename ${id} ".sam").bam ${id} &); done
# ls *.sam | while read id; do (nohup samtools sort -n -@ 2 -o ${id%%.*}.Nsort.bam $id &);done
rm *.sam
ls *.bam | xargs -i samtools index {}

(9) 统计bam文件的信息
$ cd .. && mkdir stat && cd stat
$ ls ../*.bam | while read id; do (nohup samtools flagstat -@ 2 $id > $(basename ${id} ".bam").flagstat &); done

(10) counts
$ cd ~/project/RNA/work1/clean

## root共享的gtf：/home/share/gtf/hg38.gtf
$ featureCounts -T 2 -t exon -g gene_id -a /home/share/gtf/hg38.gtf -o geneCounts.txt *.bam





$ htseq-count -f bam -r pos ./SRR1039516.subjunc.bam ~/reference/gtf/
