1、分析流程
(1) 获取sra数据

(2) 质控
    fastqc、multiqc、trim_galore
(3) alignment
    bowtie2、bwa
(4) peak calling
    macs2
(5) sam2bw
    deeptools
(6) motif
    homer、meme
(7) peak annotation
    ChIPseeker
(8) peaks差异分析、peaks基因富集分析


2、workflow
(1) 下载sra数据
$ mkdir  -p ~/project/ChIP-seq/{sra,fastq,fastqc_out,multiqc_out,clean,align,peaks,motif,stat,tss}
$ cd /home1/jxiao/project/ChIP-seq/sra
$ cat sra_id.txt | while read id; do ( nohup wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-2/$id/$id.2 -O $id.sra & ); done

# 由sra.table得到config文件，vim修改
$ cat sra.table | cut -d ':' -f 2 > config

(2) sra2fq
$ cd /home1/jxiao/project/ChIP-seq/sra
# 使用sra编号作为输出文件名
# $ ls *.sra | while read id; do ( nohup fastq-dump --gzip --split-3 -O ../fastq/ $id & ); done
# 使用sample编号作为输出文件名
$ cat config | while read id; do
arr=($id)
srr=${arr[1]}
sample=${arr[0]}
nohup fastq-dump -A $sample --gzip --split-3 -O ../fastq/ $srr.sra &
done

(3) 对原始fq文件使用fastqc进行qc
$ cd /home1/jxiao/project/ChIP-seq/fastq
$ ls *.gz | while read id; do ( nohup fastqc $id -o ../fastqc_out/raw/ & ); done
$ multiqc ../fastqc_out/raw/ -o ../multiqc_out/raw/

(4) 过滤低质量reads、去除adaptor等
$ cd /home1/jxiao/project/ChIP-seq/fastq
$ ls *gz | while read id; do ( nohup trim_galore -q 25 --phred33 --length 25 -e 0.1 --stringency 4 -o ../clean/ $id & ); done

(5) 对过滤后的fq文件再次进行qc
$ cd /home1/jxiao/project/ChIP-seq/clean
$ ls *.gz | while read id; do ( nohup fastqc $id -o ../fastqc_out/clean/ & ); done
$ multiqc ../fastqc_out/clean/ -o ../multiqc_out/clean/

(6) 将过滤后的fq比对至参考基因组
$ cd /home1/jxiao/project/ChIP-seq/align
$ ls ../clean/*gz | while read id; do
file=$(basename $id)
sample=${file%%.*}
nohup bowtie2 -p 2 -x /home/share/bowtie2_index/mm10 -U $id -S ${sample}.sam &
done

(7) sam2bam
$ cd /home1/jxiao/project/ChIP-seq/align
$ ls *.sam | while read id; do ( nohup samtools sort -O bam -@ 2 -o $(basename ${id} ".sam").bam ${id} &); done
$ rm *.sam
$ ls *.bam | xargs -i samtools index {}
## 以上(6)(7)也可一步进行
# $ cd /home1/jxiao/project/ChIP-seq/align
# $ ls ../clean/*gz | while read id; do
# file=$(basename $id)
# sample=${file%%.*}
# nohup bowtie2 -p 2 -x /home/share/bowtie2_index/mm10 -U $id | samtools sort -O bam -@ 2 -o - > ${sample}.bam &
# done

(8) 对bam文件进行qc
$ cd /home1/jxiao/project/ChIP-seq/stat
$ ls ../align/*.bam | while read id; do (nohup samtools flagstat -@ 2 $id > $(basename ${id} ".bam").flagstat &); done

(9) 如果一个样品分成多个lane进行测序，在peaks calling时，需要合并bam
$ mkdir -p /home1/jxiao/project/ChIP-seq/mergeBam
$ cd /home1/jxiao/project/ChIP-seq/align
$ ls *trimmed.bam | sed 's/[0-9]_trimmed.bam//g' | sort -u | while read id; do ( nohup samtools merge ../mergeBam/${id}merge.bam ${id}*trimmed.bam & ); done

(10) 对合并后的bam文件进行qc
$ cd /home1/jxiao/project/ChIP-seq/mergeBam
$ ls *.bam | xargs -i samtools index {}
$ cd /home1/jxiao/project/ChIP-seq/stat
$ ls ../mergeBam/*.bam | while read id; do (nohup samtools flagstat -@ 2 $id > $(basename ${id} ".bam").flagstat &); done


(11) （如果存在PCR重复）去除PCR重复（samtools、picard）
 $ cd /home1/jxiao/project/ChIP-seq/mergeBam
 $ ls *.bam | while read id ; do ( nohup samtools markdup -@ 2 -r $id $(basename $id ".bam").rmdup.bam &);done
 $ ls *.rmdup.bam | xargs -i samtools index {}

(12) 对去除PCR重复的.rmdup.bam进行qc
 $ cd /home1/jxiao/project/ChIP-seq/stat
 $ ls ../mergeBam/*.rmdup.bam | while read id; do (nohup samtools flagstat -@ 2 $id > $(basename ${id} ".bam").flagstat &); done

(13) 使用macs2进行peak calling
$ cd /home1/jxiao/project/ChIP-seq/mergeBam
# -t    实验组
# -c    对照组
# -f    指定-t、-c文件的格式。不提供则自动检测
# -g    基因组大小，bp为单位。默认提供hs、mm、ce、dm等
# -n    输出文件的前缀名
# -B    输出bedgraph文件
# -q    阈值：q值.默认为0.05
# -p    p值。指定p值就不再使用q值

$ ls *merge.bam | while read id; do
sample=$(basename $id "_merge.bam")
nohup macs2 callpeak -c Control_merge.bam -t ${sample}_merge.bam -f BAM -g mm -n $sample -B --outdir ../peaks/raw/ 2>../peaks/raw/$sample.log &
done

$ ls *merge.rmdup.bam | while read id; do
sample=$(basename $id "_merge.rmdup.bam")
nohup macs2 callpeak -c Control_merge.rmdup.bam -t ${sample}_merge.rmdup.bam -f BAM -g mm -n ${sample}_rmdup -B --outdir ../peaks/rmdup/ 2>../peaks/rmdup/${sample}_rmdup.log &
done

(14) 使用deeptools进行可视化
## step1：bam2bw
$ cd /home1/jxiao/project/ChIP-seq/mergeBam
$ ls *.bam |while read id; do ( nohup bamCoverage --normalizeUsing CPM -b $id -o $(basename $id ".bam").bw & ); done

## step2：使用IGV查看

(15) 查看TSS附近的信号强度
$ cd /home1/jxiao/project/ChIP-seq/tss

# --referencePoint  指定参考点（TSS,TES,center）。默认为TSS
# -p    指定线程
# -b    指定上游到参考点的距离
# -a    指定下游到参考点的距离
# -R    指定绘制区域的.bed 或 .gtf文件
# -S    要绘制的.bw文件
##首先对单一样本绘图：
$ computeMatrix reference-point --referencePoint TSS -p 2 \
-b 10000 -a 10000 \
-R /home1/jxiao/project/ChIP-seq/annotation/mm10/ucsc.mm10.gene.bed \
-S ../mergeBam/H2Aub1_merge.bw \
--skipZeros \
-o matrix_test_TSS.gz \
--outFileSortedRegions regions_test_genes.bed

##both plotHeatmap and plotProfile will use the output from computeMatrix
$ plotHeatmap -m matrix_test_TSS.gz -out test_Heatmap.png
$ plotHeatmap -m matrix_test_TSS.gz -out test_Heatmap.pdf --plotFileFormat pdf --dpi 720
$ plotProfile -m matrix_test_TSS.gz -out test_Profile.png
$ plotProfile -m matrix_test_TSS.gz -out test_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720


##批处理
bed="/home1/jxiao/project/ChIP-seq/annotation/mm10/ucsc.mm10.gene.bed"

# TSS上下游2K
for id in /home1/jxiao/project/ChIP-seq/mergeBam/*bw;
do 
    sample=$(basename $id ".bw")
    computeMatrix reference-point --referencePoint TSS -p 5 \
    -b 2000 -a 2000 -R $bed -S $id \
    --skipZeros -o matrix_${sample}_TSS_2K.gz \
    --outFileSortedRegions regions_${sample}_genes_2K.bed
    plotHeatmap -m matrix_${sample}_TSS_2K.gz -out ${sample}_Heatmap_2K.png
    plotHeatmap -m matrix_${sample}_TSS_2K.gz -out ${sample}_Heatmap_2K.pdf --plotFileFormat pdf --dpi 720
    plotProfile -m matrix_${sample}_TSS_2K.gz -out ${sample}_Profile_2K.png
    plotProfile -m matrix_${sample}_TSS_2K.gz -out ${sample}_Profile_2K.pdf --plotFileFormat pdf --perGroup --dpi 720
done

# TSS上下游10K
for id in /home1/jxiao/project/ChIP-seq/mergeBam/*bw;
do 
    sample=$(basename $id ".bw")
    computeMatrix reference-point --referencePoint TSS -p 2 \
    -b 10000 -a 10000 -R $bed -S $id \
    --skipZeros -o matrix_${sample}_TSS_10K.gz \
    --outFileSortedRegions regions_${sample}_genes_10K.bed
    plotHeatmap -m matrix_${sample}_TSS_10K.gz -out ${sample}_Heatmap_10K.png
    plotHeatmap -m matrix_${sample}_TSS_10K.gz -out ${sample}_Heatmap_10K.pdf --plotFileFormat pdf --dpi 720
    plotProfile -m matrix_${sample}_TSS_10K.gz -out ${sample}_Profile_10K.png
    plotProfile -m matrix_${sample}_TSS_10K.gz -out ${sample}_Profile_10K.pdf --plotFileFormat pdf --perGroup --dpi 720
done

(16) peak注释、基因富集分析
## 方法一：ChIPseeker包
# 输入：/home1/jxiao/project/ChIP-seq/peaks/raw/和/home1/jxiao/project/ChIP-seq/peaks/rmdup/下的.bed文件

## 方法二：homer寻找motif、peak注释
$ cd /home1/jxiao/project/ChIP-seq/motif/raw/
for id in /home1/jxiao/project/ChIP-seq/peaks/raw/*bed;
do
    file=$(basename $id)
    sample=${file%%.*}
    awk '{ print $4"\t"$1"\t"$2"\t"$3"\t+" }' $id > homer_peaks.tmp
    findMotifsGenome.pl homer_peaks.tmp mm10 ${sample}_motifDir -len 8,10,12
    annotatePeaks.pl homer_peaks.tmp mm10 1>${sample}.peakAnn.x1s 2> ${sample}.annLog.txt
done

$ cd /home1/jxiao/project/ChIP-seq/motif/rmdup/
for id in /home1/jxiao/project/ChIP-seq/peaks/rmdup/*bed;
do
    file=$(basename $id)
    sample=${file%%.*}
    echo $sample
    awk '{ print $4"\t"$1"\t"$2"\t"$3"\t+" }' $id > homer_peaks.tmp
    findMotifsGenome.pl homer_peaks.tmp mm10 ${sample}_motifDir -len 8,10,12
    annotatePeaks.pl homer_peaks.tmp mm10 1>${sample}.peakAnn.x1s 2> ${sample}.annLog.txt
done

## 方法三：meme