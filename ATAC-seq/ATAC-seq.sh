########################################
# ATAC-seq
# date: 2020.04.29 -- 2020.08.31
# author: Jing Xiao
########################################

#########################################################################
# 2020.04.29
#########################################################################
dir=/home1/jxiao/project/ATAC-seq/SRP055881
####################################################################(1) 下载sra数据
mkdir 0.0sra_data
cd $dir/0.0sra_data
cat srr_list.txt | while read id;
do
    nohup wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-2/$id/$id.1 -O $id.sra &
done

#########################################################################
# 2020.04.30
#########################################################################
####################################################################(2) sra2fq
mkdir -p $dir/1.0raw_data
cd $dir/0.0sra_data
# -A    使用sample编号作为输出文件名
cat config.txt | while read id;
do
    arr=($id)
    sample=${arr[0]}
    srr=${arr[1]}
    nohup fastq-dump -A $sample --gzip --split-3 -O $dir/1.0raw_data $srr.sra &
done

####################################################################(3) 过滤低质量reads、去除adaptor等
mkdir -p $dir/2.0clean_data
cd $dir/2.0clean_data
# -q                q值
# --phred33         使用ASCIl+33质量值作为Phred得分
# --length          去除长度小于参数值的reads
# -e                允许的最大误差
# --stringency      设置与接头重叠的序列
# --paired          对于双端测序文件，正反链质控通过才保留
# -o                输出目录

ls $dir/1.0raw_data/*_1.fastq.gz > raw_fq1.txt
ls $dir/1.0raw_data/*_2.fastq.gz > raw_fq2.txt
paste raw_fq1.txt raw_fq2.txt > raw_fq.txt
cat raw_fq.txt | while read id;
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    nohup trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4  --paired -o ./ $fq1 $fq2 &
done

####################################################################(4) 过滤前后的质控
mkdir -p $dir/3.0qc/{raw,clean}
cd $dir/3.0qc/raw
ls $dir/1.0raw_data/*.gz | while read id; do ( nohup fastqc -t 5 $id -o ./ & ); done
multiqc ./ -o ./ -n raw_multiqc_out

cd $dir/3.0qc/clean
ls $dir/2.0clean_data/*.gz | while read id; do ( nohup fastqc -t 5 $id -o ./ & ); done
multiqc ./ -o ./ -n clean_multiqc_out

####################################################################(5) 比对至参考基因组、sam2bam
mkdir -p $dir/4.0align
cd $dir/4.0align
ls $dir/2.0clean_data/*_1.fq.gz > clean_fq1.txt
ls $dir/2.0clean_data/*_2.fq.gz > clean_fq2.txt
paste clean_fq1.txt clean_fq2.txt > clean_fq.txt
bowtie2_index=/home/share/bowtie2_index/mm10
cat clean_fq.txt | while read id;
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    sample=$(basename $fq1 "_1_val_1.fq.gz")
    nohup bowtie2 -p 3 -X 1000 -x $bowtie2_index -1 $fq1 -2 $fq2 -S ${sample}.sam &
done

# sam2bam
ls *.sam | while read id; do (nohup samtools sort -O bam -@ 3 -o $(basename ${id} ".sam").bam ${id} &); done
ls *.bam | xargs -i samtools index {}

#########################################################################
# 2020.05.01
#########################################################################
####################################################################(6) 对bam文件进行qc
mkdir -p $dir/5.0stat
cd $dir/5.0stat
ls $dir/4.0align/*.bam | while read id; do (nohup samtools flagstat -@ 3 $id > $(basename ${id} ".bam").flagstat &); done

####################################################################(7) 去除PCR重复
cd $dir/4.0align
######## samtools markdup对PE数据去重需要按照以下步骤；而SE数据则不需要
# 使用samtools markdup去除，在这里会报错（[markdup] error: no MC tag. Please run samtools fixmate on file first.
# [markdup] error: unable to assign pair hash key.）
# 解决方法：https://github.com/samtools/samtools/issues/765
#ls *.bam | while read id ; do ( nohup samtools markdup -@ 2 -r $id $(basename $id ".bam").rmdup.bam &);done

# 按reads名排序（或许可以在第(5)步sam2bam过程中就添加-n选项）
ls *.bam | while read id; do
    sample=$(basename $id ".bam")
    nohup samtools sort -@ 3 -n -o ${sample}.sort.bam $id &
done

# Add ms and MC tags for markdup to use later
ls *.sort.bam | while read id; do
    sample=$(basename $id ".sort.bam")
    nohup samtools fixmate -@ 3 -m $id ${sample}.fixmate.bam &
done

# Markdup needs position orde
ls *.fixmate.bam | while read id; do
    sample=$(basename $id ".fixmate.bam")
    nohup samtools sort -@ 3 -o ${sample}.positionsort.bam $id &
done

# Finally mark duplicates and then remove
ls *.positionsort.bam | while read id; do
    sample=$(basename $id ".positionsort.bam")
    nohup samtools markdup -@ 2 -r  $id ${sample}.rmdup.bam &
done

rm *sort*
rm *fixmate*

# ls *.bam | while read id ; do ( nohup sambamba markdup -r $id $(basename $id ".bam").rmdup.bam &);done
# 循环使用sambamba markdup去重时会报错（sambamba-markdup: Cannot open file `/tmp/sambamba-pid155459-markdup-xbfj/PairedEndsInfoijnq2' in mode `w+' (Too many open files)）
# 解决方法：https://github.com/biod/sambamba/issues/229；https://github.com/biod/sambamba/issues/177
# 每个数据单独跑一遍，但如果数据过大进程会休眠。无法解决问题

ls *.rmdup.bam | xargs -i samtools index {}

(8) 对去重后的bam文件进行qc
cd $dir/5.0stat
ls $dir/4.0align/*.rmdup.bam | while read id; do (nohup samtools flagstat -@ 3 $id > $(basename ${id} ".bam").flagstat &); done

####################################################################
#(9) 过滤低质量reads、去除线粒体上的reads、保证两条reads比对到同一条染色体
cd $dir/4.0align
ls *.rmdup.bam | while read id; do
    sample=$(basename $id ".rmdup.bam")
    nohup samtools view -hf 2 -q 30 $id | grep -v chrM | samtools view -b -o ${sample}.tmp.bam - &
done

ls *.tmp.bam | while read id; do
    sample=$(basename $id ".tmp.bam")
    nohup samtools sort -O bam -@ 2 -o ${sample}.last.bam  $id &
done

# ！！！在这一步会报错
# [E::sam_parse1] missing SAM header
# [W::sam_read1] Parse error at line 1
# samtools sort: truncated file. Aborting
######## 已解决，显示缺少sam header内容，因此在samtools view要使用-h参数；
# 另外，输出的是bam文件便于排序使用
ls *.last.bam | xargs -i samtools index {}


####################################################################(10) 对再次过滤后的bam文件进行qc
cd $dir/5.0stat
ls $dir/4.0align/*.last.bam | while read id; do (nohup samtools flagstat -@ 3 $id > $(basename ${id} ".bam").flagstat &); done

#########################################################################
# 2020.08.12
#########################################################################
####################################################################(11) bam2bed
cd $dir/4.0align

### *.raw.bed用于计算FRiP值
cat rawbam.txt | while read id;
do
    nohup bedtools bamtobed -i $id > ${id%%.*}.raw.bed &
done

### *.last.bed用于peak calling
ls *.last.bam | while read id;
do
   nohup bedtools bamtobed -i $id > ${id%%.*}.last.bed &
done

#########################################################################
# 2020.08.29
#########################################################################
####################################################################(12) macs2进行peak calling
mkdir -p $dir/6.0macs2_peaks
cd $dir/6.0macs2_peaks
ls $dir/4.0align/*.last.bed | while read id;
do
   sample=$(basename $id ".last.bed")
   nohup macs2 callpeak -t $id -g mm --nomodel --shift -100 --extsize 200 \
   -n $sample --outdir ./ > ${sample}.macs2.log 2>&1 &
done

####################################################################(13) 计算插入片段长度
cd $dir/4.0align
ls *.last.bam | while read id;
do
  nohup samtools view $id | cut -f 9 > ${id%%.*}.insert_size.txt &
done
# 随后再使用reads_len_freq.R绘图

####################################################################(14) 计算FRiP值
cd $dir/6.0macs2_peaks
ls *narrowPeak | while read id;
do
    bed=../align/$(basename $id "_peaks.narrowPeak").raw.bed
    echo $(basename $id "_peaks.narrowPeak")
    Reads=$(bedtools intersect -a $bed -b $id | wc -l)
    totalReads=$(wc -l $bed | awk '{ print $1 }')
    echo $Reads $totalReads
    echo '==> FRiP value:' $(bc <<< "scale=2;100*$Reads/$totalReads")'%'
done

####################################################################(15) IDR计算技术重复性
cd $dir/6.0macs2_peaks
idr --samples 2-cell-1_peaks.narrowPeak 2-cell-2_peaks.narrowPeak --plot

####################################################################(16) 使用deeptools进行可视化
# step1：bam2bw
cd $dir/4.0align
ls *.last.bam | while read id;
do
    nohup bamCoverage -p 5--normalizeUsing CPM -b $id -o ${id%%.bam}.bw &
done
# step2：使用IGV查看

####################################################################(17) 查看TSS附近的信号强度
mkdir -p $dir/7.0tss_signal
cd $dir/7.0tss_signal

# --referencePoint  指定参考点（TSS,TES,center）。默认为TSS
# -p    指定线程
# -b    指定上游到参考点的距离
# -a    指定下游到参考点的距离
# -R    指定绘制区域的.bed 或 .gtf文件
# -S    要绘制的.bw文件
##首先对单一样本绘图：
# computeMatrix reference-point --referencePoint TSS -p 2 \
# -b 10000 -a 10000 \
# -R /home1/jxiao/project/ChIP-seq/annotation/mm10/ucsc.mm10.gene.bed \
# -S $dir/4.0align/2-cell-1.last.bw \
# --skipZeros \
# -o matrix_test_TSS.gz \
# --outFileSortedRegions regions_test_genes.bed

##both plotHeatmap and plotProfile will use the output from computeMatrix
#plotHeatmap -m matrix_test_TSS.gz -out test_Heatmap.png
#plotHeatmap -m matrix_test_TSS.gz -out test_Heatmap.pdf --plotFileFormat pdf --dpi 720
#plotProfile -m matrix_test_TSS.gz -out test_Profile.png
#plotProfile -m matrix_test_TSS.gz -out test_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720

## 也可以将多个样本一起画
nohup computeMatrix reference-point --referencePoint TSS -p 5 \
-b 10000 -a 10000 \
-R /home1/jxiao/project/ChIP-seq/annotation/mm10/ucsc.mm10.gene.bed \
-S $dir/4.0align/*.last.bw \
--skipZeros \
-o matrix_all_TSS.gz \
--outFileSortedRegions regions_all_genes.bed &

##both plotHeatmap and plotProfile will use the output from computeMatrix
plotHeatmap -m matrix_all_TSS.gz -out all_Heatmap.png
plotHeatmap -m matrix_all_TSS.gz -out all_Heatmap.pdf --plotFileFormat pdf --dpi 720
plotProfile -m matrix_all_TSS.gz -out all_Profile.png
plotProfile -m matrix_all_TSS.gz -out all_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720

##批处理
bed=/home1/jxiao/project/ChIP-seq/annotation/mm10/ucsc.mm10.gene.bed
# TSS上下游2K
# 还可以画10K等
for id in $dir/4.0align/*.last.bw;
do
    sample=$(basename $id ".last.bw")
    computeMatrix reference-point --referencePoint TSS -p 5 \
    -b 2000 -a 2000 -R $bed -S $id \
    --skipZeros -o matrix_${sample}_TSS_2K.gz \
    --outFileSortedRegions regions_${sample}_genes_2K.bed
    plotHeatmap -m matrix_${sample}_TSS_2K.gz -out ${sample}_Heatmap_2K.png
    plotHeatmap -m matrix_${sample}_TSS_2K.gz -out ${sample}_Heatmap_2K.pdf --plotFileFormat pdf --dpi 720
    plotProfile -m matrix_${sample}_TSS_2K.gz -out ${sample}_Profile_2K.png
    plotProfile -m matrix_${sample}_TSS_2K.gz -out ${sample}_Profile_2K.pdf --plotFileFormat pdf --perGroup --dpi 720
done

####################################################################(18) 查看gene_body附近的信号强度
mkdir -p $dir/7.1gene_body_signal
cd $dir/7.1gene_body_signal
nohup computeMatrix  scale-regions --regionBodyLength 2000 -p 5 \
-b 10000 -a 10000 -bs 10 \
-R /home1/jxiao/project/ChIP-seq/annotation/mm10/ucsc.mm10.gene.bed \
-S $dir/4.0align/*.last.bw \
--skipZeros \
-o matrix_all_body.gz &

##both plotHeatmap and plotProfile will use the output from computeMatrix
plotHeatmap -m matrix_all_body.gz -out all_Heatmap.png
plotHeatmap -m matrix_all_body.gz -out all_Heatmap.pdf --plotFileFormat pdf --dpi 720
plotProfile -m matrix_all_body.gz -out all_Profile.png
plotProfile -m matrix_all_body.gz -out all_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720

#########################################################################
# 2020.08.30
#########################################################################
####################################################################(19) peak注释
## method_1, 使用R包ChIPseeker, peakAnno.R
mkdir -p $dir/8.0peakAnno_ChIPseeker
cd $dir/8.0peakAnno_ChIPseeker
touch peakAnno.R

## method_2, 使用homer
# 下面下载mm10的genome
# per1 /usr/local/software/homer/configureHomer.pl -install mm10
mkdir -p $dir/8.1peakAnno_homer
cd $dir/8.1peakAnno_homer
ls $dir/6.0macs2_peaks/*.narrowPeak | while read id;
do
    sample=$(basename $id "_peaks.narrowPeak")
    awk '{ print $4"\t"$1"\t"$2"\t"$3"\t+" }' $id > ${sample}.homer_peaks.tmp
    nohup annotatePeaks.pl ${sample}.homer_peaks.tmp mm10 \
    1>${sample}.peakAnn.xls \
    2>${sample}.annLog.txt &
done

####################################################################(20) 寻找motif
## method_1, 使用homer
# -bg <peak/BED file> 自定义背景
# -size: 用于motif寻找得片段大小, default: 200。越大需要得计算资源越多
#        -size <#>
#        -size <#>,<#>  -size -300,100：peak上游100bp，下游300bp区域
#        -size given 设置片段大小为目标序列长度
# -len：motif大小设置，default 8,10,12。越大需要得计算资源越多
#       -len <#> or -len <#>,<#>,...
# -S <#>：结果输出motifs的个数, default 25
# -mis <#>：motif错配碱基数，default 2bp
#           允许错配可以提升灵敏度，如果寻找12-15 bp Motif ，可以设置3-4bp的错配
# -nlen <#>：消除短寡聚核苷酸引入的不平衡。default 3, "-nlen 0" to disable
# -norevopp：只在+链搜索motif
# -rna: 输出RNA motif，使用RNA motif数据库
# -h：使用超几何检验代替二项式分布
# -N：用于motif寻找得背景序列数目，default=max(50k, 2x input)；耗内存参数
mkdir -p $dir/9.0motif_homer
cd $dir/9.0motif_homer
ls $dir/8.1peakAnno_homer/*.homer_peaks.tmp | while read id;
do
    sample=$(basename $id ".homer_peaks.tmp")
    nohup findMotifsGenome.pl $id mm10  ${sample}_motifDir -size 200 -len 8,10,12 &
done

#########################################################################
# 2020.08.31
#########################################################################
## method_1, meme
# 获取序列：https://github.com/jmzeng1314/NGS-pipeline/blob/master/CHIPseq/stepZ-peaks2sequence.R
## method_1, R包，比如motifmatchr包
# https:/bioconductor.org/packages/release/b ioc/html/motifmatchr.htmlr

####################################################################(19) 差异peaks
## method_1, R包，DiffBind