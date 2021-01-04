########################################
# WGBS
# date: 2020.08.26 -- 2020.09.02
# author: Jing Xiao
# ref:
# https://www.bilibili.com/video/BV1Ub41187CW?p=15
# https://www.bilibili.com/video/BV1Ub41187CW?p=16
########################################

#########################################################################
# 2020.08.26
#########################################################################
dir=/home1/jxiao/project/WGBS
####################################################################(1) 下载sra数据
mkdir 0.0sra_data
cd $dir/0.0sra_data
# sra_info.txt的内容：
# EV1 SRR1916129
# 3bstrain1 SRR1916134
# set1rep1 SRR1916142
cat sra_info.txt | while read id;
do
    arr=($id)
    srr=${arr[1]}
    nohup wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/$srr/${srr}.1 -O ${srr}.sra &
done

####################################################################(2) fastq-dump, sra2fq
mkdir -p $dir/1.0raw_data
cd $dir/0.0sra_data
# 使用sra编号作为输出文件名
# $ ls *.sra | while read id; do ( nohup fastq-dump --gzip --split-3 -O ../fastq/ $id & ); done
# 使用sample编号作为输出文件名
cat sra_info.txt | while read id;
do
    arr=($id)
    srr=${arr[1]}
    sample=${arr[0]}
    nohup fastq-dump -A $sample --gzip --split-3 -O $dir/1.0raw_data $srr.sra &
done

####################################################################(3) fastqc, qc
cd $dir/1.0raw_data
mkdir -p $dir/2.0qc/{raw,clean}
ls *.gz | while read id; do ( nohup fastqc $id -o $dir/2.0qc/raw & ); done
multiqc $dir/2.0qc/raw -o $dir/2.0qc/raw -n raw_multiqc_report

####################################################################(4) trim_galore, 过滤低质量reads
cd $dir/1.0raw_data
mkdir -p $dir/3.0clean_data
ls *gz | while read id;
do
    nohup trim_galore -q 25 --length 20 -o $dir/3.0clean_data $id \
    > $dir/3.0clean_data/$(basename $id ".fastq.gz").log 2>&1 &
done

####################################################################(5) fastqc, 过滤后qc
cd $dir/3.0clean_data
ls *.gz | while read id; do ( nohup fastqc $id -o $dir/2.0qc/clean & ); done
multiqc $dir/2.0qc/clean -o $dir/2.0qc/clean -n clean_multiqc_report

####################################################################(6) 练习的是酵母，服务器没有酵母的bismark_index，手动构建
mkdir $dir/bismark_index
cd $dir/bismark_index
wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
cd $dir
bismark_genome_preparation bismark_index/

####################################################################(6) 将过滤后的fq比对至参考基因组
# 默认使用bowtie2进行比对
cd $dir/3.0clean_data
mkdir -p $dir/4.0bismark_align
ls *gz | while read id;
do
    nohup bismark -p 2 -N 1 -o $dir/4.0bismark_align $dir/bismark_index $id &
done

####################################################################(7) deduplicate_bismark, 去重
cd $dir/4.0bismark_align
mkdir -p $dir/5.0bismark_dedup
ls *.bam |while read id;
do
    nohup deduplicate_bismark -s --output_dir $dir/5.0bismark_dedup --bam $id \
    > $dir/5.0bismark_dedup/$(basename $id "_trimmed_bismark_bt2.bam").dedup.log &
done

####################################################################(8) bismark_methylation_extractor, 统计甲基化数据
# 以3bstrain1_trimmed_bismark_bt2.deduplicated.bam为例
cd $dir/5.0bismark_dedup
mkdir -p $dir/6.0met_extractor

nohup bismark_methylation_extractor -s --gzip --parallel 2 \
-o $dir/6.0met_extractor 3bstrain1_trimmed_bismark_bt2.deduplicated.bam \
> $dir/6.0met_extractor/3bstrain1.met_ext.log &

#########################################################################
# 2020.09.02
#########################################################################
mkdir -p $dir/6.1met_extractor
# --comprehensive 生成CHG、CHH、CpG
# --merge_non_CpG 在--comprehensive模式下，将CHG、CHH合并，生成Non_CpG
# --cytosine_report生成cytosine methylation report
nohup bismark_methylation_extractor -s --gzip --parallel 5 \
--comprehensive --merge_non_CpG --bedGraph --buffer_size 10G \
--cytosine_report --genome_folder $dir/bismark_index \
-o $dir/met_extractor_2 3bstrain1_trimmed_bismark_bt2.deduplicated.bam \
> $dir/met_extractor_2/3bstrain1.met_ext.log &

### coverage2cytosine, 生成cytosine methylation report
mkdir coverage2cytosine
ls *bismark.cov.gz | while read id;
do
    sample=${id%%_*}
    nohup coverage2cytosine --dir ./coverage2cytosine --merge_CpG --gzip \
    --genome_folder $dir/bismark_index -o $sample $id > ./coverage2cytosine/${sample}.cov2cyto.log &
done