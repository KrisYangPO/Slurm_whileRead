#!/bin/bash
# htseq human genome

# parameters
sampleID=$1
inputpath=$2
outputpath=$3
core=$4
species=$5


# index:
# 用判斷式判斷 speci 裡面有沒有 mouse/human 字串。
if [[ $species == *'mouse'* ]]; then
  gtf=${gtf_mm10}

elif [[ $species == *'human'* ]]; then
  gtf=${gtf_chm13}

fi


# 抓取所有 bam 檔案，直接執行所有 htseq-count
cd ${inputpath}
sample=$(ls ${sampleID}_*.out.bam)
samples=($(ls *.out.bam))


# Run HTseq-count
# files 要放在 gtf file 前面！！
htseq-count -m intersection-nonempty --nonunique all \
 -f bam -s reverse \
 -t exon \
 --idattr gene_name \
 -n ${core} \
 ${sample} \
 ${gtf} \
 > ${outputpath}/HTseq_${sampleID}.txt


# Run HTseq-count and collect all samples into a table
htseq-count -m intersection-nonempty --nonunique all \
 -f bam -s reverse \
 -t exon \
 --idattr gene_name \
 -n ${core} \
 ${samples[@]} \
 ${gtf} \
 > ${outputpath}/HTseq_ALL.txt


echo "Sample: "${sampleID}
echo "Species: "${species}
echo "GTF file: "$(basename ${gtf})
