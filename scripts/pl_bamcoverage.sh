#!/bin/bash

# parameters
sampleID=$1
inputpath=$2
outputpath=$3
outputformat=$4
core=$5

: '
1. STAR 會 output 出 bam 檔，可以用 bamcoverage 轉換成 bigwig 檔案，
   但是要先用 samtools 製作 bam.bai index。

2. STAR output 出的檔案名稱比較複雜：_Aligned.sortedByCoord，
   所以也不能用 ${sampleID}_Aligned.sortedByCoord.out.bam 寫入 scripts，
   這樣會僵化使用 (當今天後綴不再是"Aligned.sortedByCoord" 會爆掉)，
   因此要用抓取檔案的方式：$(ls *out.bam) 整個抓出來用。
'

# 抓取檔案
# 裡當應該只會對應到一個 ${sampleID}_*.out.bam
# 經過修改，sample 就會儲存 sampleID_Aligned.sortedByCoord
cd ${inputpath}
sample=$(ls ${sampleID}*.out.bam)
sample=${sample/.out.bam/}


# samtools build index
# 在原地製作 .bai 檔所以路徑也是用 ${inputpath}
samtools index \
 -@ ${core} ${inputpath}/${sample}.out.bam \
 ${inputpath}/${sample}.out.bam.bai


# program
bamCoverage \
 --normalizeUsing CPM \
 --outFileFormat ${outputformat} \
 --centerReads \
 --extendReads 150 \
 -p ${core} \
 -b ${inputpath}/${sample}.out.bam \
 -o ${outputpath}/${sample}.${outputformat}

