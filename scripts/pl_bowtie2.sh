#!/bin/bash

# 基本 input
sampleID=$1
inputpath=$2
out=$3
core=$4
species=$5

# index:
# 用判斷式判斷 speci 裡面有沒有 mm10/chm13/hg38 字串。
if [[ $species == *'mm10'* ]]; then
  index=${bowtie2idx_mm10}

elif [[ $species == *'chm13'* ]]; then
  index=${bowtie2idx_chm13}

elif [[ $species == *'hg38'* ]]; then
  index=${bowtie2idx_hg38}

fi


# 抓取資料，移動到 input 所在位置
# 根據 sample ID 找尋 (ls) 名稱，之後建立陣列：array=($(command))
cd ${inputpath}
samples=($(ls ${sampleID}_*.fq.gz))

# 計數陣列：有幾個 samples 被抓 (${#array[@]})
samples_len=${#samples[@]}


# report
echo "Target files: "${samples[@]}
echo "Number of samples: "${samples_len}


# 判斷陣列長度是否等於 2
# 2: 執行 paired end
# 1: 執行 single end

if [ "$samples_len" -eq 2 ];then
  bowtie2 -x ${index} -p ${core} -1 ${samples[0]} -2 ${samples[1]} -S ${out}/${sampleID}.sam

elif [ "$samples_len" -eq 1 ];then
  bowtie2 -x ${index} -p ${core} -U ${samples[0]} -S ${out}/${sampleID}.sam

else
  echo -e "receive more than 2 or less than 1 file\n file missing!"

fi
