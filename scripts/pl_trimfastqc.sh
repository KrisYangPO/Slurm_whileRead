#!/bin/bash

# input parameters
sampleID=$1
inputpath=$2
out=$3
core=$4

# 抓取資料，移動到 input 所在位置
# 根據 sample ID 找尋 (ls) 名稱，之後建立陣列：array=($(command))
cd ${inputpath}
samples=($(ls ${sampleID}_*fastq*))

# 計數陣列長度：檔案有幾個 (${#array[@]})
samples_len=${#samples[@]}

echo "target samples: "${samples[@]}
echo "sample number: "${samples_len}


# 判斷陣列長度是否等於 2
# 2: 執行 paired end
# 1: 執行 single end

if [ "$samples_len" -eq 2 ];then
  echo "perform paired-end trimming"
  trim_galore -j ${core} --illumina --gzip --paired --fastqc -q 30 -o ${out} ${samples[@]}

elif [ "$samples_len" -eq 1 ];then
  trim_galore -j ${core} --illumina --gzip --fastqc -q 30 -o ${out} ${samples[@]}

else
  echo "Files MISSING or Files ERROR"

fi
