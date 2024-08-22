#!/bin/bash

# input parameters
sampleID=$1
inputpath=$2
out=$3
core=$4

# grep files from sampleID at input folder.
# collect files based on sampleID, then create an array using outter "()"
cd ${inputpath}
samples=($(ls ${sampleID}_*fastq*))

# count the length of array: (${#array[@]})
samples_len=${#samples[@]}

echo "target samples: "${samples[@]}
echo "sample number: "${samples_len}


# determinee the length of the array
# 2: perform paired end
# 1: perform single end

if [ "$samples_len" -eq 2 ];then
  echo "perform paired-end trimming"
  trim_galore -j ${core} --illumina --gzip --paired --fastqc -q 30 -o ${out} ${samples[@]}

elif [ "$samples_len" -eq 1 ];then
  trim_galore -j ${core} --illumina --gzip --fastqc -q 30 -o ${out} ${samples[@]}

else
  echo "Files MISSING or Files ERROR"

fi
