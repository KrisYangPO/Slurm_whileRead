#!/bin/bash

# 基本 input
sampleID=$1
inputpath=$2
out=$3
core=$4
species=$5

# index:
# Use conditions to determines whether species contains string: mouse/human.
if [[ $species == *'mm10'* ]]; then
  index=${bowtie2idx_mm10}

elif [[ $species == *'chm13'* ]]; then
  index=${bowtie2idx_chm13}

elif [[ $species == *'hg38'* ]]; then
  index=${bowtie2idx_hg38}

fi


# grep files from sampleID at input folder.
# collect files based on sampleID, then create an array using outter "()"
cd ${inputpath}
samples=($(ls ${sampleID}_*.fq.gz))

# count the length of array: (${#array[@]})
samples_len=${#samples[@]}


# report
echo "Target files: "${samples[@]}
echo "Number of samples: "${samples_len}
echo "Genome index file: "${index}


# determinee the length of the array
# 2: perform paired end
# 1: perform single end

if [ "$samples_len" -eq 2 ];then
  bowtie2 -x ${index} -p ${core} -1 ${samples[0]} -2 ${samples[1]} -S ${out}/${sampleID}.sam

elif [ "$samples_len" -eq 1 ];then
  bowtie2 -x ${index} -p ${core} -U ${samples[0]} -S ${out}/${sampleID}.sam

else
  echo -e "receive more than 2 or less than 1 file\n file missing!"

fi
