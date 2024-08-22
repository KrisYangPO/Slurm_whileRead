#!/bin/bash

# input parameters:
sampleID=$1
inputpath=$2
outputpath=$3
core=$4
mapq=30
region=/staging/biology/ls807terra/0_bedfiles/hTERRA/CHM13_TERRA_region_v7.bed

# grep input iles:
cd ${inputpath}
sample=$(ls ${sampleID}_*.out.bam)

# create folder for each sample:
output_sample=${outputpath}/${sampleID}
mkdir ${output_sample}


# calculate count in whole region
samtools view \
 --fetch-pairs \
 -q ${mapq} \
 -@ ${core} \
 -h -b -M \
 -L ${region} \
 ${sample} \
 -o ${sample/.out.bam/.selectedRegions.bam}

# indexing
samtools index \
 -@ ${core} \
 ${sample/.out.bam/.selectedRegions.bam}
mv ${sample/.out.bam/.selectedRegions.bam*} ${output_sample}


# Individule region
mkdir ${output_sample}/search_each_Region_Aligned
rl=$(awk '{printf "%s:%d-%d ",$1,$2,$3}' ${region})

for each in ${rl[@]}
do
 samtools view -h --fetch-pairs -q ${mapq} -@ ${core} ${sample} ${each} > ${output_sample}/search_each_Region_Aligned/${each/:/_}.sam
done


# Collecting read ID from selected alignment
mkdir ${output_sample}/search_each_Region_Aligned_readID

for each in $(ls ${output_sample}/search_each_Region_Aligned)
do
 samtools view -@ ${core} -q ${mapq} ${output_sample}/search_each_Region_Aligned/${each} | \
 awk '{print $1}' > ${output_sample}/search_each_Region_Aligned_readID/${each/.sam/.readID.txt}
done


# Count chromosome-specific TERRA read counts
for i in $(ls ${output_sample}/search_each_Region_Aligned_readID)
do
 echo ${i/.readID.txt/},$(cat ${output_sample}/search_each_Region_Aligned_readID/${i} | sort | uniq | wc -l | awk '{print $1}')
done | sed 's/,/\t/g' > ${output_sample}/Summary_of_eachRegion_read_count.txt


# options:
# -M: This avoids re-reading the same regions of files
# --fetch-pairs: Retrieve pairs even when the mate is outside of the requested region.

