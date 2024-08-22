#!/bin/bash

# parameters
sampleID=$1
inputpath=$2
core=$3

# get a ${sampleID}_*.out.bam
# output name: sampleID_Aligned.sortedByCoord
cd ${inputpath}
sample=$(ls ${sampleID}_*.out.bam)
sample=${sample/.out.bam/}


# samtool builds index at current folder.
samtools index \
 -@ ${core} ${inputpath}/${sample}.out.bam \
 ${inputpath}/${sample}.out.bam.bai

