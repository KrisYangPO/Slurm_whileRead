#!/bin/bash

# parameters
sampleID=$1
inputpath=$2
outputpath=$3
outputformat=$4
core=$5

# get a ${sampleID}_*.out.bam
# output name: sampleID_Aligned.sortedByCoord
cd ${inputpath}
sample=$(ls ${sampleID}_*.out.bam)
sample=${sample/.out.bam/}


# samtool builds index at current folder.
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

