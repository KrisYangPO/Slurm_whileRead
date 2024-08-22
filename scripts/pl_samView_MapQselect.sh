#!/bin/bash

# input parameters:
sampleID=$1
inputpath=$2
outputpath=$3
core=$4
mapq=$5

# grep input iles:
cd ${inputpath}
sample=$(ls ${sampleID}_*.out.bam)

# calculate count in whole region
samtools view \
 -q ${mapq} \
 -@ ${core} \
 -b -h \
 ${sample} \
 -o ${outputpath}/${sample/.out.bam/.mapQ${mapq}.out.bam}

# indexing
samtools index \
 -@ ${core} \
 ${outputpath}/${sample/.out.bam/.mapQ${mapq}.out.bam} \
 ${outputpath}/${sample/.out.bam/.mapQ${mapq}.out.bam.bai}


# options:
# -M: This avoids re-reading the same regions of files
# --fetch-pairs: Retrieve pairs even when the mate is outside of the requested region.

