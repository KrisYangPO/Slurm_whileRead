#!/bin/bash

sampleID=$1
Path_main=$2
core=$3

#samtobam sorttobam rmdup build index
samtools view \
 --threads ${core} -b ${Path_main}/${sampleID}.sam > ${Path_main}/${sampleID}.bam

# sort bam files
samtools sort \
 --threads ${core} -O bam ${Path_main}/${sampleID}.bam -o ${Path_main}/${sampleID}.sorted.bam

# remove duplicates
samtools rmdup \
 ${Path_main}/${sampleID}.sorted.bam ${Path_main}/${sampleID}.rmdup.bam

# build index file
samtools index \
 -@ ${core} ${Path_main}/${sampleID}.rmdup.bam ${Path_main}/${sampleID}.rmdup.bam.bai

# mapping ratio
samtools flagstat \
 --threads ${core} ${Path_main}/${sampleID}.sorted.bam > ${Path_main}/${sampleID}.sorted.flagstat

# mapping ratio
samtools flagstat \
 --threads ${core} ${Path_main}/${sampleID}.rmdup.bam > ${Path_main}/${sampleID}.rmdup.flagstat

