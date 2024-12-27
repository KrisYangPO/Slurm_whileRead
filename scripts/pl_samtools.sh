#!/bin/bash

sampleID=$1
inputpath=$2
outputpath=$3
core=$4

#samtobam sorttobam rmdup build index
samtools view \
 --threads ${core} -b ${inputpath}/${sampleID}.sam > ${outputpath}/${sampleID}.bam

# sort bam files
samtools sort \
 --threads ${core} -O bam ${outputpath}/${sampleID}.bam -o ${outputpath}/${sampleID}.sorted.bam

# remove duplicates
samtools rmdup \
 ${outputpath}/${sampleID}.sorted.bam ${outputpath}/${sampleID}.rmdup.bam

# build index file
samtools index \
 -@ ${core} ${outputpath}/${sampleID}.rmdup.bam ${outputpath}/${sampleID}.rmdup.bam.bai

# mapping ratio
samtools flagstat \
 --threads ${core} ${outputpath}/${sampleID}.sorted.bam > ${outputpath}/${sampleID}.sorted.flagstat

# mapping ratio
samtools flagstat \
 --threads ${core} ${outputpath}/${sampleID}.rmdup.bam > ${outputpath}/${sampleID}.rmdup.flagstat

