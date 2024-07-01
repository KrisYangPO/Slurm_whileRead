#!/bin/bash

# parameters
sampleID=$1
inputpath=$2
output=$3
outputformat=$4
core=$5

# program
bamCoverage \
 --normalizeUsing CPM \
 --outFileFormat ${outputformat} \
 --centerReads \
 --extendReads 150 \
 -p ${core} \
 -b ${inputpath}/${sampleID}.rmdup.bam \
 -o ${output}/${sampleID}.${outputformat}

