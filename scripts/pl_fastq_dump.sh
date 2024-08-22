#!/bin/bash

sampleID=$1
inputpath=$2
outputpath=$3

# grep sample by their sampleID.
cd ${inputpath}
sample=$(ls ${sampleID}*.sra)

fastq-dump ${sample} --split-files --gzip -O ${outputpath}

