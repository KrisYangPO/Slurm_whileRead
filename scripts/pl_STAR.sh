#!/bin/bash

# human CHM13 genome
# basic input
sampleID=$1
inputpath=$2
outputpath=$3
core=$4
species=$5

# index
# Use conditions to determines whether species contains string: mouse/human.
if [[ $species == *'mm10'* ]]; then
  index=${staridx_mm10}

elif [[ $species == *'chm13'* ]]; then
  index=${staridx_chm13}

elif [[ $species == *'hg38'* ]]; then
  index=${staridx_hg38}

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

# There is no need to clarify whether input is paried-end or single-end,
# directly subject the sampleID array (${samples[@]}) to STAR.


# program
STAR \
 --genomeDir ${index} \
 --runMode alignReads \
 --runThreadN ${core} \
 --readFilesCommand zcat \
 --readFilesIn ${samples[@]} \
 --outSAMmapqUnique 255 \
 --outSAMprimaryFlag OneBestScore \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMattributes All \
 --outWigType wiggle \
 --outWigStrand Stranded \
 --outWigNorm RPM \
 --outFileNamePrefix ${outputpath}/${sampleID}_


