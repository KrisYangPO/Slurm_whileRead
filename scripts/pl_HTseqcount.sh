#!/bin/bash
# htseq human genome

# parameters
sampleID=$1
inputpath=$2
outputpath=$3
core=$4
species=$5


# index:
# Use conditions to determines whether species contains string: mouse/human.
if [[ $species == *'mm10'* ]]; then
  gtf=${gtf_mm10}

elif [[ $species == *'chm13'* ]]; then
  gtf=${gtf_chm13}

elif [[ $species == *'hg38'* ]]; then
  gtf=${gtf_hg38}

fi


# grep all the bam files and perform htseq-count
cd ${inputpath}
sample=$(ls ${sampleID}_*.out.bam)
samples=($(ls *.out.bam))


# Run HTseq-count
htseq-count -m intersection-nonempty --nonunique all \
 -f bam -s reverse \
 -t exon \
 --idattr gene_name \
 -n ${core} \
 ${sample} \
 ${gtf} \
 > ${outputpath}/HTseq_${sampleID}.txt


# Run HTseq-count and collect all samples into a table
htseq-count -m intersection-nonempty --nonunique all \
 -f bam -s reverse \
 -t exon \
 --idattr gene_name \
 -n ${core} \
 ${samples[@]} \
 ${gtf} \
 > ${outputpath}/HTseq_ALL.txt


echo "Sample: "${sampleID}
echo "Species: "${species}
echo "GTF file: "$(basename ${gtf})
