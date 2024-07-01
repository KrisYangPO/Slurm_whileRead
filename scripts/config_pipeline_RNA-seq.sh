#!/bin/bash
# version 7: RNA-seq pipeline (STAR alignment) from readSampleTable_framework
# pipeline analysis covers: trimgalore, STAR, bamcoverage, and HTseq-count.



# Input parameters:
# ==============================================================================
sampleID=$1
Path_main=$2
Path_fastq=$3
SB_proj=$4
SB_part=$5
SB_core=$6
SB_mem=$7
genome=$8
pipeline_scripts=(pl_trimfastqc.sh pl_STAR.sh pl_bamcoverage.sh pl_HTseqcount.sh)



# functinos:
# ==============================================================================
# A function to create folder
function mkFolder(){
name=$1
if [ ! -d ${Path_main}/${name} ]
then
mkdir ${Path_main}/${name}
fi
}



# report parameters
# ==============================================================================
echo -e "\n"
echo "Processing samples:   "${sampleID}
echo "Pipeline scripts has: "${pipeline_scripts[@]}
echo "SBATCH project:       "${SB_proj}
echo "SBATCH partition:     ""p"${SB_part}"G"
echo "SBATCH core:          "${SB_core}
echo "SBATCH memory:        "${SB_mem}"G"
echo "Species:              "${genome}



# Pipeline
# Step1: trimgalore
# ==============================================================================
# annotation
: '
script: pl_trimfastqc.sh
input:
  1. A sampleID.
  2. fastqPath orientates fastq.
  3. output directory.
  4. core number.'

echo "execute script: ""${pipeline_scripts[0]}"
mkFolder report
mkFolder Step1_output

# submit SBATCH job
A_JID=$(\
 sbatch \
  -A ${SB_proj} \
  -p ngs${SB_part}G \
  -c ${SB_core} \
  --mem=${SB_mem}g \
  -J Step1_${sampleID} \
  -o ${Path_main}/report/Step1_${sampleID}.o.txt \
  -e ${Path_main}/report/Step1_${sampleID}.e.txt \
  ${Path_main}/scripts/${pipeline_scripts[0]} \
  ${sampleID} \
  ${Path_fastq} \
  ${Path_main}/Step1_output \
  ${SB_core})

# prune Job.ID
A_JID=${A_JID/"Submitted batch job "/}
echo "Job ID: "${A_JID}" submitted"



# Step2: STAR
# ==============================================================================
# annotation
: '
script: pl_STAR.sh
input:
  1. A sampleID.
  2. input directory (trimgalore: Step1_output).
  3. output directory.
  4. core number.
  5. genome.'

echo "execute script: ""${pipeline_scripts[1]}"
mkFolder Step2_output

B_JID=$(\
 sbatch \
  -A ${SB_proj} \
  -p ngs${SB_part}G \
  -c ${SB_core} \
  --mem=${SB_mem}g \
  -J Step2_${sampleID} \
  --dependency=afterok:${A_JID} \
  -o ${Path_main}/report/Step2_${sampleID}.o.txt \
  -e ${Path_main}/report/Step2_${sampleID}.e.txt \
  ${Path_main}/scripts/${pipeline_scripts[1]} \
  ${sampleID} \
  ${Path_main}/Step1_output \
  ${Path_main}/Step2_output \
  ${SB_core} \
  ${genome})

# prune Job.ID
B_JID=${B_JID/"Submitted batch job "/}
echo "Job ID: "${B_JID}" submitted"



# Step3: bamcoverage
# ==============================================================================
# annotation
: '
script: pl_bamcoverage.sh
input:
  1. A sampleID
  2. input directory：Step2 (sam2bam output)
  3. output directory：Step2：Step4 (bigwig)
  4. output format: bigwig.
  5. core number.

Use *.rmdup.bam as input, only *.rmdup.bam has *.bai files'

echo "execute script: ""${pipeline_scripts[2]}"
mkFolder Step3_output

C_JID=$(\
 sbatch \
  -A ${SB_proj} \
  -p ngs${SB_part}G \
  -c ${SB_core} \
  --mem=${SB_mem}g \
  -J Step3_${sampleID} \
  --dependency=afterok:${B_JID} \
  -o ${Path_main}/report/Step3_${sampleID}.o.txt \
  -e ${Path_main}/report/Step3_${sampleID}.e.txt \
  ${Path_main}/scripts/${pipeline_scripts[2]} \
  ${sampleID} \
  ${Path_main}/Step2_output \
  ${Path_main}/Step3_output \
  bigwig \
  ${SB_core})

# prune Job.ID
C_JID=${C_JID/"Submitted batch job "/}
echo "Job ID: "${C_JID}" submitted"



# Step4: HTSeq-count
# ==============================================================================
# annotation
: '
script: pl_HTseqcount.sh
input:
  1. input directory (STAR output)
  2. output directory.
  3. core number.
  4. genome'

echo "execute script: ""${pipeline_scripts[3]}"
mkFolder Step4_output

D_JID=$(\
 sbatch \
  -A ${SB_proj} \
  -p ngs${SB_part}G \
  -c ${SB_core} \
  --mem=${SB_mem}g \
  -J Step4_${sampleID} \
  --dependency=afterok:${C_JID} \
  -o ${Path_main}/report/Step4_${sampleID}.o.txt \
  -e ${Path_main}/report/Step4_${sampleID}.e.txt \
  ${Path_main}/scripts/${pipeline_scripts[3]} \
  ${sampleID} \
  ${Path_main}/Step2_output \
  ${Path_main}/Step4_output \
  ${SB_core} \
  ${genome})

# prune Job.ID
D_JID=${D_JID/"Submitted batch job "/}
echo "Job ID: "${D_JID}" submitted"



#
