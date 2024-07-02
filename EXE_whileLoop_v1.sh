#!/bin/bash
# version 1


# Input description
# ==============================================================================
: '
1. Create *project_name* folder under main path (Path_main)
2. Scripts used in pipeline would be stored in folder "scripts".
3. sample table files should be stored in folder "input".

variables description:
   Path_main=main path for this analysis project.
   Path_fastq=directory storing FASTQ files.
   sampleTable=directory storing sample Table file (ex: sampleTable.xlsx).'

Path_main=/staging/biology/ls807terra/analysis/*project_name*
Path_fastq=/staging/biology/ls807terra/0_fastq/*project_name*
SampleTable=${Path_main}/input/$(basename ${Path_main}/input/*.xlsx)
SB_prj=MST109178
SB_part=186
SB_core=28
SB_mem=186



# source export genome files and Tools PATH
# ==============================================================================
# Export tools PAHT
source \
 ${Path_main}/scripts/config_PATHtools_exprSource.sh \
 ${Path_main}/scripts/config_PATHtools.sh

# Export variables storing genome files
source ${Path_main}/scripts/config_PATHgenome.sh

# Transfer xlsx to csv
sh ${Path_main}/scripts/config_xlsx2csv.sh ${SampleTable}
grep -v "#" ${SampleTable/.xlsx/.bed}


# Execute pipeline following sample table's information.
# ==============================================================================
grep -v "#" ${SampleTable/.xlsx/.bed} | while IFS=$'\t' read -r c1 c2 c3 c4
do
  # informations should be identical to columns of sample table.
  sampleID="${c1}"
  strategy="${c2}"
  layout="${c3}"
  genome="${c4}"
  pipeline=config_pipeline_${strategy}.sh

  echo -e "\n\n"
  echo "Input Summary:"
  echo "------------------------------------------"
  echo "Sample   : ${sampleID}"
  echo "Strategy : ${strategy}"
  echo "Layout   : ${layout}"
  echo "Genome   : ${genome}"
  echo "pipeline : ${pipeline}"
  echo "------------------------------------------"

  # execute pipeline script:
  sh ${Path_main}/scripts/${pipeline} \
   ${sampleID} \
   ${Path_main} \
   ${Path_fastq} \
   ${SB_prj} \
   ${SB_part} \
   ${SB_core} \
   ${SB_mem} \
   ${genome}

   sleep 0.5s

done
