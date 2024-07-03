#!/bin/bash
# version 3: optimize sample table reading


# Input description
# ==============================================================================
: '
1. Path_main depneds on the location of your current project folder.
2. Every informations that pipeline required will
   be encompassed in sample_table.xlsx.
3. '


# Setting configurations
# ==============================================================================
# Export tools PAHT

Path_main=$(pwd)
echo -e "Current working directory is: \n"${Path_main}

source \
 ${Path_main}/scripts/config_PATHtools_exprSource.sh \
 ${Path_main}/scripts/config_PATHtools.sh

# Export variables storing genome files
source ${Path_main}/scripts/config_PATHgenome.sh


# accessing prerequisites from sample_table.xlsx
# ==============================================================================
# finding sample_table.xlsx:
SampleTable=${Path_main}/input/$(basename ${Path_main}/input/*.xlsx)

# transfer it from xlsx to bed file format:
sh ${Path_main}/scripts/config_xlsx2csv.sh ${SampleTable}


# refine informations from sample_table.xlsx
# ==============================================================================
Path_fastq=$(sed -n 's/#Path_fastq:\t//1p' ${SampleTable/.xlsx/.bed})
SB_prj=$(sed -n 's/#SB_Project:\t//1p' ${SampleTable/.xlsx/.bed})
SB_part=$(sed -n 's/#SB_part:\t//1p' ${SampleTable/.xlsx/.bed})
SB_core=$(sed -n 's/#SB_core:\t//1p' ${SampleTable/.xlsx/.bed})
SB_mem=$(sed -n 's/SB_mem:\t//1p' ${SampleTable/.xlsx/.bed})

# report summary:
echo -e "\n"
echo "Summary of prerequisites and parameters:"
echo "------------------------------------------"
echo "Fastq path  : ${Path_fastq}"
echo "Project ID  : ${SB_prj}"
echo "Partition   : ${SB_part}"
echo "Core num    : ${SB_core}"
echo "Memory used : ${SB_mem}"
echo "------------------------------------------"
echo "Sample Table content:"
grep -v "#" ${SampleTable/.xlsx/.bed}

sleep 2s


# Execute pipeline by following sample table's informations.
# ==============================================================================
grep -v "#" ${SampleTable/.xlsx/.bed} | while IFS=$'\t' read -r c1 c2 c3 c4
do
  # informations should be identical to columns of sample table.
  sampleID="${c1}"
  strategy="${c2}"
  layout="${c3}"
  genome="${c4}"
  pipeline=config_pipeline_${strategy}.sh

  echo -e "\n"
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

