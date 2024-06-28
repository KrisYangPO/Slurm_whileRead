#!/bin/bash

inputname=sampleFormula_test.xlsx
outputname=sampleFormula_test.bed

/staging/biology/ls807terra/0_Programs/xlsx2csv-master/xlsx2csv.py \
 --delimiter 'tab' \
 ${inputname} \
 ${outputname}

