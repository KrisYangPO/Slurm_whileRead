#!/bin/bash

inputname=$1
outputname=${inputname/.xlsx/.bed}

xlsx2csv.py \
 --delimiter 'tab' \
 ${inputname} \
 ${outputname}
