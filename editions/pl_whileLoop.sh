#!/bin/bash

inputfile=$1

grep -v "#" ${inputfile} | while IFS=$'\t' read -r col1 col2 col3 col4
do
    var1="$col1"
    var2="$col2"
    var3="$col3"
    var4="$col4"

    echo "Sample  : $var1"
    echo "Strategy: $var2"
    echo "Layout  : $var3"
    echo "Genome  : $var4"

done
