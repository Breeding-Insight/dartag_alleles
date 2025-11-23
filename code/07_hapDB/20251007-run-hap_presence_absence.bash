#!/bin/bash

# Check if a file argument is provided
if [ $# -lt 2 ]; then
  echo "Usage: $0 <filename> <log_file>"
  echo "Error: At least 2 arguments are required."
  exit 1
fi

# Assign the first argument to a variable
input_file="$1"
log_file="$2"

SCRIPT='/Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/07_hapDB/util_haplotype_presence_absence_matrix.py'

while IFS= read -r line; do
    echo "Text read from file: $line"
    line_array=($line)
    echo ${stringarray[0]}
    python $SCRIPT ${line_array[0]} ${line_array[1]} --read_threshold ${line_array[2]} >> $log_file
done < "$input_file"
