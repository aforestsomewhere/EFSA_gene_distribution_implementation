#!/usr/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <working_directory> <species> <input_fasta> <seq_type> <plasmid_usage>"
  exit 1
fi

working_directory="$1"
species="$2"
input_fasta="$3"
seq_type="$4"
plasmid_usage="$5"

#copy the input fasta file to the output folder (create it if it does not exist already)
mkdir -p output
cp "$input_fasta" output/

./main_script.sh "$species" "$input_fasta" "$seq_type" "$plasmid_usage"
