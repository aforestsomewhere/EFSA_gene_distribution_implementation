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


docker run -v "$working_directory":/workspace/output famr:1.0 \
  bash -c "./main_script.sh \"$species\" \"$input_fasta\" \"$seq_type\" \"$plasmid_usage\""