#!/usr/bin/bash
#22/7/25 added completeness as cmd line arg
#23/7/25 added specification for output folder
# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
  echo "Usage: $0 <working_directory> <species> <input_fasta> <seq_type> <plasmid_usage> <completeness>"
  exit 1
fi

working_directory="$1"
species="$2"
input_fasta="$3"
seq_type="$4"
plasmid_usage="$5"
completeness="$6"

# Check compatibility
if [ "$completeness" = "incomplete" ] && [ "$plasmid_usage" != "remove" ]; then
  echo "Error: When completeness is 'incomplete', plasmid_usage must be 'remove'."
  exit 1
fi

#copy the input fasta file to the output folder (create it if it does not exist already)
mkdir -p "$working_directory"
cp "$input_fasta" "$working_directory"/

./main_script.sh "$working_directory" "$species" "$input_fasta" "$seq_type" "$plasmid_usage" "$completeness"
