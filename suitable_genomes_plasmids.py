import os
import argparse
from Bio import SeqIO

def filter_sequences(input_directory, output_directory, keyword, mode):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
    for filename in os.listdir(input_directory):
        if filename.endswith(".fna"):
            input_file = os.path.join(input_directory, filename)
            output_file = os.path.join(output_directory, filename)
            
            # only for debugging
            sequences = SeqIO.parse(input_file, "fasta")
            if any(keyword in seq.description for seq in sequences):
                print(f"Plasmid found in file: {filename}")

            sequences = SeqIO.parse(input_file, "fasta")

            if mode == "remove":
                filtered = (seq for seq in sequences if keyword not in seq.description)
            elif mode == "extract":
                filtered = (seq for seq in sequences if keyword in seq.description)
            else:
                raise ValueError("Mode must be either 'remove' or 'extract'")
            
            SeqIO.write(filtered, output_file, "fasta")
            
            if mode == "extract" and os.path.getsize(output_file) == 0:
                os.remove(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter plasmid sequences from genomes.')
    parser.add_argument('--input_directory', type=str, required=True, help='Input directory containing genome files.')
    parser.add_argument('--output_directory', type=str, required=True, help='Output directory.')
    parser.add_argument('--plasmid_setting', type=str, required=True, choices=['remove', 'extract'], help='Operation mode: "remove" or "extract" plasmids.')

    args = parser.parse_args()

    # Directories
    input_directory = args.input_directory
    output_directory = args.output_directory

    # Call the function with the appropriate mode
    filter_sequences(input_directory, output_directory, "plasmid", args.plasmid_setting)