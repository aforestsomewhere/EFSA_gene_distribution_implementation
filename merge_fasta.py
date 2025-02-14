import os
import sys
from Bio import SeqIO
from Bio.SeqIO import write

def process_genomic_sequences(genome_folder):
    # List all files in the genome folder
    genome_files = [os.path.join(genome_folder, f) for f in os.listdir(genome_folder) if os.path.isfile(os.path.join(genome_folder, f))]
    print(f"Processing files: {genome_files}")

    # Load the genomic sequences into a dictionary
    genome_seqs = {}
    for genome_file in genome_files:
        for seq in SeqIO.parse(genome_file, "fasta"):
            genome_seqs[seq.id] = seq
    
    return genome_seqs

def write_sequences_to_file(genome_seqs, output_file):
    # Write genomic sequences to the specified output file
    SeqIO.write(list(genome_seqs.values()), output_file, "fasta")
    print(f"Sequences written to {output_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_folder> <output_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    
    if not os.path.isdir(input_folder):
        print(f"Error: {input_folder} is not a valid directory")
        sys.exit(1)

    # Process genomic sequences from the input folder
    genome_seqs = process_genomic_sequences(input_folder)
    output_file = os.path.join(output_folder, "dict.fasta")
    write_sequences_to_file(genome_seqs, output_file)


if __name__ == "__main__":
    main()