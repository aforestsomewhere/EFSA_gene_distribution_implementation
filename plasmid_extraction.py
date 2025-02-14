import os
from Bio import SeqIO

def filter_sequences(input_directory, output_directory, keyword, mode):
    for filename in os.listdir(input_directory):
        if filename.endswith(".fna"):
            input_file = os.path.join(input_directory, filename)
            output_file = os.path.join(output_directory, filename)
            sequences = SeqIO.parse(input_file, "fasta")
            
            if mode == "remove":
                filtered = (seq for seq in sequences if keyword not in seq.description)
            elif mode == "extract":
                filtered = (seq for seq in sequences if keyword in seq.description)
            
            SeqIO.write(filtered, output_file, "fasta")
            
            if mode == "extract" and os.path.getsize(output_file) == 0:
                os.remove(output_file)

if __name__ == "__main__":
    # User-defined variable
    plasmid_setting = "extract"  # or "remove"

    # Directories
    input_directory = "suitable_query_genomes"
    output_directory_remove = "./noplasm"
    output_directory_extract = "./plasm"

    # Determine the appropriate output directory based on the setting
    output_directory = output_directory_remove if plasmid_setting == "remove" else output_directory_extract

    # Call the function with the appropriate mode
    filter_sequences(input_directory, output_directory, "plasmid", plasmid_setting)