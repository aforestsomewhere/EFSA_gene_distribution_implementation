import sys
from Bio import SeqIO
import os

#Added 22/07/25. For incomplete genomes, helper function that concatenates each multifasta to a single fasta sequence, aiding visualisation
if len(sys.argv) != 3:
    print("Usage: python concatenate_fasta.py <input.fna> <output.fna>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Use filename (without extension) as new FASTA header
base_name = os.path.splitext(os.path.basename(input_file))[0]

combined_seq = ""
for record in SeqIO.parse(input_file, "fasta"):
    combined_seq += str(record.seq)

with open(output_file, "w") as f:
    f.write(f">{base_name}\n")
    for i in range(0, len(combined_seq), 60):
        f.write(combined_seq[i:i+60] + "\n")
