#!/usr/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <species> <input_fasta> <seq_type> <plasmid_usage> <completeness>"
  exit 1
fi

# Assign arguments to variables
species="$1"
input_fasta="$2"
seq_type="$3"
plasmid_usage="$4"
completeness="$5"

# Download the genomes and save the compressed fasta into a folder
python3 mops_gen_download_with_error_handling.py -g bacteria -o ./output/ -t4 -n "$species"

# removing non complete genomes
python3 select_genomes.py --completeness "$completeness"
directory="output/refseq/"
id_file="output/refseq/non_complete_genomes_ids.txt"
while IFS= read -r id; do
    find "$directory" -type f -name "*$id*" -exec rm -v {} \;
done < "$id_file"


#Creating the genome lists for FastANI
python3 extract_ref.py --input_file output/refseq/metadata.tsv --output_file output/reference_accession.txt
directory=output/refseq/
output_file="output/ref_list.txt"
find "$directory" -type f ! -name "metadata.tsv" > "$output_file"
output_file="output/ref_genome.txt"
cat output/reference_accession.txt | xargs -I {} find "$directory" -type f -name "*{}*" > "$output_file"

#run fastani type strain to all
fastANI --ql output/ref_genome.txt --rl output/ref_list.txt  -o output/ani_one_to_all_sub.txt
#run fastani all to all for creation of matrix
#fastANI --ql output/ref_list.txt --rl output/ref_list.txt --matrix -o output/ani_all_strains_allVsAll.txt


# to produce graphical output of the fastani data as heatmap (all strains vs all strains)
#ANIclustermap -i ani_step/query_genomes/refseq/ -o clustermap_out --fig_width 20 --fig_height 15 --cmap_colors white,orange,red

# Keeping only the genomes with ANI to reference genome of at least 95
python3 extraction_suitable_genomes.py --input_file output/ani_one_to_all_sub.txt --source_folder output/refseq/ --output_folder output/filtered

# decompressing filtered genomes
DIRECTORY="output/filtered/"
for file in "$DIRECTORY"*.gz; do
    gunzip "$file"
done

#22/07/25 artificially concatenating incomplete genomes for easier data viz
#added script concatenate_fasta.py to repo
if [ "$completeness" = "incomplete" ]; then
    input_dir="$DIRECTORY"
    output_comb="output/filtered_combined"
    mkdir -p "$output_comb"

    for file in "$input_dir"/*.fna; do
	base=$(basename "$file" .fna)
	output="$output_comb/${base}.fna"
	python concatenate_fasta.py "$file" "$output"
    done

  # Remove original .fna files
  rm "$input_dir"/*.fna
  mv "$output_comb"/*.fna "$input_dir/"
fi

python3 suitable_genomes_plasmids.py --input_directory output/filtered --output_directory output/filtered_plasmids --plasmid_setting "$plasmid_usage"

# here to chose if plasmid or not
python3 merge_fasta.py output/filtered_plasmids output/

# Assign number of genomes to a variable
directory="output/filtered_plasmids"
genome_count=$(ls -1 "$directory" | wc -l)
echo "Summary of Variables:"
echo "Species: $species"
echo "Input FASTA: $input_fasta"
echo "Seq Type: $seq_type"
echo "Plasmid Usage: $plasmid_usage"
echo "Number of filtered genomes: $genome_count"


makeblastdb -in output/dict.fasta -out output/blastdb -dbtype nucl -blastdb_version 4

if [ "$seq_type" = "nucleotide" ]; then
    echo "Running blast for nucleotides"
    blastn -query "output/${input_fasta}" -db output/blastdb -out output/blast_output.txt -outfmt 6 -num_threads 4 -task blastn -evalue 0.05
elif [ "$seq_type" = "aminoacid" ]; then
    echo "Running blast for proteins"
    tblastn -query "output/${input_fasta}" -db output/blastdb -out output/blast_output.txt -outfmt 6 -num_threads 4 -evalue 0.05
else
    echo "Invalid seq_type. Please specify either 'nucleotide' or 'aminoacid'."
fi

python3 graphical_out.py --genome_num "$genome_count" --species_name "$species" --blast_out output/blast_output.txt --input_fasta "output/$input_fasta" --identity_cutoff 70 --coverage_cutoff 70

