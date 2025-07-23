#!/usr/bin/bash

#23/7/25 added output folder specification (working_dir)
# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
  echo "Usage: $0 <working_dir> <species> <input_fasta> <seq_type> <plasmid_usage> <completeness>"
  exit 1
fi

# Assign arguments to variables
working_dir="$1"
species="$2"
input_fasta="$3"
seq_type="$4"
plasmid_usage="$5"
completeness="$6"

# Download the genomes and save the compressed fasta into a folder
python3 mops_gen_download_with_error_handling.py -g bacteria -o ./"$working_dir"/ -t4 -n "$species"

# removing non complete genomes
python3 select_genomes.py --completeness "$completeness"
directory="$working_dir/refseq/"
id_file="$working_dir/refseq/non_complete_genomes_ids.txt"
while IFS= read -r id; do
    find "$directory" -type f -name "*$id*" -exec rm -v {} \;
done < "$id_file"


#Creating the genome lists for FastANI
python3 extract_ref.py --input_file "$working_dir"/refseq/metadata.tsv --output_file "$working_dir"/reference_accession.txt
directory="$working_dir/refseq/"
output_file="$working_dir/ref_list.txt"
find "$directory" -type f ! -name "metadata.tsv" > "$output_file"
output_file="$working_dir/ref_genome.txt"
cat "$working_dir"/reference_accession.txt | xargs -I {} find "$directory" -type f -name "*{}*" > "$output_file"

#run fastani type strain to all
fastANI --ql $working_dir/ref_genome.txt --rl $working_dir/ref_list.txt  -o $working_dir/ani_one_to_all_sub.txt
#run fastani all to all for creation of matrix
#fastANI --ql $working_dir/ref_list.txt --rl $working_dir/ref_list.txt --matrix -o $working_dir/ani_all_strains_allVsAll.txt


# to produce graphical output of the fastani data as heatmap (all strains vs all strains)
#ANIclustermap -i ani_step/query_genomes/refseq/ -o clustermap_out --fig_width 20 --fig_height 15 --cmap_colors white,orange,red

# Keeping only the genomes with ANI to reference genome of at least 95
python3 extraction_suitable_genomes.py --input_file "$working_dir"/ani_one_to_all_sub.txt --source_folder "$working_dir"/refseq/ --output_folder "$working_dir"/filtered

# decompressing filtered genomes
DIRECTORY="$working_dir/filtered/"
for file in "$DIRECTORY"*.gz; do
    gunzip "$file"
done

#artificially concatenating incomplete genomes for easier data viz
#added script concatenate_fasta.py to repo
if [ "$completeness" = "incomplete" ]; then
    input_dir="$DIRECTORY"
    output_comb="$working_dir/filtered_combined"
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


python3 suitable_genomes_plasmids.py --input_directory "$working_dir"/filtered --output_directory "$working_dir"/filtered_plasmids --plasmid_setting "$plasmid_usage"

# here to chose if plasmid or not
python3 merge_fasta.py "$working_dir"/filtered_plasmids "$working_dir"/

# Assign number of genomes to a variable
directory="$working_dir/filtered_plasmids"
genome_count=$(ls -1 "$directory" | wc -l)
echo "Summary of Variables:"
echo "Species: $species"
echo "Input FASTA: $input_fasta"
echo "Seq Type: $seq_type"
echo "Plasmid Usage: $plasmid_usage"
echo "Number of filtered genomes: $genome_count"

makeblastdb -in "$working_dir"/dict.fasta -out "$working_dir"/blastdb -dbtype nucl -blastdb_version 4

if [ "$seq_type" = "nucleotide" ]; then
    echo "Running blast for nucleotides"
    blastn -query "$working_dir/${input_fasta}" -db "$working_dir"/blastdb -out "$working_dir"/blast_output.txt -outfmt 6 -num_threads 4 -task blastn -evalue 0.05
elif [ "$seq_type" = "aminoacid" ]; then
    echo "Running blast for proteins"
    tblastn -query "$working_dir/${input_fasta}" -db "$working_dir"/blastdb -out "$working_dir"/blast_output.txt -outfmt 6 -num_threads 4 -evalue 0.05
else
    echo "Invalid seq_type. Please specify either 'nucleotide' or 'aminoacid'."
fi

python3 graphical_out.py --genome_num "$genome_count" --species_name "$species" --blast_out "$working_dir"/blast_output.txt --input_fasta "$working_dir/$input_fasta" --identity_cutoff 70 --coverage_cutoff 70

