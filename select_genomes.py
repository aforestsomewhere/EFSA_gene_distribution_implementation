import pandas as pd

# Load the TSV file
file_path = 'output/refseq/metadata.tsv'
metadata_df = pd.read_csv(file_path, sep='\t')

# Filter the dataframe to get rows where assembly_level is not 'Complete Genome'
incomplete_genomes_df = metadata_df[metadata_df['assembly_level'] != 'Complete Genome']
complete_genomes_df = metadata_df[metadata_df['assembly_level'] == 'Complete Genome']



# Write the non-complete genome IDs to a text file
output_txt_path = 'output/refseq/non_complete_genomes_ids.txt'
incomplete_genomes_df['ssembly_accession'].to_csv(output_txt_path, sep='\t', index=False, header=False)


# Save the new metadata file containing only complete genomes
new_metadata_file_path = 'output/refseq/metadata.tsv'
complete_genomes_df.to_csv(new_metadata_file_path, sep='\t', index=False)