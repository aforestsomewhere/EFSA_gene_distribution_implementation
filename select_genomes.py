import pandas as pd
#22/07/25
import os
import argparse

#rewrite as a function
def select_genomes_fun(completeness, output):
  # Load the TSV file
  file_path = '$output/refseq/metadata.tsv'
  metadata_df = pd.read_csv(file_path, sep='\t')
  os.makedirs('$output/refseq', exist_ok=True)

  # Select genomes based on completeness
  if completeness == "complete":
    selected_df = metadata_df[metadata_df['assembly_level'] == 'Complete Genome']
    incomplete_genomes_df = metadata_df[metadata_df['assembly_level'] != 'Complete Genome']
    output_txt_path = '$output/refseq/non_complete_genomes_ids.txt'
    incomplete_genomes_df['ssembly_accession'].to_csv(output_txt_path, sep='\t', index=False, header=False)
  elif completeness == "incomplete":
    selected_df = metadata_df  # keep all
    #write empty csv of "non-complete" genomes
    with open('$output/refseq/non_complete_genomes_ids.txt', 'w') as f:
      pass
  else:
      raise ValueError("completeness must be either 'complete' or 'incomplete'")

def main():
    parser = argparse.ArgumentParser(description='Select genomes to include, based on completeness.')
    parser.add_argument('--completeness', type=str, required=True, help='complete or incomplete')
    parser.add_argument('--output', type=str, required=True, help='output folder')
    args = parser.parse_args()

    select_genomes_fun(args.completeness, args.output)

if __name__ == '__main__':
    main()
