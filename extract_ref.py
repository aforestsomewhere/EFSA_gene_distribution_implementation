import pandas as pd
import argparse

def filter_tsv(input_file, output_file):
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')

    # Filter the rows where the 'refseq_category' column contains 'reference genome'
    filtered_df = df[df['refseq_category'] == 'reference genome']

    # If no 'reference genome' is found, search for 'representative genome'
    if filtered_df.empty:
       filtered_df = df[df['refseq_category'] == 'representative genome']

    # Extract only the 'assembly_accession' column
    assembly_accessions = filtered_df['ssembly_accession']

    # Save the result to a text file
    assembly_accessions.to_csv(output_file, sep='\t', index=False, header=False)

    # Print the result
    print(assembly_accessions)

def main():
    parser = argparse.ArgumentParser(description='Filter TSV file by refseq_category column and extract assembly_accession values.')
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input TSV file')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output TXT file')

    args = parser.parse_args()

    filter_tsv(args.input_file, args.output_file)

if __name__ == '__main__':
    main()
