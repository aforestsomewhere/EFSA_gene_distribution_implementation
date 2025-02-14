import shutil
import os
import argparse

def extract_suitable_genomes(input_file, source_folder, output_folder):
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Initialize the set to track unique values
    unique_values = set()

    # Open and read the input file
    with open(input_file, "r") as file:
        for line in file:
            columns = line.strip().split('\t')  # tab-separated columns

            # Ensure there are at least 3 columns in the line
            if len(columns) >= 3:
                file_id = columns[1]  # Assuming the file ID is in the second column
                value = float(columns[2])  # Assuming the quantitative value is in the third column

                # Remove duplicates from the third column
                if value in unique_values:
                    continue
                unique_values.add(value)

                # Check if the value is higher than or equal to 95
                if value >= 95:
                    # Remove the "query_genomes/" prefix from the file_id
                    file_id = file_id.replace(source_folder, '')
                    
                    # Construct the source and destination paths
                    source_path = os.path.join(source_folder, file_id)
                    destination_path = os.path.join(output_folder, file_id)

                    # Check if the source file exists
                    if os.path.exists(source_path):
                        # Copy the file to the output folder
                        shutil.copy(source_path, destination_path)
                    else:
                        print(f"File not found: {source_path}")

def main():
    parser = argparse.ArgumentParser(description='Extract suitable genomes based on % identity.')
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input file (e.g., ani_one_to_all_sub.txt)')
    parser.add_argument('--source_folder', type=str, required=True, help='Path to the folder containing the query genomes')
    parser.add_argument('--output_folder', type=str, required=True, help='Path to the folder to store suitable genomes')

    args = parser.parse_args()

    extract_suitable_genomes(args.input_file, args.source_folder, args.output_folder)

if __name__ == '__main__':
    main()