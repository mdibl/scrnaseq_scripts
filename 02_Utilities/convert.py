import sys
import os
import gzip
import csv
import re
import shutil

def process_features(input_file):
    temp_file = input_file + '.temp'
    with gzip.open(input_file, 'rt') as infile, gzip.open(temp_file, 'wt') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                if fields[0] != fields[1]:
                    outfile.write(f"{fields[0]}\t{fields[1]}::{fields[0]}\tExpression\n")
                else:
                    outfile.write(f"{fields[0]}\t{fields[1]}\tExpression\n")
    
    # Replace the original file with the processed file
    shutil.move(temp_file, input_file)
    print(f"Processed and updated features file: {input_file}")

def process_aux_list(aux_file, features_file, output_file, column_index):
    gene_list = []
    with open(aux_file, 'r') as infile:
        reader = csv.reader(infile)
        next(reader)  # Skip header
        for row in reader:
            if len(row) > column_index and row[column_index].strip():
                gene_list.append(re.escape(row[column_index].strip()))
    
    print(f"Found {len(gene_list)} genes in column {column_index} of {aux_file}")

    # Create a single regex pattern for all genes
    gene_pattern = re.compile(r'^(' + '|'.join(gene_list) + r')::')

    matching_genes = []
    with gzip.open(features_file, 'rt') as features:
        for line in features:
            fields = line.strip().split('\t')
            if len(fields) >= 2 and gene_pattern.match(fields[1]):
                matching_genes.append(fields[1])

    print(f"Found {len(matching_genes)} matching genes for {output_file}")

    with open(output_file, 'w') as outfile:
        outfile.write(f"{os.path.splitext(output_file)[0]}genes\n")
        outfile.write("\n".join(matching_genes))

def main(feature_file_path, aux_file_path):
    # Process features (this will now overwrite the original file)
    process_features(feature_file_path)

    # Process auxiliary gene lists
    for i, name in enumerate(['MT', 'G2M', 'S', 'RM']):
        process_aux_list(aux_file_path, feature_file_path, f'{name}.csv', i)

    # Combine results
    with open('AuxGeneList.csv', 'w') as outfile:
        files = ['MT.csv', 'G2M.csv', 'S.csv', 'RM.csv']
        contents = [open(f).readlines() for f in files]
        max_length = max(len(content) for content in contents)
        
        for i in range(max_length):
            line = []
            for content in contents:
                line.append(content[i].strip() if i < len(content) else '')
            outfile.write(','.join(line) + '\n')

    print("Created AuxGeneList.csv")

    # Clean up intermediate files
    for file in ['MT.csv', 'G2M.csv', 'S.csv', 'RM.csv']:
        os.remove(file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert.py <feature_file_path> <aux_file_path>")
        sys.exit(1)
    
    feature_file_path = sys.argv[1]
    aux_file_path = sys.argv[2]
    main(feature_file_path, aux_file_path)