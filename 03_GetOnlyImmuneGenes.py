import os
import pandas as pd


def filter_bed_files(input_dir, output_dir, immune_genes_file):
    """
    Filters all BED files in the input directory, keeping rows where any immune gene is found
    in the last three columns, considering semicolon-separated entries.
    The filtered files are saved to the output directory.
    """
    # Read and clean the immune genes list
    immune_genes = pd.read_csv(immune_genes_file, sep='\t', header=None)
    immune_genes = immune_genes[immune_genes[2] != '']  # Remove rows where column 2 is empty
    immune_genes = immune_genes.dropna()  # Drop rows with NaN
    immune_gene_set = set(immune_genes[2].tolist())  # Convert to a hash set for faster lookup

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Process each BED file in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.bed'):
            input_file_path = os.path.join(input_dir, file_name)
            output_file_path = os.path.join(output_dir, file_name)

            # Read the BED file
            bed_data = pd.read_csv(input_file_path, sep='\t')

            # Initialize a list to hold rows that match
            filtered_rows = []

            # Iterate through each row in the BED file
            for _, row in bed_data.iterrows():
                # Check the last three columns for any immune gene
                for col in row.iloc[-3:]:
                    if not pd.isna(col):  # Skip NaN values
                        genes = col.split(';')  # Split by semicolon
                        if any(gene in immune_gene_set for gene in genes):
                            filtered_rows.append(row)
                            break  # Move to the next row once a match is found

            # Convert the filtered rows back to a DataFrame
            filtered_data = pd.DataFrame(filtered_rows)

            # Save the filtered data
            filtered_data.to_csv(output_file_path, sep='\t', index=False)


# Paths for input directory, output directory, and immune genes file
immune_genes_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\immune_genes.tsv"  # Replace with the path to your immune genes file
input_directory = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\02_c_AlmostGASmotifs_samples_annotated"  # Replace with the path to your input BED files
output_directory = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\02_d_AlmostGASmotifs__samples_annotated_immuneGenes"  # Replace with the path for filtered files

# Run the function
filter_bed_files(input_directory, output_directory, immune_genes_file)
