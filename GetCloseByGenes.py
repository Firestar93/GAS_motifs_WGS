import os
import pandas as pd
import requests
from datetime import datetime
import urllib3
import time

def print_with_time(message):
    """
    Prints a message with the current timestamp.
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")

def annotate_bed_with_genes(input_dir, output_dir):
    # Define the Ensembl API endpoint
    ensembl_api_url = "https://rest.ensembl.org/overlap/region/human"

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Loop through all BED files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".bed"):
            input_file_path = os.path.join(input_dir, file_name)
            output_file_path = os.path.join(output_dir, f"annotated_{file_name}")

            # Load the BED file
            bed_data = pd.read_csv(input_file_path, sep='\t', header=0)

            # Create new columns for gene annotations
            inside_genes = []
            near_genes_2500 = []
            near_genes_10000 = []

            count = 1
            if count % 100 == 0:
                print_with_time("Checked " + str(count) + " SNPs.")

            for _, row in bed_data.iterrows():
                count = count + 1
                chrom = row['chrom']
                start = int(row['start'])
                end = int(row['end'])

                # Query genes overlapping this region
                response = requests.get(
                    f"{ensembl_api_url}/{chrom}:{start}-{end}?feature=gene",
                    headers={"Content-Type": "application/json"}
                )

                #print_with_time(response)

                if response.ok:
                    try:
                        genes = [gene.get('external_name', 'Unknown') for gene in response.json()]
                        inside_genes.append(";".join(genes))
                    except Exception as e:
                        print_with_time(f"Error processing inside genes for {chrom}:{start}-{end}: {e}")
                        inside_genes.append("")
                else:
                    inside_genes.append("")

                # Query genes within 2,500 bp
                response_2500 = requests.get(
                    f"{ensembl_api_url}/{chrom}:{start - 2500}-{end + 2500}?feature=gene",
                    headers={"Content-Type": "application/json"}
                )
                if response_2500.ok:
                    try:
                        genes_2500 = [gene.get('external_name', 'Unknown') for gene in response_2500.json()]
                        near_genes_2500.append(";".join(genes_2500))
                    except Exception as e:
                        print_with_time(f"Error processing near genes (2,500 bp) for {chrom}:{start}-{end}: {e}")
                        near_genes_2500.append("")
                else:
                    near_genes_2500.append("")

                # Query genes within 10,000 bp
                response_10000 = requests.get(
                    f"{ensembl_api_url}/{chrom}:{start - 10000}-{end + 10000}?feature=gene",
                    headers={"Content-Type": "application/json"}
                )
                if response_10000.ok:
                    try:
                        genes_10000 = [gene.get('external_name', 'Unknown') for gene in response_10000.json()]
                        near_genes_10000.append(";".join(genes_10000))
                    except Exception as e:
                        print_with_time(f"Error processing near genes (10,000 bp) for {chrom}:{start}-{end}: {e}")
                        near_genes_10000.append("")
                else:
                    near_genes_10000.append("")


                #if len(inside_genes) > 1 or len(genes_2500) > 1 or len(genes_10000) > 1:
                    #print_with_time("Found genes!")
                    #print_with_time(inside_genes)
                    #print_with_time(genes_2500)
                    #print_with_time(genes_10000)

                time.sleep(0.5)

            # Add the new columns to the DataFrame
            bed_data["Inside_Genes"] = inside_genes
            bed_data["Near_Genes_<2500bp"] = near_genes_2500
            bed_data["Near_Genes_2500-10000bp"] = near_genes_10000

            # Save the annotated file
            bed_data.to_csv(output_file_path, sep='\t', index=False)
            print(f"Annotated file saved: {output_file_path}")


# Specify the input and output directories
input_directory = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\GeneNameTest"  # Replace with the path to your input directory
output_directory = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\GeneNameTest_output"  # Replace with the path to your output directory

# Annotate all BED files in the input directory
annotate_bed_with_genes(input_directory, output_directory)
