import os
import pandas as pd
from intervaltree import IntervalTree

from datetime import datetime

def print_with_time(message):
    """
    Prints a message with the current timestamp.
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")

def load_gene_intervals_from_gtf(gtf_file):
    """
    Load gene intervals from a GTF file into an interval tree.
    """
    tree = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip header lines
            cols = line.strip().split('\t')
            if cols[2] != "gene":
                continue  # We only care about gene entries

            chrom = cols[0]  # Chromosome
            start = int(cols[3]) - 1  # GTF is 1-based; convert to 0-based
            end = int(cols[4])  # End position is inclusive
            attributes = cols[8]  # Attributes field

            # Extract gene_name from attributes
            gene_name = ""
            for attr in attributes.split(";"):
                if attr.strip().startswith("gene_name"):
                    gene_name = attr.split('"')[1]
                    break

            if chrom not in tree:
                tree[chrom] = IntervalTree()
            tree[chrom][start:end] = gene_name
    return tree


def find_nearby_genes(tree, chrom, start, end, distance=2500):
    """
    Find genes overlapping and near a genomic interval.
    """
    nearby_genes = set()
    if chrom in tree:
        # Find overlapping genes
        for interval in tree[chrom][start:end]:
            nearby_genes.add(interval.data)

        # Find genes within distance
        for interval in tree[chrom][start - distance:end + distance]:
            nearby_genes.add(interval.data)
    return nearby_genes


def annotate_bed_with_local_data(input_dir, output_dir, gene_tree):
    os.makedirs(output_dir, exist_ok=True)

    for file_name in os.listdir(input_dir):
        if file_name.endswith(".bed"):
            input_file_path = os.path.join(input_dir, file_name)
            output_file_path = os.path.join(output_dir, f"annotated_{file_name}")

            bed_data = pd.read_csv(input_file_path, sep='\t', header=0)
            inside_genes = []
            near_genes_2500 = []
            near_genes_10000 = []

            counter = 0

            for _, row in bed_data.iterrows():
                chrom = row['chrom']
                start = int(row['start'])
                end = int(row['end'])

                counter = counter + 1
                if counter % 100 == 0:
                    print_with_time("Processed " + str(counter) + " SNPs.")

                # Genes overlapping the SNP
                inside = find_nearby_genes(gene_tree, chrom, start, end, distance=0)
                inside_genes.append(";".join(inside))

                # Genes within 2,500 bp
                near_2500 = find_nearby_genes(gene_tree, chrom, start, end, distance=2500)
                near_genes_2500.append(";".join(near_2500 - inside))

                # Genes within 10,000 bp
                near_10000 = find_nearby_genes(gene_tree, chrom, start, end, distance=10000)
                near_genes_10000.append(";".join(near_10000 - near_2500))

            # Add new columns
            bed_data["Inside_Genes"] = inside_genes
            bed_data["Near_Genes_<2500bp"] = near_genes_2500
            bed_data["Near_Genes_2500-10000bp"] = near_genes_10000

            bed_data.to_csv(output_file_path, sep='\t', index=False)
            print(f"Annotated file saved: {output_file_path}")


# Example usage
gtf_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\hg38.refGene"
input_directory = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\GeneNameTest"
output_directory = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\GeneNameTest_output"

print_with_time("Start building tree.")
# Load the gene annotation data into an interval tree
gene_tree = load_gene_intervals_from_gtf(gtf_file)
print_with_time("Finished building tree.")

print_with_time("Start searching for genes.")
# Annotate BED files using the interval tree
annotate_bed_with_local_data(input_directory, output_directory, gene_tree)
