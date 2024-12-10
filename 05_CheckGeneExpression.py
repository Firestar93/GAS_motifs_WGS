import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# Function to create boxplots and perform statistical tests
def create_boxplots_and_tests(deseq2_matrix, bed_file, output_path):
    results = []  # To store statistical test results

    # Ensure output path exists
    os.makedirs(output_path, exist_ok=True)

    # Iterate over rows in the BED file
    for _, row in bed_file.iterrows():
        genes = [
            row["Inside_Genes"], row["Near_Genes_<2500bp"], row["Near_Genes_2500-10000bp"]
        ]
        vcf_samples = row["vcf_files"].split(";")
        vcf_samples = [sample.split("_")[1].split(".")[0] for sample in vcf_samples]

        # Separate LUDx samples into with and without GAS motif groups

        ##################################################################
        #### THIS NEEDS TO BE EXCHANGED THE OTHER WAY ROUND ALL THE TIME!!
        ##################################################################
        without_gas = [col for col in deseq2_matrix.columns if any(f"LUD{x}" in col for x in vcf_samples)]
        with_gas = [col for col in deseq2_matrix.columns if col not in without_gas]

        if len(without_gas) == 0:
            continue

        # Generate boxplots for each gene
        for gene in genes:
            if pd.isna(gene):  # Skip if gene column is empty
                continue

            if gene not in deseq2_matrix.index:
                continue

            gene_data = deseq2_matrix.loc[gene]
            with_gas_data = gene_data[with_gas].values if with_gas else []
            without_gas_data = gene_data[without_gas].values if without_gas else []

            # Perform statistical test
            if len(with_gas_data) > 0 and len(without_gas_data) > 0:
                stat, p_value = mannwhitneyu(with_gas_data, without_gas_data, alternative='two-sided')
            else:
                stat, p_value = None, None

            # Plotting
            plt.figure(figsize=(8, 6))
            sns.boxplot(data=[with_gas_data, without_gas_data])
            plt.xticks([0, 1], ["With GAS motif", "Without GAS motif"])
            title = (
                f"hg38 position {row['chrom']}:{row['start']}-{row['end']}\n"
                f"Ref.: {row['ref_allele']}; Alt.: {row['alt_allele']} at nucleotide position: {row['position_nucleotide']} and motif position: {row['position_in_GAS']}\n"
                f"{row['sequence']} → {row['modified_sequence']}\n"
                f"Gene: {gene} (Source: {'Inside_Genes' if gene == row['Inside_Genes'] else 'Near_Genes_<2500bp' if gene == row['Near_Genes_<2500bp'] else 'Near_Genes_2500-10000bp'})\n"
                f"p-value: {'NA' if p_value is None else f'{p_value:.3e}'}"
            )
            plt.title(title)
            plt.ylabel("Expression")
            plt.xlabel("Group")
            plt.tight_layout()

            # Save plot
            plot_filename = f"{gene}_Source-{'Inside_Genes' if gene == row['Inside_Genes'] else 'Near_Genes_2500bp' if gene == row['Near_Genes_<2500bp'] else 'Near_Genes_2500-10000bp'}.png"
            plt.savefig(os.path.join(output_path, plot_filename))
            plt.close()

            # Store result
            results.append({
                "Gene": gene,
                "Sequence": row["sequence"],
                "With_GAS_Size": len(with_gas_data),
                "Without_GAS_Size": len(without_gas_data),
                "p-value": p_value
            })

    return pd.DataFrame(results)

# Load the DESeq2 normalized matrix and the BED file
deseq2_matrix_path = 'C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\03_RNA-seq\\Others\\Others_DESeq2_normalized.tsv'
bed_file_path = 'C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\01_d_AllGASmotif_samples_annotated_immuneGenes\\annotated_others.bed'
output_file_path = 'C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\01_e_PLOTS_AllGASmotif_samples_annotated_immuneGenes\\Others\\plots'

# Read the DESeq2 normalized matrix and the BED file
deseq2_matrix = pd.read_csv(deseq2_matrix_path, sep="\t", index_col=0)
bed_file = pd.read_csv(bed_file_path, sep="\t")

# Rename BED columns for easier handling
bed_file.columns = [
    "chrom", "start", "end", "sequence", "sequence_type", "vcf_files",
    "modified_sequence", "position_chromosome", "position_nucleotide",
    "ref_allele", "alt_allele", "position_in_GAS", "Inside_Genes",
    "Near_Genes_<2500bp", "Near_Genes_2500-10000bp"
]

# Execute the function and store results
results = create_boxplots_and_tests(deseq2_matrix, bed_file, output_file_path)

# Save results to a CSV file
results.to_csv(os.path.join(output_file_path, "boxplot_statistical_results.csv"), index=False)
