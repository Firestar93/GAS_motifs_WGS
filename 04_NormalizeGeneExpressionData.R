# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
library(DESeq2)

# Input directory containing gene expression files
input_dir <- "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\03_RNA-seq\\Others\\Raw"  # Replace with your input directory
output_file <- "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\03_RNA-seq\\Others\\Others_DESeq2_normalized.tsv"  # Replace with your desired output file path

# List all gene expression files in the input directory
file_list <- list.files(input_dir, full.names = TRUE, pattern = "\\.txt$")

# Function to read and prepare data from a file
read_gene_expression <- function(file) {
  data <- read.table(file, header = TRUE, row.names = 1, sep = "\t")
  colnames(data) <- basename(file)  # Use the filename as the column name
  return(data)
}

# Read all files and combine them into a single data frame
expression_data <- lapply(file_list, read_gene_expression)
combined_counts <- do.call(cbind, expression_data)

# Create metadata for DESeq2 (dummy data if no experimental design is provided)
sample_names <- colnames(combined_counts)
col_data <- data.frame(
  sample = sample_names,
  condition = "unknown"  # Placeholder; modify as needed for real experimental data
)
rownames(col_data) <- sample_names

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                              colData = col_data,
                              design = ~1)  # No design formula since we're only normalizing

# Perform normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Save normalized counts to the output file
write.table(normalized_counts, file = output_file, sep = "\t", quote = FALSE)

cat("Normalized gene expression data has been saved to:", output_file, "\n")
