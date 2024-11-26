import pandas as pd
import os
import glob
import numpy as np
from collections import defaultdict

def build_bed_intervals(bed_file):
    """
    Reads the BED file and creates a dictionary of intervals for fast lookup.
    Uses NumPy arrays for faster iteration.
    """
    # Load the BED file
    bed_df = pd.read_csv(
        bed_file, sep="\t", header=None, names=["chrom", "start", "end", "sequence"]
    )
    # Process the 'sequence' column
    bed_df["sequence_type"] = bed_df["sequence"]
    bed_df["sequence"] = bed_df["sequence"].str.split('-').str[1]
    # Initialize the dictionary
    bed_dict = defaultdict(list)

    # Convert DataFrame columns to NumPy arrays
    chrom_array = bed_df["chrom"].astype(str).to_numpy()
    start_array = bed_df["start"].to_numpy()
    end_array = bed_df["end"].to_numpy()
    sequence_array = bed_df["sequence"].to_numpy()
    sequence_type_array = bed_df["sequence_type"].to_numpy()
    index_array = bed_df.index.to_numpy()

    # Iterate over the NumPy arrays
    for chrom, start, end, sequence, sequence_type, idx in zip(
        chrom_array, start_array, end_array, sequence_array, sequence_type_array, index_array
    ):
        bed_dict[chrom].append(
            {
                "start": start,
                "end": end,
                "sequence": sequence,
                "sequence_type" : sequence_type_array,
                "idx": idx,
            }
        )

    return bed_df, bed_dict

def parse_vcf_once(vcf_file, bed_dict, bed_df):
    """
    Processes a VCF file and updates BED entries with matching VCF file names.
    """
    counter = 0
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip header lines
            fields = line.strip().split('\t')
            vcf_chrom, vcf_pos, ref, alts = (
                str(fields[0]),
                int(fields[1]),
                fields[3],
                fields[4],
            )
            vcf_chrom = vcf_chrom.replace('chr', '')

            counter = counter + 1

            if counter % 1000 == 0:
                print("Checked " + str(counter))

            if vcf_chrom in bed_dict:
                for interval in bed_dict[vcf_chrom]:
                    if interval["start"] <= vcf_pos <= interval["end"]:
                        pos = vcf_pos - interval["start"] -1 # Position within the sequence
                        if pos < 0:
                            continue
                        for alt in alts.split(","):
                            if len(alt)>1:
                                continue
                            modified_sequence = (
                                interval["sequence"][:pos]
                                + alt
                                + interval["sequence"][pos + 1 :]
                            )
                            if check_gas_motif(modified_sequence):

                                print("Found GAS motif creating SNP! " + str(vcf_chrom) + ";" + str(vcf_pos) +  ";" + ref + ";" + alt + " in GAS motif modified sequence: " + modified_sequence + " at " + str(interval["start"]) + "-" + str(interval["end"]))
                                print(f"vcf_pos: {vcf_pos}, interval_start: {interval['start']}, pos: {pos}")
                                print(f"Original sequence: {interval['sequence']}")

                                current_vcfs = bed_df.at[interval["idx"], "vcf_files"]
                                if pd.isna(current_vcfs) or current_vcfs == "":
                                    bed_df.at[
                                        interval["idx"], "vcf_files"
                                    ] = os.path.basename(vcf_file)
                                else:
                                    bed_df.at[
                                        interval["idx"], "vcf_files"
                                    ] += f";{os.path.basename(vcf_file)}"

def check_gas_motif(sequence):
    """
    Checks if the sequence matches GAS motif patterns TTCxxxGAA or TTCxxxxGAA.
    """
    return sequence[:3] == "TTC" and sequence[-3:] == "GAA"

def process_bed_and_vcf(bed_file, vcf_folder, output_bed):
    """
    Reads the BED file and processes VCF files to determine GAS motif creation.
    Writes a new BED file with an extra column containing VCF file names.
    """
    # Build intervals from BED file
    bed_df, bed_dict = build_bed_intervals(bed_file)

    # Initialize an empty column for VCF file names
    bed_df["vcf_files"] = ""

    # Get all VCF files from the folder
    vcf_files = glob.glob(os.path.join(vcf_folder, "*.vcf"))
    # Process each VCF file once
    for vcf_file in vcf_files:
        parse_vcf_once(vcf_file, bed_dict, bed_df)

    # Clean up VCF file lists and remove trailing semicolons
    bed_df["vcf_files"] = bed_df["vcf_files"].str.rstrip(";")

    # Write the updated BED file
    bed_df.to_csv(output_bed, sep="\t", index=False, header=False)

    # Print directories for logging
    print(f"Input BED file: {os.path.abspath(bed_file)}")
    print(f"VCF folder: {os.path.abspath(vcf_folder)}")
    print(f"Output BED file: {os.path.abspath(output_bed)}")

if __name__ == "__main__":
    bed_file = (
        "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\AlmostGASmotifs\\all_motifs.sorted.bed"
    )
    vcf_folder = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\Sample2\\"
    output = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\test_output.bed"

    # Process BED and VCF files
    process_bed_and_vcf(bed_file, vcf_folder, output)
