import pandas as pd
import os
import glob
import numpy as np
from intervaltree import IntervalTree
from collections import defaultdict
from datetime import datetime

def print_with_time(message):
    """
    Prints a message with the current timestamp.
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")

def build_bed_intervals(bed_file):
    """
    Reads the BED file and creates interval trees for fast lookup.
    """
    # Load the BED file
    print_with_time("Load the bedfile")
    bed_df = pd.read_csv(
        bed_file, sep="\t", header=None, names=["chrom", "start", "end", "sequence"]
    )
    # Process the 'sequence' column
    bed_df["sequence_type"] = bed_df["sequence"]
    bed_df["sequence"] = bed_df["sequence"].str.split('-').str[1]

    print_with_time("Finished loading the bedfile")
    print_with_time("Populate the IntervalTree")

    # Initialize a dictionary of interval trees
    bed_trees = defaultdict(IntervalTree)

    # Populate the interval trees
    for idx, row in bed_df.iterrows():
        chrom = str(row["chrom"])
        start = row["start"]
        end = row["end"]
        bed_trees[chrom][start:end] = {
            "sequence": row["sequence"],
            "sequence_type": row["sequence_type"],
            "idx": idx,
        }

    print_with_time("Finished populating the IntervalTree")


    return bed_df, bed_trees

def parse_vcf_once(vcf_file, bed_trees, bed_df):
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

            counter += 1
            if counter % 10000 == 0:
                print_with_time("Checked " + str(counter) + " SNPs in " + vcf_file)

            if vcf_chrom in bed_trees:
                # Query the interval tree for overlapping intervals
                overlapping_intervals = bed_trees[vcf_chrom][vcf_pos]
                for interval in overlapping_intervals:
                    data = interval.data
                    pos = vcf_pos - interval.begin - 1  # Position within the sequence
                    if pos < 0:
                        continue
                    for alt in alts.split(","):
                        if len(alt) > 1:
                            continue
                        modified_sequence = (
                            data["sequence"][:pos]
                            + alt
                            + data["sequence"][pos + 1:]
                        )
                        if check_gas_motif(modified_sequence):
                            print_with_time("Found GAS motif creating SNP! " + str(vcf_chrom) + ";"+ str(vcf_pos) +";" + ref +";" + alt)
                            print_with_time("in GAS motif modified sequence: " + modified_sequence + " from original sequence:" + data["sequence"] + " at position " + str(pos))
                            print_with_time("at " + str(interval.begin)+"-"+str(interval.end))

                            current_vcfs = bed_df.at[data["idx"], "vcf_files"]
                            if pd.isna(current_vcfs) or current_vcfs == "":
                                bed_df.at[data["idx"], "vcf_files"] = os.path.basename(vcf_file)
                                bed_df.at[data["idx"], "modified_sequence"] = modified_sequence
                                bed_df.at[data["idx"], "position_chromosome"] = str(vcf_chrom)
                                bed_df.at[data["idx"], "position_nucleotide"] = str(vcf_pos)
                                bed_df.at[data["idx"], "ref_allele"] = ref
                                bed_df.at[data["idx"], "alt_allele"] = alt
                                bed_df.at[data["idx"], "position_in_GAS"] = str(pos)
                            else:
                                bed_df.at[data["idx"], "vcf_files"] += f";{os.path.basename(vcf_file)}"
                                bed_df.at[data["idx"], "modified_sequence"] += f";{modified_sequence}"
                                bed_df.at[data["idx"], "position_chromosome"] += f";{vcf_chrom}"
                                bed_df.at[data["idx"], "position_nucleotide"] += f";{vcf_pos}"
                                bed_df.at[data["idx"], "ref_allele"] += f";{ref}"
                                bed_df.at[data["idx"], "alt_allele"] += f";{alt}"
                                bed_df.at[data["idx"], "position_in_GAS"] += f";{pos}"



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
    bed_df, bed_trees = build_bed_intervals(bed_file)

    # Initialize an empty column for VCF file names
    bed_df["vcf_files"] = ""
    bed_df["modified_sequence"] = ""
    bed_df["position_chromosome"] = ""
    bed_df["position_nucleotide"] = ""
    bed_df["ref_allele"] = ""
    bed_df["alt_allele"] = ""
    bed_df["alt_allele"] = "position_in_GAS"


    # Get all VCF files from the folder
    vcf_files = glob.glob(os.path.join(vcf_folder, "*.vcf"))
    # Process each VCF file once
    for vcf_file in vcf_files:
        parse_vcf_once(vcf_file, bed_trees, bed_df)

    # Clean up VCF file lists and remove trailing semicolons
    bed_df["vcf_files"] = bed_df["vcf_files"].str.rstrip(";")

    # Filter out rows where 'vcf_files' is empty
    bed_df = bed_df[bed_df["vcf_files"].notna() & (bed_df["vcf_files"] != "")]

    # Write the updated BED file
    bed_df.to_csv(output_bed, sep="\t", index=False, header=True)

    # Print directories for logging
    print_with_time(f"Input BED file: {os.path.abspath(bed_file)}")
    print_with_time(f"VCF folder: {os.path.abspath(vcf_folder)}")
    print_with_time(f"Output BED file: {os.path.abspath(output_bed)}")

if __name__ == "__main__":
    bed_file = (
        "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\AlmostGASmotifs\\all_motifs.sorted.bed"
    )
    vcf_folder = "\\\\shares2.dkisilon2.niddk.nih.gov\\LGPGenomics\\shared\\Austria-WGS\\Others\\"
    output = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs_WGS\\others.bed"

    if not os.path.exists(vcf_folder):
        print_with_time(f"Network folder not found: {vcf_folder}")
    else:
        print_with_time(f"Network folder found: {vcf_folder}")

    # Process BED and VCF files
    process_bed_and_vcf(bed_file, vcf_folder, output)
