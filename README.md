#Overview of the scripts

01_AllGASmotifs.py and 01_AllGASmotifs.py take VCF files and the files with all GAS motifs or all almost GAS motifs on the genome (created by FIMO like in Hoffmann et al. (2024) "Data-driven projections of candidate enhancer-activating SNPs in immune regulation") and find individuals that either create or destroy a GAS motif


02_GetCloseByGenes.py takes the output of 01 Python scripts and annotates gene names (either inside, up to 2.5 kb or up to 10 kb)
