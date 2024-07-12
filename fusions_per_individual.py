#!/usr/bin/env python

################################################################################
# Script: fusions_per_individual.py
# Description: For each of the 717 fusion transcripts, find in how many of the 6
# individuals that fusion was found. This analysis is performed on the fusions 
# identified in each individual by minimap2.
# Usage: sbatch fusions_per_individual.py
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100GB
#SBATCH --time=7-12
#SBATCH --export=ALL

import pandas as pd
import pandas as pd
import csv 
from Bio import SeqIO

print("For each fusion transcript, determine the number of individuals where the fusion was identified after the minimap2 alignment")

# List of files containing CSV fusions identified by minimap2 per each individual
csv_individuals = 'csv_individuals.txt'

# Define the path to your FASTA file
fasta = "/path/to/fastq_flnc/FT_sequences_86bp.fasta"


# Read the individual files
print("Reading the list of csv_individuals files...")
with open(csv_individuals, "r") as f:
    dataframe_files = [line.strip() for line in f]
    print(dataframe_files)
print(f"Found {len(dataframe_files)} individuals files to process.")


def fusion_frequency_after_minimap(fasta_file,dataframes):

    # Initialize the dictionary with fusion transcript names as keys and zeros as values
    fusion_transcripts = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fusion_name = record.id
        fusion_transcripts[fusion_name] = 0
    
    # Check across all dataframes and update the dictionary
    for df in dataframes:
        df = pd.read_csv(df)
        for fusion in df['Transcript']:
            if fusion in fusion_transcripts:
                fusion_transcripts[fusion] += 1

    # Display the resulting dictionary
    print(fusion_transcripts)
    fusion_freq = pd.DataFrame(fusion_transcripts.items(), columns=['Fusion_name', 'Number_individuals'])
    fusions_with_individuals = fusion_freq[fusion_freq['Number_individuals'] > 0]

    # Count the number of such fusions
    count_fusions_with_individuals = fusions_with_individuals.shape[0]

    print(f"Number of fusions detected in at least one individual is: {count_fusions_with_individuals}")
    fusion_freq.to_csv("fusion_frequency.csv",index=False)

        
fusion_frequency_after_minimap(fasta,dataframe_files)
