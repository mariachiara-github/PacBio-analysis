#!/usr/bin/env python

################################################################################
# Script: fusions_per_individual_pbfusion.py
# Description: For each of the 717 fusion transcripts, find how many of the 6
# individuals that fusion was found. This analysis is performed on the fusions 
# the minimap2 fusions matched to pbfusion.
# Usage: sbatch fusions_per_individual_pbfusion.py
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

print("For each fusion transcript, determine the number of individuals where the fusion was detected in the pbfsuion file")

# List of files containing pbfusion fusions per each of the 12 individuals
csv_individuals_pbfusion = 'csv_individuals_pbfusion.txt'
# Define the path to your FASTA file
fasta = "/path/to/Colette/Maria/fastq_flnc/FT_sequences_86bp.fasta"

# Read the pbfusion individual files 

print("Reading the list of csv_individuals_pbfusion files...")
with open(csv_individuals_pbfusion, "r") as f:
    dataframe_files_pbfusion = [line.strip() for line in f]
    print(dataframe_files_pbfusion)
print(f"Found {len(dataframe_files_pbfusion)} csv_individuals_pbfusion files to process.")


# Step 1: Filter files by prefix
def individuals(files, prefix):
    f = [file for file in files if file.startswith(prefix)]
    print(f"The files to process for individual {prefix} are: {f}")
    return f
 

def fusions_per_individual_pbfusion(files, individual):
    transcript_list = []
    
    for file in files:
        df = pd.read_csv(file)
        transcript_list.append(df['Transcript'])
    
    # Combine all transcript columns into a single series, then remove duplicates
    combined_transcripts = pd.concat(transcript_list).drop_duplicates().reset_index(drop=True)
    # Convert to DataFrame
    combined_df = pd.DataFrame(combined_transcripts, columns=['Transcript'])
    # Save to CSV
    combined_df.to_csv(f"sample{individual}_pbfusion.csv", index=False)
    print("CSV file saved successfully")
    return combined_df

def fusion_frequency_after_minimap(fasta_file,dataframes):

    # Initialize the dictionary with fusion transcript names as keys and zeros as values
    fusion_transcripts = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fusion_name = record.id
        fusion_transcripts[fusion_name] = 0
    
    # Check across all dataframes and update the dictionary
    for df in dataframes:
        print(df)
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
    fusion_freq.to_csv("fusion_frequency_pbfusion.csv",index=False)

        

N1_individual = fusions_per_individual_pbfusion(individuals(dataframe_files_pbfusion, "N1"), "N1")
N28_individual = fusions_per_individual_pbfusion(individuals(dataframe_files_pbfusion, "N28"), "N28")
R2_individual = fusions_per_individual_pbfusion(individuals(dataframe_files_pbfusion, "R2"), "R2")
R8_individual = fusions_per_individual_pbfusion(individuals(dataframe_files_pbfusion, "R8"), "R8")
print(f"The files to process for individual C5_6 are: sample5_transcripts.csv,sample6_transcripts.csv")
fetal_cerebellum_C5_6 = fusions_per_individual_pbfusion(["sample5_pbfusion.csv","sample6_pbfusion.csv"], "_C5_6")

sample8_individual = "sample8_pbfusion.csv"
sample8_individual_dataframe = pd.read_csv(sample8_individual)
sample8_individual_finale = sample8_individual_dataframe.drop_duplicates(subset='Transcript')
dataframe_files = [N1_individual,N28_individual,R2_individual,R8_individual,fetal_cerebellum_C5_6,sample8_individual_finale]

fusion_frequency_after_minimap(fasta,dataframe_files)
