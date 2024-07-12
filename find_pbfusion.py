#!/usr/bin/env python

###################################################################################
# Script: find_pbfusions.py
# Description: This script processes a list of CSV files (output of 
# minimap2_alignment_pacbio.py) and a list of pbfusion BED files. For each CSV and 
# pbfusion file (of the same sample) this script finds matches between the 2 files 
# based on the PacBio read name. The matches are stored in a pandas dataframe and
# are then saved in a separate directory as CSV files for each sample. 
# Moreover, this script checks how many interchromosomal fusions are found in each 
# final dataframe and if the gene name of the short-read fusion transcript 
# matches the gene name of the long-read fusion (pbfusion name). 
# The number of samples processed by this script is 15 (First database: 3 samples 
# = 2 individuals, Second Database: 12 samples = 4 individuals).
# This version uses parallel processing for each sample.
# Install pysam and pandas in a new environment and activate before use.
# Usage: sbatch minimap2_alignment_pacbio.py
###################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100GB
#SBATCH --time=7-12
#SBATCH --export=ALL

import os
import re
import pysam
from multiprocessing import Pool
import math
import pandas as pd 

# Paths to files and tools
csv_files = "csv_aligned.txt"                                # Text file with list of the matched_minimap2.csv files
pbfusion_list_file = "dataframes_pbfusion_all_samples.txt"   # Text file with list of the pbfusion_output.bed files directories
interchromosomal_fusions_total = "/path/to/FusionCatcher_717_fusions/interchromosomal_fusions_total.csv" # CSV file of short reads inter-chromosomal fusions
dataframes_pbfusion = os.path.join("dataframes_pbfusion")    

# Create output directories if they don't exist
os.makedirs(dataframes_pbfusion, exist_ok=True)

print("Starting the script...")

# Read the list of matched_minimap2.csv files
print("Reading the list of csv_aligned files...")
with open(csv_files, "r") as f:
    dataframe_files = [line.strip() for line in f]
    print(dataframe_files)
print(f"Found {len(dataframe_files)} csv_aligned files to process.")

# Read the list of pbfusion_output.bed files
print("Reading the list of pbfusion_output.bed files...")
with open(pbfusion_list_file, "r") as f:
    pbfusion_files = [line.strip() for line in f]
print(f"Found {len(pbfusion_files)} pbfusion_output.bed files to process.")


def read_dataframes(sample,pbfusion,interchromosomal_fusions_total):
     
    print(f"Processing {sample} with {pbfusion}...")
    # Step 1: Create a sorted dataframe from the CSV file of read names
    print("Writing a dataframe of the CSV file of matched fusion transcripts...")
    print(f"Opening sorted BAM file for {sample}...")
    df_transcripts = pd.read_csv(sample)
    df_transcripts_sorted = df_transcripts.sort_values(by='Transcript').reset_index(drop=True)
    read_count_sample = df_transcripts_sorted['Read Count']
    all_reads= 0
    for count in read_count_sample:
        all_reads += count 
    print(f"Fusion Transcripts DataFrame created and sorted. The total number of reads that matched the 717 fusion in {sample} is: {all_reads}")
     
    # Step 2: Create a dataframe form the pbfusion file 
    print("Writing a dataframe of the pbfusion file...")
    column_names = ['chr1',	'start1', 'end1', 'chr2', 'start2', 'end2','id', 'score', 'strand1', 'strand2', 'info',	'extra']     
    df_pbfusion = pd.read_csv(pbfusion, sep='\t', comment='#', names=column_names)
    
    #Extract the arguments of info and extra in different columns

    df_pbfusion[['RC','MD', 'GN', 'GI','GC','CL','MQ','MI','SA','BP','EX','IN','AC']] = df_pbfusion['info'].str.split(';', expand=True)
    df_pbfusion[['RN','ON','CB']] = df_pbfusion['extra'].str.split(';', expand=True)

    #Delete the extra and info columns (not needed anymore)

    df_pbfusion = df_pbfusion.drop('extra', axis=1)
    df_pbfusion = df_pbfusion.drop('info', axis=1)

    #Clean the columns, ex: RC=1 --> 1
    columns_to_process = ['RC','MD', 'GN', 'GI','GC','CL','MQ','MI','SA','BP','EX','IN','AC','RN','ON','CB']
    for col in columns_to_process:
        df_pbfusion[col] = df_pbfusion[col].str.extract(r'=(.*)')

    print(df_pbfusion)
    print("pbfusion DataFrame created")
    

    # Step 3: Create a list with all the read names of that sample
    print("Create a list with all the reads to match...")
    reads = df_transcripts_sorted['Reads']
    reads_list = []
    for read in reads:
        if ',' in read:
            split_read = read.split(',')
            for r in split_read:
                reads_list.append(str(r).lstrip())
        else:
            reads_list.append(str(read).lstrip())

    print("List created")

    
    # Step 4: Find reads in read_name_pbfusion and create a new dataframe
    print("Search the reads names of the list in the pbfusion dataframe...")

    df_pbfusion_find = pd.DataFrame()
    series_transcript = pd.Series(dtype=str)
    reads_found = 0 

    # Iterate over each unique read in the reads list
    for unique_read in reads_list:
        df_grep = df_pbfusion[df_pbfusion['RN'].str.contains(unique_read, na=False)]
        if not df_grep.empty:
            reads_found += 1
            for _, row in df_grep.iterrows():
                transcript = df_transcripts_sorted[df_transcripts_sorted['Reads'].str.contains(unique_read, na=False)]['Transcript']
                if len(transcript) > 1: 
                    for name in transcript:
                        # Create a new row for each transcript found
                        new_row = row.to_frame().T
                        new_row['Transcript'] = name
                        df_pbfusion_find = pd.concat([df_pbfusion_find, new_row], ignore_index=True)
                        
                        series_transcript = pd.concat([series_transcript, pd.Series([name])], ignore_index=True)
                else:
                    new_row = row.to_frame().T
                    new_row['Transcript'] = transcript.values[0]
                    df_pbfusion_find = pd.concat([df_pbfusion_find, new_row], ignore_index=True)
                    
                    series_transcript = pd.concat([series_transcript, transcript], ignore_index=True)

    # Ensure 'Transcript' is the first column
    df_pbfusion_find = df_pbfusion_find[['Transcript'] + [col for col in df_pbfusion_find.columns if col != 'Transcript']]
    
    df_pbfusion_find = df_pbfusion_find.drop_duplicates()
    fusion_matched = len(df_pbfusion_find['Transcript'].drop_duplicates())   
    
    # Step 5: call the function inter_fusions to find the number of possible possible interchromosomal fusion transcripts in each sample
    interch = inter_fusions(interchromosomal_fusions_total,df_pbfusion_find,sample)
    message1 = (f"The number of possible interchromosomal fusion transcripts identified in {sample} is: {interch}")
    
    # Step 6: call the function check_gene_matches to check for each short read gene name, if it correspond to the long read gene name (check if it has 0,1 or 2 names in common)
    df_pbfusion_find['GeneMatches'] = df_pbfusion_find.apply(lambda row: check_gene_matches(row['Transcript'], row['GN']), axis=1)
    # Count the occurrences of 0, 1, and 2 matches
    counts = df_pbfusion_find['GeneMatches'].value_counts().sort_index()
    for matches in range(3):
        count = counts.get(matches, 0)
        print(f"There are {count} FusionIDs with {matches} genes in common with the gene names.")

   
    # Step 7: create a csv file for the df_pbfusion_find
    dataframes_pbfusion_path = os.path.join(dataframes_pbfusion, f"{sample.split('_transcripts')[0]}_pbfusion.csv")
    print(f"Writing pbfusion file to {dataframes_pbfusion_path}")
    df_pbfusion_find.to_csv(dataframes_pbfusion_path,index=False)
    message2 = (f"The number of reads matched for {sample} is: {reads_found}/{all_reads}. The number of fusion transcripts found in pbfusion is: {fusion_matched}. A DataFrame is created.")
    with open('output.txt', "a") as out_f:
        print(message1,file=out_f)
        print(message2,file=out_f)


def inter_fusions(fusions,df_pbfusion_inter,s):
    inter = pd.DataFrame()
    l = []
    df_f = pd.read_csv(fusions)
    for f in df_f["FusionID"]:
        l.append(str(f))

    for inter_f in l: 
        inter_f = inter_f.lstrip()
        df = df_pbfusion_inter.loc[df_pbfusion_inter['Transcript'] == inter_f]
        #print(df_pbfusion_inter)
        if not df.empty:
            print("There is/are inter-chromosomal fusion/s")
            different_chr = df[df['chr1'] != df['chr2']]
            #print(different_chr)
            inter = pd.concat([inter,different_chr], ignore_index=True)
            #print(inter)
        else: 
            print("There are no inter-chromosomal fusions")
    new_name = f"Transcript_{s.split('_transcripts')[0]}"
    unique_inter_transcript = inter['Transcript'].drop_duplicates()

    # Add the string to the end of each entry in the 'Transcript' column
    inter.rename(columns = {'Transcript': new_name}, inplace = True)
    inter.to_csv('interchromosomal_fusions.csv', mode='a', index=False)  
    return len(unique_inter_transcript)


def check_gene_matches(fusion_id, gene_names):
    # Extract gene names from the FusionID
    fusion_parts = fusion_id.split('_')
    genes_in_fusion = [fusion_parts[0].upper()] if len(fusion_parts) < 3 else [gene.upper() for gene in fusion_parts[:2]]
    # Handle NaN values in gene_names
    if pd.isna(gene_names):
        gene_list = []
    else:
        gene_list = [gene.upper() for gene in gene_names.split(',')]
    
    # Count matches
    matches = len(set(genes_in_fusion) & set(gene_list))
 
    return matches


# Prepare the arguments for starmap
args = [(df_file, pb_file, interchromosomal_fusions_total) for df_file, pb_file in zip(dataframe_files, pbfusion_files)]
# Parallel processing of samples
print("Starting parallel processing of samples...")
with Pool(32) as pool:
    print(args)
    pool.starmap(read_dataframes,args)

print("All samples processed successfully.")


print("Identification of a unique list of possible interchromosomal fusion transcripts...")


def unique_list_interchromosomal_f(file_inter,inter_total):
    
    # Step 1: read the 'interchromosomal_fusions.csv'file and 
    df_inter = pd.read_csv(file_inter)
    col_name = (df_inter.columns.tolist())
    
    inter_total_df = pd.read_csv(inter_total)
    # Step 2 : determine the unique list of interchromosomal fusions and determine which ones correspond to the short read version
    unique_transcripts = df_inter[col_name[0]].drop_duplicates()
    unique_transcripts_list = unique_transcripts.tolist()
    filtered_transcripts_list = sorted([name for name in unique_transcripts_list if not name.startswith('Transcript')])
    message = (f"The possibile interchromosomal fusion transcripts identified among all samples are {len(filtered_transcripts_list)} and are the following: {filtered_transcripts_list}")
 
    # Step 3: Determine which of the short inter-chromosomal fusions are found to be between the same chromosomes in long reads
    matching_rows = []
    df2_grouped = df_inter.groupby(col_name[0])

    # Iterate through each row in the first DataFrame
    for index, row in inter_total_df.iterrows():
        fusion_id = row['FusionID']
        chr1 = row['Chromosome_1.5end_partner']
        chr2 = row['Chromosome_2.3end_partner']
        
        if fusion_id in df2_grouped.groups:
            for _, matched_row in df2_grouped.get_group(fusion_id).iterrows():
                if matched_row['chr1'] == chr1 and matched_row['chr2'] == chr2:
                    matching_rows.append(fusion_id)
                    break
    
    message2 = (f"The possibile interchromosomal fusion transcripts found to be between the same chromosomes in long reads as short reads are {len(matching_rows)} and are the following: {matching_rows}")
 
    # Step 4: write the output to the output file.
    with open('output.txt', "a") as out_f:
        print(message,file=out_f)
        print(message2,file=out_f)
    print("Identification finished")   
    return 


unique_list_interchromosomal_f('interchromosomal_fusions.csv',interchromosomal_fusions_total)


################################################################################
# The number of total samples is 15, but the number of individuals is 6.
# Before the frequency analysis, run the following functions 2 
# functions for the 3 cerebellum samples (2 individuals) and the 12 samples 
# (4 individuals).
# Description: The functions are used to identify a unique set of fusions per
# individual. The results are saved in a CSV file for each of the 6 individuals
################################################################################


# Step 1: Filter files by prefix
def individuals(files, prefix):
    f = [file for file in files if file.startswith(prefix)]
    print(f"The files to process for individual {prefix} are: {f}")
    return f
    
# Step 2: Open and process the filtered files
def fusions_per_individual(files,individual):
    dataframe_list = []
    for file in files:
        df = pd.read_csv(file)
        dataframe_list.append(df)
    merged_df = pd.concat(dataframe_list, ignore_index=True) 
    merged_df = merged_df.drop(['Reads'], axis=1)
    merged_df = merged_df.groupby('Transcript', as_index=False).sum()
    print(merged_df)
    merged_df.to_csv(f"sample{individual}_transcripts.csv",index=False)
    print("CSV files saved successfully")
    
N1_individual = fusions_per_individual(individuals(dataframe_files, "N1"), "N1")
N28_individual = fusions_per_individual(individuals(dataframe_files, "N28"), "N28")
R2_individual = fusions_per_individual(individuals(dataframe_files, "R2"), "R2")
R8_individual = fusions_per_individual(individuals(dataframe_files, "R8"), "R8")
print(f"The files to process for individual C5_6 are: sample5_transcripts.csv,sample6_transcripts.csv")
fetal_cerebellum = fusions_per_individual(["sample5_transcripts.csv","sample6_transcripts.csv"], "_C5_6")


    
