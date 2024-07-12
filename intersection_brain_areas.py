#!/usr/bin/env python

################################################################################
# Script: intersection_brain_areas.py
# Description: This script takes as input CSV files of the 12 samples of dataset 
# 2, generated from the minimap2_alignment_pacio.py, and returns an output file 
# containing the number of fusions shared among the 3 different brain regions 
# (cerebellum,hypotalamus,temporal cortex) for each individual.
# Usage: sbatch intersection_brain_areas.py
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100GB
#SBATCH --time=7-12
#SBATCH --export=ALL

import pandas as pd 

csv_files = "csv_aligned.txt" 

# Read the individual files
print("Reading the list of csv_aligned files...")
with open(csv_files, "r") as f:
    dataframe_files = [line.strip() for line in f]
    print(dataframe_files)
print(f"Found {len(dataframe_files)} csv_aligned files to process.")

# Step 1: Filter files by prefix
def individuals(files, prefix):
    f = [file for file in files if file.startswith(prefix)]
    print(f"The files to process for individual {prefix} are: {f}")
    return f
    
# Step 2: Open and process the filtered files
def intersection_brain_areas(files,individual):
    dataframe_list = []
    for file in files:
        df = pd.read_csv(file)
        dataframe_list.append(df)
    df_C = dataframe_list[0]
    df_H = dataframe_list[1]
    df_T = dataframe_list[2]
    
    # Extract the first columns as sets
    fusions_C = set(df_C.iloc[:, 0])
    fusions_H = set(df_H.iloc[:, 0])
    fusions_T = set(df_T.iloc[:, 0])

    # Find intersections
    common_C_H = fusions_C.intersection(fusions_H)
    common_C_T = fusions_C.intersection(fusions_T)
    common_H_T = fusions_H.intersection(fusions_T)
    
    common_C_H_sorted = sorted(list(common_C_H))
    common_C_T_sorted = sorted(list(common_C_T))
    common_H_T_sorted = sorted(list(common_H_T))
    
    # Print results in an output file 
    message1 = (f"The common fusions between the cerebellum and hypotalamus of individual {individual} are {len(common_C_H_sorted)}: {common_C_H_sorted}")
    message2 = (f"The common fusions between the cerebellum and temporal cortex of individual {individual} are {len(common_C_T_sorted)}: {common_C_T_sorted}")
    message3 = (f"The common fusions between the hypotalamus and temporal cortex of individual {individual} are {len(common_H_T_sorted)}: {common_H_T_sorted}")

  
    with open('output_intersection_brain_areas.txt', "a") as out_f:
        print(message1,file=out_f)
        print(message2,file=out_f)
        print(message3,file=out_f)
  

N1_individual = intersection_brain_areas(individuals(dataframe_files, "N1"), "N1")
N28_individual = intersection_brain_areas(individuals(dataframe_files, "N28"), "N28")
R2_individual = intersection_brain_areas(individuals(dataframe_files, "R2"), "R2")
R8_individual = intersection_brain_areas(individuals(dataframe_files, "R8"), "R8")
