#!/usr/bin/env python

################################################################################
# Script: get_long_reads.py
# Description: This script was used to retrieve the long reads sequences among 
# the 29 samples, of KANSL1_ARL17A_17, KANSL1_ARL17B_17, KANSL1_LRRC37A3,
# NAIP_OCLN and NSF_LRRC37A3 fusion transcripts
# Usage: sbatch get_long_reads.py
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100GB
#SBATCH --time=7-12
#SBATCH --export=ALL

import os
import glob

dir_fusions = 'directories.txt'      # txt file with sample directories, containing the FASTA files for each FTs (each sample directory is found in the SAMPLE directory created with the minimap2_alignment_pacbio.py script)

# Read the individual files
print("Reading the list of csv_aligned files...")
with open(dir_fusions, "r") as f:
    directories = [line.strip() for line in f]
    print(directories)
print(f"Found {len(directories)} directories to process.")


# List of base FASTA names to search for (without the _N1, _N2, etc.)
fasta_base_names = ['KANSL1_ARL17A_17:46094560:-_17:46517233:-', 
                    'KANSL1_ARL17A_17:46094560:-_17:46570869:-',
                    'KANSL1_ARL17A_17:46170855:-_17:46517233:-',
                    'KANSL1_ARL17A_17:46170855:-_17:46570869:-',
                    'KANSL1_ARL17B_17:46094560:-_17:46299665:-',
                    'KANSL1_ARL17B_17:46094560:-_17:46352930:-',
                    'KANSL1_ARL17B_17:46170855:-_17:46299665:-',
                    'KANSL1_ARL17B_17:46170855:-_17:46352930:-',
                    'KANSL1_LRRC37A3_17:46152904:-_17:64868536:-',
                    'KANSL1_LRRC37A3_17:46152904:-_17:64869166:-',
                    'KANSL1_LRRC37A3_17:46152904:-_17:64892600:-',
                    'KANSL1_LRRC37A3_17:46170855:-_17:64868536:-',
                    'KANSL1_LRRC37A3_17:46170855:-_17:64869166:-',
                    'KANSL1_LRRC37A3_17:46170855:-_17:64892600:-',
                    'NAIP_OCLN_5:70974129:-_5:69534694:+',
                    'NAIP_OCLN_5:70979869:-_5:69534694:+',
                    'NAIP_OCLN_5:70983775:-_5:69534694:+',
                    'NAIP_OCLN_5:70998724:-_5:69534694:+',
                    'NSF_LRRC37A3_17:46704854:+_17:64860973:-',
                    'NSF_LRRC37A3_17:46704854:+_17:64869166:-',
                    'NSF_LRRC37A3_17:46704854:+_17:64892600:-',
                    'NSF_LRRC37A3_17:46704854:+_17:64897654:-'
                    ]

# Function to concatenate FASTA files
def concatenate_fasta_files(fasta_base_name, directories):
    output_file_path = f'/zfs/jacobs/Colette/Maria/fastq_flnc/long_reads_fusions/{fasta_base_name}_combined.fasta'
    content_written = False
    with open(output_file_path, 'w') as outfile:
        for directory in directories:
            # Search for files matching the pattern in the current directory
            pattern = os.path.join(directory, f'{fasta_base_name}_*.fasta')
            matching_files = glob.glob(pattern)
            for file_path in matching_files:
                try:
                    with open(file_path, 'r') as infile:
                        content = infile.read()
                        outfile.write(content)
                        content_written = True
                except Exception as e:
                    print(f"An error occurred while processing {file_path}: {e}")
    
    # Remove the file if no content was written
    if not content_written:
        os.remove(output_file_path)
        print(f"No reads found for {fasta_base_name}, no file created.")
    else:
        print(f"Concatenated content written to {output_file_path}")

# Process each FASTA base name
for fasta_base_name in fasta_base_names:
    concatenate_fasta_files(fasta_base_name, directories)
