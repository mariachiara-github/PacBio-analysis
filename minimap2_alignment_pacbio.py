#!/usr/bin/env python

################################################################################
# Script: minimap2_alignment_pacbio.py
# Description: This script processes a list of PacBio HiFi reads, aligns the reads 
# to the indexed transcripts using Minimap2, converts SAM to BAM, sorts and indexes 
# the BAM files, and extracts reads for each transcript. The results are saved 
# in separate directories for SAM/BAM files and FASTA files for each fusion 
# transcript and sample combination. This version uses parallel processing for
# each sample.
# Install pysam, biopython and minimap2 in a new environment and activate before use.
# Usage: sbatch minimap2_alignment_pacbio.py
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100GB
#SBATCH --time=7-12
#SBATCH --export=ALL

import os
import csv
import subprocess
import pysam
import gzip
from Bio import SeqIO
from multiprocessing import Pool
import math
import pandas as pd 


# Paths to files and tools
transcripts_fasta = "FT_sequences_86bp.fasta"
fastq_list_file = "fastq_files.txt"  # a list of files containing RNAseq reads
output_dir = "minimap2_aligned_reads" 
minimap2_index = "FT_sequences_86bp_index.mmi"  
samples_dir = os.path.join(output_dir, "samples")
sam_bam_dir = os.path.join(output_dir, "sam_bam")
csv_dir = os.path.join(output_dir, "csv")


# Create output directories if they don't exist
os.makedirs(samples_dir, exist_ok=True)
os.makedirs(sam_bam_dir, exist_ok=True)
os.makedirs(csv_dir, exist_ok=True)

print("Starting the script...")

# Step 1: Index the transcripts using Bowtie2
print("Indexing transcripts with minimap2...")
subprocess.run(["minimap2","-xmap-hifi", "-d", minimap2_index, transcripts_fasta], check=True)
print("Indexing completed.")

# Read the list of fastq.gz files
print("Reading the list of fastq.gz files...")
with open(fastq_list_file, "r") as f:
    fastq_files = [line.strip() for line in f]
    print(fastq_files)
print(f"Found {len(fastq_files)} fastq.gz files to process.")

# Load the transcripts
print("Loading transcripts from the FASTA file...")
transcripts = {record.id: str(record.seq) for record in SeqIO.parse(transcripts_fasta, "fasta")}
print(f"Loaded {len(transcripts)} transcripts.")

# Create subdirectories for each sample in the samples directory
print("Creating subdirectories for each sample in the samples directory...")
for i in fastq_files:
    names = i.replace('.fastq.gz', '')
    os.makedirs(os.path.join(samples_dir, names), exist_ok=True)
print("Subdirectories created.")

def process_sample(rna_seq_fastq):
    # Extract the sample name from the fastq names 
    sample_name = os.path.basename(rna_seq_fastq).replace('.fastq.gz', '')
    print(f"Processing sample: {sample_name}...")

    # Generate paths for SAM, BAM, and sorted BAM files
    aligned_sam = os.path.join(sam_bam_dir, f"{sample_name}_aligned_reads.sam")
    aligned_bam = os.path.join(sam_bam_dir, f"{sample_name}_aligned_reads.bam")
    sorted_bam = os.path.join(sam_bam_dir, f"{sample_name}_aligned_reads.sorted.bam")

    # Step 2: Align the RNA-seq reads to the indexed transcripts
    print(f"Aligning reads for {rna_seq_fastq} with Minimap2...")
    subprocess.run(["minimap2","-a","-xmap-hifi","--eqx","-uf","--secondary=no","--sam-hit-only","-s 70", "-O 20,26", transcripts_fasta,rna_seq_fastq,"-o", aligned_sam], check=True)
    print(f"Alignment completed for {rna_seq_fastq}.")

    # Step 3: Convert SAM to BAM, sort and index BAM
    print(f"Converting SAM to BAM for {sample_name}...")
    subprocess.run(["samtools", "view", "-bS", aligned_sam, "-o", aligned_bam], check=True)
    print(f"SAM to BAM conversion completed for {sample_name}.")

    print(f"Sorting BAM file for {sample_name}...")
    subprocess.run(["samtools", "sort", aligned_bam, "-o", sorted_bam], check=True)
    print(f"BAM sorting completed for {sample_name}.")

    print(f"Indexing sorted BAM file for {sample_name}...")
    subprocess.run(["samtools", "index", sorted_bam], check=True)
    print(f"BAM indexing completed for {sample_name}.")

    # Open the sorted BAM file
    print(f"Opening sorted BAM file for {sample_name}...")
    bamfile = pysam.AlignmentFile(sorted_bam, "rb")

    # Dictionary to store reads for each transcript
    transcript_reads = {transcript: [] for transcript in transcripts}

    # Step 4: Extract reads for each transcript
    print(f"Extracting reads for each transcript for {sample_name}...")
    for read in bamfile.fetch():
        if not read.is_unmapped:
            transcript = bamfile.get_reference_name(read.reference_id)
            transcript_reads[transcript].append(read)
    print(f"Read extraction completed for {sample_name}.")

    # Step 5: Write each transcript and its reads to a separate FASTA file
    print(f"Writing FASTA files for {sample_name}...")
    count_transcripts = 0 
    csv_lines = []
    for transcript, reads in transcript_reads.items():
        if reads: 
            transcript_fasta_path = os.path.join(samples_dir, sample_name, f"{transcript}_{sample_name}.fasta")
            count_transcripts += 1
            transcript_read_ids = []
            count_reads = 0 
            with open(transcript_fasta_path, "w") as f:
                for read in reads:
                    count_reads += 1 
                    read_seq = read.query_sequence
                    read_id = read.query_name
                    f.write(f">{transcript}_{sample_name}_{read_id}\n{read_seq}\n")
                    transcript_read_ids.append(read_id)
                read_ids_str = ', '.join(transcript_read_ids)
                csv_lines.append([transcript,read_ids_str,count_reads])         
    print(f"FASTA files written for {sample_name}. The number of transcripts matched is: {count_transcripts}.") #count number of tot reads
    print("Call the CSV_file function to create a CSV file with the Transcipts, Reads and Read Count for each sample")
    CSV_file(sample_name,csv_lines)
    bamfile.close()
    print(f"Processing completed successfully for {rna_seq_fastq}.")

def CSV_file(s,l):
    # Step 6: Write CSV file for each sample
    transcripts_csv_path = os.path.join(csv_dir, f"{s}_transcripts.csv")
    print(f"Writing CSV file to {transcripts_csv_path}")
    with open(transcripts_csv_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Transcript', 'Reads', 'Read Count'])  # Write header
        csv_writer.writerows(l)
    print(f"CSV file saved successfully to {transcripts_csv_path}")
    print("Save the CSV file names in a list...")
    
    # Step 7: Write a txt file with the CSV files name 
    print("Writing the CSV file name in the csv_aligned.txt file")
    txt_file_path = os.path.join(csv_dir, "csv_aligned.txt")
    with open(txt_file_path, 'a') as txt_file:
        txt_file.write(f"{s}_transcripts.csv\n")
    print(f"txt file saved successfully to {txt_file_path}")
    

def sort_file_contents(csv_dir):
    txt_file_path = os.path.join(csv_dir, "csv_aligned.txt")
    
    # Read the contents of the file
    with open(txt_file_path, 'r') as txt_file:
        lines = txt_file.readlines()
        print(lines)
    
    # Sort the lines alphabetically
    sorted_lines = sorted(set(lines))
    
    # Write the sorted lines back to the file
    with open(txt_file_path, 'w') as txt_file:
        txt_file.writelines(sorted_lines)
    
    print(f"Sorted the contents of {txt_file_path}")



# Parallel processing of samples
print("Starting parallel processing of samples...")
with Pool(32) as pool:
    pool.map(process_sample, fastq_files)

#Sort the CSV file 
sort_file_contents(csv_dir)

print("All samples processed successfully.")


