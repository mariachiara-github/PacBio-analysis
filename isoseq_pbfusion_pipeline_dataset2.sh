#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --time=15-12
#SBATCH --export=ALL


# SET ALL THE VARIABLES NECESSARY FOR THIS SCRIPT 
# Set variables for the pbfusions workflow
h38mmi="/path/to/GRCh38.p13/GRCh38.p13.genome.mmi"
gtfv38="/pathto/GRCh38.p13/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.bin.xz"
# Set the variable for the primer FASTA file
FASTA="/path/to/SRA_PRJDB15555/isoseq_primers.fasta"


# CREATE ALL THE NECESSARY DIRECTORIES
# Create a directory to store the flnc.fastq.gz files that will be used in further analyses
mkdir -p "/path/to/fastq_flnc/"
fastq_flnc="/path/to/Maria/fastq_flnc/"


# LISTING THE VERSION OF EACH TOOL USED IN THE PIPELINE
echo "extracthifi version"
extracthifi --version
echo "lima version"
lima --version
echo "isoseq version (used for refine)"
isoseq refine --version
echo "bam2fastq version"
bam2fastq --version
echo "pbmm2 version"
pbmm2 --version 
echo "pbfusion discover version"
pbfusion discover --version

# ARRAY OF SAMPLES NAMES
samples=("N1_C" "N1_H" "N1_T" "N28_C" "N28_H" "N28_T" "R2_C" "R2_H" "R2_T" "R8_C" "R8_H" "R8_T") 

# RUN THE ISO-SEQ AND PBFUSION WORKFLOW FOR EACH SAMPLE
for sample in "${samples[@]}"; do
    echo "Creating directory for sample: $sample"
    mkdir -p "./${sample}_isoseq"
    sample_dir="./${sample}_isoseq"

    echo "Processing sample: $sample"

    # RUNNING THE ISO-SEQ WORKFLOW
    echo "1. Extracting HiFi reads from sample_ccs.bam file"
    extracthifi Analysis.${sample}.reads.bam "${sample_dir}/${sample}.hifi.bam"
    echo "extracthifi finished"
    
    #Change directory to the sample directory
    cd "$sample_dir"
    
    echo "2. Count the number of HiFi reads"
    samtools view ${sample}.hifi.bam | wc -l > hifi_reads_${sample}.txt
    echo "count finished"

    echo "3. Running lima for primer removal and demultiplexing"
    lima ${sample}.hifi.bam ${FASTA} ${sample}.fl.bam --isoseq --dump-clips --peek-guess
    echo "lima finished"

    echo "4. Running isoseq refine for trimming PolyA tails and concatemer removal"
    isoseq refine ${sample}.fl.5p--3p.bam ${FASTA} ${sample}.flnc.bam --require-polya
    echo "isoseq refine finished"
    
    echo "5. Count the number of flnc reads"
    samtools view ${sample}.flnc.bam | wc -l > flnc_${sample}.txt
    echo "count finished"
    
    # RUNNING THE pbfusion WORKFLOW
    #Set the variable for the flnc file
    export flnc="${sample}.flnc.bam"
    
    # Create the pbfusion directory 
    mkdir -p "./pbfusion_${sample}"
    pbfusion_dir="./pbfusion_${sample}"
    
    echo "6. Running pbmm2 for flnc.bam reads"
    echo "Use index to align reads"
    pbmm2 align ${h38mmi} ${flnc} "${pbfusion_dir}/${sample}.aligned.bam" --preset ISOSEQ --sort --log-level INFO --log-file "${pbfusion_dir}/align_${sample}.txt"
    echo "pbmm2 flnc.bam finished"
    
    cd ${pbfusion_dir}
    
    echo "7. Running pbfusion discover"
    pbfusion discover \
        --gtf ${gtfv38} \
        --output-prefix isoseq_${sample} \
        --threads 8 \
        --min-fusion-quality LOW \
        --emit-readthrough \
        --emit-novel-exons \
        -v \
        ${sample}.aligned.bam
    echo "pbfusion discover finished"

    # Change directory to the previous one (sample_dir)
    cd ..
    
    # We want to use the flnc.bam file in further analyses, but in the flnc.fastq.gz format 
    echo "8. Converting flnc.bam file into a fastq.gz file for further analyses"
    bam2fastq -o "${fastq_flnc}/${sample}" ${flnc}
    echo "fastq.gz finished"
    
    # Change the directory to the working directory to process the other samples
    cd ..

done

echo "All samples processed successfully."


