#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=300GB
#SBATCH --time=15-12
#SBATCH --export=ALL


# SET ALL THE VARIABLES NECESSARY FOR THIS SCRIPT 
# Set variables for the pbfusions workflow
h38mmi="/path/to/GRCh38.p13/GRCh38.p13.genome.mmi"
gtfv38="/path/to/GRCh38.p13/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.bin.xz"

# CREATE ALL THE NECESSARY DIRECTORIES
# Create a directory to store the flnc.fastq.gz files that will be used in further analyses
mkdir -p "/path/to/fastq_flnc/"
fastq_flnc="/path/to/fastq_flnc/"

# LISTING THE VERSION OF EACH TOOL USED IN THE PIPELINE
echo "bam2fastq version"
bam2fastq --version
echo "pbmm2 version"
pbmm2 --version
echo "pbfusion discover version"
pbfusion discover --version

#DOWNLOADING THE FLNC.BAM FILES FOR EACH OF THE 3 SAMPLES

echo "Downloading flnc_sample5.bam and flnc_sample5.bam.pbi files..."
wget -O sample5.flnc.bam https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/DATA-Revio-SCRI-Sample5-Cerebellum/2-FLNC/flnc.bam
wget -O sample5.flnc.bam.pbi https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/DATA-Revio-SCRI-Sample5-Cerebellum/2-FLNC/flnc.bam.pbi
echo "Downloading flnc_sample6.bam and flnc_sample6.bam.pbi files..."
wget -O sample6.flnc.bam https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/DATA-Revio-SCRI-Sample6-Cerebellum/2-FLNC/flnc.bam
wget -O sample6.flnc.bam.pbi https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/DATA-Revio-SCRI-Sample6-Cerebellum/2-FLNC/flnc.bam.pbi
echo "Downloading flnc_sample8.bam and flnc_sample8.bam.pbi files..."
wget -O sample8.flnc.bam https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/DATA-Revio-SCRI-Sample8-Cerebellum/2-FLNC/flnc.bam
wget -O sample8.flnc.bam.pbi https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/DATA-Revio-SCRI-Sample8-Cerebellum/2-FLNC/flnc.bam.pbi

# ARRAY OF SAMPLES NAMES
samples=("sample5" "sample6" "sample8") 

# RUN THE PBFUSION WORKFLOW FOR EACH SAMPLE
for sample in "${samples[@]}"; do

    echo "Processing sample: $sample"
    
    echo "1. Count the number of flnc reads"
    samtools view ${sample}.flnc.bam | wc -l > flnc_${sample}.txt
    echo "count finished"
   
    # RUNNING THE pbfsion WORKFLOW
    #Set the varaible for the flnc file
    export flnc="${sample}.flnc.bam"
    
    # Create the pbfusion directory 
    mkdir -p "./pbfusion_${sample}"
    pbfusion_dir="./pbfusion_${sample}"
    
    echo "2. Running pbmm2 for flnc.bam reads"
    echo "Use index to align reads"
    pbmm2 align ${h38mmi} ${flnc} "${pbfusion_dir}/${sample}.aligned.bam" --preset ISOSEQ --sort --log-level INFO --log-file "${pbfusion_dir}/align_${sample}.txt"
    echo "pbmm2 flnc.bam finished"
    
    cd ${pbfusion_dir}
    
    echo "3.. Running pbfusion discover"
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
    

    # We want to use the flnc.bam file in further anlyses, but in the flnc.fastq.gz format 
    echo "4. Converting flnc.bam file into a fastq.gz file for further analyses"
    bam2fastq -o "${fastq_flnc}/${sample}" ${flnc}
    echo "fastq.gz finished"
    
done

echo "All samples processed successfully."
