#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --time=15-12
#SBATCH --export=ALL

############################################################
This dataset contains samples with old and new chemistry 
(see paper (PMID: 34788620)), thefore the script can be run 
with one set of sample or the otehr, making sure to use the 
right ccs version.
############################################################

# SET ALL THE VARIABLES NECESSARY FOR THIS SCRIPT 
# Set variables for Iso-Seq workflow
# Depending on the set of samples, use either ccs_3 (old chemistry samples) or ccs_8 (new chemistry samples)
export ccs_3="/path/to/smrtlink/install/smrtlink-release_7.0.1.66975/bundles/smrttools/install/smrttools-release_7.0.1.66768/smrtcmds/bin/ccs"
export ccs_8="/path/to/SMRT_ROOT/install/smrtlink-release_13.0.0.214433/smrtcmds/bin/ccs" 
# Set variables for the pbfusions workflow
h38mmi="/path/to/GRCh38.p13/GRCh38.p13.genome.mmi"
gtfv38="/path/to/GRCh38.p13/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.bin.xz"
# Set the variable for the primer FASTA file
FASTA="/path/to/SRA_PRJNA664117/primer.fasta"

# CREATE ALL THE NECESSARY DIRECTORIES
# Create a directory to store the flnc.fastq.gz files that will be used in further analyses
mkdir -p "/path/to/fastq_flnc/"
fastq_flnc="/path/to/fastq_flnc/"


# LISTING THE VERSION OF EACH TOOL USED IN THE PIPELINE
# Use the appropriate ccs version
echo "ccs version"
$ccs_3 --version/$ccs_8 --version 
echo "lima version"
lima --version
echo "isoseq version (used for refine)"
isoseq refine --version
echo "bam2fastq version"
bam2fastq --version
echo "pbfusion discover version"
pbfusion discover --version

# Use the set of samples depending on the version used
# ARRAY OF SAMPLES NAMES --> ccs version 3 (old chemistry)
samples=("FetalF" "FetalE" "FetalD" "FetalC1S" "FetalB1S" "FetalGS") 
# ARRAY OF SAMPLES NAMES --> ccs version 8 (new chemistry)
samples=("AdultA" "AdultB" "AdultC1" "AdultC2" "AdultD" "FetalA" "FetalB2" "FetalC2") 

# RUN THE ISO-SEQ AND PBFUSION WORKFLOW FOR EACH SAMPLE
for sample in "${samples[@]}"; do
    echo "Creating directory for sample: $sample"
    mkdir -p "./${sample}_isoseq"
    sample_dir="./${sample}_isoseq"

    echo "Processing sample: $sample"
    
    # RUNNING THE ISO-SEQ WORKFLOW

    echo "1. Running ccs to get HiFi reads...
     #--min-passes=3 --min-rq=0.99 are the defults of ccs to get HiFi reads (default for ccs_8)
    $ccs_3 --reportFile "${sample_dir}/${sample}_ccs_report.txt" --minPasses=3 --minPredictedAccuracy=0.99 --logFile "mylog_${sample_dir}/${sample}.txt" ${sample}_subreads.bam "${sample_dir}/${sample}.ccs.ba
    $ccs_8 ${sample}_subreads.bam "${sample_dir}/${sample}.ccs.bam" --report-file "${sample_dir}/${sample}_ccs_report.txt" --logFile mylog_"${sample_dir}/${sample}.txt"
    echo "ccs finished"
    
    #Change directory to the sample directory
    cd "$sample_dir"
    
    echo "2. Count the number of HiFi reads..."
    samtools view ${sample}.ccs.bam | wc -l > hifi_reads_${sample}.txt
    echo "count finished"

    echo "3. Running lima"
    lima ${sample}.ccs.bam $FASTA ${sample}.fl.bam --isoseq --dump-clips --dump-removed
    echo "lima finished"

    echo "4. Running isoseq refine for trimming PolyA tails and concatemer removal..."
    isoseq refine ${sample}.fl.primer_5p--primer_3p.bam $FASTA ${sample}.flnc.bam --require-polya
    echo "isoseq refine finished"

    echo "5. Count the number of flnc reads..."
    samtools view ${sample}.flnc.bam | wc -l > flnc_${sample}.txt
    echo "count finished"

    
    # RUNNING THE pbfsion WORKFLOW
    #Set the varaible for the flnc file
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
    
    # We want to use the flnc.bam file in further anlyses, but in the flnc.fastq.gz format 
    echo "8. Converting flnc.bam file into a fastq.gz file for further analyses"
    bam2fastq -o "${fastq_flnc}/${sample}" ${flnc}
    echo "fastq.gz finished"
    
    # Change directory to the working directory to process the other samples
    cd ..

done

echo "All samples processed successfully."
