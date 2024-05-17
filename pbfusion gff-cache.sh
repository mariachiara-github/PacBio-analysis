#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=256GB
#SBATCH --time=5-12
#SBATCH --export=ALL
#
#GENCODE (https://www.gencodegenes.org/human/release_38.html)
#GTF file --> Comprehensive gene annotation, regions = ALL (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz)
#FASTA file --> (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz)

echo "1.pbfusion gff-cache version "
pbfusion gff-cache --version 

echo "2.Serializing the input gtf file"
pbfusion gff-cache \
    --gtf gencode.v38.chr_patch_hapl_scaff.annotation.gtf \
    --gtf-out gencode.v38.chr_patch_hapl_scaff.annotation.gtf.bin
       
echo "3.Compressing the annotation bin file"
xz -c -9 --extreme gencode.v38.chr_patch_hapl_scaff.annotation.gtf.bin > gencode.v38.chr_patch_hapl_scaff.annotation.gtf.bin.xz
#ls -oh gencode.v38.chr_patch_hapl_scaff.annotation.gtf.bin* | awk '{print $4, $NF}'

echo "4.Done"
