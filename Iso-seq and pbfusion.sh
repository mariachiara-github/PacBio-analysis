#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=256GB
#SBATCH --time=5-12
#SBATCH --export=ALL
#The Iso-Seq workflow code was adapted from the original one: https://github.com/SziKayLeung/Whole_Transcriptome_Paper/tree/f3190b3d6edcbf9703f69acb58fa6f210968f5bb
#The dataset used for this particular Iso-Seq and pbfusion pipeline was taken from this study: https://doi.org/10.1016/j.celrep.2021.110022
#
#ISO-SEQ and pbfusion WORKFLOW 
#1. ccs--> produce HiFi reads from subreads. Use version SMRT LINK V7.0 (v7) or SMRT LINK V13.0 (v13) depending on the CHEMISTRY of the samples
  #Use ccs_v7 for samples with OLD CHEMISTRY (see IsoSeq_Processing/Human/All_Human_RawData_OldChem.txt from the GitHub mentioned above)
  #Use ccs_v13 for samples with NEW CHEMISTRY (see IsoSeq_Processing/Human/All_Human_RawData_NewChem.txt from the GitHub mentioned above)
  # <ccs movie.subreads.bam movie.ccs.bam>
#2. lima --> primer removal and demultiplexing.
  # <lima movie.ccs.bam barcoded_primers.fasta movieX.fl.bam --isoseq --dump-clips --dump-removed>
#3. isoseq refine --> trimming of poly(A) tail and rapid concatemer identification and removal.
  # <isoseq refine movieX.5p--3p.fl.bam movieX.flnc.bam --require-polya>
#4. isoseq cluster2 --> FLNC reads are clustered by their sequencing similarity to produce isoform consensus sequences. This step is the last step of Iso-Seq analysis if no genome is provided
  # <isoseq cluster2 flnc.bam clustered.bam --singletons --log-file cluster2log.txt>
#5. pbmm2 --> align the Iso-seq HiFi reads (flnc.bam file, full-length non-concatemer reads)
#6. pbfusion gtf-cache

#7. pbfusions discover


#Export variables for Iso-Seq workflow
export ccs_v7="/home/mbocchi/smrtlink/install/smrtlink-release_7.0.1.66975/bundles/smrttools/install/smrttools-release_7.0.1.66768/smrtcmds/bin/ccs"
export ccs_v13="/home/mbocchi/SMRT_ROOT/install/smrtlink-release_13.0.0.214433/smrtcmds/bin/ccs"
export lima_v7="/home/mbocchi/smrtlink/install/smrtlink-release_7.0.1.66975/bundles/smrttools/install/smrttools-release_7.0.1.66768/smrtcmds/bin/lima"
export lima_v13="/home/mbocchi/SMRT_ROOT/install/smrtlink-release_13.0.0.214433/smrtcmds/bin/lima"
export isoseq3_v7="/home/mbocchi/smrtlink/install/smrtlink-release_7.0.1.66975/bundles/smrttools/install/smrttools-release_7.0.1.66768/smrtcmds/bin/isoseq3"
export isoseq_v13="/home/mbocchi/SMRT_ROOT/install/smrtlink-release_13.0.0.214433/smrtcmds/bin/isoseq"

#Export variables for the pbfusions workflow
export pbmm2_v13="/home/mbocchi/SMRT_ROOT/install/smrtlink-release_13.0.0.214433/smrtcmds/bin/pbmm2" 
export h38mmi="/zfs/jacobs/Colette/Maria/SRA_PRJNA664117/GRCh38.p13/GRCh38.p13.genome.mmi"
export gtfv38="/zfs/jacobs/Colette/Maria/SRA_PRJNA664117/GRCh38.p13/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.bin.xz"

#Path from the GitHub: Utilities/primer.fasta
FASTA="/zfs/jacobs/Colette/Maria/SRA_PRJNA664117/primer.fasta"
cat $FASTA

#List the versions for each tool used (use the right version depending on the CHEMISTRY of the samples)
echo "ccs version"
$ccsv_7 / $ccsv_13 --version
echo "lima version"
$lima_v7 / $lima_v13 --version
echo "isoseq3 version (used for refine)"
$isoseq3_v7 / $isoseq3_v13 --version
echo "isoseq cluster2 version"
$isoseq_v13 cluster2 --version

#For each step of the workflow use the correct version (either one or the other)

echo "1. Running ccs to get HiFi reads"
$ccsv_7 --reportFile --minPasses=3 --minPredictedAccuracy=0.99 --logFile mylog.txt subreads.bam ccs.bam
$ccsv_13 subreads.bam ccs.bam
echo "ccs finished"

echo "2. Running lima"
$lima_v7 ccs.bam $FASTA fl.bam --isoseq --dump-clips --dump-removed
echo "lima finished"

echo "3. Running isoseq refine"
$isoseq3_v7 refine fl.primer_5p--primer_3p.bam $FASTA flnc.bam --require-polya
$isoseq_v13
echo "isoseq refine finished"

echo "4. Running isoseq cluster2"
$isoseq_v13 cluster2 flnc.bam clustered_2.bam --singletons --log-file cluster2log.txt
echo "isoseq cluster2 finished"

echo "5. Running pbmm2 for flnc.bam reads"
echo "pbmm2 flnc.bam finished"

echo "6. Running pbmm2 for clustered.bam reads"
echo "pbmm2 clustered.bam finished"

echo "7. Running pbmm2 for clustered.bam reads"
echo "pbmm2 clustered.bam finished"

echo "8. Running pbfusion gtf-cache"
echo "pbfusion gtf-cache finished"

echo "9. Running pbfusion discover"
echo "pbfusion discover finished"







