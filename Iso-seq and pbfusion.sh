#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=256GB
#SBATCH --time=5-12
#SBATCH --export=ALL
#
#The following code was adapted from the original one: https://github.com/SziKayLeung/Whole_Transcriptome_Paper/tree/f3190b3d6edcbf9703f69acb58fa6f210968f5bb
#The dataset used for this particular Iso-Seq and pbfusion pipeline was taken from this study: https://doi.org/10.1016/j.celrep.2021.110022
#ISO-SEQ WORKFLOW
#CSS analysis workflow --> produce HiFi reads from subreads. Use version v7 or v13 depending on the CHEMISTRY of the samples
#Use ccs7 for samples with OLD CHEMISTRY 
#Use ccs13 for samples with NEW CHEMISTRY

export ccsv_7="/home/mbocchi/smrtlink/install/smrtlink-release_7.0.1.66975/bundles/smrttools/install/smrttools-release_7.0.1.66768/smrtcmds/bin/ccs"
export lima_v7="/home/mbocchi/smrtlink/install/smrtlink-release_7.0.1.66975/bundles/smrttools/install/smrttools-release_7.0.1.66768/smrtcmds/bin/lima"
export isoseq3_v7="/home/mbocchi/smrtlink/install/smrtlink-release_7.0.1.66975/bundles/smrttools/install/smrttools-release_7.0.1.66768/smrtcmds/bin/isoseq3"
export isoseq_v13="/home/mbocchi/SMRT_ROOT/install/smrtlink-release_13.0.0.214433/smrtcmds/bin/isoseq"

FASTA="/zfs/jacobs/Colette/Maria/SRA_PRJNA664117/primer.fasta"
cat $FASTA

#List the versions for each tool used (use the right version depending on the CHEMISTRY of the samples)
echo "1.ccs version"
$ccsv_7 / $ccsv_13 --version

echo "2.lima version"
$lima_v7 / $lima_v13 --version

echo "3.isoseq3 version (used for refine)"
$isoseq3_v7 / $isoseq3_v13 --version

echo "4.isoseq cluster2 version"
$isoseq_v13 cluster2 --version



$ccs7 --reportFile --minPasses=3 --minPredictedAccuracy=0.99 --logFile mylog.txt subreads.bam ccs.bam




#Primer removal and demultiplexing

Version used : paper




echo "5.Running ccs"

$ccsv_7 --reportFile --minPasses=3 --minPredictedAccuracy=0.99 --logFile mylog.txt FetalE_subreads.bam FetalE_ccs.bam

echo "6.Finish ccs"

echo "7.Running lima"

$lima_v7 FetalE_ccs.bam $FASTA FetalE.fl.bam --isoseq --dump-clips --dump-removed

echo "8.Finish lima"

echo "9.Running isoseq refine"

$isoseq3_v7 refine FetalE.fl.primer_5p--primer_3p.bam $FASTA FetalE.flnc.bam --require-polya

echo "10.Finish isoseq refine"

echo "11.Running isoseq cluster2"

$isoseq_v13 cluster2 FetalE.flnc.bam FetalE.clustered.bam --singletons --log-file clusterlog.txt

echo "12.Finish isoseq cluster2"


Isoform clustering
FLNC reads are clustered by their sequencing similarity to produce isoform consensus sequences. This step is the last step of 
Iso-Seq analysis if no genome is provided


