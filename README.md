# PacBio Analyses for fusion transcripts validation in RNA sequencing HiFi data

Two different analyses were developed to validate a list of 717 short-read fusion transcripts (FTs), detected by FusionCatcher, in long-read PacBio data (Iso-Seq data). 

1. The first one aimed to validate FTs in long reads using an alignment tool,```minimap2```. The ```__minimap2_alignment_pacbio.py__```Python code was employed for this purpose.

4. The second analysis aimed to perform a _de novo_ analysis for finding fusion transcripts in long read data employing the __pbfusion__ tool [pbfusion](https://github.com/PacificBiosciences/pbfusion/tree/master?tab=readme-ov-file).

## Selected datasets to perform the two analyses

1. [First PacBio dataset](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/). Three cerebellum samples: the flnc.bam files(full-length non-concatemer reads) of each sample were used to perform the analyses (the flnc.bam files are provided in the database).

2. Second dataset [Bioproject 975746] (https://www.ncbi.nlm.nih.gov/bioproject/975746). RNA was extracted from 4 individuals in the temporal cortex, hypothalamus, and cerebellum, 12 samples in total. The CCS files were disclosed, so the Iso-Seq pipeline was carried out to get the flnc.bam files.

## FIRST ANALYSIS: find matches of short read fusions in long reads

The Python code __minimap2_alignment_pacbio.py__ was implemented to align long-reads data to a reference (the 717 short-reads fusions). For each alignment found, a FASTA file is created for that transcript so that all the matches (long-reads) found for each fusion transcript are stored in a single FASTA file.

Parameters used for minimap2 (https://github.com/lh3/minimap2, https://lh3.github.io/minimap2/minimap2.html): 
* indexing of the target sequence
  __"minimap2","-xmap-hifi", "-d", minimap2_index, transcripts_fasta__

## SECOND ANALYSIS: pbfusion analysis

The PacBio tool pbfusion was used to look for fusion transcripts in long reads Iso-Seq (PacBio) data. Differently from the analysis performed before, which looks for matches of short read fusions in long read data, this tool allows the discovery of novel fusion transcripts (https://github.com/PacificBiosciences/pbfusion/tree/master?tab=readme-ov-file). 
The __pipeline_pbfusion__ script was used to find fusion transcripts in the samples (the general pipeline used can be found in the GitHub repository of pbfusion); in particular a BED file is generated with all the information about the detected fusion. 

Since we want to visualize the fusions in genome browser, the BED file is cleaned and adapted to the interact track format of Genome Browser (https://genome.ucsc.edu/goldenPath/help/interact.html). This track format is used to displays pairwise interactions, therefore the BED file was cleaned from interactions between more than 2 genes (fusions which resulted from the joining of exons of 2 or more genes). The __bed_gb_interact.R__ script was used for this purpose

If we want to visualize the fusions with more than 2 genes detected by pbfusion, the __bed_complexfusions.R__ script can be used to create a BED track which can be uploaded on the Genome Browser. This script is used just for visualization purposes, it aims only to visualize all at once the reads classified as fusions between more than 2 genes.

## Iso-seq pipeline for Dataset 2

For the second dataset the ccs.bam files were disclosed. Therfore the first step is to generate the HiFi reads from the ccs.bam, and then perforom the Iso-Seq pipeline up to the flnc.bam generation. This file is then used to look for the fusions.

# CCS and Iso-Seq pipeline for Dataset 3

The first step is to generate HiFi PacBio long-read sequences from the subreads.bam files (sample files) provided by the dataset. The SMRT LINK V7.0 (https://downloads.pacbcloud.com/public/software/installers/smrtlink_7.0.1.66975.zip) version was used to perform the analysis (version compatible with the data of the third dataset). For more information about SMRTLink refer to the following pdf: https://www.pacb.com/wp-content/uploads/SMRT_Tools_Reference_Guide_v700.pdf.  

The CCS pipeline was performed with the default parameters of the tool (version: ccs 3.4.1). In particular, to generate PacBio HiFi data, the minimum number of passes was set to 3 and the minimum predicted accuracy for a read was set to 0.99 (only reads expected to be 99% accurate are emitted). 

## Citation

If you use minimap2 in your work, please cite:

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191

and/or:

Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574. doi:10.1093/bioinformatics/btab705

