# PacBio Analyses for fusion transcripts detection in RNAseq data

Two different analyses were developed to detected fusion transcripts in long read PacBio data (Iso-Seq data). 

1. The first one aimed to find MATCHES of previously identified short reads fusion transcripts (660 filtered fusions transcripts detected by FusionCatcher tool) in long reads. The __long_read_seq_analysis.R__ R code was employed for this purpose.

2. The second anaylsis aimed to perform a _de novo_ analysis for finding fusion transcripts in long read data employing the __pbfusion__ tool (https://github.com/PacificBiosciences/pbfusion/tree/master?tab=readme-ov-file).

# Selected datasets to perform the two analyses

1. First PacBio dataset: https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/. Three cerebellum samples were used for 'toy' analysis. Therefore these samples were used to perform the two analyses previously explained in order to develop and test the codes for the real analysis. The flnc.bam files(full-length non-concatemer reads) of each sample were used to perform the analyses (the flnc.bam files are provided in the database, so I did not perform the Iso-Seq pipeline on these samples).

2. Second dataset [Bioproject 975746] (https://www.ncbi.nlm.nih.gov/bioproject/975746). RNA was extracted from 4 individuals in the temporal cortex, hypothalamus, and cerebellum. The 12 samples were used to perform the real analyses. 

3. Third dataset [Bioproject PRJNA664117] (https://doi.org/10.1016/j.celrep.2021.110022). RNA was extracted from 11 samples. For each sample the subreads.bam files were provided, therefore the Iso-Seq pipeline was employed (https://isoseq.how/clustering/cli-workflow.html) to generate HiFi reads that can then be used to carry out the 2 analyses.


# FIRST ANALYSIS: find matches of short read fusions in long reads

The R code found in the repository (long_read_seq_analysis.R) was implemented to look for MATCHES (with a maximum mismatch of 2) of each short read fusion transcript in long reads RNA seq data, employing the "Biostrings" library in R. For each match found, a new FASTA file is created, so that all the matches (long-reads) found for each fusion transcripts are stored in a single FASTA file.

# SECOND ANALYSIS: pbfusion analysis

The PacBio tool pbfusion was used to look for fusion transcripts in long reads Iso-Seq (PacBio) data. Differently from the analysis performed before, which looks for matches of short read fusions in long read data, this tool allows the discovery of novel fusion transcripts (https://github.com/PacificBiosciences/pbfusion/tree/master?tab=readme-ov-file). 
The __pipeline_pbfusion__ script was used to find fusion transcripts in the samples (the general pipeline used can be found in the GitHub repository of pbfusion); in particular a BED file is generated with all the information about the detected fusion. 

Since we want to visualize the fusions in genome browser, the BED file is cleaned and adapted to the interact track format of Genome Browser (https://genome.ucsc.edu/goldenPath/help/interact.html). This track format is used to displays pairwise interactions, therefore the BED file was cleaned from interactions between more than 2 genes (fusions which resulted from the joining of exons of 2 or more genes). The __bed_gb_interact.R__ script was used for this purpose

If we want to visualize the fusions with more than 2 genes detected by pbfusion, the __bed_complexfusions.R__ script can be used to create a BED track which can be uploaded on the Genome Browser. This script is used just for visualization purposes, it aims only to visualize all at once the reads classified as fusions between more than 2 genes.

# Iso-seq pipeline for Dataset 2

For the second dataset the ccs.bam files were disclosed. Therfore the first step is to generate the HiFi reads from the ccs.bam, and then perforom the Iso-Seq pipeline up to the flnc.bam generation. This file is then used to look for the fusions.

# CCS and Iso-Seq pipeline for Dataset 3

The first step is to generate HiFi PacBio long-read sequences from the subreads.bam files (sample files) provided by the dataset. The SMRT LINK V7.0 (https://downloads.pacbcloud.com/public/software/installers/smrtlink_7.0.1.66975.zip) version was used to perform the analysis (version compatible with the data of the third dataset). For more information about SMRTLink refer to the following pdf: https://www.pacb.com/wp-content/uploads/SMRT_Tools_Reference_Guide_v700.pdf.  

The CCS pipeline was performed with the default parameters of the tool (version: ccs 3.4.1). In particular, to generate PacBio HiFi data, the minimum number of passes was set to 3 and the minimum predicted accuracy for a read was set to 0.99 (only reads expected to be 99% accurate are emitted). 

