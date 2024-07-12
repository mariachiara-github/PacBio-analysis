# PacBio Analyses for fusion transcripts validation in RNA sequencing HiFi data

We collected RNA sequencing brain data from three publicly available datasets to validate the 717 fusion transcripts (FTs), detected by FusionCatcher, in PacBio long reads. The flnc.bam(full-length, non-concatemer reads)files for each sample were used for this purpose. 

In the first dataset, PacBio Kinnex full-length RNA (Homo sapiens – Cerebellum) [1], the FLNC files of the 3 cerebellum samples were directly available for download. In the second dataset (DDBJ: PRJDB15555)[2], the CCS reads for each of the 12 samples were disclosed, and the FLNC reads were generated following the Iso-Seq pipeline [3]. HiFi reads (predicted accuracy ≥Q20) were first extracted from the CCS reads of each sample using extracthifi (v3.1.1) [4] PacBio tool, and with the subsequent Iso-Seq workflow, primers were removed using lima (v2.9.0) tool and isoseq refine (v4.0.0) was used for trimming PolyA tails and concatemer removal. From the third dataset [5] the subreads.bam files were disclosed; therefore, the ccs (v3.4.1 and v8.0.0) tool [6] was first used to extract HiFi reads, and FLNC reads were then generated using lima and isoseq refine Iso-Seq tools as previously described.                                                                          After we obtained all the flnc.bam files, minimap2 (v2.28-r1209) [7] was used to first index the 7171 FTs, and then to map the FLNC long reads on the FT sequences, allowing for a maximum of 3 mismatches in a sequence of 86bp. 
Validation of FTs Detected by minimap2 in pbfusion 
A de novo 

### Fusion Transcripts Detections in Long Reads PacBio HiFi Data by minimap2
First, we aimed to validate FTs in long reads using an alignment tool,```minimap2(v2.28-r1209)```, that we used to first index the 717 FTs, and then to map the FLNC long reads on the FT sequences, allowing for a maximum of 3 mismatches in a sequence of 86bp. The ```minimap2_alignment_pacbio.py```Python code was employed for this purpose. The results of the minimap2 were recorded for the samples of the first and the second dataset. The samples of the third dataset were solely used to search for KANSL1_ARL17A/B, KANSL1_LRRC37A3, NAIP_OCLN and NSF_LRRC37A3 FTs.

### Validation of FTs Detected by minimap2 in pbfusion 
A _de novo_ FTs detection in Iso-Seq HiFi data was performed on the samples of the first and second dataset using the  [```pbfusion (v0.4.1) ```](https://github.com/PacificBiosciences/pbfusion/tree/master?tab=readme-ov-file). The Iso-Seq HiFi reads (FLNC) were first aligned to the human genome (GRCh38.p13.genome.fa) by pbmm2(v1.14.99)[10]. Then pbfusion was run on aligned reads, allowing the emission of both LOW and MEDIUM fusions. A python(v3.11.8) script was next developed to look for matches of FTs, identified by minimap2, in FTs detected by pbfusion. The matching was based on the read name of each FTs and using pandas library.

2. The second analysis aimed to perform a _de novo_ analysis for finding fusion transcripts in long read data employing the __pbfusion__ tool [pbfusion](https://github.com/PacificBiosciences/pbfusion/tree/master?tab=readme-ov-file).

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

