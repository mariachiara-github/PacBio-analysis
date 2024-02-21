# PacBio-analysis
Analysis of long-read RNAseq data for the identification of fusion transcripts

The R code found in the repository (long_read_seq_analysis.R) was used to detect fusion transcripts, previously identified by the FusionCatcher tool in short-read RNAseq data, in long-read RNAseq data. 
The code was implemented to look for matches (maximum mistmatch=2) of each fusion transcript in long reads RNAseq data, and every time a match is found, it is saved in a new FASTA file. 
