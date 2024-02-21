# PacBio-analysis
Analysis of long-read RNAseq data for the identification of fusion transcripts

The R code found in the repository (long_read_seq_analysis.R) was used to detect fusion transcripts, previously identified by the FusionCatcher tool in short-read RNAseq data, in long-read RNAseq data. 
The code was implemented to look for matches (with a maximum mismatch of 2) of each fusion transcript in long reads RNAseq data, employing the "Biostring" library in R. For each match found, a new FASTA file is created, so that all the matches (long-reads) found for each fusion transcripts are stored in a single FASTA file.
