# PacBio Analysis for fusion transcripts detection

SELECTED DATASETS

1. First PacBio dataset: https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/ , the three cerebellum samples were used for a the toy analysis. 

2. Second PacBio dataset (taken from a Bioproject): https://www.ncbi.nlm.nih.gov/bioproject/975746 . In the Bioproject 12 samples were employed; in particular from each sample, the long-read sequence data was extracted (Iso-seq, Pacific Biosciences (PacBio)) using RNA extracted from three regions of the human postmortem brain: temporal cortex, hypothalamus, and cerebellum. ((PACBIO human RNA seq brain) AND bioproject_sra[filter] NOT bioproject_gap[filter])






Analysis of long-read RNAseq data for the identification of fusion transcripts

The R code found in the repository (long_read_seq_analysis.R) was used to detect fusion transcripts, previously identified by the FusionCatcher tool in short-read RNAseq data, in long-read RNAseq data. 
The code was implemented to look for matches (with a maximum mismatch of 2) of each fusion transcript in long reads RNAseq data, employing the "Biostrings" library in R. For each match found, a new FASTA file is created, so that all the matches (long-reads) found for each fusion transcripts are stored in a single FASTA file.
