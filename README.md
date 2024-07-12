## PacBio Analyses for fusion transcripts validation in RNA sequencing HiFi data

We collected RNA sequencing brain data from three publicly available datasets to validate the 717 fusion transcripts (FTs), detected by FusionCatcher, in PacBio long reads. The ```flnc.bam``` (full-length, non-concatemer reads) files for each sample were used for this purpose. 

### Fusion Transcripts Detections in Long Reads PacBio HiFi Data by minimap2
First, we aimed to validate FTs in long reads using an alignment tool,```minimap2(v2.28-r1209)```, that we initially used to index the 717 FTs, and then to map the FLNC long reads on the FT sequences, allowing for a maximum of 3 mismatches in a sequence of 86bp. The ```minimap2_alignment_pacbio.py``` Python code was employed for this purpose. The results of the ```minimap2``` were recorded for the first and second datasets' samples. The samples of the third dataset were solely used to search for __KANSL1_ARL17A/B__, __KANSL1_LRRC37A3__, __NAIP_OCLN__ and __NSF_LRRC37A3__ FTs.

### Validation of FTs Detected by minimap2 in pbfusion 
A _de novo_ FTs detection in Iso-Seq HiFi data was performed on the samples of the first and second dataset using the  [```pbfusion (v0.4.1) ```](https://github.com/PacificBiosciences/pbfusion/tree/master?tab=readme-ov-file). The Iso-Seq HiFi reads (FLNC) were first aligned to the human genome (GRCh38.p13.genome.fa) by [```pbmm2(v1.14.99)```](https://github.com/PacificBiosciences/pbmm2). Then ```pbfusion``` was run on aligned reads, allowing the emission of both LOW and MEDIUM fusions. A [```python(v3.11.8)```](https://www.python.org/) script was next developed to look for matches of FTs, identified by minimap2, in FTs detected by pbfusion. The matching was based on the read name of each FTs and using ```pandas``` library.

### Selected datasets 

1. [PacBio Kinnex full-length RNA (Homo sapiens – Cerebellum)](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/): the FLNC files of the 3 cerebellum samples were directly available for download. ```isoseq_pbfusion_pipeline_dataset1.sh``` was used to 

2. Second dataset [Bioproject 975746](https://www.ncbi.nlm.nih.gov/bioproject/975746). RNA was extracted from 4 individuals in the temporal cortex, hypothalamus, and cerebellum for a pool of 12 samples in total. CCS reads for each of the 12 samples were disclosed, and the FLNC reads were generated following the Iso-Seq pipeline [3]. HiFi reads (predicted accuracy ≥Q20) were first extracted from the CCS reads of each sample using extracthifi (v3.1.1) [4] PacBio tool, and with the subsequent Iso-Seq workflow, primers were removed using lima (v2.9.0) tool and isoseq refine (v4.0.0) was used for trimming PolyA tails and concatemer removal.

3. subreads.bam files were disclosed; therefore, the ccs (v3.4.1 and v8.0.0) tool [6] was first used to extract HiFi reads, and FLNC reads were then generated using lima and isoseq refine Iso-Seq tools as previously described.



Parameters used for minimap2 (https://github.com/lh3/minimap2, https://lh3.github.io/minimap2/minimap2.html): 
* indexing of the target sequence
  __"minimap2","-xmap-hifi", "-d", minimap2_index, transcripts_fasta__



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

