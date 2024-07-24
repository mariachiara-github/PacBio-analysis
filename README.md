# PacBio Analyses for fusion transcripts validation in RNA sequencing HiFi data
We collected RNA sequencing brain data from three publicly available datasets to validate the 717 fusion transcripts (FTs), detected by FusionCatcher, in PacBio long reads. The ```flnc.bam``` (full-length, non-concatemer reads) files for each sample were used for this purpose. 

### Fusion Transcripts Detections in Long Reads PacBio HiFi Data by minimap2
First, we aimed to validate FTs in long reads using an alignment tool,[```minimap2(v2.28-r1209)```](https://github.com/lh3/minimap2), that we initially used to index the 717 FTs, and then to map the FLNC long reads on the FT sequences, allowing for a maximum of 3 mismatches in a sequence of 86bp. The ```minimap2_alignment_pacbio.py``` Python code was employed for this purpose. The results of the ```minimap2``` were recorded for the first and second datasets' samples. The samples of the third dataset were solely used to search for __KANSL1_ARL17A/B__, __KANSL1_LRRC37A3__, __NAIP_OCLN__ and __NSF_LRRC37A3__ FTs.

### Validation of FTs Detected by minimap2 in pbfusion 
A _de novo_ FTs detection in Iso-Seq HiFi data was performed on the samples of the first and second dataset using the  [```pbfusion (v0.4.1) ```](https://github.com/PacificBiosciences/pbfusion/tree/master?tab=readme-ov-file). The Iso-Seq HiFi reads (FLNC) were first aligned to the human genome (GRCh38.p13.genome.fa) by [```pbmm2(v1.14.99)```](https://github.com/PacificBiosciences/pbmm2). Then ```pbfusion``` was run on aligned reads, allowing the emission of both LOW and MEDIUM fusions. A [```python(v3.11.8)```](https://www.python.org/) script was next developed to look for matches of FTs, identified by minimap2, in FTs detected by pbfusion. The matching was based on the read name of each FTs and using ```pandas``` library.

### Selected datasets 

1. [PacBio Cerebellum)](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/): the FLNC files of the 3 cerebellum samples were directly available for download. ```isoseq_pbfusion_pipeline_dataset1.sh``` was used to download and process the ```flnc.bam``` samples for the ```pbfusion``` pipeline.

2. [Bioproject PRJDB15555](https://www.ncbi.nlm.nih.gov/sra/?term=PRJDB15555). RNA was extracted from 4 individuals in the temporal cortex, hypothalamus, and cerebellum for a pool of 12 samples in total. CCS reads for each of the 12 samples were disclosed, and the FLNC reads were generated following the [Iso-Seq pipeline](https://isoseq.how/). HiFi reads (predicted accuracy ≥Q20) were first extracted from the CCS reads of each sample using [```extracthifi (v3.1.1)```](https://github.com/PacificBiosciences/extracthifi?tab=readme-ov-file) PacBio tool, and with the subsequent Iso-Seq workflow, primers were removed using [```lima (v2.9.0)```](https://lima.how/) tool and ```isoseq refine (v4.0.0)``` was used for trimming PolyA tails and concatemer removal.

3. [Bioproject PRJNA664117](https://www.ncbi.nlm.nih.gov/sra?term=PRJNA664117&cmd=DetailsSearch): the ```subreads.bam``` files were disclosed in this dataset; therefore, the [```ccs (v3.4.1 and v8.0.0)```](https://ccs.how/) tool was first used to extract HiFi reads, and FLNC reads were then generated using ```lima``` and ```isoseq refine``` Iso-Seq tools as previously described. In particular, the version of the  ```ccs ``` tool used was different for some samples (old chemistry samples), which required an older version of the tool (v3.4.1). Moreover, to generate PacBio HiFi data, the parameters of the ccs tools were set such that the minimum number of passes was set to 3 and the minimum predicted accuracy for a read was set to 0.99 (only reads expected to be 99% accurate are emitted). 


*Note: for each of the 3 datasets, a different ```isoseq_pbfusion_pipeline``` was used to better deal with the different parameters of the ```Iso-Seq pipeline``` required for each specific dataset. However, the ```pbfusion``` pipeline and parameters were the same for all the samples of the 3 datasets.

### Codes of this repository

#### ```minimap2_alignment_pacbio.py```
The first script is used to validate short-read FTs in long-reads HiFi (FLNC) data, which takes as input ```flnc.fastq.gz``` files for each sample and returns the fusion transcripts matched between short and long reads (FASTA files for each fusion found in each sample). Moreover, it returns a CSV file for each sample which includes the transcript name and the read names associated with each transcript.

#### ```find_pbfusion.py```

The second script processes a list of CSV files (output of minimap2_alignment_pacbio.py) and a list of ```pbfusion``` BED files (one of the outputs of the isoseseq_pbufusion pipeline). For each CSV and pbfusion file (of the same sample) this script finds matches between the 2 files based on the PacBio read name. The matches are stored as CSV files for each sample. 
Moreover, this script checks how many interchromosomal fusions are found in each final dataframe and if the gene name of the short-read fusion transcript matches the gene name of the long-read fusion (pbfusion name). 

#### ```fusions_per_indvidual.py``` and ```fusions_per_indvidual_pbfusion.py```

These two codes were developed to check for each of the 717 fusion transcripts how many of 6 individuals had it. This analysis was performed both on the fusions matched by the ```minimap2_alignment_pacbio.py``` script and on the fusions matched from the ```find_pbfusion.py```



#### ```get_long_reads.py```

#### ```intersection_brain_areas.py``` 
