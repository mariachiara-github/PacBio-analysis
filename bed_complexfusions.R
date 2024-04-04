#Create a data frame with the name of the fusion and the reads associated with it 
#The following columns are taken from the bed_gb_interact.R code

name <- name_column
name_ID <- name_column_ID
read_name <- bed_clean_1$V12             #RN (Names of reads contributing support)

bed_reads <- data.frame(name,name_ID,read_name)
write.table(bed_reads, "C:/.../bed_reads.txt", row.names=FALSE, quote=FALSE)  #directory where you want to save the file

#Create a subset from the bed_reads containing just fusions between more than 2 genes --> COMPLEX FUSIONS

bed_morethan2_reads <- data.frame()

for (i in 1:nrow(bed_reads)){
  print(i)
  n <- bed_reads$name[i]
  d <- data.frame(strsplit(n ,"_"))
  row_n <- nrow(d)
  if (row_n > 2){
    bed_morethan2_reads <- rbind(bed_morethan2_reads,bed_reads[i,])
  }
  
}

reads<- bed_morethan2_reads$read_name
write.table(reads, "C:/.../reads.txt", row.names=FALSE, quote=FALSE)
#reads <- read.table("C:/.../reads.txt")   #used to upload the file on R without creating it every time

#Extract just the RN name from the read_name column
read_col <- c()
for (i in reads){
  rn <- gsub("RN=","",i)
  rn_name <- gsub(";ON=.;CB=.","",rn)
  read_col <- append(read_col,rn_name)
}

#Create a file where each line is the name of ONE read name 
sep <- c()
reads_BED <- c()
for (i in read_col){
  if (grepl(",",i)){
    sep <- str_split(i, ",")
    for (j in sep){
      reads_BED <- append(reads_BED,j)
    }
  }
  if (grepl(",",i) == FALSE){
    reads_BED <- append(reads_BED,i)
  }
}

write.table(reads_BED, "C:/.../reads_BED.txt", row.names=FALSE, quote=FALSE)   

#####
#SBATCH CODE 
#To create a BED file of the reads_BED reads, use the following sbatch script
#The aligned.bed file is the BED file derived from the aligned.bam file used in pipeline_pbfusion script

sort -u reads_BED.txt > unique_reads_BED.txt
cat aligned.bed | grep -wf unique_reads_BED.txt > bed_complexfusions.bed

#The fusions_not_pairwise.bed can be directly uploaded on the custom track of Genome Browser

