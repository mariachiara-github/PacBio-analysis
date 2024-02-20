library("Biostrings")
library(seqRFLP)     #to tranform a dataframe in FASTA file

fasta <- readDNAStringSet("SRR23409825_1.fasta")  #FASTA file
head(fasta)
fasta.df <- data.frame(as.character(fasta))


fusions <- readDNAStringSet("short_reads_fusion.txt")
seq_name = names(fusions)
sequence = paste(fusions)
df.fusions <- data.frame(seq_name, sequence)


for (i in 1:length(sequence)){
  print(i)
  seq <- (sequence[i])
  print(sequence[i])
  fusion_name <- seq_name[i]
  match <- vmatchPattern(seq,fasta,max.mismatch=2)
  #print(match)
  match.df <- data.frame(as.character(unlist(match)))
  seq.matched <-subset(fasta.df,rownames(fasta.df) %in% rownames(match.df))
  if (nrow(seq.matched) != 0) {
    fasta.match<- dataframe2fas(seq.matched, file = paste(fusion_name,"mismatch",".fasta"))
  }
  
}


