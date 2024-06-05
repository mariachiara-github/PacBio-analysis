#Match short reads fusions in long reads 

library("Biostrings")
library(seqRFLP)         #to tranform a dataframe in FASTA file


#Read and transform the FASTA file in dataframe
fasta <- readDNAStringSet("/zfs/jacobs/Colette/Maria/SRA_PRJNA664117/SRR12660778/AdultB_fasta_flnc.fasta")     #name or path of the FASTA file (long reads)
head(fasta)

#Create a dataframe from the long-reads FASTA file
fasta.df <- data.frame(sequences=as.character(fasta))


#Read the FASTA file with the intereseting fusions
fusions <- readDNAStringSet("/zfs/jacobs/Colette/Maria/SRA_PRJDB15555/KANSL1_NAIP_NSF.fasta")
#Extract the sequence name
seq_name <- names(fusions)
#Extract the sequence
sequence <- as.character(fusions)

setwd("/zfs/jacobs/Colette/Maria/SRA_PRJNA664117/SRR12660778/FindMatches_AdultB/")

#For loop to look for matches of each filtered fusion transcript in each long read
for (i in 1:length(sequence)){                                              
  print(i)
  seq <- (sequence[i])
  print(seq)                                  #sequence 
  fusion_name <- seq_name[i]  
  match <- vmatchPattern(seq,fasta,max.mismatch=2)    #find the matchings with a mx mistmatch between the 2 sequences of 2 
  match.df <- data.frame(unlist(match))               #create a dataframe with the matches
  print(match.df)
  seq.matched <-subset(fasta.df,rownames(fasta.df) %in% match.df$names) 
  #print(seq.matched)
  if (nrow(seq.matched) != 0) {                       #if a match is found create a FASTA file, and name it as the fusion transcript name 
    fasta.match<- dataframe2fas(seq.matched, file = paste(fusion_name,".fasta", sep =''))
  }
}

