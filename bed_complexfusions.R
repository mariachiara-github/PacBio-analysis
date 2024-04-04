#Create a dataframe for BAM --> just get the reads and the name 

name <- name_column
name_ID <- name_column_ID
read_name <- bed$V12

bed_reads <- data.frame(name,name_ID,read_name)
write.table(bed_reads, "C:/Users/Maria Chiara/Downloads/bed_reads.txt", row.names=FALSE, quote=FALSE)

#Write a dataset with all the fusions between more than 2 genes 


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
write.table(reads, "C:/Users/Maria Chiara/Downloads/reads.txt", row.names=FALSE, quote=FALSE)
reads <- read.table("C:/Users/Maria Chiara/Downloads/reads.txt")

read_col <- c()
for (i in reads){
  rn <- gsub("RN=","",i)
  rn_name <- gsub(";ON=.;CB=.","",rn)
  read_col <- append(read_col,rn_name)
}

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

write.table(reads_BED, "C:/Users/Maria Chiara/Downloads/reads_BED.txt", row.names=FALSE, quote=FALSE)

