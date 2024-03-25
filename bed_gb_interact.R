#Read BED file of pbfusion
bed <- read.table("C:/Users/Maria Chiara/Downloads/isoseq.breakpoints.groups.bed")
col <- bed$V11

##### CREATE THE COLUMNS FOR THE INTERACT TRACK

#Create a new column with the NAME of the fusion
new_column <- c()
for (i in col){
  d <- data.frame(strsplit(i ,";"))
  v <- d[3,]
  new_column <- append(new_column,v)
}

name_column <-c()
for (i in new_column){
  new <- gsub("GN=","",i)
  name <- gsub(",","-",new)
  name_column <- append(name_column,name)
}

#Create a new column with the ID NAME of the gene

new_column_ID <- c()
for (i in col){
  d_ID <- data.frame(strsplit(i ,";"))
  v_ID <- d_ID[4,]
  new_column_ID <- append(new_column_ID,v_ID)
}

name_column_ID <-c()
for (i in new_column_ID){
  new_ID <- gsub("GI=","",i)
  name_ID <- gsub(","," ",new_ID)
  name_column_ID <- append(name_column_ID,name_ID)
}

#Create a column for the COLOR 

color <- c()
for (i in 1:nrow(bed)){
  print(i)
  f <- bed$V9[i]
  r <- bed$V10[i]
  if (f == r){
    color <- append(color,"#7A67EE")
  }
  if (f != r){
    color <- append(color,"#F73BB9")
  }
}

#Create a column for the SCORE --> read count column

score <- c()
for (i in col){
  d <- data.frame(strsplit(i ,";"))
  v <- d[1,]
  score <- append(score,v)
}

score_column <-c()
for (i in score){
  s <- gsub("RC=","",i)
  score_column <- append(score_column,s)
}

#Create a column for the VALUE --> LOW = 0, MEDIUM = 1

value <- bed$V8

value_column <- c()
for(i in value){
  if(i == "LOW"){
    val<- gsub("LOW",0,i)
    value_column <- append(value_column,val)
  }
  if (i == "MEDIUM"){
    val<- gsub("MEDIUM",1,i)
    value_column <- append(value_column,val)
  }
}

#Create a column for EXP --> PacBio

exp_col <- rep("PacBio",177614)


#Create a new dataframe in BED format to upload to GenomeBrowser
#COLUMN: chrom	chromStart	chromEnd	name	score	value	exp	color	sourceChrom	sourceStart	sourceEnd	sourceName	sourceStrand	targetChrom	targetStart	targetEnd	targetName	targetStrand
#First we create a dataframe without the sourceName and targetName --> do this after the cleaning. These two columns will be created after cleaning the dataset from fusions between more than 2 genes

chrm <- bed$V1
chromStart <- bed$V2
chromEnd <- bed$V3
name <- name_column
name_ID <- name_column_ID
score <- score_column
value <- value_column
exp <- exp_col
color <- color
sourceChrom <- bed$V1
sourceStart <- bed$V2
sourceEnd <- bed$V3
#sourceName <- source_name
sourceStrand <- bed$V9
targetChrom <- bed$V4
targetStart <- bed$V5
targetEnd <- bed$V6
#targetName <- target_name
targetStrand <- bed$V10

bed_genome_browser <- data.frame(chrm,chromStart,chromEnd,name,name_ID,score,value, exp, color, sourceChrom,sourceStart,sourceEnd,sourceStrand,targetChrom,targetStart,targetEnd,targetStrand)

#Remove fusions between 3 or more genes

bed_gb <- data.frame()
bed_gb_clean <- data.frame()

for (i in 1:nrow(bed_genome_browser)){
  print(i)
  n <- bed_genome_browser$name[i]
  d <- data.frame(strsplit(n ,"-"))
  row_n <- nrow(d)
  if (row_n <= 2){
    bed_gb <- bed_genome_browser[i,]
    bed_gb_clean <- rbind(bed_gb_clean,bed_gb)
  }
}


#Create a column for the source and the target names 

name_col_clean <- bed_gb_clean$name
name_ID_clean <- bed_gb_clean$name_ID
source_name <- c()
target_name <- c()
for (i in 1:nrow(bed_gb_clean)){
  gene_n<- bed_gb_clean$name[i]
  gene_name <- data.frame(strsplit(gene_n ,"-"))
  ID_n <- bed_gb_clean$name_ID[i]
  ID_name <- data.frame(strsplit(ID_n," "))
  source_name <- append(source_name, paste(gene_name[1,],"(",ID_name[1,],")",sep ="" ))
  target_name <- append(target_name, paste(gene_name[2,],"(",ID_name[2,],")",sep ="" ))
  }

#CREATE FINAL DATASET
chrm <- bed$V1
chromStart <- bed$V2
chromEnd <- bed$V3
name <- name_column
name_ID <- name_column_ID
score <- score_column
value <- value_column
exp <- exp_col
color <- color
sourceChrom <- bed$V1
sourceStart <- bed$V2
sourceEnd <- bed$V3
sourceName <- source_name
sourceStrand <- bed$V9
targetChrom <- bed$V4
targetStart <- bed$V5
targetEnd <- bed$V6
targetName <- target_name
targetStrand <- bed$V10


bed_genome_browser_clean <- cbind(bed_gb_clean[,1:4],bed_gb_clean[,6:12],source_name,bed_gb_clean[,13:16],target_name, bed_gb_clean[,17])
colnames(bed_genome_browser_clean)[18] <- "targetStrand"

#Create a txt file --> upload it on Genome Browser. The dataset is also cleaned from the fusions found in alt haplotypes (for now, will fix this)
write.table(bed_genome_browser_clean, "C:/Users/Maria Chiara/Downloads//bed_genome_browser_clean.txt", row.names=FALSE, quote=FALSE)

bed_gb_subset <- bed_genome_browser_clean[6583:143315,]
write.table(bed_gb_subset, "C:/Users/Maria Chiara/Downloads//bed_bg_subset.txt", row.names=FALSE, quote=FALSE)
