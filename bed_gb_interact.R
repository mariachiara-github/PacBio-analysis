#####
#Read BED file of pbfusion (isoseq.breakpoints.groups.bed)
bed <- read.table("C:/Users/../isoseq.breakpoints.groups.bed")  #directory of the isoseq.breakpoints.groups. bed file created by the pbfusion pipeline

#Filter the BED file 
#1ST FILTER: remove all the fusions inter-chromosomal

bed_clean_1 <- data.frame()
chr1 <- bed$V1
chr2 <- bed$V4

for (i in 1:nrow(bed)){
  print(i)
  if (chr1[i] == chr2[i]){
    b <- bed[i,]
    bed_clean_1 <- rbind(bed_clean_1,b)
  }
}

write.table(bed_clean_1, "C:/.../bed_clean_1.txt", row.names=FALSE, quote=FALSE)

##### 
#CREATE THE COLUMNS FOR THE INTERACT TRACK

col <- bed_clean_1$V11             #column with info about the fusions detected

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
  name <- gsub(",","_",new)
  name_column <- append(name_column,name)
}

#Create a new column with the ID NAME of the gene (gene name from reference GTF)

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

#Create a column for the COLOR (the fusion will be purple if strand 1 and strand 2 are the same; + + or - -; the fusion will be pink otherwise)

color <- c()
for (i in 1:nrow(bed_clean_1)){
  print(i)
  f <- bed_clean_1$V9[i]
  r <- bed_clean_1$V10[i]
  if (f == r){
    color <- append(color,"#7A67EE")
  }
  if (f != r){
    color <- append(color,"#F73BB9")
  }
}

#Create a column for the SCORE --> RC (Number of reads supporting breakpoints) 

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

value <- bed_clean_1$V8

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

#Create a column for EXP --> PacBio (name of the experiment)

exp_col <- rep("PacBio",177614)


#Create a new dataframe in BED format to upload to GenomeBrowser
#COLUMN: chrom	chromStart	chromEnd	name	score	value	exp	color	sourceChrom	sourceStart	sourceEnd	sourceName	sourceStrand	targetChrom	targetStart	targetEnd	targetName	targetStrand
#First we create a dataframe without the sourceName and targetName --> do this after the cleaning. These two columns will be created after cleaning the dataset from fusions between more than 2 genes

chrm <- bed_clean_1$V1
chromStart <- bed_clean_1$V2
chromEnd <- bed_clean_1$V3
name <- name_column
name_ID <- name_column_ID
score <- score_column
value <- value_column
exp <- exp_col
color <- color
sourceChrom <- bed_clean_1$V1
sourceStart <- bed_clean_1$V2
sourceEnd <- bed_clean_1$V3
#sourceName <- source_name
sourceStrand <- bed_clean_1$V9
targetChrom <- bed_clean_1$V4
targetStart <- bed_clean_1$V5
targetEnd <- bed_clean_1$V6
#targetName <- target_name
targetStrand <- bed_clean_1$V10

bed_genome_browser <- data.frame(chrm,chromStart,chromEnd,name,name_ID,score,value, exp, color, sourceChrom,sourceStart,sourceEnd,sourceStrand,targetChrom,targetStart,targetEnd,targetStrand)

#Remove fusions between 2 or more genes and write in a new data frame all the interactions between more than 2 genes (bed_morethan2_clean dataframe)
#Create a data frame just with fusions between 2 genes(bed_gb_clean data frame)

bed_gb <- data.frame()
bed_gb_clean <- data.frame()
bed_morethan2_clean <- data.frame()

for (i in 1:nrow(bed_genome_browser)){
  print(i)
  n <- bed_genome_browser$name[i]
  d <- data.frame(strsplit(n ,"_"))
  row_n <- nrow(d)
  if (row_n <= 2){
    bed_gb <- bed_genome_browser[i,]
    bed_gb_clean <- rbind(bed_gb_clean,bed_gb)
  }
   if (row_n > 2){
    bed_morethan2_clean <- rbind(bed_morethan2_clean,bed_genome_browser[i,])
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
write.table(bed_genome_browser_clean, "C:/.../bed_genome_browser_clean.txt", row.names=FALSE, quote=FALSE)

bed_gb_subset <- bed_genome_browser_clean[6583:143315,]
write.table(bed_gb_subset, "C:/.../bed_bg_subset.txt", row.names=FALSE, quote=FALSE)
