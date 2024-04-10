#####
#Read BED file of pbfusion (isoseq.breakpoints.groups.bed)
bed <- read.table("C:/Users/../isoseq.breakpoints.groups.bed")  #directory of the isoseq.breakpoints.groups. bed file created by the pbfusion pipeline

#Filter the BED file 
#1ST FILTER: remove all the fusions inter-chromosomal and write them in a different data frame

bed_clean_1 <- data.frame()
bed_inter <- data.frame()
chr1 <- bed$V1
chr2 <- bed$V4

for (i in 1:nrow(bed)){
  print(i)
  if (chr1[i] == chr2[i]){
    b <- bed[i,]
    bed_clean_1 <- rbind(bed_clean_1,b)
  }
  if (chr1[i] != chr2[i]){
    inter <- bed[i,]
    bed_inter <- rbind(bed_inter,inter)
  }
}

write.table(bed_clean_1, "C:/.../bed_clean_1.txt", row.names=FALSE, quote=FALSE)
write.table(bed_clean_1, "C:/.../bed_inter.txt", row.names=FALSE, quote=FALSE)

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
  gene_name <- data.frame(strsplit(gene_n ,"_"))
  ID_n <- bed_gb_clean$name_ID[i]
  ID_name <- data.frame(strsplit(ID_n," "))
  source_name <- append(source_name, paste(gene_name[1,],"(",ID_name[1,],")",sep ="" ))
  target_name <- append(target_name, paste(gene_name[2,],"(",ID_name[2,],")",sep ="" ))
  }


#CREATE FINAL DATASET
chrm <- bed_gb_clean$chrm
chromStart <- bed_gb_clean$chromStart
chromEnd <- bed_gb_clean$chromEnd
name <- bed_gb_clean$name
score <- bed_gb_clean$score
value <- bed_gb_clean$value
exp <- bed_gb_clean$exp
color <- bed_gb_clean$color
sourceChrom <- bed_gb_clean$sourceChrom
sourceStart <- bed_gb_clean$sourceStart
sourceEnd <- bed_gb_clean$sourceEnd
sourceName <- source_name
sourceStrand <- bed_gb_clean$sourceStrand
targetChrom <- bed_gb_clean$targetChrom
targetStart <- bed_gb_clean$targetStart
targetEnd <- bed_gb_clean$targetEnd
targetName <- target_name
targetStrand <- bed_gb_clean$targetStrand

bed_genome_browser_clean <- data.frame(chrm,chromStart,chromEnd,name,score,value,exp,color,sourceChrom,sourceStart,sourceEnd,sourceName,sourceStrand,targetChrom,targetStart,targetEnd,targetName,targetStrand)

#Substitute the name of the chr with the one used by Genome Browser for the costumed tracks

#Get all the chr names (unique list)
c <- unique(bed_genome_browser_clean$chrm)

# Replace the name in the chr,sourcechr,and targetchr
library(dplyr)

bed1 <- bed_genome_browser_clean %>%
  mutate(chrm = recode(chrm, 'GL383534.2' = 'chr7_GL383534v2_alt', 'GL383545.1' = 'chr10_GL383545v1_alt', 'GL383574.1' = 'chr19_GL383574v1_alt','GL383575.2'= 'chr19_GL383575v2_alt', 'GL383580.2'='chr21_GL383580v2_alt','GL383582.2'='chr22_GL383582v2_alt','GL000250.2'='chr6_GL000250v2_alt',
                       'GL000251.2'= 'chr6_GL000251v2_alt','GL000252.2' = 'chr6_GL000252v2_alt','GL000195.1'='chrUn_GL000195v1','GL000219.1' ='chrUn_GL000219v1','GL000256.2'='chr6_GL000256v2_alt','GL000258.2' = 'chr17_GL000258v2_alt','GL949749.2'= 'chr19_GL949749v2_alt','GL949750.2'= 'chr19_GL949750v2_alt',
                       'GL949752.1' = 'chr19_GL949752v1_alt','GL949753.2'= 'chr19_GL949753v2_alt', 'KB663609.1' = 'chr22_KB663609v1_alt', 'KI270713.1' ='chr1_KI270713v1_random	', 'GL339449.2' = 'chr5_GL339449v2_alt', 'KI270772.1' = 'chr2_KI270772v1_alt', 'KI270790.1' = 'chr4_KI270790v1_alt', 
                       'KI270795.1' = 'chr5_KI270795v1_alt', 'KI270803.1' = 'chr7_KI270803v1_alt', 'KI270805.1' = 'chr7_KI270805v1_alt', 'KI270759.1' = 'chr1_KI270759v1_alt', 'KI270763.1' = 'chr1_KI270763v1_alt', 'KI270862.1' = 'chr17_KI270862v1_alt', 'KI270876.1' = 'chr22_KI270876v1_alt', 
                       'KI270877.1' = 'chr22_KI270877v1_alt', 'KI270878.1' = 'chr22_KI270878v1_alt', 'KI270879.1' = 'chr22_KI270879v1_alt', 'KI270897.1' = 'chr5_KI270897v1_alt', 'GL000253.2' = 'chr6_GL000253v2_alt', 'GL000254.2' = 'chr6_GL000254v2_alt', 'GL000255.2' = 'chr6_GL000255v2_alt', 
                       'KI270847.1' = 'chr14_KI270847v1_alt', 'KI270850.1' = 'chr15_KI270850v1_alt', 'KI270853.1' = 'chr16_KI270853v1_alt', 'KI270904.1' = 'chr12_KI270904v1_alt','KI270905.1' = 'chr15_KI270905v1_alt', 'KI270816.1' = 'chr8_KI270816v1_alt', 'KI270821.1' = 'chr8_KI270821v1_alt', 
                       'KI270830.1' = 'chr11_KI270830v1_alt', 'KI270832.1' = 'chr11_KI270832v1_alt', 'KI270844.1' = 'chr14_KI270844v1_alt', 'KN196475.1' = 'chr3_KN196475v1_fix', 'KN196479.1' = 'chr9_KN196479v1_fix', 'KN196481.1' = 'chr11_KN196481v1_fix', 'KN196484.1' = 'chr19_KN196484v1_fix', 
                       'KN538369.1' = 'chr12_KN538369v1_fix', 'KI270908.1' = 'chr17_KI270908v1_alt', 'KI270924.1' = 'chr3_KI270924v1_alt','KI270934.1' = 'chr3_KI270934v1_alt', 'KI270938.1' = 'chr19_KI270938v1_alt', 'KQ031389.1' = 'chr15_KQ031389v1_alt','KQ458383.1' = 'chr1_KQ458383v1_alt',
                       'KQ458386.1' = 'chr19_KQ458386v1_fix', 'KI270855.1' = 'chr16_KI270855v1_alt', 'KI270857.1' = 'chr17_KI270857v1_alt', 'KZ208911.1' = 'chr6_KZ208911v1_fix', 'KZ208912.1' = 'chr7_KZ208912v1_fix', 'KZ208913.1' = 'chr7_KZ208913v1_alt', 'KZ208914.1' = 'chr8_KZ208914v1_fix',
                       'KZ208915.1' = 'chr8_KZ208915v1_fix', 'KV766194.1' = 'chr6_KV766194v1_fix', 'KV766198.1' = 'chr17_KV766198v1_alt', 'KV880765.1' = 'chr7_KV880765v1_fix', 'KV880768.1' = 'chr16_KV880768v1_fix', 'KZ208906.1' = 'chr1_KZ208906v1_fix' , 'ML143371.1' = 'chr15_ML143371v1_fix',  
                       'KQ983256.1' = 'chr2_KQ983256v1_alt', 'KV575243.1' = 'chr5_KV575243v1_alt', 'KV575244.1' = 'chr5_KV575244v1_fix', 'KZ208916.1' = 'chr12_KZ208916v1_fix', 'KZ208920.1' = 'chr14_KZ208920v1_fix', 'KZ208922.1' = 'chr18_KZ208922v1_fix', 'KZ559109.1' = 'chr11_KZ559109v1_fix',
                       'ML143355.1' = 'chr10_ML143355v1_fix', 'ML143361.1' = 'chr12_ML143361v1_fix' , 'ML143362.1' = 'chr12_ML143362v1_fix', 'ML143367.1' = 'chr14_ML143367v1_fix', 'ML143370.1' = 'chr15_ML143370v1_fix'))


bed2 <- bed1 %>%
  mutate(sourceChrom = recode(sourceChrom, 'GL383534.2' = 'chr7_GL383534v2_alt', 'GL383545.1' = 'chr10_GL383545v1_alt', 'GL383574.1' = 'chr19_GL383574v1_alt','GL383575.2'= 'chr19_GL383575v2_alt', 'GL383580.2'='chr21_GL383580v2_alt','GL383582.2'='chr22_GL383582v2_alt','GL000250.2'='chr6_GL000250v2_alt',
                       'GL000251.2'= 'chr6_GL000251v2_alt','GL000252.2' = 'chr6_GL000252v2_alt','GL000195.1'='chrUn_GL000195v1','GL000219.1' ='chrUn_GL000219v1','GL000256.2'='chr6_GL000256v2_alt','GL000258.2' = 'chr17_GL000258v2_alt','GL949749.2'= 'chr19_GL949749v2_alt','GL949750.2'= 'chr19_GL949750v2_alt',
                       'GL949752.1' = 'chr19_GL949752v1_alt','GL949753.2'= 'chr19_GL949753v2_alt', 'KB663609.1' = 'chr22_KB663609v1_alt', 'KI270713.1' ='chr1_KI270713v1_random	', 'GL339449.2' = 'chr5_GL339449v2_alt', 'KI270772.1' = 'chr2_KI270772v1_alt', 'KI270790.1' = 'chr4_KI270790v1_alt', 
                       'KI270795.1' = 'chr5_KI270795v1_alt', 'KI270803.1' = 'chr7_KI270803v1_alt', 'KI270805.1' = 'chr7_KI270805v1_alt', 'KI270759.1' = 'chr1_KI270759v1_alt', 'KI270763.1' = 'chr1_KI270763v1_alt', 'KI270862.1' = 'chr17_KI270862v1_alt', 'KI270876.1' = 'chr22_KI270876v1_alt', 
                       'KI270877.1' = 'chr22_KI270877v1_alt', 'KI270878.1' = 'chr22_KI270878v1_alt', 'KI270879.1' = 'chr22_KI270879v1_alt', 'KI270897.1' = 'chr5_KI270897v1_alt', 'GL000253.2' = 'chr6_GL000253v2_alt', 'GL000254.2' = 'chr6_GL000254v2_alt', 'GL000255.2' = 'chr6_GL000255v2_alt', 
                       'KI270847.1' = 'chr14_KI270847v1_alt', 'KI270850.1' = 'chr15_KI270850v1_alt', 'KI270853.1' = 'chr16_KI270853v1_alt', 'KI270904.1' = 'chr12_KI270904v1_alt','KI270905.1' = 'chr15_KI270905v1_alt', 'KI270816.1' = 'chr8_KI270816v1_alt', 'KI270821.1' = 'chr8_KI270821v1_alt', 
                       'KI270830.1' = 'chr11_KI270830v1_alt', 'KI270832.1' = 'chr11_KI270832v1_alt', 'KI270844.1' = 'chr14_KI270844v1_alt', 'KN196475.1' = 'chr3_KN196475v1_fix', 'KN196479.1' = 'chr9_KN196479v1_fix', 'KN196481.1' = 'chr11_KN196481v1_fix', 'KN196484.1' = 'chr19_KN196484v1_fix', 
                       'KN538369.1' = 'chr12_KN538369v1_fix', 'KI270908.1' = 'chr17_KI270908v1_alt', 'KI270924.1' = 'chr3_KI270924v1_alt','KI270934.1' = 'chr3_KI270934v1_alt', 'KI270938.1' = 'chr19_KI270938v1_alt', 'KQ031389.1' = 'chr15_KQ031389v1_alt','KQ458383.1' = 'chr1_KQ458383v1_alt',
                       'KQ458386.1' = 'chr19_KQ458386v1_fix', 'KI270855.1' = 'chr16_KI270855v1_alt', 'KI270857.1' = 'chr17_KI270857v1_alt', 'KZ208911.1' = 'chr6_KZ208911v1_fix', 'KZ208912.1' = 'chr7_KZ208912v1_fix', 'KZ208913.1' = 'chr7_KZ208913v1_alt', 'KZ208914.1' = 'chr8_KZ208914v1_fix',
                       'KZ208915.1' = 'chr8_KZ208915v1_fix', 'KV766194.1' = 'chr6_KV766194v1_fix', 'KV766198.1' = 'chr17_KV766198v1_alt', 'KV880765.1' = 'chr7_KV880765v1_fix', 'KV880768.1' = 'chr16_KV880768v1_fix', 'KZ208906.1' = 'chr1_KZ208906v1_fix' , 'ML143371.1' = 'chr15_ML143371v1_fix',  
                       'KQ983256.1' = 'chr2_KQ983256v1_alt', 'KV575243.1' = 'chr5_KV575243v1_alt', 'KV575244.1' = 'chr5_KV575244v1_fix', 'KZ208916.1' = 'chr12_KZ208916v1_fix', 'KZ208920.1' = 'chr14_KZ208920v1_fix', 'KZ208922.1' = 'chr18_KZ208922v1_fix', 'KZ559109.1' = 'chr11_KZ559109v1_fix',
                       'ML143355.1' = 'chr10_ML143355v1_fix', 'ML143361.1' = 'chr12_ML143361v1_fix' , 'ML143362.1' = 'chr12_ML143362v1_fix', 'ML143367.1' = 'chr14_ML143367v1_fix', 'ML143370.1' = 'chr15_ML143370v1_fix'))

bed_genome_browser_clean_final <- bed2 %>%
  mutate(targetChrom = recode(targetChrom, 'GL383534.2' = 'chr7_GL383534v2_alt', 'GL383545.1' = 'chr10_GL383545v1_alt', 'GL383574.1' = 'chr19_GL383574v1_alt','GL383575.2'= 'chr19_GL383575v2_alt', 'GL383580.2'='chr21_GL383580v2_alt','GL383582.2'='chr22_GL383582v2_alt','GL000250.2'='chr6_GL000250v2_alt',
                       'GL000251.2'= 'chr6_GL000251v2_alt','GL000252.2' = 'chr6_GL000252v2_alt','GL000195.1'='chrUn_GL000195v1','GL000219.1' ='chrUn_GL000219v1','GL000256.2'='chr6_GL000256v2_alt','GL000258.2' = 'chr17_GL000258v2_alt','GL949749.2'= 'chr19_GL949749v2_alt','GL949750.2'= 'chr19_GL949750v2_alt',
                       'GL949752.1' = 'chr19_GL949752v1_alt','GL949753.2'= 'chr19_GL949753v2_alt', 'KB663609.1' = 'chr22_KB663609v1_alt', 'KI270713.1' ='chr1_KI270713v1_random	', 'GL339449.2' = 'chr5_GL339449v2_alt', 'KI270772.1' = 'chr2_KI270772v1_alt', 'KI270790.1' = 'chr4_KI270790v1_alt', 
                       'KI270795.1' = 'chr5_KI270795v1_alt', 'KI270803.1' = 'chr7_KI270803v1_alt', 'KI270805.1' = 'chr7_KI270805v1_alt', 'KI270759.1' = 'chr1_KI270759v1_alt', 'KI270763.1' = 'chr1_KI270763v1_alt', 'KI270862.1' = 'chr17_KI270862v1_alt', 'KI270876.1' = 'chr22_KI270876v1_alt', 
                       'KI270877.1' = 'chr22_KI270877v1_alt', 'KI270878.1' = 'chr22_KI270878v1_alt', 'KI270879.1' = 'chr22_KI270879v1_alt', 'KI270897.1' = 'chr5_KI270897v1_alt', 'GL000253.2' = 'chr6_GL000253v2_alt', 'GL000254.2' = 'chr6_GL000254v2_alt', 'GL000255.2' = 'chr6_GL000255v2_alt', 
                       'KI270847.1' = 'chr14_KI270847v1_alt', 'KI270850.1' = 'chr15_KI270850v1_alt', 'KI270853.1' = 'chr16_KI270853v1_alt', 'KI270904.1' = 'chr12_KI270904v1_alt','KI270905.1' = 'chr15_KI270905v1_alt', 'KI270816.1' = 'chr8_KI270816v1_alt', 'KI270821.1' = 'chr8_KI270821v1_alt', 
                       'KI270830.1' = 'chr11_KI270830v1_alt', 'KI270832.1' = 'chr11_KI270832v1_alt', 'KI270844.1' = 'chr14_KI270844v1_alt', 'KN196475.1' = 'chr3_KN196475v1_fix', 'KN196479.1' = 'chr9_KN196479v1_fix', 'KN196481.1' = 'chr11_KN196481v1_fix', 'KN196484.1' = 'chr19_KN196484v1_fix', 
                       'KN538369.1' = 'chr12_KN538369v1_fix', 'KI270908.1' = 'chr17_KI270908v1_alt', 'KI270924.1' = 'chr3_KI270924v1_alt','KI270934.1' = 'chr3_KI270934v1_alt', 'KI270938.1' = 'chr19_KI270938v1_alt', 'KQ031389.1' = 'chr15_KQ031389v1_alt','KQ458383.1' = 'chr1_KQ458383v1_alt',
                       'KQ458386.1' = 'chr19_KQ458386v1_fix', 'KI270855.1' = 'chr16_KI270855v1_alt', 'KI270857.1' = 'chr17_KI270857v1_alt', 'KZ208911.1' = 'chr6_KZ208911v1_fix', 'KZ208912.1' = 'chr7_KZ208912v1_fix', 'KZ208913.1' = 'chr7_KZ208913v1_alt', 'KZ208914.1' = 'chr8_KZ208914v1_fix',
                       'KZ208915.1' = 'chr8_KZ208915v1_fix', 'KV766194.1' = 'chr6_KV766194v1_fix', 'KV766198.1' = 'chr17_KV766198v1_alt', 'KV880765.1' = 'chr7_KV880765v1_fix', 'KV880768.1' = 'chr16_KV880768v1_fix', 'KZ208906.1' = 'chr1_KZ208906v1_fix' , 'ML143371.1' = 'chr15_ML143371v1_fix',  
                       'KQ983256.1' = 'chr2_KQ983256v1_alt', 'KV575243.1' = 'chr5_KV575243v1_alt', 'KV575244.1' = 'chr5_KV575244v1_fix', 'KZ208916.1' = 'chr12_KZ208916v1_fix', 'KZ208920.1' = 'chr14_KZ208920v1_fix', 'KZ208922.1' = 'chr18_KZ208922v1_fix', 'KZ559109.1' = 'chr11_KZ559109v1_fix',
                       'ML143355.1' = 'chr10_ML143355v1_fix', 'ML143361.1' = 'chr12_ML143361v1_fix' , 'ML143362.1' = 'chr12_ML143362v1_fix', 'ML143367.1' = 'chr14_ML143367v1_fix', 'ML143370.1' = 'chr15_ML143370v1_fix'))


write.table(bed_genome_browser_clean_final, "C:/.../bed_genome_browser_clean_final.txt", row.names=FALSE, quote=FALSE)


