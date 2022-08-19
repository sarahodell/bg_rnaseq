#!/usr/bin/env Rscript

library('xlsx')
library('tidyverse')
library('data.table')

setwd('~/projects/biogemma/expression/raw_reads')

allfiles=list.files('.','*.fastq.gz',recursive=T,full.names=F)
samples=sapply(seq(1,length(allfiles)),function(x) strsplit(allfiles[x],"[.]")[[1]][1])
sample_names=data.frame(sample=samples,stringsAsFactors=F)
fwrite(sample_names,'sample_names.txt',row.names=F,col.names=F,quote=F)
#filenames=list.files('../batch_1','*.fastq.gz',full.names=F)
#samples=sapply(seq(1,length(filenames)),function(x) strsplit(filenames[x],"[.]")[[1]][1])

sample_table=fread('BG_WD_OPT_sample_submission.csv',data.table=F)
seq_table=fread('../metadata/BG_sequencing_log.csv',data.table=F)

#samples_names=data.frame(sample=samples,stringsAsFactors=F)
#fwrite(samples_names,'sample_names.txt',row.names=F,col.names=F,quote=F)

seq_table$plate_num=sapply(seq(1,dim(seq_table)[1]),function(x) as.integer(strsplit(seq_table$plate,"P")[[x]][2]))

sample_names=fread('sample_names.txt',data.table = F,header=F)
names(sample_names)=c("seq_name")
split0=strsplit(sample_names$seq_name,"/")
sample_names$batch=sapply(seq(1,dim(sample_names)[1]),function(x) split0[[x]][1])

unique(sample_names$batch)
#[1] "18048-85-01-03"       "18048-85-04-11-Extra" "18048-85-04-11"
#[4] "18048-85-lane12-16"   "batch_1"              "batch_2"

#batch_1 pattern: "batch_1/18048FL-06-01-01_S1_L001_R1_001"
#batch_2 pattern:  "batch_2/well1_S259_L007_R1_001"
#18048-85-01-03 pattern: "18048-85-01-03/BG-P12-well-1_S2_L006_R1_001"
#18048-85-04-11-Extra pattern: 18048-85-04-11-Extra/Extra-GATGACAA_S1_L001_R1_001"
#18048-85-04-11 pattern: "18048-85-04-11/BG-P15-well-1_S1_L001_R1_001"
#18048-85-lane12-16 pattern: "18048-85-lane12-16/BG-P23-well-1_S1_L001_R1_001"
pattern1=c('18048-85-01-03','18048-85-04-11')
pattern2=c('batch_1','18048-85-04-11-Extra')
pattern3=c('batch_2')
pattern4=c('18048-85-lane12-16')


sample=c()
batch_number=c()
lane=c()
read=c()
plate=c()
for(p in seq(1,length(split0))){
  if(split0[[p]][1] %in% pattern1){
    splita=strsplit(split0[[p]][2],'_')
    lane=c(lane,as.integer(strsplit(splita[[1]][3],"L")[[1]][2]))
    read=c(read,as.integer(strsplit(splita[[1]][4],"R")[[1]][2]))
    sample=c(sample,as.integer(strsplit(splita[[1]][2],"S")[[1]][2]))
    splitb=strsplit(splita[[1]][1],'-')
    batch_number=c(batch_number,as.integer(splitb[[1]][4]))
    plate=c(plate,as.integer(strsplit(splitb[[1]][2],'P')[[1]][2]))
  }
  else if(split0[[p]][1] %in% pattern2){
    splita=strsplit(split0[[p]][2],'_')
    lane=c(lane,as.integer(strsplit(splita[[1]][3],"L")[[1]][2]))
    read=c(read,as.integer(strsplit(splita[[1]][4],"R")[[1]][2]))
    sample=c(sample,as.integer(strsplit(splita[[1]][2],"S")[[1]][2]))
    splitb=strsplit(splita[[1]][1],'-')
    batch_number=c(batch_number,as.integer(splitb[[1]][4]))
    if(split0[[p]][1] == "batch_1"){
      plate=c(plate,ifelse(lane[p]==1,4,ifelse(lane[p]==2,7,8)))
    }else{
      plate=c(plate,NA)
    }
  }
  else if(split0[[p]][1] %in% pattern4){
    splita=strsplit(split0[[p]][2],'_')
    lane=c(lane,as.integer(strsplit(splita[[1]][3],"L")[[1]][2]))
    read=c(read,as.integer(strsplit(splita[[1]][4],"R")[[1]][2]))
    sample=c(sample,as.integer(strsplit(splita[[1]][2],"S")[[1]][2]))
    splitb=strsplit(splita[[1]][1],'-')
    if(length(grep('Phos',splitb[[1]][1]))==1){
      batch_number=c(batch_number,as.integer(splitb[[1]][3]))
      plate=c(plate,29)
    }else{
      batch_number=c(batch_number,as.integer(splitb[[1]][4]))
      plate=c(plate,as.integer(strsplit(splitb[[1]][2],'P')[[1]][2]))
    }
  }
  else{
    splita=strsplit(split0[[p]][2],'_')
    lane=c(lane,as.integer(strsplit(splita[[1]][3],"L")[[1]][2]))
    read=c(read,as.integer(strsplit(splita[[1]][4],"R")[[1]][2]))
    sample=c(sample,as.integer(strsplit(splita[[1]][2],"S")[[1]][2]))
    batch_number=c(batch_number,as.integer(strsplit(splita[[1]][1],'well')[[1]][2]))
    plate=c(plate,ifelse(lane[p]==7,10,11))
  }
}

sample_names$sample=sample
sample_names$batch_number=batch_number
sample_names$read=read
sample_names$lane=lane
sample_names$plate=plate

bar75=sample_table[sample_table$barcode=="ILLSINHA75",]
sample_names[is.na(sample_names$plate),]$batch_number=70
fwrite(sample_names,'sample_seq_info.txt',quote=F,sep='\t',row.names=F)

#sample_names$sample_number=sapply(seq(1,dim(sample_names)[1]),function(x) as.integer(strsplit(split[[x]][2],"S")[[1]][2]))


let=c("A","B","C","D","E","F","G","H")
well_code=c()
for(i in let){
  for(j in seq(1,12)){
    well_code=c(well_code,paste0(i,j))
  }
}
sample_table$well_no=match(sample_table$well,well_code)

sample_names$well=sample_table[match(sample_names$batch_number,sample_table$well_no),]$well
#sample_names$plate=seq_table[match(sample_names$lane,seq_table$Lane),]$plate_num

sample_name=c()
genotype=c()
exp=c()
for(i in seq(1,dim(sample_names)[1])){
  line=sample_names[1,]
  m=sample_table[sample_table$plate==line$plate & sample_table$well==line$well,]$sample
}

sample_merge=left_join(sample_names,sample_table,by=c('plate','well'))
fwrite(sample_merge,'BG_sequencing_sample_conversion_table.txt',row.names=F,quote=F,sep='\t')


sample_complete=sample_merge[complete.cases(sample_merge),]
sample_complete %>% group_by(experiment) %>% count()
## A tibble: 7 x 2
# Groups:   experiment [7]
#  experiment     n
#  <chr>      <int>
#1 OPT_0711     146
#2 OPT_0728     338
#3 WD_0712      338
#4 WD_0718      530
#5 WD_0720      721
#6 WD_0725      722
#7 WD_0727      721

# divide number by 2 for number of genotypes

geno_count=sample_complete %>% group_by(genotype) %>% count()
summary(geno_count$n)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  6.000   8.000  10.000   9.794  10.000  20.000

# divide number by 2 for number of genotypes
