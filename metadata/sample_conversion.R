#!/usr/bin/env Rscript

library('xlsx')
library('tidyverse')
library('data.table')

#setwd('~/projects/biogemma/expression/raw_reads')

#allfiles=list.files('.','*.fastq.gz',recursive=T,full.names=F)
#samples=sapply(seq(1,length(allfiles)),function(x) strsplit(allfiles[x],"[.]")[[1]][1])
#sample_names=data.frame(sample=samples,stringsAsFactors=F)
#fwrite(sample_names,'sample_names.txt',row.names=F,col.names=F,quote=F)


#filenames=list.files('../batch_1','*.fastq.gz',full.names=F)
#samples=sapply(seq(1,length(filenames)),function(x) strsplit(filenames[x],"[.]")[[1]][1])
sample_table=fread('metadata/BG_WD_OPT_sample_submission_FIXED_WELL.csv',data.table=F)
#sample_table=fread('metadata/BG_WD_OPT_sample_submission.csv',data.table=F)
seq_table=fread('metadata/BG_sequencing_log.csv',data.table=F)

#samples_names=data.frame(sample=samples,stringsAsFactors=F)
#fwrite(samples_names,'sample_names.txt',row.names=F,col.names=F,quote=F)

seq_table$plate_num=sapply(seq(1,dim(seq_table)[1]),function(x) as.integer(strsplit(seq_table$plate,"P")[[x]][2]))

sample_names=fread('metadata/sample_names.txt',data.table = F,header=F)
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

#let=c("H","G","F","E","D","C","B","A")
# A1 - A12
let=c("A","B","C","D","E","F","G","H")
well_code=c()
for(i in let){
  for(j in seq(1,12)){
    well_code=c(well_code,paste0(i,j))
  }
}
# A1-H1

#let=c("A","B","C","D","E","F","G","H")
#well_code2=c()
#for(i in seq(1,12)){
#  for(j in let){
#    well_code2=c(well_code2,paste0(j,i))
#  }
#}

# H1-H12
#let=c("H","G","F","E","D","C","B","A")
#well_code3=c()
#for(i in let){
#  for(j in seq(1,12)){
#    well_code3=c(well_code3,paste0(i,j))
#  }
#}
# H1-A1
#let=c("H","G","F","E","D","C","B","A")
#well_code4=c()
#for(i in seq(1,12)){
#  for(j in let){
#    well_code4=c(well_code4,paste0(j,i))
#  }
#}

oldv=1:96
newv=c(t(matrix(1:96,nc=8))[8:1,])

# Old way - A1 starts
sample_table$old_well_no=match(sample_table$well,well_code)

# Option 1 - H1 starts
new_well_code=well_code[newv]
sample_table$new_well_no=match(sample_table$well,new_well_code)
#sample_table$well_check=match(sample_table$new_well,well_code)
# Option 2 
#new_well_code2=well_code[match(well_code,new_well_code)]
#sample_table$new_well_no2=match(sample_table$new_well,new_well_code2)

# Option 3

#sample_table$new_well_no3=match(sample_table$well,well_code3)

#sample_table$old_well_no=match(sample_table$well,well_code)

#sample_names$old_well=sample_table[match(sample_names$batch_number,sample_table$old_well_no),]$well
#sample_names$new_well=sample_table[match(sample_names$batch_number,sample_table$new_well_no),]$new_well
#sample_names$plate=seq_table[match(sample_names$lane,seq_table$Lane),]$plate_num

#sample_name=c()
#genotype=c()
#exp=c()
#for(i in seq(1,dim(sample_names)[1])){
#  line=sample_names[1,]
#  m=sample_table[sample_table$plate==line$plate & sample_table$new_well==line$new_well,]$sample
#}
sample_names$plate_well=paste0(sample_names$plate,'_',sample_names$batch_number)
sample_table$plate_well=paste0(sample_table$plate,'_',sample_table$new_well_no)
sample_merge=merge(sample_names,sample_table,by=c('plate_well'))
fwrite(sample_merge,'metadata/BG_sequencing_sample_conversion_table_FIXED.txt',row.names=F,quote=F,sep='\t')

#sample_complete=sample_merge[!is.na(sample_merge$genotype) & sample_merge$genotype!="",]
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

metadata=fread('metadata/BG_sequencing_sample_conversion_table_FIXED.txt',data.table=F)
linenames=fread('metadata/Lines_Names_BALANCE.txt',data.table=F)
#rownames(K)=gsub("-",".",rownames(K))
#etadata$dh_genotype=linenames[,match]

meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
meta$plate_well_read = paste0(meta$plate,'_',meta$batch_number,'_',meta$read)

metadata$plate_well_read = paste0(metadata$plate.x,'_',metadata$batch_number,'_',metadata$read)
metadata=metadata[metadata$plate_well_read %in% meta$plate_well_read,]
metadata$dh_genotype=linenames[match(metadata$genotype,linenames$CODE_HYBRIDE),]$BGA_ID
#metadata=metadata[!is.na(metadata$dh_genotype),]
metadata$fastq_name=sapply(seq(1,nrow(metadata)),function(x) strsplit(metadata[x,]$seq_name,'/')[[1]][2])
metadata$sample_name=sapply(seq(1,nrow(metadata)),function(x) strsplit(metadata[x,]$fastq_name,'_R')[[1]][1])

fwrite(metadata,'metadata/BG_completed_sample_list_FIXED.txt',row.names=F,quote=F,sep='\t')


# What samples did I not analyze?
more=metadata[!(metadata$seq_name %in% og$seq_name),]
more=more[!(more$experiment %in% c("","Blanco_lab")),]
more=more[more$dh_genotype!="",]
fwrite(more,'metadata/unanalyzed_samples.csv',row.names=F,quote=F,sep=',')
# test

metadata=fread('metadata/BG_completed_sample_list_FIXED.txt',data.table=F)
#metadata$plate_well_read = paste0(metadata$plate,'_',metadata$batch_number,'_',metadata$read)

og=fread('metadata/BG_completed_sample_list.txt',data.table=F)
og$plate_well_read = paste0(og$plate,'_',og$batch_number,'_',og$read)
check=fread('metadata/match_between_expression_pairs.csv',data.table=F)
correct=0
correct_i=c()

subcheck=check %>% group_by(dh_1,expt,expt2) %>% top_n(1, cor)
for(i in 1:nrow(subcheck)){
	row=subcheck[i,]
	dh1=row$dh_1
	exp1=row$expt
	dh2=row$dh_2
	exp2=row$expt2
	old=og[og$dh_genotype==dh1 & og$experiment==exp1 & og$read==1,]
	old2=og[og$dh_genotype==dh2 & og$experiment==exp2 & og$read==1,]
	old_well1=paste0(old$plate,'_',old$batch_number,'_',old$read)
	old_well2=paste0(old2$plate,'_',old2$batch_number,'_',old2$read)

	new_geno1=metadata[metadata$plate_well_read==old_well1,]$dh_genotype
	new_geno2=metadata[metadata$plate_well_read==old_well2,]$dh_genotype
	#if()
	if(new_geno1==new_geno2){
		correct=correct+1
		correct_i=c(correct_i,i)
	}

}
print(correct)
# old well 1000 are correct
# What plates are in the ones that are correct? - In all of these, they are the same 
# individual in both time points
c1=check[correct_i,]
c1$plate1=metadata[metadata$experiment==c1$expt & metadata$dh_genotype==c1$dh_1,]$plate.x

metadata=metadata[,complete.cases(metadata),]
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
metadata=metadata[metadata$experiment %in% times,]
plate1=c()
plate2=c()
for(i in 1:nrow(c1)){
	row=c1[i,]
	p1=metadata[metadata$experiment==row$expt & metadata$dh_genotype==row$dh_1 & metadata$read==1,]$plate.x
	p2=metadata[metadata$experiment==row$expt2 & metadata$dh_genotype==row$dh_2 & metadata$read==1,]$plate.x
	plate1=c(plate1,p1)
	plate2=c(plate2,p2)
}
c1$plate1=plate1
c1$plate2=plate2

#new_well_no 361 correct
c2=check[correct_i,]

c2=check[correct_i,]

metadata=metadata[complete.cases(metadata),]
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
metadata=metadata[metadata$experiment %in% times,]
plate1=c()
plate2=c()
for(i in 1:nrow(c2)){
	row=c2[i,]
	p1=metadata[metadata$experiment==row$expt & metadata$dh_genotype==row$dh_1 & metadata$read==1,]$plate.x
	p2=metadata[metadata$experiment==row$expt2 & metadata$dh_genotype==row$dh_2 & metadata$read==1,]$plate.x
	plate1=c(plate1,p1)
	plate2=c(plate2,p2)
}
c2$plate1=plate1
c2$plate2=plate2
#new_well_no2 24 correct

#new_well_no 10

#well_check 40 correct


files=fread('/group/runciegrp2/Projects/Runcie/biogemma/updated_file_genotype_info_FIXED_genotypes.csv',data.table=F)
newfiles=files[files$fq_size>0 & files$bam_size==0,]

split0=strsplit(newfiles$fq_file_R1,"/")
newfiles$batch=sapply(seq(1,nrow(newfiles)),function(x) strsplit(newfiles[x,]$fq_file_R1,'/')[[1]][8])
newfiles=newfiles[!is.na(newfiles$fq_file_R1),]
newfiles[newfiles$batch=="18048-69",]$batch="batch_2"
#sample_names$batch=sapply(seq(1,dim(sample_names)[1]),function(x) split0[[x]][1])
drops=c('BG-P24-well-86_S182_L002',"BG-P12-well-46_S47_L006",'BG-P13-well-72_S170_L007')
newfiles=newfiles[!(newfiles$sample_name %in% drops),]
fwrite(newfiles,'metadata/updated_file_genotype_info_FIXED_genotypes_fastq.csv',row.names=F,quote=F,sep=',')

unique(newfiles$batch)


newfiles$fastq_name=sapply(seq(1,nrow(newfiles)),function(x) strsplit(newfiles[x,]$fq_file_R1,'/')[[1]][9])


#newfiles$fastq_name=sapply(seq(1,nrow(newfiles)),function(x) strsplit(newfiles[x,]$fq_file_R1,'/')[[1]][9])
#newfiles$sample_name=sapply(seq(1,nrow(newfiles)),function(x) strsplit(newfiles[x,]$fastq_name,'_R')[[1]][1])
pattern1=c('18048-85-01-03','18048-85-04-11','redos')
pattern2=c('batch_1','18048-85-04-11-Extra')
pattern3=c('batch_2')
pattern4=c('18048-85-lane12-16')

pattern5=c("18048-69") # "well89_S347_L007_R1_001.fastq.gz"

sample=c()
well_no=c()
lane=c()
read=c()
plate=c()
for(p in seq(1,length(split0))){
  if(split0[[p]][8] %in% pattern1){
    splita=strsplit(split0[[p]][9],'_')
    lane=c(lane,as.integer(strsplit(splita[[1]][3],"L")[[1]][2]))
    read=c(read,as.integer(strsplit(splita[[1]][4],"R")[[1]][2]))
    sample=c(sample,as.integer(strsplit(splita[[1]][2],"S")[[1]][2]))
    splitb=strsplit(splita[[1]][1],'-')
    well_no=c(well_no,as.integer(splitb[[1]][4]))
    plate=c(plate,as.integer(strsplit(splitb[[1]][2],'P')[[1]][2]))
  }
  else if(split0[[p]][8] %in% pattern2){
    splita=strsplit(split0[[p]][9],'_')
    lane=c(lane,as.integer(strsplit(splita[[1]][3],"L")[[1]][2]))
    read=c(read,as.integer(strsplit(splita[[1]][4],"R")[[1]][2]))
    sample=c(sample,as.integer(strsplit(splita[[1]][2],"S")[[1]][2]))
    splitb=strsplit(splita[[1]][1],'-')
    well_no=c(well_no,as.integer(splitb[[1]][4]))
    if(split0[[p]][8] == "batch_1"){
      plate=c(plate,ifelse(lane[p]==1,4,ifelse(lane[p]==2,7,8)))
    }else{
      plate=c(plate,NA)
    }
  }
  else if(split0[[p]][8] %in% pattern4){
    splita=strsplit(split0[[p]][9],'_')
    lane=c(lane,as.integer(strsplit(splita[[1]][3],"L")[[1]][2]))
    read=c(read,as.integer(strsplit(splita[[1]][4],"R")[[1]][2]))
    sample=c(sample,as.integer(strsplit(splita[[1]][2],"S")[[1]][2]))
    splitb=strsplit(splita[[1]][1],'-')
    if(length(grep('Phos',splitb[[1]][1]))==1){
      well_no=c(well_no,as.integer(splitb[[1]][3]))
      plate=c(plate,29)
    }else{
      well_no=c(well_no,as.integer(splitb[[1]][4]))
      plate=c(plate,as.integer(strsplit(splitb[[1]][2],'P')[[1]][2]))
    }
  }
  else{
    splita=strsplit(split0[[p]][9],'_')
    lane=c(lane,as.integer(strsplit(splita[[1]][3],"L")[[1]][2]))
    read=c(read,as.integer(strsplit(splita[[1]][4],"R")[[1]][2]))
    sample=c(sample,as.integer(strsplit(splita[[1]][2],"S")[[1]][2]))
    well_no=c(well_no,as.integer(strsplit(splita[[1]][1],'well')[[1]][2]))
    plate=c(plate,ifelse(lane[p]==7,10,11))
  }
}

newfiles$sample=sample
newfiles$well_no=well_no
newfiles$read=read
newfiles$lane=lane
newfiles$plate=plate

bar75=sample_table[sample_table$barcode=="ILLSINHA75",]
sample_names[is.na(sample_names$plate),]$batch_number=70
fwrite(sample_names,'sample_seq_info.txt',quote=F,sep='\t',row.names=F)
