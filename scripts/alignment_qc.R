#!/usr/bin/env Rscript

library('xlsx')
library('tidyverse')
library('data.table')
library('ggplot2')
library('stringr')

# Need to find
# mapping rates

# # of gene counts of at least X length (X = 20,50,100,200,300)
#Per timepoint how many reads per gene, how many counts per sample (mean)
#How many samples have at least a million reads
#How many genes that have at least 10 reads in at least X samples
#Histogram of number of samples that have at least 10 (across genes)
#allfiles=list.files('star/batch_2','bg_*_data/multiqc_fastqc.txt',recursive=T)

# change to whatever file is 
meta=fread('metadata/updated_file_genotype_info_FIXED_genotypes_2.csv',data.table=F)
dropexp=c('Blanco_lab','blanco_lab')
meta=meta[!(meta$experiment %in% dropexp),]

# check that samples from Blanco lab are really not ours
dropexp=c('Blanco_lab','blanco_lab')
unique(meta[meta$experiment %in% dropexp,]$plate)

empty=meta[meta$experiment=="",]
# 12 16 20 24 28  4  8

# check duplicates
#dups=meta[meta$match_type=="possible_duplicate",]
#dups=dups[!is.na(dups$final_genotype),]
#eg=interaction(dups$experiment,dups$final_genotype)

#new_plants=c()
#for(i in eg){
#	split=strsplit(as.character(i),'\\.')
#	e=split[[1]][1]
#	g=split[[1]][2]
#	sub=meta[meta$experiment==e & meta$final_genotype==g,]
#	sub=sub[!is.na(sub$final_genotype),]
#	if(length(unique(sub$final_X_Y))>1){
#		sdup=dups[dups$experiment==e & dups$final_genotype==g,]
#		xys=c(sub[sub$match_type!="duplicate",]$final_X_Y)
#		for(j in 1:nrow(sdup)){
#			row=sdup[j,]
#			xy=row$final_X_Y
#			if(!(xy %in% xys)){
#				new_plants=c(new_plants,row$sample_name)
#				xys=c(xys)
#			}
#		}
#	}
#}

keep=c('good','switch')
meta=meta[meta$match_type %in% keep,]
meta=meta[!is.na(meta$final_genotype),]


#meta=fread('metadata/BG_completed_sample_list_FIXED.txt',data.table=F)

# which ones that are missin pass2 are needed? from our experiment?
star=fread('metadata/BG_completed_sample_list_FIXED_2.txt',data.table=F)
dirs=c('18048-85-lane12-16','batch_2','18048-85-04-11','18048-85-01-03','batch_1','batch_2')

missing=c()
allsamples=c()
for(d in dirs){
	pass2=Sys.glob(sprintf('star/%s/*_pass2',d))
	pass2_s1=sapply(seq(1,length(pass2)),function(x) strsplit(pass2[x],'/')[[1]][[3]])
  	pass2_s2=sapply(seq(1,length(pass2)),function(x) strsplit(pass2_s1[x],'_pass2')[[1]][[1]])
	pass2u=Sys.glob(sprintf('star/update/%s/*-pass2',d))
  	pass2u_s1=sapply(seq(1,length(pass2u)),function(x) strsplit(pass2u[x],'/')[[1]][[4]])
  	pass2u_s2=sapply(seq(1,length(pass2u)),function(x) strsplit(pass2u_s1[x],'-pass2')[[1]][[1]])

	sample=union(pass2_s2,pass2u_s2)
	allsamples=c(allsamples,sample)
}
allsamples=unique(allsamples)
missing=star$sample_name[!(star$sample_name %in% allsamples)]

miss=star[star$sample_name %in% missing,]
miss=miss[miss$experiment!="",]
qual=fread('metadata/sample_names_qual.txt',data.table=F)
miss$batch_sample=paste0(miss$batch,'/',miss$fastq_name)
miss$qual=qual[match(miss$batch_sample,qual$V1),]$V2
miss=miss[miss$qual==TRUE,]

miss=fread('metadata/BG_completed_sample_list_FIXED_2.txt',data.table=F)
miss$qual=qual[match(miss$batch_sample,qual$V1),]$V2
miss=miss[miss$qual==TRUE,]
dropexp=c('Blanco_lab','blanco_lab','',NA)
miss=miss[!(miss$experiment %in% dropexp),]

miss=miss[!(miss$sample_name %in% missing),]
fwrite(miss,'metadata/BG_completed_sample_list_FIXED_2.txt',row.names=F,quote=F,sep='\t')

#meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
meta=fread('metadata/updated_file_genotype_info_FIXED_genotypes_2.csv',data.table=F)


#dirs=c('18048-85-lane12-16','batch_2','18048-85-04-11','18048-85-01-03','batch_1')

#d=dirs[2] #batch_2 only
  #sample=unique(samples[samples$batch==d,]$sample_name)

# Combine read counts into one file
dirs=c('18048-85-lane12-16','batch_2','18048-85-04-11','18048-85-01-03','batch_1','batch_2')

all_files=c()
for(d in dirs){
	if(d=="redos"){
		pass2u=Sys.glob(sprintf('star/update/%s/*-pass2/ReadsPerGene.out.tab',d))
  		pass2u_s1=sapply(seq(1,length(pass2u)),function(x) strsplit(pass2u[x],'/')[[1]][[4]])
  		pass2u_s2=sapply(seq(1,length(pass2u)),function(x) strsplit(pass2u_s1[x],'-pass2')[[1]][[1]])
		sample=unique(pass2u_s2)
	}else{
		pass2=Sys.glob(sprintf('star/%s/*_pass2/ReadsPerGene.out.tab',d))
		pass2_s1=sapply(seq(1,length(pass2)),function(x) strsplit(pass2[x],'/')[[1]][[3]])
  		pass2_s2=sapply(seq(1,length(pass2)),function(x) strsplit(pass2_s1[x],'_pass2')[[1]][[1]])
		pass2u=Sys.glob(sprintf('star/update/%s/*-pass2/ReadsPerGene.out.tab',d))
  		pass2u_s1=sapply(seq(1,length(pass2u)),function(x) strsplit(pass2u[x],'/')[[1]][[4]])
  		pass2u_s2=sapply(seq(1,length(pass2u)),function(x) strsplit(pass2u_s1[x],'-pass2')[[1]][[1]])
		sample=union(pass2_s2,pass2u_s2)
	}
	sample=sample[sample %in% meta$sample_name]
	sample=unique(sample)
  	for(s in sample){
  		if(s %in% pass2u_s2){
  			path=sprintf('star/update/%s/%s-pass2/ReadsPerGene.out.tab',d,s)
  		}else{
  			path=sprintf('star/%s/%s_pass2/ReadsPerGene.out.tab',d,s)
  		}
  		file=fread(path,data.table=F)
  		if(s==sample[1] & d==dirs[1]){
      		#column 3 is the first strand gene counts, dont know if this is what we want
      		file=file[,c(1,3)]
      		names(file)=c('Gene_ID',s)
      		files=file
      		rownames(files)=files$Gene_ID
      	}else{
      		#file=file[,c(3),drop=F]
      		file=file[,c(1,3)]
      		names(file)=c('V1',s)
      		#files=cbind(files,file)
      		files[,s]=file[match(files$Gene_ID,file$V1),s]
      	}
    }
    if(d==dirs[1]){
    	all_files=files
    }else{
    	all_files=cbind(all_files,files[,-1])
    }
}

fwrite(all_files,'star/update_all_Reads_per_Gene.txt',row.names=F,quote=F,sep='\t')

# of these, how many of them are my experiment and good quality?
# need to run qualimap

x=dim(all_files)[1]
y=dim(all_files)[2]
files=all_files[5:x,]

f=sapply(seq(2:y),function(x) is.numeric(files[1,x]))
drop=which(f==F)
files=files[,-drop]
y=dim(files)[2]

# For originals
mapping=fread('qc/star_multiqc_data/multiqc_star.txt',data.table=F)
passed=mapping[mapping$uniquely_mapped>=1e6,]
#passed=mapping[mapping$uniquely_mapped_percent>=50,]
passed$sample_name=sapply(seq(1,nrow(passed)),function(x) strsplit(passed$Sample[x],'_pass')[[1]][[1]])


#for updated
mappingu=fread('qc/star_update_multiqc_data/multiqc_star.txt',data.table=F)
passedu=mappingu[mappingu$uniquely_mapped>=1e6,]
#passed=mapping[mapping$uniquely_mapped_percent>=50,]
passedu$sample_name=sapply(seq(1,nrow(passedu)),function(x) strsplit(passedu$Sample[x],'-pass')[[1]][[1]])

passedsamps=union(passed$sample_name,passedu$sample_name)

meta=fread('metadata/updated_file_genotype_info_FIXED_genotypes_2.csv',data.table=F)
meta=meta[meta$sample_name %in% passedsamps,]

keep=c('good','switch')
meta=meta[meta$match_type %in% keep,]
meta %>% group_by(experiment) %>% count()
# A tibble: 7 × 2
# Groups:   experiment [7]
#  experiment     n
#  <chr>      <int>
#1 OPT_0711      36
#2 OPT_0728      40
#3 WD_0712       83
#4 WD_0718      144
#5 WD_0720      221
#6 WD_0725      110
#7 WD_0727      195

meta=meta[!is.na(meta$final_genotype),]

meta$batch=sapply(seq(1,nrow(meta)),function(x) strsplit(meta$fq_file_R1[x],'/')[[1]][8])
meta[meta$batch=="18048-69",]$batch="batch_2"
#meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
linenames=fread('metadata/Lines_Names_BALANCE.txt',data.table=F)
meta$dh_genotype=linenames[match(meta$final_genotype,linenames$CODE_HYBRIDE),]$BGA_ID


fwrite(meta,'metadata/samples_passed_genotype_check.txt',row.names=F,quote=F,sep='\t')

# look at qualimap output

meta=fread('metadata/samples_passed_genotype_check.txt',data.table=F)
qual=fread(;m)
#all_files=fread('star/update_all_Reads_per_Gene.txt',data.table=F)
timep=unique(meta$experiment)
for(t in timep){
	samples=unique(meta[meta$experiment==t,]$sample_name)
	submeta=meta[meta$experiment==t,]
	# get path to read counts file
	#files=c()
	for(s in samples){
		d=submeta[submeta$sample_name==s,'batch']
		path1=sprintf('star/%s/%s_pass2/ReadsPerGene.out.tab',d,s)
		path2=sprintf('star/update/%s/%s-pass2/ReadsPerGene.out.tab',d,s)
		if(file.exists(path1) & file.exists(path2)){
			path=path2
		}else if(file.exists(path1)){
			path=path1
		}else if(file.exists(path2)){
			path=path2
		}else{
			next
		}
		file=fread(path,data.table=F)
		if(s==samples[1]){
			#column 3 is the first strand gene counts, dont know if this is what we want
			file=file[,c(1,3)]
			names(file)=c('Gene_ID',s)
			files=file
			rownames(files)=files$Gene_ID
		}else{
			#file=file[,c(3),drop=F]
			file=file[,c(1,3)]
			names(file)=c('V1',s)
			#files=cbind(files,file)
			files[,s]=file[match(files$Gene_ID,file$V1),s]
		}
	}
	print(dim(files))
	fwrite(files,sprintf('star/%s_updated_gene_counts.txt',t),row.names=F,quote=F,sep='\t')
}


  #samples=samples[samples %in% unique(passed$sample_name)]
  #sub_files=files[,colnames(files) %in% samples]
  #names=colnames(sub_files)

  y=dim(sub_files)[2]
  totals=data.frame(colSums(sub_files[,2:y]))
  names(totals)='total'
  totals$sample=rownames(totals)


  pers10=data.frame(apply(sub_files[,2:y],MARGIN=2,function(x) sum(x>10)))
  names(pers10)='pers10'
  pers10$sample=rownames(pers10)

  pers100=data.frame(apply(sub_files[,2:y],MARGIN=2,function(x) sum(x>100)))
  names(pers100)='pers100'
  pers100$sample=rownames(pers100)

  perg10=data.frame(apply(sub_files[,2:y],MARGIN=1,function(x) sum(x>10)))
  names(perg10)='perg10'
  perg10$gene=all_files$Gene_ID[5:x]

  perg100=data.frame(apply(sub_files[,2:y],MARGIN=1,function(x) sum(x>100)))
  names(perg100)='perg100'
  perg100$gene=all_files$Gene_ID[5:x]

  sub_files$Gene_ID=all_files$Gene_ID[5:x]
  sub_files=sub_files[,c('Gene_ID',names)]
  fwrite(sub_files,sprintf('star/%s_gene_counts.txt',t),row.names=F,quote=F)

  line0=sprintf("Number of Samples: %.0f, Number of Genes: %.0f",y,x)
  line1=sprintf("Avg. Number of Total Reads per Sample: %.0f, SD %.0f",mean(totals$total),sd(totals$total))
  line2=sprintf("Avg. Number of Samples with >10 Reads per Gene: %.1f, SD %.1f",mean(perg10$perg10),sd(perg10$perg10))
  line3=sprintf("Number of Genes with >10 Reads in >10 Samples: %.0f",sum(perg10$perg10>10))
  line4=sprintf("Number of Genes with >10 Reads in >100 Samples: %.0f",sum(perg10$perg10>100))
  line5=sprintf("Number of Genes with >100 Reads in >100 Samples: %.0f",sum(perg100$perg100>100))
  #line6=sprintf("Avg. Mapping Rate (%%): %.1f, SD: %.1f",mean(mapping$uniquely_mapped_percent),sd(mapping$uniquely_mapped_percent))
  #line7=sprintf("Avg. Mapped Read Length: %.0f, SD: %.1f",mean(mapping$avg_mapped_read_length),sd(mapping$avg_mapped_read_length))
  timepoint_lines=c(timepoint_lines,c(sprintf("Time Point: %s",t),line0,line1,line2,line3,line4,line5,"\n"))
}

totals=data.frame(colSums(files[,2:y]))
names(totals)='total'
totals$sample=rownames(totals)

pers10=data.frame(apply(files[,2:y],MARGIN=2,function(x) sum(x>10)))
names(pers10)='pers10'
pers10$sample=rownames(pers10)

pers100=data.frame(apply(files[,2:y],MARGIN=2,function(x) sum(x>100)))
names(pers100)='pers100'
pers100$sample=rownames(pers100)

perg10=data.frame(apply(files[,2:y],MARGIN=1,function(x) sum(x>10)))
names(perg10)='perg10'
perg10$gene=all_files$Gene_ID[5:x]

perg100=data.frame(apply(files[,2:y],MARGIN=1,function(x) sum(x>100)))
names(perg100)='perg100'
perg100$gene=all_files$Gene_ID[5:x]

y=dim(files)[2]
line0=sprintf("Number of Samples: %.0f, Number of Genes: %.0f",y,x)
line1=sprintf("Avg. Number of Total Reads per Sample: %.0f, SD %.0f",mean(totals$total),sd(totals$total))
line2=sprintf("Avg. Number of Samples with >10 Reads per Gene: %.1f, SD %.1f",mean(perg10$perg10),sd(perg10$perg10))
line3=sprintf("Number of Genes with >10 Reads in >10 Samples: %.0f",sum(perg10$perg10>10))
line4=sprintf("Number of Genes with >10 Reads in >100 Samples: %.0f",sum(perg10$perg10>100))
line5=sprintf("Number of Genes with >100 Reads in >100 Samples: %.0f",sum(perg100$perg100>100))
line6=sprintf("Avg. Mapping Rate (%%): %.1f, SD: %.1f",mean(mapping$uniquely_mapped_percent),sd(mapping$uniquely_mapped_percent))
line7=sprintf("Avg. Mapped Read Length: %.0f, SD: %.1f",mean(mapping$avg_mapped_read_length),sd(mapping$avg_mapped_read_length))
line8=sprintf("Avg. Mapped Read Length: %.0f, SD: %.1f",mean(mapping$uniquely_mapped),sd(mapping$uniquely_mapped))
fileConn<-file("star/gene_count_summary.txt")
writeLines(c(line0,line1,line2,line3,line4,line5,line6,line7,'\n',timepoint_lines), fileConn)
close(fileConn)



p1=ggplot(pers10,aes(x=pers10)) + geom_histogram() + ggtitle('# of Genes per Sample with >10 Reads')
png('images/all_samples_10genecounts.png')
print(p1)
dev.off()

p2=ggplot(pers100,aes(x=pers100)) + geom_histogram() + ggtitle('# of Genes per Sample with >100 Reads')
png('images/all_samples_100genecounts.png')
print(p2)
dev.off()

p3=ggplot(totals,aes(x=total/1e6)) + geom_histogram() + ggtitle('Total # Reads Mapped to Genes') + xlab("Million Reads")
png('images/all_reads_mapped_to_genes.png')
print(p3)
dev.off()

mapping2=mapping[grep('pass2',mapping$Sample),]
p4=ggplot(mapping2,aes(y=log10(uniquely_mapped),x=uniquely_mapped_percent)) + geom_point() + ylab('log10(# Uniquely Mapped)') + xlab("% Uniquely Mapped")

png('images/all_mapping.png')
print(p4)
dev.off()


#  path=sprintf('star/%s/%s_pass2/ReadsReadsPerGene.out.tab')
#  file=fread(sprintf('%s/ReadsPerGene.out.tab',d),data.table=F)
#  file$group=d
#  files=rbind(files,file)
