#!/usr/bin/env Rscript

library('xlsx')
library('tidyverse')
library('data.table')
library('ggplot2')

# Need to find
# mapping rates

# # of gene counts of at least X length (X = 20,50,100,200,300)
#Per timepoint how many reads per gene, how many counts per sample (mean)
#How many samples have at least a million reads
#How many genes that have at least 10 reads in at least X samples
#Histogram of number of samples that have at least 10 (across genes)
#allfiles=list.files('star/batch_2','bg_*_data/multiqc_fastqc.txt',recursive=T)

# change to whatever file is 
meta=fread('metadata/BG_completed_sample_list_FIXED.txt',data.table=F)

#meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
dirs=c('18048-85-lane12-16','batch_2','18048-85-04-11','18048-85-01-03','batch_1','batch_2','redos')

#dirs=c('18048-85-lane12-16','batch_2','18048-85-04-11','18048-85-01-03','batch_1')
all_files=c()
#d=dirs[2] #batch_2 only
  #sample=unique(samples[samples$batch==d,]$sample_name)
for(d in dirs){
  pass2=Sys.glob(sprintf('star/update/%s/*-pass2',d))
  pass2_s1=sapply(seq(1,length(pass2)),function(x) strsplit(pass2[x],'/')[[1]][[3]])
  pass2_s2=sapply(seq(1,length(pass2)),function(x) strsplit(pass2_s1[x],'-pass2')[[1]][[1]])
  sample=unique(pass2_s2)
  for(s in sample){
    path=sprintf('star/update/%s/%s_pass2/ReadsPerGene.out.tab',d,s)
    file=fread(path,data.table=F)
    if(s==sample[1] & d==dirs[1]){
      #column 3 is the first strand gene counts, dont know if this is what we want
      file=file[,c(1,3)]
      names(file)=c('Gene_ID',s)
      files=file
    }else{
      file=file[,c(3),drop=F]
      names(file)=c(s)
      files=cbind(files,file)
    }
  }
  if(d==dirs[1]){
    all_files=files
  }else{
    all_files=cbind(all_files,files)
  }
}

all_files=fread('star/update/all_Reads_per_Gene.txt',data.table=F)

x=dim(all_files)[1]
y=dim(all_files)[2]
files=all_files[5:x,]

f=sapply(seq(2:y),function(x) is.numeric(files[1,x]))
drop=which(f==F)
files=files[,-drop]
y=dim(files)[2]


mapping=fread('multiqc_data/multiqc_star.txt',data.table=F)
passed=mapping[mapping$uniquely_mapped>=10e5,]
#passed=mapping[mapping$uniquely_mapped_percent>=50,]

passed$sample_name=sapply(seq(1,nrow(passed)),function(x) strsplit(passed$Sample[x],'_pass')[[1]][[1]])
meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
timep=unique(meta$experiment)
timepoint_lines=c()
for(t in timep){
  samples=unique(meta[meta$experiment==t,]$sample_name)
  samples=samples[samples %in% unique(passed$sample_name)]
  sub_files=files[,colnames(files) %in% samples]
  names=colnames(sub_files)

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
