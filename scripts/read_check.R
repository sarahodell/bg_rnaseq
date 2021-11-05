#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('ggplot2')

sizes=fread('raw_reads/batch_2_file_sizes.txt',data.table=F)
names(sizes)=c('file','size')
tsizes=fread('raw_reads/batch_2_trimmed_sizes.txt',data.table=F)

mapping=fread('salmon_quant/mapping_rates.txt',data.table=F)
names(mapping)=c('sample','mapping_perc')

tmp=sapply(seq(1,nrow(sizes)),function(x) strsplit(sizes[x,1],'.fastq.gz')[[1]][1])
sample=sapply(seq(1,nrow(sizes)),function(x) strsplit(tmp[x],'_R')[[1]][1])
read=sapply(seq(1,nrow(sizes)),function(x) strsplit(tmp[x],'_')[[1]][4])
well=sapply(seq(1,nrow(sizes)),function(x) strsplit(tmp[x],'_')[[1]][1])
well2=sapply(seq(1,nrow(sizes)),function(x) as.numeric(strsplit(well[x],'well')[[1]][2]))

lane=sapply(seq(1,nrow(sizes)),function(x) strsplit(tmp[x],'_')[[1]][3])
lane2=sapply(seq(1,nrow(sizes)),function(x) as.numeric(strsplit(lane[x],'L00')[[1]][2]))

num=sapply(seq(1,nrow(sizes)),function(x) strsplit(tmp[x],'_')[[1]][2])
num2=sapply(seq(1,nrow(sizes)),function(x) as.numeric(strsplit(num[x],'S')[[1]][2]))

sample2=sapply(seq(1,nrow(tsizes)),function(x) strsplit(tsizes[x,1],'_001.pe.qc.fastq.gz')[[1]][1])
sample3=sapply(seq(1,nrow(tsizes)),function(x) strsplit(sample2[x],'_R')[[1]][1])
read=sapply(seq(1,nrow(tsizes)),function(x) strsplit(sample2[x],'_')[[1]][4])

tsizes$sample=sample3
sizes$trimmed_size=tsizes[match(sizes$sample,tsizes$sample),]$V2

sizes$size_diff=sizes$size-sizes$trimmed_size

sizes$sample=sample
sizes$read=read
sizes$well=well2
sizes$lane=lane2
sizes$num=num2

sizes$mapping_perc=mapping[match(sizes$sample,mapping$sample),]$mapping_perc

avg_size=sizes %>% group_by(sample) %>% summarize(size=mean(size),size_diff=mean(size_diff),trimmed_size=mean(trimmed_size))
avg_size=as.data.frame(avg_size,stringsAsFactors=F)
avg_size$size_diff=as.numeric(avg_size$size_diff)
avg_size$size=as.numeric(avg_size$size)
avg_size$trimmed_size=as.numeric(avg_size$trimmed_size)
mapping$size=avg_size[match(avg_size$sample,mapping$sample),]$size
mapping$size_diff=avg_size[match(avg_size$sample,mapping$sample),]$size_diff
mapping$trimmed_size=avg_size[match(avg_size$sample,mapping$sample),]$trimmed_size
cor(mapping$size,mapping$mapping_perc)
#correlation between file size and mapping rate 0.1575612
p10_water=fread('raw_reads/p10_amnt_water.csv',data.table=F)
p11_water=fread('raw_reads/p11_amnt_water.csv',data.table=F)

sizes=sizes[order(sizes$num),]
rownames(sizes)=seq(1,nrow(sizes))

num=sapply(seq(1,nrow(mapping)),function(x) strsplit(mapping$sample[x],'_')[[1]][2])
num2=sapply(seq(1,nrow(mapping)),function(x) as.numeric(strsplit(num[x],'S')[[1]][2]))
mapping$num=num2
mapping=mapping[order(mapping$num),]
rownames(mapping)=seq(1,nrow(mapping))
mapping$num2=mapping$num-28

lane=sapply(seq(1,nrow(mapping)),function(x) strsplit(mapping$sample[x],'_')[[1]][3])
lane2=sapply(seq(1,nrow(mapping)),function(x) as.numeric(strsplit(lane[x],'L00')[[1]][2]))
mapping$lane=lane2

p10=mapping[mapping$lane==8,]
p10$conc=p10_water$conc
cor(p10$conc,p10$size)
# -0.06444644
p11=mapping[mapping$lane==7,]
p11$conc=p11_water$conc
cor(p11$conc,p11$size)
# 0.04422894

p1=ggplot(mapping,aes(x=size/1e6,y=mapping_perc)) + geom_point() + xlab("Size (Megabytes)") + ylab("Mapping Rate (%)")
png('raw_reads/batch2_size_by_mappingrate.png')
print(p1)
dev.off()

p2=ggplot(p10,aes(x=size/1e6,y=conc)) + geom_point() + xlab("Size (Megabytes)") + ylab("Library Concentration (ng/ml)")
png('raw_reads/batch2_plate10_size_by_conc.png')
print(p2)
dev.off()

p3=ggplot(p11,aes(x=size/1e6,y=conc)) + geom_point() + xlab("Size (Megabytes)") + ylab("Library Concentration (ng/ml)")
png('raw_reads/batch2_plate11_size_by_conc.png')
print(p3)
dev.off()


cor(p10$conc,p10$mapping_perc)
#0.108069
cor(p11$conc,p11$mapping_perc)
#0.1222891

p2=ggplot(p10,aes(x=mapping_perc ,y=conc)) + geom_point() + xlab("Size (Megabytes)") + ylab("Library Concentration (ng/ml)")
png('raw_reads/batch2_plate10_maping_by_conc.png')
print(p2)
dev.off()

p3=ggplot(p11,aes(x=size/1e6,y=conc)) + geom_point() + xlab("Size (Megabytes)") + ylab("Library Concentration (ng/ml)")
png('raw_reads/batch2_plate11_size_by_conc.png')
print(p3)
dev.off()

allp=rbind(p10,p11)
sub=allp[allp$conc>=1.8,]
sub=sub[sub$size<3e8,]

p4=ggplot(allp,aes(x=size_diff/1e6,y=conc)) + geom_point() + xlab("Post-Trimming File Size Difference") + ylab("Library Conc. (ng/ml)")
png('raw_reads/trimmed_by_conc.png')
print(p4)
dev.off()

p5=ggplot(allp,aes(x=trimmed_size/1e6,y=conc)) + geom_point() + xlab("Trimmed File Size (Mb)") + ylab("Library Conc. (ng/ml)")
png('raw_reads/trimmed_by_conc2.png')
print(p5)
dev.off()


mapping=fread('salmon_quant/batch_1_mapping_rates.txt',data.table=F)
names(mapping)=c('sample','mapping_perc')

lane=sapply(seq(1,nrow(mapping)),function(x) strsplit(mapping$sample[x],'_')[[1]][3])
lane2=sapply(seq(1,nrow(mapping)),function(x) as.numeric(strsplit(lane[x],'L00')[[1]][2]))
mapping$lane=lane2
p4=mapping[mapping$lane==1,]
