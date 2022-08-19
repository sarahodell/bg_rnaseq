#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)

i=as.numeric(args[[1]])

library('tidyverse')
library('data.table')
#library('DESeq2')
library('limma')
library('xlsx')

#mat=fread(sprintf('star/%s_gene_counts.txt',time),data.table=F)
#rownames(mat)=mat$Gene_ID
#mat=mat[,-1]
#meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
#snames=colnames(mat)
#gnames=meta[match(snames,meta$sample_name),]$genotype
#mat=as.matrix(mat)
#coldata=data.frame(sample=snames,genotype=gnames)

#dds <- DESeqDataSetFromMatrix(countData = mat,
#                              colData = coldata,
#                              design = ~ 1)

#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

#dds2 = DESeq(dds)

#dds <- estimateSizeFactors(dds)
#normalized_counts <- counts(dds, normalized=TRUE)

#vsd <- vst(dds, blind=TRUE)


#g=500
#trans_cts=data.frame()
times=c("WD_0712","WD_0718","WD_0720","WD_0727")

#for(time in times){
#  print(time)
#  cts = fread(sprintf('star/vst_%s_gene_counts_results2.csv',time),data.table=F)
  #cts_melt=melt(cts,'V1')
#  rownames(cts)=cts$V1
#  cts=cts[,-1]
#  colnames(cts)=paste0(colnames(cts),'_',time)
#  if(ncol(trans_cts)==0){
#    trans_cts=cts[1:g,]
#  }
#  else{
#    i=intersect(rownames(trans_cts),rownames(cts))
#    i=i[1:g]
#    cts=cts[i,]
#    trans_cts=trans_cts[i,]
#    trans_cts=cbind(trans_cts,cts)
#  }
#}
#drop=which(rowSums(is.na(trans_cts))>0)
#trans_cts=trans_cts[-drop,]
#g=500

#trans_cts=trans_cts[1:g,]
#cnames=colnames(trans_cts)
#snames=sapply(seq(1,length(cnames)), function(x) strsplit(cnames[x],'_WD_')[[1]][1])
#tnames=sapply(seq(1,length(cnames)), function(x) strsplit(cnames[x],'_WD_')[[1]][2])
meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)

#gnames=meta[match(snames,meta$sample_name),]$genotype
#trans_cts=trans_cts[,-1]
#trans_cts=as.data.frame(lapply(trans_cts,as.numeric),stringsAsFactors=F)

#trans_cts=as.matrix(trans_cts)
#trans_cts=t(trans_cts)
#coldata=data.frame(cnames=cnames,sample=snames,genotype=gnames,time=tnames,stringsAsFactors=F)

#pca <- prcomp(trans_cts, center = TRUE)

#pcs=as.data.frame(pca$x,stringsAsFactors=F)
#pcs$time=coldata[match(rownames(trans_cts),coldata$cnames),]$time

#pca$x %>%
#  as.data.frame %>%
#  rownames_to_column("sample") %>%
#  separate(sample,c("_WD_")) %>%


#var_explained <- pca$sdev^2/sum(pca$sdev^2)
#png('top_500_PCA.png',width=800,height=800)
#print(ggplot(aes(x=PC1,y=PC2),data=pcs) + geom_point(aes(color=time),size=4) +
#  theme_bw(base_size=24) +
#  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#  theme(legend.position="top"))
#dev.off()


file='Listes_prèlèvements_transcriptome_GS17.xlsx'
ordering=read.xlsx(file, 2, header=TRUE, colClasses=NA)
date_times=c('2017-07-12','2017-07-18','2017-07-20','2017-07-27')

# Each timepoint by itself
#for(i in 1:length(times)){
  time=times[i]
  cts = fread(sprintf('star/vst_%s_gene_counts_results2.csv',time),data.table=F)
  #cts_melt=melt(cts,'V1')
  rownames(cts)=cts$V1
  cts=cts[,-1]
  #cts=t(cts)
  #pca <- prcomp(cts, center = TRUE)

  #pcs=as.data.frame(pca$x,stringsAsFactors=F)

  #var_explained <- pca$sdev^2/sum(pca$sdev^2)
  #png(sprintf('%s_top_500_PCA.png',time),width=800,height=800)
  #print(ggplot(aes(x=PC1,y=PC2),data=pcs) + geom_point(size=2) +
#    theme_bw(base_size=24) +
#    labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#         y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#    theme(legend.position="top"))
#  dev.off()

  #MDS
#character, "pairwise" to choose the top genes separately
#for each pairwise comparison between the samples or
#"common" to select the same genes for all comparisons.
#If gene.selection is "common", then the top genes are those
#with the largest standard deviations between samples.
#If gene.selection is "pairwise", then a different set of top genes
# is selected for each pair of samples. The pairwise feature selection
# may be appropriate for microarray data when different molecular pathways
# are relevant for distinguishing different pairs of samples.
  top=500
  #cts=cts[1:10,1:10]
  mds_pairwise=plotMDS(cts,top=top,gene.selection='pairwise')
  png(sprintf('%s_MDS_pairwise.png',time),width=800,height=800)
  print(plot(mds_pairwise))
  dev.off()


  mds_common=plotMDS(cts,top=top,gene.selection='common')
  png(sprintf('%s_MDS_common.png',time),width=800,height=800)
  print(plot(mds_common))
  dev.off()

  #sample_names=rownames(pcs)
  #order=ordering[ordering$DATE==date_times[i],]
  #genos=meta[match(sample_names,meta$sample_name),]$genotype
  #coldata=data.frame(sample_names=sample_names,genotype=genos)
  #ord=order[match(coldata$genotype,order$GENOTYPE),]$ORDER
  #coldata$order=ord
  #coldata=coldata[order(coldata$order),]
  #coldata$scale_order=seq(1,nrow(coldata))
  #print(time)
  #print("PC1")
  #print(cor(coldata$order,pcs$PC1))
  #print("PC2")
  #print(cor(coldata$order,pcs$PC2))

#}
