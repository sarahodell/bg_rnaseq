#!/usr/bin/env Rscript

library('data.table')

gtf=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46.gtf',data.table=F)
#names(gtf)=c("CHROM",)
keep=c("exon")
gtf=gtf[gtf$V3 %in% keep,]
bed=data.frame(chr=gtf[,1],start=gtf$V4,end=gtf$V5,stringsAsFactors=F)

gene1=sapply(seq(1,nrow(gtf)),function(x) strsplit(gtf$V9[x],';')[[1]][1])
gene2=sapply(seq(1,nrow(gtf)),function(x) strsplit(gene1[x],'\"')[[1]][2])
tran1=sapply(seq(1,nrow(gtf)),function(x) strsplit(gtf$V9[x],';')[[1]][2])
tran2=sapply(seq(1,nrow(gtf)),function(x) strsplit(tran1[x],'\"')[[1]][2])

exon1=sapply(seq(1,nrow(gtf)),function(x) strsplit(gtf$V9[x],';')[[1]][8])
exon2=sapply(seq(1,nrow(gtf)),function(x) strsplit(exon1[x],'\"')[[1]][2])
#exon1=sapply(seq(1,nrow(gtf)),function(x) strsplit(gtf$V9[x],';')[[1]][2])
#tran2=sapply(seq(1,nrow(gtf)),function(x) strsplit(tran1[x],'\"')[[1]][2])
bed$name=paste0(exon2,'.',gtf$V7)
bed$score='.'
bed$strand=gtf$V7
bed[,7:12]='.'
#bed$gene=gene2
#bed$transcript=tran2
bed=bed[bed$end-bed$start>0,]

fwrite(gtf,'ASE/data/Zea_mays.B73_RefGen_v4.46_exons.gtf',row.names=F,quote=F,col.names=F,sep='\t')
fwrite(bed,'ASE/data/Zea_mays.B73_RefGen_v4.46_exons.bed',row.names=F,col.names=F,quote=F,sep='\t')

# For later, it's probably worh making a bed12 file for the transcripts
# Will be useful
