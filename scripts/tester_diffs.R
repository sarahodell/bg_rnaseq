#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
c=as.numeric(args[[1]])

library('abind')
library('data.table')
library('dplyr')
#c=10

bed=fread('ASE/data/Zea_mays.B73_RefGen_v4.46_mRNA_regions.bed',data.table=F)
names(bed)=c("CHROM","START","END","TYPE","GENE_ID","TRANSCRIPT_ID")
bed$CHROM=as.numeric(bed$CHROM)
bed=bed[bed$CHROM==c,]
bed=bed[complete.cases(bed),]
rownames(bed)=seq(1,nrow(bed))
genos=fread(sprintf('ASE/Biogemma_WGS_transcripts_chr%.0f.csv',c),data.table=F)

inds=genos$ind
f_init=c("A632_usa","A654_inra","B73_inra","B96","C103_inra","CO255_inra","D105_inra",
"DK63","EP1_inra","F492","FV252_inra","FV2_inra","MBS847","ND245","OH43_inra",
"VA85","W117_inra")
genos$ind=f_init
rownames(genos)=genos$ind
founders=c("MBS847","B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
genos=genos[founders,]
genos=as.data.frame(t(genos),stringsAsFactors=F)
names(genos)=genos[1,]
genos=genos[-1,]
genos$SNP=rownames(genos)

chr1=sapply(seq(1,nrow(genos)),function(x) strsplit(genos$SNP[x],'_')[[1]][1])
chr=sapply(seq(1,nrow(genos)),function(x) strsplit(chr1[x],'r')[[1]][2])

pos=sapply(seq(1,nrow(genos)),function(x) strsplit(genos$SNP[x],'_')[[1]][2])
genos$CHROM=as.numeric(chr)
genos$POS=as.numeric(pos)


env1=genos
env1$END=env1$POS+1
env1=as.data.table(env1)
env2=as.data.table(bed)
#env2$end=env2$end-1
#setkey(env1,POS,END)
setkey(env2,START,END)
comparison=foverlaps(env1,env2,by.x=c('POS','END'),by.y=c('START','END'),nomatch=NULL)
comparison=as.data.frame(comparison,stringsAsFactors=F)
comparison=comparison[!is.na(comparison$MBS847),]

fwrite(comparison,sprintf('ASE/Biogemma_WGS_mRNA_overlap_chr%.0f.txt',c),row.names=F,quote=F,sep='\t')
#overlap


# for each gene, how many SNPs within the transcript are different from the tester for each founder?
#matrix with rows as genes and columns as founders, value is the number of SNPs that are different from MBS847
genes=unique(comparison$GENE_ID)
f16=founders[-1]
founder_counts=data.frame(matrix(0,nrow=length(genes),ncol=length(f16)),stringsAsFactors=F)
names(founder_counts)=f16
rownames(founder_counts)=genes
for(j in 1:length(genes)){
  gene=genes[j]
  t=comparison[comparison$GENE_ID==gene,]
  snps=unique(t$SNP)
  for(i in snps){
    st=t[t$SNP==i,]
    st=st[1,]
    diff_f=f16[which(st[,c(f16)]!=st$MBS847 & !is.na(st[,c(f16)]))]
    founder_counts[gene,diff_f]=founder_counts[gene,diff_f]+1
  }
}

fwrite(founder_counts,sprintf('ASE/Biogemma_WGS_founder_variant_counts_chr%.0f.csv',c),row.names=T,quote=F,sep='\t')

#total number of SNPs different between each founder and the tester
colSums(founder_counts)
# total number of genes where a SNP is present between the founder and the tester
apply(founder_counts,MARGIN=2,function(x) sum(x!=0))
