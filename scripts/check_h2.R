#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
chr=as.character(args[[2]])
#cores=as.numeric(args[[3]])


library('data.table')
library('xlsx')
library('lme4qtl')
library('lmtest')
library('edgeR')
library('dplyr')
library('lme4')
library('ggplot2')
#1) Check the difference in residuals from lme4qtl
# with and without K for 100 genes. Highly correlated?

#2) What is the heritability of 100 genes with lme4qtl? Is it close to 1?

#3) Check the PC and MDS plots of the top 500 genes


#4) Look at qqplots of GridLMM/lme4qtl results



K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)
genetable=fread('eqtl/data/Zea_mays.B73_v4_generanges.txt',data.table=F)
genetable=genetable[genetable$TXCHROM==chr,]
genes=unique(genetable$Gene_ID)
phenotypes=fread(sprintf('star/vst_%s_gene_counts_results2.csv',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
samples=colnames(phenotypes)[-1]
genos=metadata[match(samples,metadata$sample_name),]$genotype
testsnps=readRDS(sprintf('eqtl/gene_focal_snps_c%s.rds',chr))
founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
genes=intersect(genes,phenotypes$V1)
genes=genes[1:100]

correlation=data.frame(matrix(ncol=4,nrow=length(genes)))
names(correlation)=c('Gene','K_LRT_pvalue','noK_LRT_pvalue','K_heritability')
#heritability=c()
for(g in 1:length(genes)){
  gene=genes[g]
  snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
  if(length(snp)!=0){
    data=data.frame(sample=colnames(phenotypes)[-1],y=unlist(phenotypes[phenotypes[,1]==gene,][-1]),stringsAsFactors=F)
    data=data[!is.na(data$y),]
    data$ID=metadata[match(data$sample,metadata$sample_name),]$dh_genotype
    rownames(data)=seq(1,nrow(data))
    data=data[!is.na(data$ID),]
    X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
    inds=rownames(X_list[[1]])

    i=intersect(data$ID,inds)

    founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C
103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
    K=K[i,i]
    if(length(unique(data$ID))<nrow(data)){
      data=data%>%group_by(ID)%>%summarize(y=mean(y))
    }
    data=as.data.frame(data,stringsAsFactors=F)
    rownames(data)=data$ID
    data=data[i,]
    data$ID2=data$ID
    data=data[,c('ID','ID2','y')]


    if(length(snp)>1){
      #betas=unlist(unname(gwas[2,6:21]))
      X = do.call(cbind,lapply(X_list,function(x) x[,snp[2]]))
      colnames(X) = founders
      rownames(X) = dimnames(X_list[[1]])[[1]]
      snp=snp[2]
      X=X[i,]
    }else{
      #betas=unlist(unname(gwas[6:21]))
      X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
      colnames(X) = founders
      rownames(X) = dimnames(X_list[[1]])[[1]]
      X=X[i,]
    }
    frep2=apply(X,MARGIN=2,function(x) sum(x[x>0.8]))
    fkeep=founders[frep2>2]
    X=X[,fkeep]
    m4_K = relmatLmer(y ~ 1 + X + (1|ID),data=data,relmat = list(ID=K))
    m4_K_null = relmatLmer(y ~ 1 + (1|ID),data=data,relmat = list(ID=K))
    resK=lrtest(m4_K,m4_K_null)
    m4 = lm(y ~ 1 + X,data=data)
    m4_null=lm(y ~ 1, data=data)
    res=lrtest(m4,m4_null)
    #Is kinship changing significance of test?
    corK=resK[2,5]
    cor=res[2,5]
    #corr=cor(summary(m4)$coef[,'Estimate'],summary(m4_K)$coef[,'Estimate'])
    h2=VarProp(m4_K)[1,'prop']
#    grp        var1 var2       vcov     sdcor      prop
#1       ID (Intercept) <NA> 0.04472859 0.2114914 0.2974863
#2 Residual        <NA> <NA> 0.10562654 0.3250024 0.7025137
    correlation[g,]=c(gene,corK,cor,h2)
    #heritability=c(heritability,h2)
  }
}
correlation$K_LRT_pvalue=as.numeric(correlation$K_LRT_pvalue)
correlation$noK_LRT_pvalue=as.numeric(correlation$noK_LRT_pvalue)
correlation$K_heritability=as.numeric(correlation$K_heritability)
correlation=correlation[complete.cases(correlation),]
cor(correlation$K_LRT_pvalue,correlation$noK_LRT_pvalue)


p=ggplot(aes(x=K_heritability),data=correlation) + geom_histogram(bins=10)
png(sprintf('images/%s_chr%s_h2.png',time,chr))
print(p)
dev.off()

summary(correlation)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9730  0.9981  1.0000  0.9980  1.0000  1.0000

summary(heritability)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.0000  0.0000  0.0778  0.1385  0.4441

summary(1-heritability)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.5559  0.8615  1.0000  0.9222  1.0000  1.0000

all_res=as.data.frame(all_res,stringsAsFactors=F)
genenames=names(all_res)
all_res$ID=i
#names(all_res)=c(genes,'ID')
all_res=all_res[,c('ID',genenames)]


all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
#all_betas=as.data.frame(all_betas,stringsAsFactors=F)
#names(all_betas)=c('Gene_ID','SNP',paste0('beta.',seq(1,16)))

fwrite(all_res,sprintf('eqtl/results/eQTL_%s_c%s_residuals.txt',time,chr),row.names=F,quote=F,sep='\t')
fwrite(all_gwas,sprintf('eqtl/results/eQTL_%s_c%s_results.txt',time,chr),row.names=F,quote=F,sep='\t')

meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)

file='Listes_prèlèvements_transcriptome_GS17.xlsx'
ordering=read.xlsx(file, 2, header=TRUE, colClasses=NA)
date_times=c('2017-07-12','2017-07-18','2017-07-20','2017-07-27')
order=ordering[ordering$DATE==date_times[i],]

cts = fread(sprintf('star/vst_%s_gene_counts_results2.csv',time),data.table=F)
rownames(cts)=cts$V1
cts=cts[,-1]
cts=t(cts)

pca <- prcomp(cts, center = TRUE)
pcs=as.data.frame(pca$x,stringsAsFactors=F)
var_explained <- pca$sdev^2/sum(pca$sdev^2)

sample_names=rownames(pcs)
order=ordering[ordering$DATE==date_times[i],]
genos=meta[match(sample_names,meta$sample_name),]$genotype
coldata=data.frame(sample_names=sample_names,genotype=genos)
ord=order[match(coldata$genotype,order$GENOTYPE),]$ORDER
coldata$order=ord
coldata$pc1=pcs$PC1
coldata$pc2=pcs$PC2
coldata$pc3=pcs$PC3
coldata$pc4=pcs$PC4
ent=order[match(coldata$genotype,order$GENOTYPE),]$ENTRY
coldata$entry=ent
#tube=order[match(coldata$genotype,order$GENOTYPE),]$TUBE
#coldata$tube=tube
cle=order[match(coldata$genotype,order$GENOTYPE),]$CLE
coldata$cle=cle
x=order[match(coldata$genotype,order$GENOTYPE),]$X
coldata$x=x

tmp=coldata[,c('entry','pc1','pc2','pc3','pc4')]
tmpm=melt(tmp,'entry')

p1=ggplot(aes(x=entry,y=value),data=tmpm) + geom_point() + facet_grid(.~variable)+ labs(y=paste0("PC1: ",round(var_explained[1]*100,1),"%"))
png(sprintf('images/%s_pcs_by_entry.png',time),width=1000,height=800)
print(p1)
dev.off()

tmp=coldata[,c('order','pc1','pc2','pc3','pc4')]
tmpm=melt(tmp,'order')

p2=ggplot(aes(x=order,y=value),data=tmpm) + geom_point() + facet_grid(.~variable)+ labs(y=paste0("PC1: ",round(var_explained[1]*100,1),"%"))
png(sprintf('images/%s_pcs_by_order.png',time),width=1000,height=800)
print(p2)
dev.off()

tmp=coldata[,c('cle','pc1','pc2','pc3','pc4')]
tmpm=melt(tmp,'cle')

p3=ggplot(aes(x=cle,y=value),data=tmpm) + geom_point() + facet_grid(.~variable)+ labs(y=paste0("PC1: ",round(var_explained[1]*100,1),"%"))
png(sprintf('images/%s_pcs_by_cle.png',time),width=1000,height=800)
print(p3)
dev.off()

batch=meta[match(coldata$genotype,meta$genotype),]$batch
coldata$batch=batch
plate=meta[match(coldata$genotype,meta$genotype),]$plate
coldata$plate=plate
coldata$plate=as.factor(coldata$plate)

p4=ggplot(aes(x=pc1,y=pc2),data=coldata) + geom_point(aes(color=plate))
png(sprintf('images/%s_pcs_plate.png',timep),width=1000,height=800)
print(p4)
dev.off()

p4=ggplot(aes(x=pc1,y=pc2),data=coldata) + geom_point(aes(color=batch))
png(sprintf('images/%s_pcs_batch.png',timep),width=1000,height=800)
print(p4)
dev.off()
#p4=ggplot(aes(x=pc1,y=pc2))

#coldata=coldata[order(coldata$order),]
#coldata$scale_order=seq(1,nrow(coldata))
#print(time)
#print("PC1")
#print(cor(coldata$order,pcs$PC1))
#print("PC2")
#print(cor(coldata$order,pcs$PC2))
