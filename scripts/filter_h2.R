#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
cores=as.numeric(args[[2]])

library('data.table')
library('xlsx')
library('lme4qtl')
library('lmtest')
library('edgeR')
library('dplyr')
library('lme4')
library('ggplot2')
library('parallel')
library('MASS')

exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
#sub=fread(sprintf('test_samples/%s_test_%ssamples.txt',time,n),data.table=F)
key=exp[,c('V1'),drop=F]
rownames(exp)=exp$V1
Y = exp[,-1]
#Y=Y[,sub$sub]

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
genes=unique(genetable$Gene_ID)

K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

inter=intersect(rownames(K),rownames(Y))
K=K[inter,inter]
Y=Y[inter,]

pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates_2.txt',time),data.table=F)

allweights=fread(sprintf('eqtl/normalized/%s_voom_weights_2.txt',time),data.table=F)
allweights=allweights[,c('V1',inter)]


#h2s=c()
#allgenes=c()
#chroms=c()
#for(chr in 1:10){
#  subgene=genetable[genetable$CHROM==chr,]
#  intergene=intersect(subgene$Gene_ID,colnames(Y))
#  for(g in intergene){
#      data=data.frame(ID=rownames(Y),y=Y[,g],stringsAsFactors=F)
#      data=data[!is.na(data$y),]
#      gweights=unlist(allweights[allweights$V1==g,inter])
#      founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
#      "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra",
#      "B96","DK63","F492","ND245","VA85")

#      rownames(data)=data$ID
#      data$ID2=data$ID
#      data=data[,c('ID','ID2','y')]

#     m4_K = relmatLmer(y ~ 1 + (1|ID),data=data,relmat = list(ID=K),weights=gweights)
#      h2=VarProp(m4_K)[1,'prop']
#      h2s=c(h2s,h2)
#      allgenes=c(allgenes,g)
#      chroms=c(chroms,chr)
#  }
#}

#herit=data.frame(chr=chroms,gene=allgenes,h2=h2s,stringsAsFactors=F)
#fwrite(herit,sprintf('eqtl/data/lme4qtl_weighted_%s_h2s.txt',time),row.names=F,quote=F,sep='\t')


#p1=ggplot(aes(x=h2),data=herit) + geom_histogram() + xlab('Gene h2')
#png(sprintf('images/gene_weighted_h2_%s_hist.png',time))
#print(p1)
#dev.off()

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
     	 "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra",
      	"B96","DK63","F492","ND245","VA85")

#subgene=genetable[genetable$CHROM==chr,]
intergene=intersect(genetable$Gene_ID,colnames(Y))

calc_h2=function(rep){
	g=intergene[rep]
	chr=genetable[genetable$Gene_ID==g,]$CHROM
	data=data.frame(ID=rownames(Y),y=Y[,g],stringsAsFactors=F)
	data=data[!is.na(data$y),]
	data$PC1=pcs[match(data$ID,pcs$V1),]$PC1
    data$PC2=pcs[match(data$ID,pcs$V1),]$PC2
    data$PC3=pcs[match(data$ID,pcs$V1),]$PC3
	m1=lm(y ~ PC1 + PC2 + PC3,data)
	resid=summary(m1)$residuals
	data$resid=resid
    rownames(data)=data$ID
    data$ID2=data$ID
    data=data[,c('ID','ID2','resid')]
    gweights=unlist(allweights[allweights$V1==g,inter])
	m4_K = relmatLmer(resid ~ 1 + (1|ID),data=data,relmat = list(ID=K),weights=gweights)
    h2=VarProp(m4_K)[1,'prop']
    line=data.frame(chr=chr,gene=g,h2=h2,stringsAsFactors=F)
	return(line)
}

#for(chr in 1:10){
#	subgene=genetable[genetable$CHROM==chr,]
#	intergene=intersect(subgene$Gene_ID,colnames(Y))
#	for(g in intergene){
#		data=data.frame(ID=rownames(Y),y=Y[,g],stringsAsFactors=F)
#		data=data[!is.na(data$y),]
#		gweights=unlist(allweights[allweights$V1==g,inter])
#		founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
#    	 "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra",
#      	"B96","DK63","F492","ND245","VA85")
#      	data$PC1=pc1[match(data$ID,pcs$V1),]$PC1
#      	data$PC2=pc1[match(data$ID,pcs$V1),]$PC2
#      	data$PC3=pc1[match(data$ID,pcs$V1),]$PC3
#		m1=lm(y ~ PC1 + PC2 + PC3,d1)
#		resid=summary(m1)$residuals
#		data$resid=resid
#     	rownames(data)=data$ID
#      	data$ID2=data$ID
#      	data=data[,c('ID','ID2','resid')]

#      	m4_K = relmatLmer(resid ~ 1 + (1|ID),data=data,relmat = list(ID=K),weights=gweights)
#      	h2=VarProp(m4_K)[1,'prop']
#      	h2s=c(h2s,h2)
#      	allgenes=c(allgenes,g)
#      	chroms=c(chroms,chr)
#      }
#}

#herit=data.frame(chr=chroms,gene=allgenes,h2=h2s,stringsAsFactors=F)
#fwrite(herit,sprintf('eqtl/data/lme4qtl_resid_%s_h2s.txt',time),row.names=F,quote=F,sep='\t')


#p1=ggplot(aes(x=h2),data=herit) + geom_histogram() + xlab('Gene h2')
#png(sprintf('images/gene_resid_h2_%s_hist.png',time))
#print(p1)
#dev.off()


n_reps=1:length(intergene)

print(system.time({
results=mclapply(n_reps,calc_h2,mc.cores=cores)
}))
d=rbindlist(results)
d=as.data.frame(d)
#saveRDS(results,sprintf('eqtl/data/founder_ld_c%.0f_c%.0f.rds',c1,c2))
fwrite(d,sprintf('eqtl/data/%s_SNP_resid_h2_2.txt',time),row.names=F,quote=F,sep='\t')

#intermed=herit[herit$h2 > 0,]
#rownames(intermed)=seq(1,nrow(intermed))
#draw=sample(seq(1,nrow(intermed)),1000)


#sub=intermed[draw,]
#sub=sub[,'gene',drop=F]
#fwrite(sub,sprintf('test_samples/%s_test_1000samples.txt',time),row.names=F,quote=F,sep='\t')
