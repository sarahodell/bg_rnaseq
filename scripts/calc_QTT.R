#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')
library('qqman')

find_qtts=function(x,e,p){
	xexp=exp[,x,drop=F]
	
	pheno=phenotypes[phenotypes$Loc.Year.Treat==e,c('Genotype_code',p)]
	pheno$exp=xexp[match(pheno$Genotype_code,rownames(xexp)),x]
	test=cor.test(pheno[,p],pheno$exp,use="complete.obs")
	return(test)
}

# Test all genes as QTTs for phenotypes in ST_PAUL_2017_WD
all_qtts=data.frame(matrix(ncol=6,nrow=0))
names(all_qtts)=c('time','pheno','env','gene','r','pvalue')


times=c("WD_0712","WD_0718","WD_0720","WD_0727")
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:8]
e="EXP_STPAUL_2017_WD"

for(time in times){
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
	rownames(exp)=exp$V1
	exp=exp[,-1]
	metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
	metadata=metadata[metadata$experiment==time,]
	metadata=metadata[metadata$read==1,]
	#geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
	#kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
	#exp=exp[,kept_genes]

	#eqtl=fread('eqtl/results/all_eQTL_hits.txt',data.table=F)
	#eqtl=eqtl[eqtl$class!="factor",]
	#eqtl=eqtl[eqtl$time==time,]
	#exp=exp[,unique(eqtl$Trait)]
	genes=names(exp)

	for(p in phenos){
		correlations=sapply(genes,function(x) find_qtts(x,e,p))
		line=data.frame(time=time,pheno=p,env=e,gene=genes,r=unlist(correlations['estimate',]),pvalue=unlist(correlations['p.value',]),stringsAsFactors=F)
		all_qtts=rbind(all_qtts,line)
	}
}

all_qtts=all_qtts[order(all_qtts$pvalue),]
rownames(all_qtts)=seq(1,nrow(all_qtts))
all_qtts$padjust=p.adjust(all_qtts$pvalue,method='fdr')
sig=all_qtts[all_qtts$padjust<=0.05,]

png('QTT/allgene_QTT_qqplot.png')
print(qqman::qq(all_qtts$pvalue))
dev.off()

#genes=genes[1:5]



fwrite(all_qtts,'eqtl/results/all_cis_trans_QTTs.txt',row.names=F,quote=F,sep='\t')


ntests=nrow(all_qtts)
#388224 WD_0712 6 eQTL 288
#451536 WD_0718 3 eQTL 144
#431184 WD_0720 2 eQTL 96
#259248 WD_0727 3 eQTL 144

#ntests=3474240
threshold=0.05/ntests

#WD_0718
# Zm00001d042747 and flowering time -0.216
# Zm00001d042747 and tkw_15 0.17

#WD_0720
# Zm00001d025017 and flowering time -0.172, grain_yield_15, -0.15

# 5% FDR

    #no significant QTTs for WD_0712

#sig=qtl_genes[sum(qtl_genes[,c(27,29,31,33)]<=threshold) > 0,]
#which.main
#fwrite(qtl_genes,'eqtl/QTL_GeneExp_correlations.txt',row.names=F,quote=F,sep='\t')

#times=c("WD_0712","WD_0718",'WD_0720','WD_0727')
#all_qtts=c()
#for(t in times){
#	qtt=fread(sprintf('eqtl/results/%s_QTTs.txt',t),data.table=F)
#	qtt$time=t
#	all_qtts=rbind(all_qtts,qtt)
#}
sig=c()

comp=fread('QTT/cis_eQTL_STPAUL_QTL_overlaps.txt',data.table=F)
for(i in 1:nrow(comp)){
	test=comp[i,]
	gene=test$Trait
	time=test$time
	chr=test$CHR
	env="EXP_STPAUL_2017_WD"
	pheno=test$phenotype
	qsnp=test$SNP
	xid=test$X_ID
	null=fread(sprintf('QTT/permute/%s_%s_%s_random_permutation.txt',gene,pheno,env))
	cutoff=quantile(null$pvalue,0.05)
	pvalue=test$pvalue
	if(pvalue<=cutoff){
		sig=c(sig,i)
	}

}
compsig=comp[sig,]

