#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])

library('data.table')
library('ggplot2')
library('dplyr')

# How correlated with phenotype effect sizes are genes within the support interval?
#times=c("WD_0712","WD_0718",'WD_0720','WD_0727')
time="WD_0712"

#qtl_genes=fread('metadata/QTL_support_interval_genes.txt',data.table=F)
#qtl_genes$pheno_env=paste0(qtl_genes$Phenotype,'-',qtl_genes$Environment)


# Combine all phenotypes together into one file
#blups=fread('phenotypes/phenotypes_BLUPs.csv',data.table=F)
#phenotypes=fread('phenotypes/EXP_STPAUL_2017_WD_phenotypes.csv',data.table=F)
#phenotypes_full=fread('phenotypes/phenotypes_asi.csv',data.table=F)
#phenotypes_full$Genotype_code=gsub("-",".",phenotypes_full$Genotype_code)
#phenotypes_full=phenotypes_full[,c('Loc.Year.Treat','Genotype_code','female_flowering_d6','male_flowering_d6','total_plant_height','harvest_grain_moisture','grain_yield_15','tkw_15','asi')]

#phenotypes$Loc.Year.Treat="EXP_STPAUL_2017_WD"
#names(phenotypes)[1]="Genotype_code"
#phenotypes=phenotypes[,c('Loc.Year.Treat','Genotype_code','female_flowering_d6','male_flowering_d6','total_plant_height','harvest_grain_moisture','grain_yield_15','tkw_15','asi')]


#phenotypes_full=rbind(phenotypes_full,phenotypes)

#blups$Loc.Year.Treat="ALL"
#names(blups)=c("Genotype_code",'female_flowering_d6','male_flowering_d6','total_plant_height','harvest_grain_moisture','grain_yield_15','tkw_15','asi','Loc.Year.Treat')
#blups=blups[,c('Loc.Year.Treat','Genotype_code','female_flowering_d6','male_flowering_d6','total_plant_height','harvest_grain_moisture','grain_yield_15','tkw_15','asi')]

#phenotypes_full=rbind(phenotypes_full,blups)

#fwrite(phenotypes_full,'phenotypes/phenotypes_all.csv',row.names=F,quote=F,sep='\t')

phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
rownames(exp)=exp$V1
exp=exp[,-1]
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)

metadata=metadata[metadata$experiment==time,]


geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
exp=exp[,kept_genes]

eqtl=fread('eqtl/results/cis_trans_eQTL_hits.txt',data.table=F)
eqtl=eqtl[eqtl$time==time,]
exp=exp[,unique(eqtl$Trait)]

genes=names(exp)

#genes=genes[1:5]
find_qtts=function(x,e,p){
	xexp=exp[,x,drop=F]
	
	pheno=phenotypes[phenotypes$Loc.Year.Treat==e,c('Genotype_code',p)]
	pheno$exp=xexp[match(pheno$Genotype_code,rownames(xexp)),x]
	test=cor.test(pheno[,p],pheno$exp,use="complete.obs")
	return(test)
}

envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:8]

all_qtts=data.frame(matrix(ncol=5,nrow=0))
names(all_qtts)=c('pheno','env','gene','r','pvalue')

for(p in phenos){
	for(e in envs){
		correlations=sapply(genes,function(x) find_qtts(x,e,p))
		line=data.frame(pheno=p,env=e,gene=genes,r=unlist(correlations['estimate',]),pvalue=unlist(correlations['p.value',]),stringsAsFactors=F)
		all_qtts=rbind(all_qtts,line)
	}
}

fwrite(all_qtts,sprintf('eqtl/results/%s_QTTs.txt',time),row.names=F,quote=F,sep='\t')


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
all_qtts=all_qtts[order(all_qtts$pvalue),]
all_qtts$padjust=p.adjust(all_qtts$pvalue,method='fdr')
sig=all_qtts[all_qtts$padjust<=0.05,]
    #no significant QTTs for WD_0712

#sig=qtl_genes[sum(qtl_genes[,c(27,29,31,33)]<=threshold) > 0,]
#which.main
#fwrite(qtl_genes,'eqtl/QTL_GeneExp_correlations.txt',row.names=F,quote=F,sep='\t')

times=c("WD_0712","WD_0718",'WD_0720','WD_0727')
all_qtts=c()
for(t in times){
	qtt=fread(sprintf('eqtl/results/%s_QTTs.txt',t),data.table=F)
	qtt$time=t
	all_qtts=rbind(all_qtts,qtt)
}

