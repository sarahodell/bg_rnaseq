#!/usr/bin/env Rscript
library('data.table')
library('dplyr')
library('ggplot2')

options(scipen=20)
thresh=0.01
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:9]

full_reps=as.data.frame(matrix(nrow=0,ncol=4))
names(full_reps)=c('phenotype','environment','method','threshold')
for(p in phenos){
	for(e in envs){
		chr_reps=c()
		for(c in 1:10){
  			df=fread(sprintf('QTL/permute/updated/chr%.0f_%s_x_%s_fp_updated_1000rep_max_pvalues.txt',c,p,e),data.table=F)
   			df=df[!is.na(df$pval),]
   			chr_reps=rbind(chr_reps,df)
		}
		minp = chr_reps %>% group_by(replicate) %>% summarize(pval=min(pval))
		minp=as.data.frame(minp)
		threshold=quantile(minp$pval,thresh,lower.tail=T)
		method="founder_probs"
		line=data.table(phenotype=p,environment=e,method=method,threshold=-log10(threshold))
		full_reps=rbind(full_reps,line)
	}
}


#print(threshold)
#print(-log10(threshold))

fwrite(full_reps,file=sprintf("QTL/updated_threshold_%.2f_table.txt",thresh),sep=',',row.names=F,quote=F)


#png(sprintf('GridLMM_founderprobs/permute/images/%s_x_%s_perm_1000_pval_dist.png',pheno,env))
#print(ggplot(minp,aes(x=pval)) + geom_histogram() + geom_vline(xintercept=threshold))
#dev.off()

#png(sprintf('GridLMM_founderprobs/permute/images/%s_x_%s_%.2f_perm_1000_log10pval_dist.png',pheno,env,thresh))
#print(ggplot(minp,aes(x=-log10(pval))) + geom_histogram() + geom_vline(xintercept=-log10(threshold)))
#dev.off()
thresh=0.1


full_reps=as.data.frame(matrix(nrow=0,ncol=4))
names(full_reps)=c('phenotype','environment','method','threshold')
for(p in phenos){
	for(e in envs){
		chr_reps=c()
		for(c in 1:10){
  			df=fread(sprintf('QTL/permute/adjusted/chr%.0f_%s_x_%s_fp_adjusted_1000rep_max_pvalues.txt',c,p,e),data.table=F)
   			df=df[!is.na(df$pval),]
   			chr_reps=rbind(chr_reps,df)
		}
		minp = chr_reps %>% group_by(replicate) %>% summarize(pval=min(pval))
		minp=as.data.frame(minp)
		threshold=quantile(minp$pval,thresh,lower.tail=T)
		method="founder_probs"
		line=data.table(phenotype=p,environment=e,method=method,threshold=-log10(threshold))
		full_reps=rbind(full_reps,line)
	}
}


#print(threshold)
#print(-log10(threshold))

fwrite(full_reps,file=sprintf("QTL/adjusted_threshold_%.2f_table.txt",thresh),sep=',',row.names=F,quote=F)
