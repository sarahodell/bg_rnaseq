#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('ggplot2')

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

alldf=c()
adj_chr=c(5,9)
for(chr in 1:10){
	#pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
	if(chr %in% adj_chr){
		X_list=readRDS(sprintf('phenotypes/bg%.0f_adjusted_genoprobs.rds',chr))
	}else{
		X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',chr))
	}
	snps=dimnames(X_list[[1]])[[2]]
	frep=sapply(snps,function(x) apply(do.call(cbind,lapply(X_list,function(y) y[,x])),MARGIN=2,function(z) round(sum(z>0.75))))
	fvar=apply(frep,MARGIN=2,function(x) var(x[x>3]))
	n_est=apply(frep,MARGIN=2,function(x) length(x[x>3]))
	df=data.frame(chr=chr,snp=snps,fvar=fvar,n_est=n_est,stringsAsFactors=F)
	alldf=rbind(alldf,df)
}

fwrite(alldf,'eqtl/data/founder_variance.txt',row.names=F,quote=F,sep='\t')

## For QTL

qtl=fread('QTL/all_adjusted_QTL_peaks_trimmed.txt',data.table=F)
qtlmerge=merge(qtl,alldf,by.x="SNP",by.y="snp")

p1=ggplot(qtlmerge,aes(x=prop_var,y=fvar)) + geom_point() + xlab("Proportion Variance Explained by QTL") +
ylab("Variance in Founder Representation at QTL")

png('images/QTL_prop_var_x_fvar.png')
print(p1)
dev.off()

p2=ggplot(qtlmerge,aes(x=value,y=fvar)) + geom_point() + xlab("QTL log10(pvalue)") +
ylab("Variance in Founder Representation at QTL")

png('images/QTL_log10_x_fvar.png')
print(p2)
dev.off()

# What about ciseqtl 
all_prop=c()
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	for(chr in 1:10){
		pvar=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_weights_prop_var_FIXED.txt',time,chr),data.table=F)
		all_prop=rbind(all_prop,pvar)
	}
}
fwrite(all_prop,'eqtl/results/cis_eqtl_prop_variance.txt',row.names=F,quote=F,sep='\t')


all_prop$gene_time_snp=paste0(all_prop$gene,'-',all_prop$time,'-',all_prop$snp)
eqtl$gene_time_snp=paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)

eqtlmerge=merge(eqtl,all_prop,by="gene_time_snp")
eqtlmerge=merge(eqtlmerge,alldf,by.x="X_ID",by.y="snp")

p3=ggplot(eqtlmerge,aes(x=prop_var,y=fvar)) + geom_point()+ xlab("Proportion Variance Explained by eQTL") +
ylab("Variance in Founder Representation at QTL")

png('images/eQTL_prop_var_x_fvar.png')
print(p3)
dev.off()

p4=ggplot(eqtlmerge,aes(x=value,y=fvar)) + geom_point()+ xlab("eQTL log10(qvalue)") +
ylab("Variance in Founder Representation at eQTL")

png('images/eQTL_log10_x_fvar.png')
print(p4)
dev.off()



### QTL prop_var and top r value in interval

comparison2=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
pheno_env_ids=unique(comparison2$pheno_env_ID)

n=10
top_r=c()
for(pei in pheno_env_ids){
	subdf=comparison2[comparison2$pheno_env_ID==pei,]
	subdf=subdf[order(subdf$i.value),]
	rownames(subdf)=seq(1,nrow(subdf))
	subdf$rank=seq(nrow(subdf),1)
	max_loc=which.max(abs(subdf$r))
	max_r=max(abs(subdf$r))
	top10=subdf %>% slice_max(i.value, n = n)
	top10=as.data.frame(top10)
	#top10_r=top10_r[order(top10_r$i.value),]
	#top10_r$rank=seq(1,10)
	top10_r_loc=unlist(which.max(abs(top10$r)))
	top_rank=top10[top10_r_loc,]$rank
	top10_r=abs(top10[top10_r_loc,]$r)
	newline=data.frame(pheno_env_ID=pei,max_r=max_r,top10_r=top10_r,max_rank=top_rank,stringsAsFactors=F)
	top_r=rbind(top_r,newline)
}
qtlmerge$pheno_env_ID=paste0(qtlmerge$phenotype,'-',qtlmerge$environment,'-',qtlmerge$ID)
top_r=merge(top_r,qtlmerge,by="pheno_env_ID")

p5=ggplot(top_r,aes(x=prop_var,y=max_r)) + geom_point() + xlab("Proportion Variance Explained by QTL") +
ylab("Max Correlation of overlapping eQTL (|r|)")

png('images/QTL_prop_var_by_r.png')
print(p5)
dev.off()

p6=ggplot(top_r,aes(x=value,y=max_r)) + geom_point() + xlab("QTL log10(pvalue)") +
ylab("Max Correlation of overlapping eQTL (|r|)")

png('images/QTL_log10_by_r.png')
print(p6)
dev.off()