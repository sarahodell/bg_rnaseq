#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=F)
pei=as.character(args[[1]])
cores=as.numeric(args[[2]])

library('ggplot2')
library('data.table')
library('dplyr')
library('stringr')
library('parallel')
library('MASS')

#qtl=fread('QTL/all_adjusted_QTL_all_methods.txt',data.table=F)
#qtl$pheno_env_ID=paste0(qtl$phenotype,'-',qtl$environment,'-',qtl$ID)
#qtl=qtl[qtl$method=="Founder_probs",]
#qtl2=fread('QTL/all_adjusted_QTL_peaks_trimmed.txt',data.table=F)
#qtl2$pheno_env_ID=paste0(qtl2$phenotype,'-',qtl2$environment,'-',qtl2$ID)

#qtl$prop_var=qtl2[match(qtl$pheno_env_ID,qtl2$pheno_env_ID),]$prop_var
#fwrite(qtl,'QTL/all_adjusted_QTL_SIs.txt',row.names=F,quote=F)
qtl=fread('QTL/all_adjusted_QTL_SIs.txt',data.table=F)


#subcomp=fread('QTT/QTL_trans_eQTL_EXPSTPAUL_overlap.txt',data.table=F)
subcomp2=fread('QTT/QTL_trans_eQTL_overlap.txt',data.table=F)
trans=fread('eqtl/results/all_trans_fdr_SIs_FIXED.txt',data.table=F)

#subcomp2$pheno_env_ID=paste0(subcomp2$phenotype,'-',subcomp2$environment,'-',subcomp2$ID)
#prop_var=c()
#for(i in 1:nrow(subcomp2)){
#	row=subcomp2[i,]
#	time=row$time
#	chr=row$CHR
#	pv=fread(sprintf('eqtl/trans/results/eQTL_%s_c%s_weights_prop_var_FIXED.txt',time,chr),data.table=F)
#	gene=row$gene
#	esnp=row$SNP
#	p=pv[pv$snp==esnp & pv$gene==gene,]$prop_var
#	prop_var=c(prop_var,p)
#}
#subcomp2$prop_var=prop_var

#subcomp2$qtl_prop_var=qtl[match(subcomp2$pheno_env_ID,qtl$pheno_env_ID),]$prop_var

#fwrite(subcomp2,'QTT/QTL_trans_eQTL_overlap.txt',row.names=F,quote=F,sep='\t')

genecount=subcomp2 %>% group_by(pheno_env_ID) %>% count()
pheno_env_ids=unique(subcomp2$pheno_env_ID)


times=c("WD_0712","WD_0718","WD_0720","WD_0727")

df=c()
for(time in times){
	print(time)
	for(c in 1:10){
		print(c)
		d=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%.0f_weights_results_filtered_FIXED.txt',time,c),data.table=F)
		#dlist=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.rds',time,c))
		#d=rbindlist(dlist)
		#d=as.data.frame(d)
		d$time=time
		#rm(dlist)
		pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
		d$CHR=c
		d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
		rm(pmap)
		#d=d[,c('Trait','X_ID','p_value_ML','time','CHR','BP')]
		df=rbind(df,d)
		rm(d)
	}
}
df=as.data.frame(df,stringsAsFactors=F)
df$gene_time_snp=paste0(df$Trait,'-',df$time,'-',df$X_ID)
df=df[df$gene_time_snp %in% trans$gene_time_snp,]

get_cor=function(row1,d2){
	row2=qtl[d2,]
	pheno=row2$phenotype
	env=row2$environment
	chr1=row1$CHR
	chr2=row2$CHR
	gene=row1$gene
	esnp=row1$SNP
	qsnp=row2$SNP
	id=row2$ID
	time1=row1$time
	gts=row1$gene_time_snp
	
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr2,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	
	#results=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.txt',time,chr1),data.table=F)
	results=df[df$gene_time_snp==gts,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	return(test$estimate)
}


random_ten=function(pei){
	d2=which(qtl$pheno_env_ID==pei)
	ndraws=genecount[genecount$pheno_env_ID==pei,]$n
	draw=sample(seq(1,nrow(trans)),ndraws)
	subeqtl=trans[draw,]
	n=min(10,ndraws)
	top10=subeqtl %>% slice_max(value, n = n)
	top10=as.data.frame(top10,stringsAsFactors=F)
	all_r=sapply(seq(1,n),function(x) get_cor(top10[x,],d2))
	highest_r=max(abs(all_r))
	return(highest_r)
}

#tl$pheno_env_ID=paste0(qtl$phenotype,'-',qtl$environment,'-',qtl$ID)

#peis=c()
#reps=c()
#for(pei in pheno_env_ids){
#	peis=c(peis,rep(pei,300))
#	reps=c(reps,seq(1,300))
#}
#repdf=data.frame(pheno_env_id=peis,rep=reps,stringsAsFactors=F)

bootstrap=funtion(rep){
	#row1=repdf[rep,]
	#pei=row1$pheno_env_id
	#n=row1$rep
	res=random_ten(pei))
	line=data.frame(pei=pei,rep=rep,max_r=res,stringsAsFactors=F)
	return(line)
}

#all_perms=c()
#for(pei in pheno_env_ids){
#	bootstrap=sapply(seq(1,300),function(x) random_ten(pei))
#	c
#	all_perms=rbind(all_perms,perms)
#}
n_reps=1:300

print(system.time({
results=mclapply(n_reps,bootstrap,mc.cores=cores)
}))
d=rbindlist(results)
d=as.data.frame(d,stringsAsFactors=F)
fwrite(d,sprintf('QTT/top10_%s_trans_eqtl_permutations.txt',pei),row.names=F,quote=F,sep='\t')
