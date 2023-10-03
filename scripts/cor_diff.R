#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
cores=as.numeric(args[[1]])

library('ggplot2')
library('data.table')
library('dplyr')
library('stringr')
library('parallel')
library('MASS')
library('DescTools')

qtl=fread('QTL/all_adjusted_QTL_SIs.txt',data.table=F)
distalcomp=fread('QTT/QTL_trans_eQTL_interval_overlap.txt',data.table=F)
localcomp=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)


trans=fread('eqtl/results/all_trans_fdr_SIs_ld_FIXED.txt',data.table=F)
eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)
eqtl$gene_time_snp=paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)

distalcount=distalcomp %>% group_by(pheno_env_ID) %>% count()
localcount=localcomp %>% group_by(pheno_env_ID) %>% count()
pheno_env_ids=unique(localcomp$pheno_env_ID)


times=c("WD_0712","WD_0718","WD_0720","WD_0727")

df=c()
for(time in times){
	print(time)
	for(c in 1:10){
		print(c)
		d=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%.0f_weights_results_filtered_FIXED.txt',time,c),data.table=F)
		d$time=time
		d$gene_time_snp=paste0(d$Trait,'-',d$time,'-',d$X_ID)
		d=d[d$gene_time_snp %in% trans$gene_time_snp,]
		#dlist=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.rds',time,c))
		#d=rbindlist(dlist)
		#d=as.data.frame(d)
		
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
#df$gene_time_snp=paste0(df$Trait,'-',df$time,'-',df$X_ID)
#df=df[df$gene_time_snp %in% trans$gene_time_snp,]

get_cor_d=function(row1,d2){
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

get_cor_l=function(row1,d2){
	row2=qtl[d2,]
	pheno=row2$phenotype
	env=row2$environment
	chr1=row1$CHR
	chr2=row2$CHR
	gene=row1$Trait
	esnp=row1$X_ID
	qsnp=row2$SNP
	id=row2$ID
	time=row1$time
	gene=row1$Trait
	
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr2,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr1),data.table=F)
	results=results[results$X_ID==esnp & results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	return(test$estimate)
}



random_tenl=function(pei){
	d2=which(qtl$pheno_env_ID==pei)
	ln=localcount[localcount$pheno_env_ID==pei,]$n
	dropgenes=unique(localcomp[localcomp$pheno_env_ID==pei,]$gene_time_SNP)
	subcis=eqtl[!(eqtl$gene_time_snp %in% dropgenes),]
	drawl=sample(seq(1,nrow(subcis)),ln)
	subeqtl=subcis[drawl,]
	#n=min(10,ndraws)
	#top10=subeqtl %>% slice_max(value, n = n)
	#top10=as.data.frame(top10,stringsAsFactors=F)
	all_r=sapply(seq(1,ln),function(x) get_cor_l(subeqtl[x,],d2))
	highest_r=max(abs(all_r))
	return(highest_r)
}

random_tend=function(pei){
	d2=which(qtl$pheno_env_ID==pei)
	dn=distalcount[distalcount$pheno_env_ID==pei,]$n
	dropgenes=unique(distalcomp[distalcomp$pheno_env_ID==pei,]$gene_time_snp)
	subtrans=trans[!(trans$gene_time_snp %in% dropgenes),]
	drawd=sample(seq(1,nrow(subtrans)),dn)
	subeqtl=subtrans[drawd,]
	#n=min(10,ndraws)
	#top10=subeqtl %>% slice_max(value, n = n)
	#top10=as.data.frame(top10,stringsAsFactors=F)
	all_r=sapply(seq(1,dn),function(x) get_cor_d(subeqtl[x,],d2))
	highest_r=max(abs(all_r))
	return(highest_r)
}


bootstrap=function(rep){
	row1=repdf[rep,]
	pei=row1$pheno_env_id
	n=row1$n
	resd=random_tend(pei)
	resl=random_tenl(pei)
	line=data.frame(pei=pei,rep=n,max_r_d=resd,max_r_l=resl,stringsAsFactors=F)
	return(line)
}

peis=c()
reps=c()
for(pei in pheno_env_ids){
	peis=c(peis,rep(pei,100))
	reps=c(reps,seq(1,100))
}
repdf=data.frame(pheno_env_id=peis,n=reps,stringsAsFactors=F)

n_reps=1:nrow(repdf)

print(system.time({
results=mclapply(n_reps,bootstrap,mc.cores=cores)
}))
all_perms=rbindlist(results)
all_perms=as.data.frame(all_perms,stringsAsFactors=F)
fwrite(all_perms,'QTT/top10_all_trans_eqtl_permutations.txt',row.names=F,quote=F,sep='\t')



#all_perms=c()
#for(pei in pheno_env_ids){
#	bootstrapl=sapply(seq(1,100),function(x) random_tenl(pei))
#	bootstrapd=sapply(seq(1,100),function(x) random_tend(pei))
#	perms=data.frame(pei=pei,rep=seq(1,100),max_r_l=bootstrapl,max_r_d=bootstrapd,stringsAsFactors=F)
#	all_perms=rbind(all_perms,perms)
#}

#all_perms=as.data.frame(all_perms,stringsAsFactors=F)
all_perms$local_z=FisherZ(abs(all_perms$max_r_l))
all_perms$distal_z=FisherZ(abs(all_perms$max_r_d))
all_perms$difference=all_perms$local_z-all_perms$distal_z
fwrite(all_perms,'QTT/all_z_permutations.txt',row.names=F,quote=F,sep='\t')

# true diff
#max_r_l= localcomp %>% group_by(pheno_env_ID) %>% slice(which.max(abs(r)))
#max_r_l=as.data.frame(max_r_l)
#max_r_l$z=FisherZ(abs(max_r_l$r))

#max_r_d= distalcomp %>% group_by(pheno_env_ID) %>% slice(which.max(abs(r)))
#max_r_d=as.data.frame(max_r_d)
#max_r_d$z=FisherZ(abs(max_r_d$r))

#truediff=data.frame(pei=pheno_env_ids,stringsAsFactors=F)
#truediff$local_z=max_r_l[match(truediff$pei,max_r_l$pheno_env_ID),]$z
#truediff$distal_z=max_r_d[match(truediff$pei,max_r_d$pheno_env_ID),]$z

#truediff$zdiff=truediff$local_z-truediff$distal_z
#fwrite(truediff,'QTT/local_distal_max_z.txt',row.names=F,quote=F,sep='\t')
