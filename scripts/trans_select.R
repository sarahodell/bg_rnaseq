#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
t_chr=as.character(args[[2]])
cores=as.numeric(args[[3]])

library('data.table')
library('dplyr')
library('parallel')
library('MASS')
library('stringr')
# Get cis fitted values
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)


trans=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_FIXED.rds',time,t_chr))
tgenes=unlist(lapply(trans,function(x) unique(x$Trait)))
genetable=genetable[genetable$Gene_ID %in% tgenes,]
genetable=genetable[genetable$CHROM==t_chr,]
cgenes=unique(genetable$Gene_ID)

samples=fread(sprintf('eqtl/data/%s_samples_FIXED.txt',time),data.table=F)
inter=samples$id

transbv=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_bvs_FIXED.rds',time,t_chr))

#c_chr=genetable[genetable$Gene_ID==gene,]$CHROM
cisbv=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_bv_FIXED.txt',time,t_chr),data.table=F)
rownames(cisbv)=cisbv$V1
cisbv=cisbv[,-1]
cisbv=cisbv[inter,]

#get_bv=function(row){
#	snp=row$X_ID
#	gene=row$Trait
#	adj_chr=c("5","9")
#	if(t_chr %in% adj_chr){
#		X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',t_chr))
#
#	}else{
#		X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',t_chr))
#	}
#	X=do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
#	colnames(X) = founders
#	rownames(X) = inter
#	X_list_ordered=lapply(X_list,function(x) x[inter,snp,drop=F])
#	frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
#	fkeep=founders[frep2>3]
#	X_list_ordered = X_list_ordered[c(fkeep)]
#	#X=X[,fkeep]
#	betas=unlist(row[,fkeep])
#	full_betas=betas[founders]
#	names(full_betas)=founders
#	full_betas[!(names(full_betas) %in% fkeep)]=0
#	wn=which(!is.na(betas))[1]
#	intercept=names(betas)[wn]
#	full_betas[names(full_betas)!=intercept]=full_betas[intercept]+full_betas[names(full_betas)!=intercept]
#	bv=X %*% full_betas
#	bv=as.data.frame(bv,stringsAsFactors=F)
#	names(bv)=snp
#	return(bv)
#}

get_r=function(rep){
	results=trans[[rep]]
	gene=unique(results$Trait)
	bvpos=which(cgenes==gene)
	tbv=transbv[[bvpos]]
	#tbvs=sapply(seq(1,nrow(results)),function(x) get_bv(results[x,]))
	#tbv=do.call(cbind,tbvs)
	#rownames(tbv)=inter
	rs=apply(tbv,MARGIN=2,function(x) cor(cisbv[,gene],x))
	df=data.frame(gene=gene,snp=names(rs),r=rs,stringsAsFactors=F)
	return(df)
}

get_kept=function(rep){
	results=trans[[rep]]
	gene=unique(results$Trait)
	if(gene %in% cgenes){
		df=get_r(rep)
		cutoff=0.5
		kept_snps=df[abs(df$r)<cutoff,]$snp
		new_results=results[results$X_ID %in% kept_snps,]
	}else{
		new_results=results
	}
	return(new_results)
}

#n_reps=match(cgenes,tgenes)
n_reps=1:length(tgenes)
#n_reps=1:5

print(system.time({
full_results=mclapply(n_reps,get_kept,mc.cores=cores)
}))
d=rbindlist(full_results)
d=as.data.frame(d)

fwrite(d,sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.txt',time,t_chr),row.names=F,quote=F,sep='\t')

#for(g in 1:length(n_reps)){
#	b[[n_reps[g]]]=results[[g]]
#}
#saveRDS(full_results,sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.rds',time,t_chr))


