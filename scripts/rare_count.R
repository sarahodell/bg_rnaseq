#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])
time1=as.character(args[[2]])
cores=as.numeric(args[[3]])

library('data.table')
library('dplyr')
library('stringr')
library('ggplot2')
library('parallel')
library('MASS')

min_founder=readRDS(sprintf('datasets/hapmap/minor_allele_founders_c%s.rds',chr))

tmpdf=fread(sprintf('eqtl/results/eQTL_%s_freq_chi_data.txt',time1),data.table=F)
tmpdf$gene_time_snp=paste0(tmpdf$Trait,'-',tmpdf$time,'-',tmpdf$X_ID)

#top5k=unique(tmpdf[tmpdf$rank<=5000,]$Trait)
top5k=unique(tmpdf[tmpdf$rank>5000,]$Trait)

ranks=fread(sprintf('eqtl/results/%s_all_exp_ranks.txt',time1),data.table=F)
#ranks=fread(,sprintf('eqtl/results/%s_top5k_exp_ranks.txt',time1),data.table=F)
rownames(ranks)=ranks$V1

upstream=fread('eqtl/data/upstream_gene_list.txt',data.table=F)

#genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
upstream=upstream[upstream$CHROM==chr,]
chr5k=top5k[top5k %in% upstream$Gene_ID]

ranks=ranks[,chr5k]

bimbam=fread(sprintf('datasets/chr%s_biogemma_rare_allele_probs.txt',chr),data.table=F)

inds=rownames(ranks)
rmelt=as.data.frame(expand.grid(chr5k,inds),stringsAsFactors=F)
names(rmelt)=c("Gene_ID","ID")
cutoff=0.75

n_reps=1:nrow(rmelt)

get_rare_counts=function(rep){
	row1=rmelt[rep,]
	gene=as.character(row1$Gene_ID)
	id=as.character(row1$ID)
	r=ranks[id,gene]
	
	rstart=upstream[upstream$Gene_ID==gene,]$UP_START
	rend=upstream[upstream$Gene_ID==gene,]$UP_END
	#rstart=rend-5000
	t=bimbam[bimbam$pos>=rstart & bimbam$pos<=rend,]
	t=t[,c('marker','alt1','ref',id)]
	t[,id]=t[,id]/2
	rcount=sum(t[,id]>=cutoff)
	line=data.frame(Gene_ID=gene,ID=id,rank=r,rare_count=rcount,stringsAsFactors=F)
	return(line)
}

print(system.time({
results=mclapply(n_reps,get_rare_counts,mc.cores=cores)
}))
d=rbindlist(results)
#saveRDS(results,sprintf('eqtl/trans/permute/chr%s_%s_%s_%.0frep_min_pvalues.rds',chr,time,factor,reps))
#all_props=do.call(rbind,lapply(results,function(x) x))
d=as.data.frame(d,stringsAsFactors=F)

#fwrite(d,sprintf('eqtl/results/%s_%s_5kb_rare_counts_exp.txt',time1,chr),row.names=T,quote=F,sep='\t')

fwrite(d,sprintf('eqtl/results/%s_%s_5kb_rare_counts_lower_exp.txt',time1,chr),row.names=F,quote=F,sep='\t')
