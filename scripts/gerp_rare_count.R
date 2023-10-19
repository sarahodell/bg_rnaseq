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


allgerp=fread('datasets/all_founder_rare_alleles_GERP_scores.txt',data.table=F)
allgerp=allgerp[allgerp$CHR.x==chr,]
allgerp=allgerp[!is.na(allgerp$POS),]
tmpdf=fread(sprintf('eqtl/results/eQTL_%s_freq_chi_data.txt',time1),data.table=F)
tmpdf$gene_time_snp=paste0(tmpdf$Trait,'-',tmpdf$time,'-',tmpdf$X_ID)

#top5k=unique(tmpdf[tmpdf$rank<=5000,]$Trait)
#top5k=unique(tmpdf[tmpdf$rank>5000,]$Trait)

ranks=fread(sprintf('eqtl/results/%s_all_exp_ranks.txt',time1),data.table=F)
#ranks=fread(,sprintf('eqtl/results/%s_top5k_exp_ranks.txt',time1),data.table=F)
rownames(ranks)=ranks$V1
ranks=ranks[,-1]
upstream=fread('eqtl/data/upstream_gene_list.txt',data.table=F)

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
genetable=genetable[genetable$CHROM==chr,]
#upstream=upstream[upstream$CHROM==chr,]
#chr5k=top5k[top5k %in% upstream$Gene_ID]
chr5k=intersect(names(ranks),upstream$Gene_ID)
#ranks=ranks[,chr5k]

bimbam=fread(sprintf('datasets/chr%s_biogemma_rare_allele_probs.txt',chr),data.table=F)



gerpdel=allgerp[allgerp$bin!=1,]
# are there sites within a 10kb window around these genes
genetable$LEFT=genetable$START-5000
genetable$RIGHT=genetable$END+5000

gerpdel$POS_END=gerpdel$POS+1
env1=gerpdel
env1=as.data.table(env1)
env2=as.data.table(genetable)

setkey(env2,CHROM,LEFT,RIGHT)
comp=foverlaps(env1,env2,by.y=c('CHROM','LEFT','RIGHT'),by.x=c('CHR.x','POS','POS_END'),nomatch=NULL)
comp$marker=paste0('S',comp$CHR.x,'_',comp$POS)

chr5k=unique(comp$Gene_ID)
chr5k=intersect(chr5k,colnames(ranks))
#chr10 799 genes

#gerp4=allgerp[allgerp$bin==4,]
#gerp3=allgerp[allgerp$bin==3,]
#gerp2=allgerp[allgerp$bin==2,]
#gerp4=merge(gerp4,bimbam,by.x='POS',by.y='pos')
#gerp3=merge(gerp3,bimbam,by.x='POS',by.y='pos')
#gerp2=merge(gerp2,bimbam,by.x='POS',by.y='pos')

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
	
	sub1=comp[comp$Gene_ID==gene,]
	#rstart=upstream[upstream$Gene_ID==gene,]$UP_START
	#rend=upstream[upstream$Gene_ID==gene,]$UP_END
	#rstart=rend-5000
	#t4=gerp4[gerp4$pos>=rstart & gerp4$pos<=rend,]
	#t3=gerp3[gerp3$pos>=rstart & gerp3$pos<=rend,]
	#t2=gerp2[gerp2$pos>=rstart & gerp2$pos<=rend,]
	# How many alt alleles of each gerp bin does this ind have
	t1=bimbam[bimbam$marker %in% sub1$marker,c('marker','alt1','ref',id)]
	t1[,id]=t1[,id]/2
	sub1$id=t1[match(sub1$marker,t1$marker),id]
	#sub1$bin_f=factor(sub1$bin,levels=c(2,3,4))
	#t4=t4
	
	#sub1 %>% group_by(bin) %>% summarize(rcount=sum(id>=cutoff))
	
	rcount2=sum(sub1[sub1$bin==2,id]>=cutoff,na.rm=T)
	rcount3=sum(sub1[sub1$bin==3,id]>=cutoff,na.rm=T)
	rcount4=sum(sub1[sub1$bin==4,id]>=cutoff,na.rm=T)

	line=data.frame(Gene_ID=gene,ID=id,rank=r,gerp4=rcount4,gerp3=rcount3,gerp2=rcount2,stringsAsFactors=F)
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

fwrite(d,sprintf('eqtl/results/%s_%s_5kb_rare_counts_high_GERP.txt',time1,chr),row.names=F,quote=F,sep='\t')
