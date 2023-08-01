#!/usr/bin/env Rscript
library('abind')
library('data.table')
library('ggplot2')
library('dplyr')

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}
eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end

qtl=fread('QTL/all_adjusted_QTL_peaks.txt',data.table=F)
qtl$block_start=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$start
qtl$block_end=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$end
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,block_start,block_end)
comparison1=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)


founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
plot_list=list()
count=1
#adj_chr=c(5,9)
cors=c()
pvals=c()
for(i in 1:nrow(comparison1)){
	row=comparison1[i,]
	pheno=row$phenotype
	env=row$environment
	chr=row$CHR
	gene=row$Trait
	snp=row$SNP
	id=row$ID
	time=row$time
	#exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	gene=row$Trait
	snp=row$SNP
	chr=row$CHR
	#if(chr %in% adj_chr){
	#	X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
	#}else{
	#	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	#
	#}
	#inds=rownames(X_list[[1]])
	#inter=intersect(exp$V1,inds)
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	#X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	#colnames(X) = founders
	#rownames(X) = inter
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	#w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	#results=results[[w]]
	results=results[results$X_ID==snp& results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	pvals=c(pvals,p)
	
	df=data.frame(founder=founders,pheno=effect_size,eqts=betas,stringsAsFactors=F)
	if(abs(r)>0.5){
		p1=ggplot(df,aes(x=eqts,y=pheno)) + geom_point(aes(color=founder)) +
		xlab("eQTL effect size (log2cpm)") + ylab("QTL effect size") +
		ggtitle(sprintf("%s %s and %s %s, r=%.2f",gene,time,id,env,r))
		plot_list[[count]]=p1
		count=count+1
	}
	
}

comparison1$r=cors
comparison1$pvalue=pvals

fwrite(comparison1,'QTT/QTL_cis_eQTL_overlap.txt',row.names=F,quote=F,sep='\t')

pdf('QTL/images/QTL_eQTL_effect_size_correlations.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()

plotlist2=list()
count=1
comparison1$env_ID=paste0(comparison1$environment,'-',comparison1$ID)
ids=unique(comparison1$env_ID)
for(i in ids){
	sub=comparison1[comparison1$env_ID==i,]
	loc=which.max(abs(sub$r))
	hval=sub[loc,]$i.value
	hr=sub[loc,]$r
	gene=sub[loc,]$gene_time
	h1=ggplot(sub,aes(x=i.value)) + geom_histogram() + geom_vline(xintercept=hval) +
	ggtitle(label=sprintf('eQTL qvalues within %s',i),subtitle=sprintf('Highest %s r=%.2f',gene,hr))
	
	plotlist2[[count]]=h1
	count=count+1
}

pdf('QTL/images/QTL_eQTL_power_correlation_hist.pdf')
for(i in 1:length(plotlist2)){
	print(plotlist2[[i]])
}
dev.off()

### Look at permutation thresholds

comp=fread('QTT/cis_eQTL_STPAUL_QTL_overlaps.txt',data.table=F)
sig=c()
for(i in 1:nrow(comp)){
	row=comp[i,]
	gene=row$Trait
	pheno=row$phenotype
	env=row$environment
	time=row$time
	pfile=sprintf('QTT/permute/%s_%s_%s_%s_random_permutation.txt',gene,time,pheno,env)
	if(file.exists(pfile)){
		perms=fread(pfile,data.table=F)
		threshold=quantile(perms$pvalue,0.05)
		if(row$pvalue<threshold){
			sig=rbind(sig,row)
		}
	}	
	
}
