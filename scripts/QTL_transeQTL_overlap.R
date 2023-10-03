library('abind')
library('data.table')
library('ggplot2')
library('dplyr')

eqtl=fread('eqtl/results/all_trans_fdr_SIs_ld_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$gene,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$gene,'-',eqtl$time,'-',eqtl$SNP)




qtl=fread('QTL/all_adjusted_QTL_SIs.txt',data.table=F)
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
fwrite(comp,'QTT/QTL_trans_eQTL_interval_overlap.txt',row.names=F,quote=F,sep='\t')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
cors=c()
pvals=c()
for(i in 1:nrow(comp)){
	row=comp[i,]
	pheno=row$phenotype
	env=row$environment
	chr=row$CHR
	gene=row$gene
	qsnp=row$SNP
	time=row$time
	esnp=row$i.SNP
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	results=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.txt',time,chr),data.table=F)
	results=results[results$X_ID==esnp & results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	pvals=c(pvals,p)
}

comp$r=cors
comp$pvalue=pvals
fwrite(comp,'QTT/QTL_trans_eQTL_interval_overlap.txt',row.names=F,quote=F,sep='\t')

max_r=comp %>% group_by(pheno_env_ID) %>% slice(which.max(r))
max_r=as.data.frame(max_r,stringsAsFactors=F)
fwrite(max_r,'QTT/local_eQTL_candidates.txt',row.names=F,quote=F,sep='\t')

