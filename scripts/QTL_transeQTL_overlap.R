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

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
comp=merge(comp,genetable,by.x='gene',by.y="Gene_ID")
comp=comp[comp$CHR==comp$CHROM,]
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
comp$r=0
comp$pvalue=1

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time1 in times){
	subcomp0=comp[comp$time==time1,]
	chroms=unique(subcomp0$CHR)
	for(chr in chroms){
		results=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.txt',time1,chr),data.table=F)
		subcomp1=subcomp0[subcomp0$CHR==chr,]
		for(i in 1:nrow(subcomp1)){
			row=subcomp1[i,]
			pheno=row$phenotype
			env=row$environment
			chr=row$CHR
			gene=row$gene
			qsnp=row$SNP
			time=row$time
			esnp=row$i.SNP
			gts=row$gene_time_snp
			effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr,pheno,env),data.table=F)
			effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
			effect_size=unlist(effect_size[,c(6:21)])
			wn=which(!is.na(effect_size))[1]
			effect_size[-wn]=effect_size[-wn]+effect_size[wn]
			#results=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.txt',time,chr),data.table=F)
			result1=results[results$X_ID==esnp & results$Trait==gene,]
			betas=unlist(result1[,c(6,10:24)])
			wn=which(!is.na(betas))[1]
			betas[-wn]=betas[-wn]+betas[wn]
			test=cor.test(effect_size,betas,use="complete.obs")
			r=test$estimate
			p=test$p.value
			comp[comp$gene_time_snp==gts,]$r=r
			comp[comp$gene_time_snp==gts,]$pvalue=p
		}
	}
}

fwrite(comp,'QTT/QTL_trans_eQTL_interval_overlap2.txt',row.names=F,quote=F,sep='\t')



#comp$r=cors
#comp$pvalue=pvals
#fwrite(comp,'QTT/QTL_trans_eQTL_interval_overlap2.txt',row.names=F,quote=F,sep='\t')

max_r=comp %>% group_by(pheno_env_ID) %>% slice(which.max(r))
max_r=as.data.frame(max_r,stringsAsFactors=F)
fwrite(max_r,'QTT/distal_eQTL_candidates.txt',row.names=F,quote=F,sep='\t')

