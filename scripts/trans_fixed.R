#!/usr/bin/env Rscript

library('data.table')
library('dplyr')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
#plot_list=list()
#count=1
#adj_chr=c(5,9)
trans=fread('eqtl/results/all_trans_fdr_SIs_ld_FIXED.txt',data.table=F)

#trans=fread('eqtl/results/all_trans_fdr_SIs_FIXED.txt',data.table=F)
qtl=fread('QTL/all_adjusted_QTL_SIs.txt',data.table=F)
env1=qtl
env1=as.data.table(env1)
env2=as.data.table(trans)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
subcomp=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


subcomp$time_chr=paste0(subcomp$time,'-',subcomp$CHR)
subcomp=subcomp[subcomp$method=="Founder_probs",]
time_chr=unique(subcomp$time_chr)
subcomp$r=0
subcomp$pvalue=1
for(tc in time_chr){
	subcomp2=subcomp[subcomp$time_chr==tc,]
	time=unique(subcomp2$time)
	chr=unique(subcomp2$CHR)
	results=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.txt',time,chr),data.table=F)
	for(i in 1:nrow(subcomp2)){
		row=subcomp2[i,]
		gts=row$gene_time_snp
		pheno=row$phenotype
		env=row$environment
		gene=row$gene
		tsnp=row$SNP
		id=row$ID
		qsnp=row$i.SNP
		effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr,pheno,env),data.table=F)
		effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
		effect_size=unlist(effect_size[,c(6:21)])
		wn=which(!is.na(effect_size))[1]
		effect_size[-wn]=effect_size[-wn]+effect_size[wn]
		result=results[results$X_ID==tsnp& results$Trait==gene,]
		betas=unlist(result[,c(6,10:24)])
		wn=which(!is.na(betas))[1]
		betas[-wn]=betas[-wn]+betas[wn]
	
		test=cor.test(effect_size,betas,use="complete.obs")
		r=test$estimate
		p=test$p.value
		subcomp[subcomp$gene_time_snp==gts & subcomp$phenotype==pheno & subcomp$environment==env & subcomp$ID==id,]$r=r
		subcomp[subcomp$gene_time_snp==gts & subcomp$phenotype==pheno & subcomp$environment==env & subcomp$ID==id,]$pvalue=p
	
	}
}

fwrite(subcomp,'QTT/QTL_trans_eQTL_interval_overlap.txt',row.names=F,quote=F,sep='\t')

max_r=subcomp %>% group_by(pheno_env_ID) %>% slice(which.max(abs(r)))
max_r=as.data.frame(max_r,stringsAsFactors=F)
fwrite(max_r,'QTT/distal_eQTL_candidates.txt',row.names=F,quote=F,sep='\t')


#trans=fread('eqtl/results/all_trans_fdr_SIs_FIXED.txt',data.table=F)
#founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
#genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

#trans=merge(trans,genetable,by.x='gene',by.y="Gene_ID")
#trans$gene_time_snp=paste0(trans$gene,'-',trans$time,'-',trans$SNP)
 #121,730 distinct trans-eQTL
  
#eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
#eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
#eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
#eqtl=as.data.frame(eqtl2)

#all_founder_blocks=c()
#for(chr in 1:10){#
#  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
#  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
#}
#eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
#eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end


###### How many cis-eQTL also are trans-eQTL? ######
#env1=eqtl
#env1=as.data.table(env1)
#env2=as.data.table(trans)
#setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
#comp3=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

# How many cis-eQTL/trans-eQTL overlaps are for genes that are loaded on a factor together?










######## How many factor-eQTL overlap with trans-eQTL? #####

##### Whole factors #####
#factoreqtl=fread('eqtl/results/all_factor_fdr_SIs_FIXED.txt',data.table=F)
#factoreqtl$factor_time=paste0(factoreqtl$time,'-',factoreqtl$factor)
#env1=factoreqtl
#env1=as.data.table(env1)
#env2=as.data.table(trans)
#setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
#comp1=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
# 3,822 instances of overlap
# 3,820 trans-eQTL with all 4 factor eQTL
#neqtl=comp1 %>% group_by(factor_time) %>% count
# A tibble: 4 × 2
# Groups:   factor_time [4]
#  factor_time          n
#  <chr>            <int>
#1 WD_0712-Factor12   644
#2 WD_0712-Factor2   1011
#3 WD_0718-Factor4   1150
#4 WD_0720-Factor17  1017
#plot_list=list()
#plot_list2=list()
#count=1

#tdf=c()
#factors=unique(comp1$factor_time)
#for(f in factors){
#	subcomp=comp1[comp1$factor_time==f,]
#	timef=unique(subcomp$i.time)
#	subcomp=subcomp[subcomp$time==timef,]
#	factor=unique(subcomp$factor)
#	prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance_FIXED.txt',timef),data.table=F)
#	prop_var2=prop_var[prop_var[,factor]>0.1,c('V1',factor)]
#	tgenes=unique(subcomp$gene)
#	ninf=sum(tgenes %in% prop_var2$V1)
#	#inf=prop_var2[prop_var2$V1 %in% tgenes,]
#	prop_var2$inf=prop_var2$V1 %in% tgenes
#	names(prop_var2)=c("V1",'factor','inf')
#	full=prop_var2$factor
#	inf=prop_var2[prop_var2$inf==TRUE,]$factor
#	test=t.test(inf,full,alternative="greater")
#	prop=round(ninf/length(tgenes)*100,2)
#	line=data.frame(factor_time=f,ninf=ninf,noverlap=length(tgenes),total=nrow(prop_var2),prop=prop,pvalue=test$p.value)
#	tdf=rbind(tdf,line)
	#p1=ggplot(prop_var2,aes(x=factor,fill=inf)) + geom_histogram(alpha=0.25) + xlab("Proportion Variance Explained") + ylab("Frequency") + 
	#ggtitle(sprintf('%s',f))
	#p2=ggplot(prop_var2,aes(x=factor,fill=inf)) + geom_density(alpha=0.25) + xlab("Proportion Variance Explained") + ylab("Density") + 
	#ggtitle(sprintf('%s',f))
	#plot_list[[count]]=p1
	#plot_list2[[count]]=p2
	#count=count+1
#	print(f)
	#print(summary(prop_var2[,factor]))
	#print(summary(inf[,factor]))
#	print(length(tgenes))
#	print(ninf)
#}

#tdf=as.data.frame(tdf,stringsAsFactors=F)
#cutoff=0.05/16
#tdf[tdf$pvalue<=cutoff,]

# For WD_0720-Factor17, for the genes 
# with distal-eQTL overlapping the factor-eQTL, the facotr explained
# significantly more variation in the gene's expression than most genes
# loaded on the factor


#pdf('eqtl/images/factor_eqtl_trans_propvar.pdf')
#for(i in 1:length(plot_list)){
#	print(plot_list[[i]])
#}
#dev.off()
#
#pdf('eqtl/images/factor_eqtl_trans_density.pdf')
#for(i in 1:length(plot_list2)){
#	print(plot_list2[[i]])
#}
#dev.off()
####### How many of these overlapping genes are in the factors they overlap with?
#
###### Residual factors #####
#resids=fread('eqtl/results/all_residual_factor_fdr_SIs_FIXED.txt',data.table=F)
#resids$factor_time=paste0(resids$time,'-',resids$factor)
#
#env1=resids
#env1=as.data.table(env1)
#env2=as.data.table(trans)
#setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
#comp2=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
## 16,511 instances of overlap
## 15106 trans-eQTL with all 14 factor eQTL
#neqtl=comp2 %>% group_by(factor_time) %>% count
#
#plot_list=list()
#plot_list2=list()
#count=1
#
#tdf=c()
#factors=unique(comp2$factor_time)
#for(f in factors){
#	subcomp=comp2[comp2$factor_time==f,]
#	timef=unique(subcomp$i.time)
#	subcomp=subcomp[subcomp$time==timef,]
#	factor=unique(subcomp$factor)
#	prop_var=fread(sprintf('MegaLMM/MegaLMM_residuals_%s_prop_variance_FIXED.txt',timef),data.table=F)
#	prop_var2=prop_var[prop_var[,factor]>0.1,c('V1',factor)]
#	tgenes=unique(subcomp$gene)
#	prop_var2$inf=prop_var2$V1 %in% tgenes
#	names(prop_var2)=c("V1",'factor','inf')
#	full=prop_var2$factor
#	ninf=sum(tgenes %in% prop_var2$V1)
#	inf=prop_var2[prop_var2$inf==TRUE,]$factor
#	test=t.test(inf,full,alternative="greater")
#	prop=round(ninf/length(tgenes)*100,2)
#	line=data.frame(factor_time=f,ninf=ninf,noverlap=length(tgenes),total=nrow(prop_var2),prop=prop,pvalue=test$p.value)
#	tdf=rbind(tdf,line)
#	#p1=ggplot(prop_var2,aes(x=factor,fill=inf)) + geom_histogram(alpha=0.25) + xlab("Proportion Variance Explained") + ylab("Frequency") + 
#	#ggtitle(sprintf('Residusal %s',f))
#	#p2=ggplot(prop_var2,aes(x=factor,fill=inf)) + geom_density(alpha=0.25) + xlab("Proportion Variance Explained") + ylab("Density") + 
#	#ggtitle(sprintf('Residual %s',f))
#	#plot_list[[count]]=p1
#	#plot_list2[[count]]=p2
#	#count=count+1
#}
#tdf=as.data.frame(tdf,stringsAsFactors=F)
#cutoff=0.05/16
#tdf[tdf$pvalue<=cutoff,]
## For 5 of the residual factor trans-eQTL, 
## WD_0712-Factor5, WD_0712-Factor23, WD_0720-Factor19, WD_0712-Factor9, WD_0720-Factor17
##for the genes 
## with distal-eQTL overlapping the factor-eQTL, the facotr explained
## significantly more variation in the gene's expression than most genes
## loaded on the factor
#
#pdf('eqtl/images/resid_factor_eqtl_trans_propvar.pdf')
#for(i in 1:length(plot_list)){
#	print(plot_list[[i]])
#}
#dev.off()
#
#pdf('eqtl/images/resid_factor_eqtl_trans_density.pdf')
#for(i in 1:length(plot_list2)){
#	print(plot_list2[[i]])
#}
#dev.off()
########### Overlap with QTL ############
#
#qtl=fread('QTL/all_adjusted_QTL_all_methods.txt',data.table=F)
#env1=qtl
#env1=as.data.table(env1)
#env2=as.data.table(trans)
#setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
#comp=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
# 90,795 instances of overlap with 31,297 trans-eQTL
# Trans-eQTL overlap with all 33 QTL
# trans-eQTL for 13,995 
# What is the correlation in effect sizes?

# subcomp
#subcomp=comp[comp$environment=="EXP_STPAUL_2017_WD",]
# 6257 instances of overlap
# 5444 distal eQTL
# 4,251 genes
# 5 QTL

#subcomp %>% group_by(ID) %>% count
# A tibble: 5 × 2
# Groups:   ID [5]
#  ID          n
#  <chr>   <int>
#1 qDTA3_2  1444
#2 qDTA8     813
#3 qDTA9     954
#4 qDTS8    1823
#5 qTKW2    1223

#subcomp %>% group_by(ID) %>% summarize(n=length(unique(gene)))
# A tibble: 5 × 2
#  ID          n
#  <chr>   <int>
#1 qDTA3_2  1169
#2 qDTA8     728
#3 qDTA9     822
#4 qDTS8    1552
#5 qTKW2    1123


#fwrite(subcomp,'QTT/QTL_trans_eQTL_EXPSTPAUL_overlap.txt',row.names=F,quote=F,sep='\t')


#df=data.frame(founder=founders,pheno=effect_size,eqts=betas,stringsAsFactors=F)
#	df=df[complete.cases(df),]
#	p1=ggplot(df,aes(x=eqts,y=pheno)) + geom_point(aes(color=founder)) +
#	xlab("eQTL effect size (log2cpm)") + ylab("QTL effect size") +
#	ggtitle(sprintf("%s %s and %s %s, r=%.2f",gene,time,pheno,env,r))


#png(sprintf('QTL/images/%s_%s_%s_cor.png',gene,env,id))
#print(p1)
#dev.off()
