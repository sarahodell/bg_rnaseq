#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')

### Look for outliers
factoreqtl=fread('eqtl/results/all_residual_factor_fdr_SIs_FIXED.txt',data.table=F)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

factoreqtl$factor_time=paste0(factoreqtl$time,'-',factoreqtl$factor)
factoreqtls=unique(factoreqtl$factor_time)

########## Are F values corrleated with overlapping phenotype QTL?
time1="WD_0720"
factor="Factor19"
pheno="grain_yield_15"
env="NERAC_2016_WD"
#-0.05
time1="WD_0727"
factor="Factor12"
pheno="total_plant_height"
env="ALL"
#0.134
Fvalues=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_all_F_means_FIXED.txt',time1),data.table=F)
rownames(Fvalues)=Fvalues$V1
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

subpheno=phenotypes[phenotypes$Loc.Year.Treat==env,]
rownames(subpheno)=subpheno$Genotype_code
inter=intersect(subpheno$Genotype_code,Fvalues$V1)
Fvalues=Fvalues[inter,]
subpheno=subpheno[inter,]
cor(subpheno[,pheno],Fvalues[,factor])

# Are the founder effect sizes correlated? I can only do this for whole gene counts
chr2=3
pheno="male_flowering_d6"
env="EXP_STPAUL_2017_WD"
qsnp="PZE-103093413"
#-0.2196711

pheno="harvest_grain_moisture"
env="BLOIS_2017_OPT"
qsnp="PZE-103093413"
#-0.452876

pheno="harvest_grain_moisture"
env="ALL"
qsnp="PZE-103084819"

time1="WD_0720"
factor1="Factor17"
esnp="AX-90826809"
#-0.6890176

effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr2,pheno,env),data.table=F)
effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
effect_size=unlist(effect_size[,c(6:21)])
wn=which(!is.na(effect_size))[1]
effect_size[-wn]=effect_size[-wn]+effect_size[wn]

results=fread(sprintf('eqtl/trans/results/%s_c%s_%s_trans_results_FIXED.txt',time1,chr2,factor1),data.table=F)
result=results[results$X_ID==esnp,]
betas=unlist(result[,founders])
wn=which(!is.na(betas))[1]
betas[-wn]=betas[-wn]+betas[wn]

cor(betas,effect_size,use="complete.obs")

# 

axisdf=fread('eqtl/data/chromosome_axis.txt',data.table=F)
cumtot=2105119857
############### 
#Plot out position by prop_var
plot_list=list()

for(f in factoreqtls){
	sub=factoreqtl[factoreqtl$factor_time==f,]
	time=unique(sub$time)
	factor=unique(sub$factor)
	
	prop_var=fread(sprintf('MegaLMM/MegaLMM_residuals_%s_prop_variance_FIXED.txt',time),data.table=F)
	df=prop_var[prop_var[,factor]>0.1,c('V1',factor)]
	names(df)=c('Gene_ID','prop_var')
	df=merge(df,genetable,by="Gene_ID")
	df$midgene=round(df$START + (df$END-df$START)/2)
	df$tot=axisdf[match(df$CHROM,axisdf$gene_chr),]$tot
	df$midgene_cum=df$midgene + df$tot 
	
	sub$tot=axisdf[match(sub$CHR,axisdf$gene_chr),]$tot
	sub$bp_cum=sub$BP + sub$tot 
 
	ngenes=nrow(df)
	print(sprintf('%.0f Genes loaded on %s %s',ngenes,time, factor))
	
	p=ggplot(df,aes(x=midgene_cum,y=prop_var)) +
	scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
	ylim(0.1,1) +
    #scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
    #geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
	geom_point(aes(color=prop_var)) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
	xlab("Position") + ylab("Proportion Variance") +
	ggtitle(sprintf("Gene Loadings on %s %s, n=%.0f",time,factor,ngenes)) +
	theme_classic() +
    theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_line(colour="black"))
 	
 	for(i in 1:nrow(sub)){
 		row=sub[i,]
 		BP_cum=row$bp_cum
 		p=p+geom_vline(xintercept=BP_cum,color="coral")
 	}
 
	
	pdf(sprintf('eqtl/images/%s_%s_residuals_plot.pdf',time,factor))
	print(p)
	dev.off()
}

#factoreqtl=fread('eqtl/results/all_residual_factor_fdr_peaks_FIXED.txt',data.table=F)
factoreqtl=fread('eqtl/results/all_factor_fdr_peaks_FIXED.txt',data.table=F)
factoreqtl$factor_time=paste0(factoreqtl$time,'-',factoreqtl$factor)
factoreqtls=unique(factoreqtl$factor_time)
for(f in factoreqtls){
	sub=factoreqtl[factoreqtl$factor_time==f,]
	time=unique(sub$time)
	factor=unique(sub$factor)
	
	prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance_FIXED.txt',time),data.table=F)
	df=prop_var[prop_var[,factor]>0.1,c('V1',factor)]
	names(df)=c('Gene_ID','prop_var')
	df=merge(df,genetable,by="Gene_ID")
	df$midgene=round(df$START + (df$END-df$START)/2)
	df$tot=axisdf[match(df$CHROM,axisdf$gene_chr),]$tot
	df$midgene_cum=df$midgene + df$tot 
	
	sub$tot=axisdf[match(sub$CHR,axisdf$gene_chr),]$tot
	sub$bp_cum=sub$BP + sub$tot 
 
	ngenes=nrow(df)
	print(sprintf('%.0f Genes loaded on %s %s',ngenes,time, factor))
	
	p=ggplot(df,aes(x=midgene_cum,y=prop_var)) +
	scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
	ylim(0.1,1) +
    #scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
    #geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
	geom_point(aes(color=prop_var)) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
	xlab("Position") + ylab("Proportion Variance") +
	ggtitle(sprintf("Gene Loadings on %s %s, n=%.0f",time,factor,ngenes)) +
	theme_classic() +
    theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_line(colour="black"))
 	
 	for(i in 1:nrow(sub)){
 		row=sub[i,]
 		BP_cum=row$bp_cum
 		p=p+geom_vline(xintercept=BP_cum,color="coral")
 	}
 
	
	pdf(sprintf('eqtl/images/%s_%s_plot.pdf',time,factor))
	print(p)
	dev.off()
}

factoreqtl=fread('eqtl/results/all_factor_fdr_SIs_FIXED.txt',data.table=F)
resid=fread('eqtl/results/all_residual_factor_fdr_SIs_FIXED.txt',data.table=F)
env1=factoreqtl
env1=as.data.table(env1)
env2=as.data.table(resid)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


qtl=fread('QTL/all_adjusted_QTL_all_methods.txt',data.table=F)
env1=qtl
env1=as.data.table(env1)
env2=as.data.table(factoreqtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp2=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

env1=qtl
env1=as.data.table(env1)
env2=as.data.table(resid)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp4=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

# Is T20 Factor 17 enriched for flowering time genes
time="WD_0727"
prop_var=fread(sprintf('MegaLMM/MegaLMM_residuals_%s_prop_variance_FIXED.txt',time),data.table=F)
ft_genes=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
names(ft_genes)=c('CHR','START','END','GENE_ID')

allft=c()

factors=names(prop_var)[-1]
for(factor in factors){
	fgenes=prop_var[prop_var[,factor]>=0.1,]$V1
	tn=sum(fgenes %in% ft_genes$GENE_ID)
	ft_geneoverlap=data.frame(factor=factor,ngenes=length(fgenes),inft=tn,stringsAsFactors=F)
	allft=rbind(allft,ft_geneoverlap)
}

allgenes=prop_var$V1

null=c()
for(i in 1:nrow(allft)){
	row=allft[i,]
	ngenes=row$ngenes
	tgenes=row$inft
	alltn=c()
	for(j in 1:1000){
		draw=sample(allgenes,ngenes)
		tn=sum(draw %in% ft_genes$GENE_ID)
		alltn=c(alltn,tn)
	}
	adjust=68
	q5=quantile(alltn,1-(0.05/adjust))
	null=c(null,q5)

}
allft$null=null

allft[allft$inft>allft$null,]

#WD_0712 Factor3 WD_0727 Factor4 0.47
#WD_0718 Factor2 WD_0720 Factor4 0.7565007
#WD_0718 Factor2 WD_0727 Factor6 0.7461929
#WD_0718 Factor4 WD_0720 Factor1 0.7645210
#WD_0718 Factor4 WD_0727 Factor4 0.4809296
#WD_0720 Factor4 WD_0727 Factor6 0.6893948

# WD_0712 Factor 3 
# WD_0718 Factor 2, Factor 4, Factor 8, Factor 10
# WD_0720 Factor 4
# WD_0727 Factor 4, Factor 6

# WD_0718 Factor 4 has a factor-eQTL on chr2

time="WD_0712"
prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance_FIXED.txt',time),data.table=F)
ft_genes=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
names(ft_genes)=c('CHR','START','END','GENE_ID')

allft=c()

factors=names(prop_var)[-1]
for(factor in factors){
	fgenes=prop_var[prop_var[,factor]>=0.1,]$V1
	tn=sum(fgenes %in% ft_genes$GENE_ID)
	ft_geneoverlap=data.frame(factor=factor,ngenes=length(fgenes),inft=tn,stringsAsFactors=F)
	allft=rbind(allft,ft_geneoverlap)
}

allgenes=prop_var$V1

null=c()
for(i in 1:nrow(allft)){
	row=allft[i,]
	ngenes=row$ngenes
	tgenes=row$inft
	alltn=c()
	for(j in 1:1000){
		draw=sample(allgenes,ngenes)
		tn=sum(draw %in% ft_genes$GENE_ID)
		alltn=c(alltn,tn)
	}
	adjust=68
	q5=quantile(alltn,1-(0.05/adjust))
	null=c(null,q5)

}
allft$null=null

allft[allft$inft>allft$null,]
#WD_0712 Factor 3, Factor 8
#WD_0718 Factor 2
#WD_0720 Factor 6, Factor 8
#WD_0727 Factor5 Factor 7

########
factoreqtl=fread('eqtl/results/all_factor_fdr_peaks.txt',data.table=F)
factordf=fread('eqtl/results/all_factor_trans_eqtl_fdr_genes.txt',data.table=F)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
axisdf=fread('eqtl/data/chromosome_axis.txt',data.table=F)
cumtot=2105119857

# For each factor
#		- How many genes are loaded on them
#		- Where are the genes located in the genome
#		- What prop variance is explained by the factor for each gene
#		- How much variation in the factor is explained by the factor eQTL
#		- Are there trans-eQTL or cis-eQTL that overlap with the factor-eQTL variant? 
#		- Are there trans-eQTL or cis-eQTL genes that load on the factor?

factoreqtl$time_factor=paste0(factoreqtl$time,'-',factoreqtl$Trait)
tf=unique(factoreqtl$time_factor)
for(f in tf){
	sub=factoreqtl[factoreqtl$time_factor==f,]
	time=unique(sub$time)
	factor=unique(sub$factor)
	df=factordf[factordf$time==time & factordf$Trait==factor,]
	df$midgene=round(df$gene_start + (df$gene_end-df$gene_start)/2)
	df$tot=axisdf[match(df$gene_chr,axisdf$gene_chr),]$tot
	df$tot.x=axisdf[match(df$CHR,axisdf$gene_chr),]$tot
	BP_cum=unique(df$BP) + unique(df$tot.x)
	df$midgene_cum=df$midgene + df$tot 
	prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance.txt',time),data.table=F)
	print(time)
	print(factor)
	fgenes=prop_var[prop_var[,factor]>0.1,]$V1
	ngenes=length(fgenes)
	print(sprintf('%.0f Genes loaded on %s',ngenes,factor))
	
	p=ggplot(df,aes(x=midgene_cum,y=prop_var)) +
	geom_vline(xintercept=BP_cum,color="coral") +
	scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
	ylim(0.1,1) +
    #scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
    #geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
	geom_point(aes(color=prop_var)) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
	xlab("Position") + ylab("Proportion Variance") +
	ggtitle(sprintf("Gene Loadings on %s %s, n=%.0f",time,factor,ngenes)) +
	theme_classic() +
    theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_line(colour="black"))
 
	
	pdf(sprintf('eqtl/images/%s_%s_plot.pdf',time,factor))
	print(p)
	dev.off()
	# Plot x-axis location
	# y-axis prop_var
	# vhline of location of factor eQTL

}




##############

time="WD_0712"
factor="Factor2"
chr="2"
snp="AX-90742183"
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-4.72816 -0.47420  0.37245  0.02671  0.81115  2.19274 

time="WD_0712"
factor="Factor5"
chr="3"
snp="AX-91393272"
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-7.01527 -0.14237  0.10774 -0.01013  0.44081  0.90041

# EB.09S.H.00417 -7.01527 to -1.3517660

time="WD_0712"
factor="Factor9"
chr="5"
snp="AX-90962437"
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.29899 -0.56094 -0.27138 -0.02123  0.43270  2.76859

time="WD_0712"
factor="Factor23"
chr="4"
snp="PZE-104044271"

time="WD_0712"
factor="Factor23"
chr="6"
snp="SYN25266"

time="WD_0712"
factor="Factor23"
chr="7"
snp="AX-91067829"
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.941511 -0.440604 -0.107716  0.006501  0.107382  6.772298 

# EB.09S.H.00005 6.7722978 to 1.9960830

time="WD_0712"
factor="Factor45"
chr="2"
snp="AX-91512352"

time="WD_0712"
factor="Factor45"
chr="4"
snp="AX-90861422"
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.22449 -0.45388 -0.11319 -0.03125  0.19684  6.38950

#EB.09S.H.00232 6.3895000 to 1.5069615

time="WD_0718"
factor="Factor24"
chr="1"
snp="AX-90602079"
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-1.909594 -0.314032 -0.070451 -0.007557  0.214186  4.876477 

# EB.09S.H.00238 4.8764768 to 2.9610413

time="WD_0718"
factor="Factor4"
chr="2"
snp="AX-91511436"
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.64221 -0.33769  0.11737  0.04057  0.48563  1.28875 

time="WD_0720"
factor="Factor19"
chr="4"
snp="AX-91628035"
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.2774905 -0.0732565  0.0059726 -0.0006669  0.0751776  0.2479977

time="WD_0718"
factor="Factor7"
chr="2"
snp="SYN14630"
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.896650 -0.226921  0.007321  0.025942  0.200922  1.227611

time="WD_0718"
factor="Factor12"
chr="5"
snp="AX-90973660"
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-2.73879 -0.61673  0.07850  0.03223  0.83934  2.89556

time="WD_0718"
factor="Factor18"
chr="2"
snp="AX-90734033"
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.8369866 -0.1598194  0.0003251  0.0030659  0.1954363  0.7724557

time="WD_0720"
factor="Factor19"
chr="4"
snp="AX-91628035"
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.2774905 -0.0732565  0.0059726 -0.0006669  0.0751776  0.2479977 

time="WD_0720"
factor="Factor2"
chr="4"
snp="AX-91221743"
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.27508 -0.27144  0.02224  0.02590  0.31406  1.42915

time="WD_0720"
factor="Factor17"
chr="7"
snp="AX-91066014"
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -2.03323 -0.54775 -0.01378  0.02044  0.53642  2.68465

time="WD_0727"
factor="Factor12"
chr="3"
snp="PZE-103115047"
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.792397 -0.300084  0.013460  0.002684  0.275017  1.046949
 
time="WD_0727"
factor="Factor17"
chr="3"
snp="AX-91587910" 
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-2.40680 -0.69669 -0.03544  0.03041  0.67546  3.19159

time="WD_0727"
factor="Factor19"
chr="8"
snp="AX-91777755" 
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-1.557740 -0.438290 -0.020812 -0.003534  0.342791  2.386980

#EB.09S.H.00482 2.386980 to 1.669264




for(i in 1:nrow(factoreqtl)){
	row=factoreqtl[i,]
	time=row$time
	factor=row$factor
	chr=row$CHR
	snp=row$X_ID
	F_values=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_all_F_means_FIXED.txt',time),data.table=F)
	summary(F_values[,factor])
	df=F_values[,c('V1',factor)]

	df=df[order(df[,factor]),]
	rownames(df)=seq(1,nrow(df))
	df$ID_f=factor(df$V1,levels=c(unique(df$V1)))
	rownames(df)=df$V1
	adj_chr=c("5","9")
	if(chr %in% adj_chr){
		X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
	}else{
		X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	}
	X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
	colnames(X) = founders
	rownames(X) = dimnames(X_list[[1]])[[1]]
	X=X[rownames(df),]
	founder_id=apply(X,MARGIN=1,function(x) names(x[which.max(x)]))
	df$founder=founder_id[match(df$V1,names(founder_id))]
	names(df)=c('V1' ,'factor','ID_f','founder')
	p=ggplot(aes(x=ID_f,y=factor),data=df) + 
    	geom_point(aes(color=founder)) +
     	xlab('Sample') +
     	ylab('F-value') + geom_hline(yintercept=0)

	pdf(sprintf('eqtl/trans/images/%s_%s_F_by_ind.pdf',time,factor))
	print(p)
	dev.off()
}


# DO they overlap with each other?
chroms=unique(factoreqtl$CHR)
for(chr in chroms){
	sub=factoreqtl[factoreqtl$CHR==chr,]
	print(sub[,c('time','factor','CHR','leftmost','alt_rightmost')])
}

#1 WD_0718  Factor4   2  11217105      14861212
#2 WD_0718 Factor18   2  12493114      14861212
####

pairs=fread('MegaLMM/MegaLMM_timepoint_residuals_factor_pairs.txt',data.table=F)

lambda20=fread('MegaLMM/MegaLMM_WD_0720_residuals_all_Lambda_means_FIXED.txt',data.table=F)
lambda12=fread('MegaLMM/MegaLMM_WD_0712_residuals_all_Lambda_means_FIXED.txt',data.table=F)
rownames(lambda12)=lambda12$V1
rownames(lambda20)=lambda20$V1
ginter=intersect(lambda12$V1,lambda20$V1)


lambda12=lambda12[ginter,]
lambda20=lambda20[ginter,]

f12=fread('MegaLMM/MegaLMM_WD_0712_residuals_all_F_means_FIXED.txt',data.table=F)
f20=fread('MegaLMM/MegaLMM_WD_0720_residuals_all_F_means_FIXED.txt',data.table=F)
inter=intersect(f12$V1,f20$V1)

###### Do they overlap with cis-eQTL?
eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

# For 10% FDR, 

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}
eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end

env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(factoreqtl)
setkey(env2,CHR,leftmost,alt_rightmost)
comp=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)

comp$Trait=paste0(comp$time,'-',comp$factor)
comp$ID=paste0(comp$Trait,'-',comp$CHR)


factoreqtl$Trait=paste0(factoreqtl$time,'-',factoreqtl$factor)
factoreqtl$ID=paste0(factoreqtl$Trait,'-',factoreqtl$CHR)

# Overlap of cis-eQTL for 15 of the 16 factor eqtl (not WD_0712-Factor23-4)

comp %>% group_by(ID) %>% count()



ngenes=comp %>% group_by(ID) %>% summarize(n=length(unique(i.Trait)))
ngenes=as.data.frame(ngenes)
ngenes
#                   ID   n
#1   WD_0712-Factor2-2  55
#2  WD_0712-Factor23-6   4
#3  WD_0712-Factor23-7  11
#4   WD_0712-Factor5-3  14
#5   WD_0712-Factor9-5   1
#6  WD_0718-Factor12-5   2
#7  WD_0718-Factor18-2  39
#8   WD_0718-Factor4-2  68
#39   WD_0718-Factor7-2   4
#10 WD_0720-Factor17-7  14
#11 WD_0720-Factor19-4 126
#12  WD_0720-Factor2-4   6
#13 WD_0727-Factor12-3  20
#14 WD_0727-Factor17-3  25
#15 WD_0727-Factor19-8   2
fwrite(comp,'eqtl/results/residual_factor_eQTL_cis_eQTL_overlap.txt',row.names=F,quote=F,sep='\t')

##### Do they overlap with QTL
zqtl=fread('QTL/all_adjusted_QTL_all_methods.txt',data.table=F)
#qtl=fread('QTL/all_adjusted_QTL_support_intervals.txt',data.table=F)

env1=qtl
env1=as.data.table(env1)
env2=as.data.table(factoreqtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp2=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

factoreqtl=fread('eqtl/results/all_factor_fdr_SIs_FIXED.txt',data.table=F)
env1=qtl
env1=as.data.table(env1)
env2=as.data.table(factoreqtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp1=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


#########################
factoreqtl=fread('eqtl/results/all_factor_trans_eqtl_fdr_hits_FIXED.txt',data.table=F)

# Plot out individuals by F-value colored by founder at SNP
#snp="AX-91833136"
#time="WD_0712"
#factor="Factor12"
#chr="10"

snp="AX-91511569"
time="WD_0718"
factor="Factor4"
chr="2"

F_values=fread(sprintf('MegaLMM/MegaLMM_%s_all_F_means_FIXED.txt',time),data.table=F)

summary(F_values[,factor])

interactors=c("Zm00001d020430","Zm00001d012584",
"Zm00001d031064",
"Zm00001d031451",
"Zm00001d039267",
"Zm00001d052890",
#GRMZM2G017586
"Zm00001d002364",
"Zm00001d002799",
"Zm00001d027900",
"Zm00001d002799",
"Zm00001d002799",
"Zm00001d005692",
"Zm00001d006551",
"Zm00001d010634",
"Zm00001d012527",
"Zm00001d013443",
"Zm00001d014858",
"Zm00001d014995",
"Zm00001d015421",
"Zm00001d044232",
"Zm00001d015549",
"Zm00001d016793",
"Zm00001d018571",
"Zm00001d020025",
"Zm00001d020409",
"Zm00001d020492",
"Zm00001d021019",
"Zm00001d021019",
"Zm00001d024324",
"Zm00001d024679",
"Zm00001d026398",
"Zm00001d026542",
"Zm00001d027846",
"Zm00001d031044",
"Zm00001d031717",
"Zm00001d034417",
"Zm00001d032024",
"Zm00001d032265",
"Zm00001d033267",
"Zm00001d034298",
"Zm00001d034601",
"Zm00001d035604",
"Zm00001d035604",
"Zm00001d037221",
"Zm00001d037221",
"Zm00001d037605",
"Zm00001d037605",
"Zm00001d037605",
"Zm00001d037605",
"Zm00001d037605",
"Zm00001d038357",
"Zm00001d038843",
"Zm00001d039260",
"Zm00001d042777",
"Zm00001d042907",
"Zm00001d042907",
"Zm00001d043950",
"Zm00001d047017",
"Zm00001d049364")
interactors=unique(interactors)
# of the 59 genes, 19 are loaded on the factor
ginter=intersect(genetable$Gene_ID,prop_var$V1)
ginter=ginter[!(ginter %in% interactors)]
totsum=c()
for(i in 1:500){
	draw=sample(ginter,49)
	sums=sum(draw %in% f4genes$V1)
	totsum=c(totsum,sums)
}
# Not more than you'd expect by chance

# Does this overlap with ciseQTL
factoreqtl=fread('eqtl/results/all_factor_fdr_peaks_FIXED.txt',data.table=F)
eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)
# Look only in time
eqtl=eqtl[eqtl$time==time,]

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}
eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end

factoreqtl$block_start=all_founder_blocks[match(factoreqtl$X_ID,all_founder_blocks$focal_snp),]$start
factoreqtl$block_end=all_founder_blocks[match(factoreqtl$X_ID,all_founder_blocks$focal_snp),]$end


env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(factoreqtl)
setkey(env2,CHR,leftmost,alt_rightmost)
comparison1=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)
# 46 cis-eQTL in this region in time WD_0718
# 36 cis-eQTL in this region in time WD_0712

# Do any of them have correlated effect sizes?
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
plot_list=list()
count=1
#adj_chr=c(5,9)
cors=c()
pvals=c()
for(i in 1:nrow(comparison1)){
	row=comparison1[i,]
	env=row$environment
	chr=row$CHR
	factor=row$Trait
	gene=row$i.Trait
	fsnp=row$X_ID
	esnp=row$i.X_ID
	etime=row$time
	
	effect_sizes=fread(sprintf('eqtl/trans/results/%s_residuals_c%.0f_%s_trans_results_FIXED.txt',time,chr,factor),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==fsnp,]
	effect_size=unlist(effect_size[,c(founders)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',etime,chr),data.table=F)
	results=results[results$X_ID==esnp & results$Trait==gene,]
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
		xlab("eQTL effect size (log2cpm)") + ylab("Factor trans-eQTL effect size") +
		ggtitle(sprintf("%s and %s, %s, r=%.2f",gene,factor,time,r))
		plot_list[[count]]=p1
		count=count+1
	}
	
}

comparison1$r=cors
comparison1$pvalue=pvals

fwrite(comparison1,'eqtl/results/factor_cis_eQTL_overlap.txt',row.names=F,quote=F,sep='\t')

pdf('eqtl/images/factor_cis_eQTL_effect_size_correlations.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()

# Does this overlap with QTL


qtl=fread('QTL/all_adjusted_QTL_peaks.txt',data.table=F)
qtl$block_start=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$start
qtl$block_end=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$end
env1=factoreqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,leftmost,alt_rightmost)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)
# No overlap




# What factor is closest to residual WD_0720 Factor 17?
lambda_r=fread('MegaLMM/MegaLMM_WD_0720_residuals_all_Lambda_means_FIXED.txt',data.table=F)
r17=lambda_r[,c('V1','Factor17')]
rownames(r17)=r17$V1

lambda=fread('MegaLMM/MegaLMM_WD_0720_all_Lambda_means_FIXED.txt',data.table=F)
rownames(lambda)=lambda$V1

ginter=intersect(r17$V1,lambda$V1)
lambda=lambda[ginter,]
r17=r17[ginter,]

for(i in 2:ncol(lambda)){
	column=lambda[,i]
	r=cor(r17$Factor17,column)
	print(names(lambda)[i])
	print(r)
}
# Factor 15 in full is closest to Factor 17 residuals r=0.808
prop_var=fread('MegaLMM/MegaLMM_WD_0720_prop_variance_FIXED.txt',data.table=F)
f15genes=prop_var[prop_var$Factor15>0.1,c('V1','Factor15')]

prop_varr=fread('MegaLMM/MegaLMM_residuals_WD_0720_prop_variance_FIXED.txt',data.table=F)
f17genes=prop_varr[prop_varr$Factor17>0.1,c('V1','Factor17')]


interactors=c(
"Zm00001d001945",
"Zm00001d053819",
"Zm00001d013707",
"Zm00001d051451",
"Zm00001d031064",
"Zm00001d014858",
"GRMZM2G017586",
"Zm00001d000184",
"Zm00001d034417",
"Zm00001d002762",
"Zm00001d024324",
"Zm00001d002799",
"Zm00001d005578",
"Zm00001d005692",
"Zm00001d006236",
"Zm00001d010634",
"Zm00001d012527",
"Zm00001d012810",
"Zm00001d031182",
"Zm00001d013399",
"Zm00001d013777",
"Zm00001d015407",
"Zm00001d015549",
"Zm00001d016793",
"Zm00001d016861",
"Zm00001d017575",
"Zm00001d020492",
"Zm00001d051439",
"Zm00001d017726",
"Zm00001d018081",
"Zm00001d020025",
"Zm00001d020267",
"Zm00001d021019",
"Zm00001d021291",
"Zm00001d024679",
"Zm00001d025141",
"Zm00001d044117",
"Zm00001d026398",
"Zm00001d034298",
"Zm00001d038357",
"Zm00001d026542",
"Zm00001d027846",
"Zm00001d029963",
"Zm00001d030678",
"Zm00001d031717",
"Zm00001d032024",
"Zm00001d032265",
"Zm00001d033050",
"Zm00001d050195",
"Zm00001d033267",
"Zm00001d033378",
"Zm00001d043950",
"Zm00001d034601",
"Zm00001d035604",
"Zm00001d037221",
"Zm00001d037605",
"Zm00001d038843",
"Zm00001d048623",
"Zm00001d039893",
"Zm00001d041576",
"Zm00001d042777",
"Zm00001d042907",
"Zm00001d044232",
"Zm00001d044409",
"Zm00001d046925",
"Zm00001d047017",
"Zm00001d048208",
"Zm00001d048647",
"Zm00001d049364",
"Zm00001d050816",
"Zm00001d051458",
"Zm00001d052229",
"Zm00001d053859")

# are any interactor genes higher loaded in the factor than expected by chance?
overlaps=f17genes[f17genes$V1 %in% interactors,]

# More than expected by chance
# If I grab 73 genes, how likely is it to get 6 genes in the top 95% quantile?
quantile(prop_varr$Factor17,0.975)
#      95% 
quant95=0.7115427 
#0.5142244
sgenes=prop_varr$V1[!(prop_varr$V1 %in% interactors)]
fsums=c()
for(i in 1:1000){
	draw=sample(sgenes,73)
	sub=prop_varr[prop_varr$V1 %in% draw,c('V1','Factor17')]
	fsum=sum(sub$Factor17>quant95)
	fsums=c(fsums,fsum)

}
overlaps[overlaps$Factor17>0.5142244,]
#                  V1  Factor17
#1763  Zm00001d032265 0.7912046
#2667  Zm00001d034601 0.6563778
#12231 Zm00001d038843 0.9661856
#12932 Zm00001d020267 0.7620315
#15449 Zm00001d012527 0.9310063
#16238 Zm00001d047017 0.7780453

# There are more highly loaded genes in this subset than expected by chance (5%)

# Are there cis-eQTL in this region on chr 7?
time="WD_0712"
eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)
# Look only in time
eqtl=eqtl[eqtl$time==time,]

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}
eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end

factoreqtl$block_start=all_founder_blocks[match(factoreqtl$X_ID,all_founder_blocks$focal_snp),]$start
factoreqtl$block_end=all_founder_blocks[match(factoreqtl$X_ID,all_founder_blocks$focal_snp),]$end

feqtl=factoreqtl[2,]
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(feqtl)
setkey(env2,CHR,leftmost,alt_rightmost)
comparison1=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)
# Factor 4 WD_0718 46 cis-eQTL in this region in time WD_0718
# Factor 4 WD_0718 36 cis-eQTL in this region in time WD_0712

# For WD_0720 Factor 17, 10 cis-eQTL overlap in WD_0720

# For WD_0712 Factor 23, No cis-eQTL overlap in WD_0712 (only a 2.7kb window, so not super suprising), but there are no cis-eQTL at that marker at all, even in other timepoints?
# Is this unusual? How many markers don't have cis-eQTL there?
# Do any of them have correlated effect sizes?
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
plot_list=list()
count=1
#adj_chr=c(5,9)
cors=c()
pvals=c()
for(i in 1:nrow(comparison1)){
	row=comparison1[i,]
	env=row$environment
	chr=row$CHR
	factor=row$Trait
	gene=row$i.Trait
	fsnp=row$X_ID
	esnp=row$i.X_ID
	etime=row$time
	
	effect_sizes=fread(sprintf('eqtl/trans/results/%s_residuals_c%.0f_%s_trans_results_FIXED.txt',time,chr,factor),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==fsnp,]
	effect_size=unlist(effect_size[,c(founders)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',etime,chr),data.table=F)
	results=results[results$X_ID==esnp & results$Trait==gene,]
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
		xlab("eQTL effect size (log2cpm)") + ylab("Factor trans-eQTL effect size") +
		ggtitle(sprintf("%s and %s, %s, r=%.2f",gene,factor,time,r))
		plot_list[[count]]=p1
		count=count+1
	}
	
}

comparison1$r=cors
comparison1$pvalue=pvals

fwrite(comparison1,'eqtl/results/WD_0720_Factor17_cis_eQTL_overlap.txt',row.names=F,quote=F,sep='\t')

pdf('eqtl/images/WD_0720_Factor17_cis_eQTL_effect_size_correlations.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()


qtl=fread('QTL/all_adjusted_QTL_peaks.txt',data.table=F)
qtl$block_start=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$start
qtl$block_end=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$end
env1=factoreqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,leftmost,alt_rightmost)
comparison2=foverlaps(env1,env2,by.x=c('CHR','leftmost','alt_rightmost'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)

qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
#qtl$block_start=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$start
#qtl$block_end=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$end
env1=factoreqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison2=foverlaps(env1,env2,by.x=c('CHR','leftmost','alt_rightmost'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
No Overlaps


####### Residuals ########
factoreqtl=fread('eqtl/results/all_residual_factor_fdr_peaks_FIXED.txt',data.table=F)
#factoreqtl=fread('eqtl/results/all_residuals_factor_trans_eqtl_fdr_hits_FIXED.txt',data.table=F)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

# Plot out individuals by F-value colored by founder at SNP
#snp="AX-91833136"
#time="WD_0712"
#factor="Factor12"
#chr="10"

snp="AX-91511569"
time="WD_0718"
factor="Factor4"
chr="2"

F_values=fread(sprintf('MegaLMM/MegaLMM_%s_all_F_means_FIXED.txt',time),data.table=F)

summary(F_values[,factor])



# Does this overlap with ciseQTL
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

factoreqtl$block_start=all_founder_blocks[match(factoreqtl$X_ID,all_founder_blocks$focal_snp),]$start
factoreqtl$block_end=all_founder_blocks[match(factoreqtl$X_ID,all_founder_blocks$focal_snp),]$end

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

# Look only in time
times=unique(factoreqtl$time)
for(time in times){
	teqtl=eqtl[eqtl$time==time,]
	feqtl=factoreqtl[factoreqtl$time==time,]
	env1=teqtl
	env1=as.data.table(env1)
	env2=as.data.table(feqtl)
	setkey(env2,CHR,leftmost,alt_rightmost)
	comparison1=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)
	if(nrow(comparison1)!=0){
		# WD_0718  37 Genes 
		# Do any of them have correlated effect sizes?
		plot_list=list()
		count=1
		#adj_chr=c(5,9)
		cors=c()
		pvals=c()
		for(i in 1:nrow(comparison1)){
			row=comparison1[i,]
			env=row$environment
			chr=row$CHR
			factor=row$Trait
			gene=row$i.Trait
			fsnp=row$X_ID
			esnp=row$i.X_ID
			effect_sizes=fread(sprintf('eqtl/trans/results/%s_residuals_c%.0f_%s_trans_results_FIXED.txt',time,chr,factor),data.table=F)
			effect_size=effect_sizes[effect_sizes$X_ID==fsnp,]
			effect_size=unlist(effect_size[,c(founders)])
			wn=which(!is.na(effect_size))[1]
			effect_size[-wn]=effect_size[-wn]+effect_size[wn]
			
			results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
			results=results[results$X_ID==esnp & results$Trait==gene,]
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
				xlab("eQTL effect size (log2cpm)") + ylab("Factor trans-eQTL effect size") +
				ggtitle(sprintf("%s and %s, %s, r=%.2f",gene,factor,time,r))
				plot_list[[count]]=p1
				count=count+1
			}
		}
		comparison1$r=cors
		comparison1$pvalue=pvals

		fwrite(comparison1,sprintf('eqtl/results/%s_resid_factor_cis_eQTL_overlap.txt',time),row.names=F,quote=F,sep='\t')
		pdf(sprintf('eqtl/images/%s_resid_factor_cis_eQTL_effect_size_correlations.pdf',time))
		for(i in 1:length(plot_list)){
			print(plot_list[[i]])
		}
		dev.off()
	}
}








# Does this overlap with QTL


qtl=fread('QTL/all_adjusted_QTL_peaks.txt',data.table=F)
qtl$block_start=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$start
qtl$block_end=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$end
env1=factoreqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,leftmost,alt_rightmost)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)
# No overlap

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


# WD_0712 Factor 12
 #   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-7.23787 -0.10174  0.09730  0.01619  0.46738  1.52755 

#WD_0718 Factor 4
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-5.4374 -1.1530  0.3931  0.1327  1.5628  4.2934 

df=F_values[,c('V1',factor)]

df=df[order(df[,factor]),]
rownames(df)=seq(1,nrow(df))
df$ID_f=factor(df$V1,levels=c(unique(df$V1)))
rownames(df)=df$V1
adj_chr=c("5","9")
if(chr %in% adj_chr){
	X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))

}else{
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
}

X = do.call(cbind,lapply(X_list,function(x) x[,snp]))

colnames(X) = founders
rownames(X) = dimnames(X_list[[1]])[[1]]
X=X[rownames(df),]
founder_id=apply(X,MARGIN=1,function(x) names(x[which.max(x)]))
df$founder=founder_id[match(df$V1,names(founder_id))]
names(df)=c('V1' ,'factor','ID_f','founder')
p=ggplot(aes(x=ID_f,y=factor),data=df) + 
    geom_point(aes(color=founder)) +
     xlab('Sample') +
     ylab('F-value') + geom_hline(yintercept=0)

pdf(sprintf('eqtl/trans/images/%s_%s_F_by_ind.pdf',time,factor))
print(p)
dev.off()


# What are the effect sizes for these Factors?
# Do they overlap with any cis-eQTL? Are their effect sizes correlated?

# If we drop out EB.09S.H.00417, is the factor trans-eQTL still significant?


##### Residuals #######



factoreqtl=fread('eqtl/results/all_residuals_factor_trans_eqtl_fdr_hits_FIXED.txt',data.table=F)

snp="AX-91833136"
time="WD_0712"
factor="Factor5"
chr="10"
# Outlier

snp="AX-91393272"
time="WD_0712"
factor="Factor5"
chr="3"
# Outlier

time="WD_0712"
factor="Factor23"
snp="PZE-104044271"
chr="4"
# Outlier

snp="AX-90733794"
time="WD_0718"
factor="Factor4"
chr="2"

time="WD_0720"
factor="Factor19"
snp="AX-91628035"
chr="4"

time="WD_0720"
factor="Factor17"
snp="AX-91066014"
chr="4"


F_values=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_all_F_means_FIXED.txt',time),data.table=F)

summary(F_values[,factor])



##################
factoreqtl=fread('eqtl/results/all_factor_fdr_peaks.txt',data.table=F)
factordf=fread('eqtl/results/all_factor_trans_eqtl_fdr_genes.txt',data.table=F)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
axisdf=fread('eqtl/data/chromosome_axis.txt',data.table=F)
cumtot=2105119857

# For each factor
#		- How many genes are loaded on them
#		- Where are the genes located in the genome
#		- What prop variance is explained by the factor for each gene
#		- How much variation in the factor is explained by the factor eQTL
#		- Are there trans-eQTL or cis-eQTL that overlap with the factor-eQTL variant? 
#		- Are there trans-eQTL or cis-eQTL genes that load on the factor?

factoreqtl$time_factor=paste0(factoreqtl$time,'-',factoreqtl$Trait)
tf=unique(factoreqtl$time_factor)
for(f in tf){
	sub=factoreqtl[factoreqtl$time_factor==f,]
	time=unique(sub$time)
	factor=unique(sub$factor)
	df=factordf[factordf$time==time & factordf$Trait==factor,]
	df$midgene=round(df$gene_start + (df$gene_end-df$gene_start)/2)
	df$tot=axisdf[match(df$gene_chr,axisdf$gene_chr),]$tot
	df$tot.x=axisdf[match(df$CHR,axisdf$gene_chr),]$tot
	BP_cum=unique(df$BP) + unique(df$tot.x)
	df$midgene_cum=df$midgene + df$tot 
	prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance.txt',time),data.table=F)
	print(time)
	print(factor)
	fgenes=prop_var[prop_var[,factor]>0.1,]$V1
	ngenes=length(fgenes)
	print(sprintf('%.0f Genes loaded on %s',ngenes,factor))
	
	p=ggplot(df,aes(x=midgene_cum,y=prop_var)) +
	geom_vline(xintercept=BP_cum,color="coral") +
	scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
	ylim(0.1,1) +
    #scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
    #geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
	geom_point(aes(color=prop_var)) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
	xlab("Position") + ylab("Proportion Variance") +
	ggtitle(sprintf("Gene Loadings on %s %s, n=%.0f",time,factor,ngenes)) +
	theme_classic() +
    theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_line(colour="black"))
 
	
	pdf(sprintf('eqtl/images/%s_%s_plot.pdf',time,factor))
	print(p)
	dev.off()
	# Plot x-axis location
	# y-axis prop_var
	# vhline of location of factor eQTL

}

#[1] "WD_0727"
#[1] "Factor22"
#[1] "6 Genes loaded on Factor22"

#[1] "WD_0712"
#[1] "Factor1"
#[1] "5030 Genes loaded on Factor1"

#[1] "WD_0718"
#[1] "Factor16"
#[1] "76 Genes loaded on Factor16"

#[1] "WD_0720"
#[1] "Factor9"
#[1] "2725 Genes loaded on Factor9"
#[1] "WD_0727"

#[1] "Factor14"
#[1] "955 Genes loaded on Factor14"

#[1] "WD_0727"
#[1] "Factor2"
#[1] "2103 Genes loaded on Factor2"

df=fread('eqtl/results/all_factor_fdr_peaks.txt',data.table=F)

qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
fqtl=qtl[qtl$Method=="Founder_probs",]
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

df$block_start=all_founder_blocks[match(df$X_ID,all_founder_blocks$focal_snp),]$start
df$block_end=all_founder_blocks[match(df$X_ID,all_founder_blocks$focal_snp),]$end

fqtl$block_start=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$start
fqtl$block_end=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$end


# Overlap of factor eQTL and QTL_F recombination blocks
env1=df
env1=as.data.table(env1)
env2=as.data.table(fqtl)
setkey(env2,Chromosome,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','block_start','block_end'),nomatch=NULL)
#### Factor 2 WD_0727 overlaps with qDTA8 in NERAC_2016_WD

# Overlap of factor eQTL and any QTL support interval
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
# WD_0727 Factor 2 overlaps with qDTA8/qDTS8 [1] "NERAC_2016_WD"     "GRANEROS_2015_OPT" "ALL"              
#[4] "STPAUL_2017_WD"    "BLOIS_2017_OPT"    "SZEGED_2017_OPT"

# Overlap of factor eQTL and any 10% QTL support interval
qtl10=fread('../GridLMM/Biogemma_10p_QTL.csv',data.table=F)
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl10)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
###Factor 2 WD_0727 eQTL hit overlaps with qTPH8 in GRANEROS_2015_OPT

# Overlap of factor eQTL and St.Paul recomb blocks
#DTA
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 4 SNPs for male flowering time in EXP_STPAUL overlap with Factor 2 WD_0727

#Just DTA peaks
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison5=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# None

#DTS
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 4 SNPs for female flowering time in EXP_STPAUL overlap with Factor 2 WD_0727

#Just DTS peaks
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison5=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 1 SNPs for female flowering time in EXP_STPAUL overlap with Factor 2 WD_0727

#tkw_15
qtl2=fread('QTL/tkw_15_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

### Overlap of factor eQTL with cis-eQTL?
ciseqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits.txt',data.table=F)
ciseqtl$block_start=all_founder_blocks[match(ciseqtl$X_ID,all_founder_blocks$focal_snp),]$start
ciseqtl$block_end=all_founder_blocks[match(ciseqtl$X_ID,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(ciseqtl)
setkey(env2,CHR,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# No overlap with cis-eQTL

# cis and trans
trans=fread('eqtl/results/all_trans_fdr_hits.txt',data.table=F)
trans$block_start=all_founder_blocks[match(trans$X_ID,all_founder_blocks$focal_snp),]$start
trans$block_end=all_founder_blocks[match(trans$X_ID,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(ciseqtl)
setkey(env2,CHR,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 8 overlaps 
# cis WD_0712 Zm00001d031961 overlaps with trans Zm00001d047592 WD_0727 and Zm00001d035217 WD_0727
# cis WD_0720 Zm00001d025017 overlaps with WD_0720 Zm00001d043918
# cis WD_0718 Zm00001d044657 overlaps with WD_0720 Zm00001d053411, WD_0718 Zm00001d017387, and WD_0718 Zm00001d019724
# cis WD_0718 Zm00001d032099 overlaps with WD_0720 Zm00001d046714 and WD_0712 Zm00001d031961

# With peaks - it is just  one cis-eQTL Zm00001d032099 (WD0718) overlaps with trans-eQTL for Zm00001d046714 (WD0720)

# trans and factor
env1=df
env1=as.data.table(env1)
env2=as.data.table(trans)
setkey(env2,CHR,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 97 overlapping SNPs for trans and factor eQTL 
# Factor 1 WD_0712 overlaps with Zm00001d013230 (cytochrome P450), loaded on the factor

# Factor 2 WD_0727 overlaps with one trans-eQTL from WD_0712, Zm00001d035087, not loaded on the factor

# Factor 14 WD_0727 overlaps 21 trans-eQTL

# Factor 16 WD_0718 overlaps with 3 trans-eqTL
# WD_0727 trans-eQTL Zm00001d023472, Zm00001d042948; WD_0720 Zm00001d016783 - none of which are loaded on the Factor

# Factor 22 in WD0727 overlaps with two trans-eQTL: Zm00001d041712 and Zm00001d009969 - both are loaded on the factor


# cis eQTL and QTL
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(fqtl)
setkey(env2,Chromosome,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','block_start','block_end'),nomatch=NULL)
# None

# Overlap of cis-eQTL and any QTL support interval
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
#None

# Overlap of cis eQTL and any 10% QTL support interval
qtl10=fread('../GridLMM/Biogemma_10p_QTL.csv',data.table=F)
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(qtl10)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
#None

# Overlap of factor eQTL and St.Paul recomb blocks
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

#DTS
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

#tkw_15
qtl2=fread('QTL/tkw_15_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

# trans eQTL and QTL highest SNP recomb block
env1=trans
env1=as.data.table(env1)
env2=as.data.table(fqtl)
setkey(env2,Chromosome,block_start,block_end)
comparison1=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','block_start','block_end'),nomatch=NULL)
# 17 genes and 10 QTL
#  [1] "qDTA3_2" "qDTS3_2" "qTPH7"   "qDTA9"   "qHGM3_1" "qDTS8"   "qDTA8"  
# [8] "qHGM3_2" "qTKW7_1" "qDTS9"

# Overlap of trans-eQTL and any QTL support interval
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
#2504 overlaps with QTL support intervals
# 114 trans eQTL overlap with 20 QTL IDs

# Overlap of trans eQTL and any 10% QTL support interval
qtl10=fread('../GridLMM/Biogemma_10p_QTL.csv',data.table=F)
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl10)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison3=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
#1-29 overlaps with 10% QTL support intervals
# 170 trans-eQTL overlap with 17 QTL IDs

# Overlap of trans eQTL and St.Paul recomb blocks
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison4=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 33 overlaps, 14 genes, 6 SNPs
#qDTA3_2 and WD_0727 genes Zm00001d038732 and WD_0712 Zm00001d053060
#qDTA8_1 WD_0712 Zm00001d027535 (8:127Mb)
#qDTA8_2 (135Mb-137Mb) WD_0712 Zm00001d023869, WD_0720 Zm00001d017387
#qDTA8_2 (142-148Mb) WD_0718 Zm00001d027184, WD_0720 Zm00001d033138
# qDTA8_3 (150-151Mb) WD_0718 Zm00001d035087 and WD_0720 Zm00001d007286
#qDTA9  and WD_0712 Zm00001d046357, WD_0712 Zm00001d034175, WD_0718 Zm00001d028161, WD_0727 Zm00001d026703, and WD_0712 Zm00001d012087

#Just DTA peaks
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison5=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#qDTA3_2 and WD_0727 genes Zm00001d038732 and WD_0712 Zm00001d053060
#qDTA9  and WD_0712 Zm00001d046357, WD_0712 Zm00001d034175, WD_0727 Zm00001d026703, and WD_0712 Zm00001d012087
#qDTA8_2 (142-148Mb) WD_0718 Zm00001d027184, WD_0720 Zm00001d033138


#DTS
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# qDTS8_2 (146Mb) WD_0720 Zm00001d033138 and WD_0718 Zm00001d027184, and WD_0720 Zm00001d007286


#Just DTS peaks
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison5=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

#tkw_15
qtl2=fread('QTL/tkw_15_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# Zm00001d041608 overlaps with qTKW2


# What about just trans peaks, not all SNPs? Do QTL and trans-eQTL have the same peak SNPs?
trans=fread('eqtl/results/all_trans_fdr_peaks.txt',data.table=F)
trans$block_start=all_founder_blocks[match(trans$X_ID,all_founder_blocks$focal_snp),]$start
trans$block_end=all_founder_blocks[match(trans$X_ID,all_founder_blocks$focal_snp),]$end

#DTA
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# WD_0712 Zm00001d053060 and qDTA3_2
# WD_0712 Zm00001d046357, WD_0712 Zm00001d034175, WD_0712 Zm00001d012087, WD_0727 Zm00001d026703 and qDTA9

#DTS
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

#TKW
qtl2=fread('QTL/tkw_15_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# None



# Calculate correlation between effect sizes for eQTL and phenotypes in overlapping regions
# First do this with ones that overlap for the peak

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


phenotypes=fread('phenotypes/EXP_STPAUL_2017_WD_phenotypes.csv',data.table=F)

# First, DTA and trans-eQTL peaks
pheno="male_flowering_d6"
cors=c()
ps=c()
for(i in 1:nrow(comparison)){
	row=comparison[i,]
	time=row$time
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
	gene=row$Trait
	snp=row$SNP
	chr=row$CHR
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])
	inter=intersect(exp$V1,inds)

	effect_sizes=fread(sprintf('QTL/Biogemma_chr%s_male_flowering_d6_x_EXP_STPAUL_2017_WD_vst_founderprobs.txt',chr),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
	
	results=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results.rds',time,chr))
	w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	results=results[[w]]
	results=results[results$X_ID==snp,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	ps=c(ps,p)
	
	df=data.frame(founder=founders,pheno=effect_size,eqts=betas,stringsAsFactors=F)
	p1=ggplot(df,aes(x=eqts,y=pheno)) + geom_point(aes(color=founder))
	png(sprintf('QTT/%s_%s_%s_%s_effect_sizes.png',time,gene,snp,pheno))
	print(p1)
	dev.off()
	#subpheno=phenotypes[phenotypes$ID %in% inter,c('ID','female_flowering_d6','male_flowering_d6','tkw_15')]
	#fdata$dta=subpheno[match(fdata$ID,subpheno$ID),]$male_flowering_d6
	#fdata$dts=subpheno[match(fdata$ID,subpheno$ID),]$female_flowering_d6

}

comparison$cor=cors
comparison$pvalue=ps

# trans Peaks and DTA
# Highest correlations are  with AX-91145110 for Zm00001d034175 (r=-0.5179,p=0.102) and Zm00001d046357 (r=-0.507,p=0.110)

# all trans hits and DTA




####### Fvalue for Factor 2 and DTA

# Factor 2 and DTS r=-0.2436067, p-value = 0.3816

cors=c()
ps=c()

for(i in 1:nrow(comparison)){
	row=comparison[i,]
	time=row$time
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
	gene=row$factor
	snp=row$SNP
	chr=row$CHR
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])
	inter=intersect(exp$V1,inds)
	
	effect_sizes=fread(sprintf('QTL/Biogemma_chr%s_female_flowering_d6_x_EXP_STPAUL_2017_WD_vst_founderprobs.txt',chr),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,founders])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
		
	results=fread(sprintf('eqtl/trans/results/%s_c%s_%s_trans_results.txt',time,chr,gene),data.table=F)
	#w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	#results=results[[w]]
	results=results[results$X_ID==snp,]
	betas=unlist(results[,founders])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	ps=c(ps,p)
}
comparison$cor=cors
comparison$pvalue=ps

# highest correlation with DTS and Factor 2 trans-eQTL is AX-91202104, r=-0.327, p.value=0.253


# DTA
cors=c()
ps=c()

for(i in 1:nrow(comparison)){
	row=comparison[i,]
	time=row$time
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
	gene=row$factor
	snp=row$SNP
	chr=row$CHR
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])
	inter=intersect(exp$V1,inds)
	
	effect_sizes=fread(sprintf('QTL/Biogemma_chr%s_male_flowering_d6_x_EXP_STPAUL_2017_WD_vst_founderprobs.txt',chr),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,founders])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
		
	results=fread(sprintf('eqtl/trans/results/%s_c%s_%s_trans_results.txt',time,chr,gene),data.table=F)
	#w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	#results=results[[w]]
	results=results[results$X_ID==snp,]
	betas=unlist(results[,founders])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	ps=c(ps,p)
}
comparison$cor=cors
comparison$pvalue=ps

# AX-9110697 r=-0.383, p-value=0.176 is strongest correlation

# How do I get a null expectation?

### Look more closely at WD_0727 Factor 2 3####
library('data.table')
library('ggplot2')

ld=fread('../stats/ld_decay/circos/ld_bundled_links_filtered.txt',data.table=F)
names(ld)=c('CHR_A','START_A','END_A','CHR_B','START_B','END_B','size')
chra=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(ld[x,]$CHR_A,'chr')[[1]][2]))
ld$CHR_A=chra
chrb=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(ld[x,]$CHR_B,'chr')[[1]][2]))
ld$CHR_B=chrb

pmap=c()
for(c in 1:10){
	p=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
	pmap=rbind(pmap,p)
}

snp="AX-91772415"
snp_pos=pmap[pmap$marker==snp,]$pos

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)


df=fread('eqtl/results/factor_transQTL_all.txt',data.table=F)
df=df[df$Factor=='Factor2',]
factor_groups=readRDS('MegaLMM/MegaLMM_WD_0727_factor_groups.rds')
fgenes=factor_groups[['Factor2']]$genes

genetable=genetable[genetable$Gene_ID %in% fgenes,]
 genetable %>% group_by(CHROM) %>% count
# A tibble: 10 × 2
# Groups:   CHROM [10]
#CHROM     n
#   <int> <int>
# 1     1   338
# 2     2   237
# 3     3   228
# 4     4   197
# 5     5   270
# 6     6   167
# 7     7   153
# 8     8   203
# 9     9   161
# 10    10   149


# grab rows where SNPs in factor eQTL are being compared to regions in genetanble

## plot out genes in factor relative to eQTL snp locations

## What regions are in high interchrom ld

## Group inds by founder and plot F-values

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

exp=fread('eqtl/normalized/WD_0727_voom_normalized_gene_counts_formatted.txt',data.table=F)
K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])

inter=intersect(exp$V1,inds)

X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
colnames(X) = founders
rownames(X) = inter




founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))
fdata=data.frame(ID=inter,founder=founder,stringsAsFactors=F)

f_all_means=fread('MegaLMM/MegaLMM_WD_0727_all_F_means.txt',data.table=F)

fdata$fvalue=f_all_means[match(fdata$ID,f_all_means$V1),]$Factor2

p1=ggplot(aes(x=founder,y=fvalue),data=fdata) + geom_boxplot() + geom_jitter() +
theme(axis.text.x = element_text(angle = 45, size=10,hjust=1))
png('images/WD_0727_Factor2_F_by_founder.png')
print(p1)
dev.off()

# Ind EB.10H.H.00025 may be a little bit of an outlier?

# Correlation with flowering time?


phenotypes=fread('phenotypes/EXP_STPAUL_2017_WD_phenotypes.csv',data.table=F)
#phenotypes$Genotype_code=gsub("-",".",phenotypes$Genotype_code)
subpheno=phenotypes[phenotypes$ID %in% inter,c('ID','female_flowering_d6','male_flowering_d6')]
fdata$dta=subpheno[match(fdata$ID,subpheno$ID),]$male_flowering_d6
fdata$dts=subpheno[match(fdata$ID,subpheno$ID),]$female_flowering_d6

cor(fdata$fvalue,fdata$dts)
#[1] -0.01771362
> cor(fdata$fvalue,fdata$dta)
#[1] -0.07511269
highest_SNPs=c("AX-91772402","AX-91772415","AX-91202104","AX-91107495")


results=fread('eqtl/trans/results/WD_0727_c8_Factor2_trans_results.txt',data.table=F)
results$log10p=-log10(results$p_value_ML)
results=results[results$X_ID %in% df$SNP,]
highest=results[which.max(results$log10p),]
betas=unlist(highest[,c(6,8:22)])
betas[-1]=betas[-1]+betas[1]

bvs=X %*% betas
fdata$bv=bvs
cor(fdata$bv,fdata$dta)
 #         [,1]
#[1,] 0.1572119

cor(fdata$dts,fdata$bv)
#          [,1]
#[1,] 0.1186966

pheno="male_flowering_d6"
env="STPAUL_2017_WD"
effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%s_%s_x_%s_founderprobs.rds',chr,pheno,env))
effect_sizes=unlist(unname(effect_sizes[effect_sizes$X_ID==snp,6:21]))
effect_sizes[-1]=effect_sizes[1] + effect_sizes[-1]

ft_bvs=X %*% effect_sizes
fdata$ft_bvs_stpaul=ft_bvs

cor(fdata$bv,fdata$ft_bvs)
# DTA BLUPS
#          [,1]
#[1,] 0.2387707

cor(fdata$bv,fdata$ft_bvs_stpaul)
#         [,1]
#[1,] 0.1786436

fdata=fdata[order(fdata$fvalue),]
rownames(fdata)=seq(1,nrow(fdata))
fdata$variable_f=factor(fdata$founder,levels=c(founders))

phenotypes=fread('phenotypes/phenotypes_BLUPS.csv',data.table=F)
fdata$dta_blup=phenotypes[match(fdata$ID,phenotypes$ID),'BLUP-male_flowering_d6']
cor(fdata$fvalue,fdata$dta_blup)
#[1] -0.0713566

phenotypes=fread('phenotypes/phenotypes_asi.csv',data.table=F)
phenotypes$Genotype_code=gsub("-",".",phenotypes$Genotype_code)
envs=c("BLOIS_2014_OPT","BLOIS_2017_OPT","GRANEROS_2015_OPT","NERAC_2016_WD","STPAUL_2017_WD","SZEGED_2017_OPT")
for(env in envs){
	subpheno=phenotypes[phenotypes$Loc.Year.Treat==env,]
	print(env)
	phenos=unlist(subpheno[match(fdata$ID,phenotypes$Genotype_code),]$male_flowering_d6)
	print(cor(fdata$fvalue,phenos,use='complete.obs'))
}
#[1] "BLOIS_2014_OPT"
#[1] -0.02306991
#[1] "BLOIS_2017_OPT"
#[1] -0.1415513
#[1] "GRANEROS_2015_OPT"
#[1] -0.02039115
#[1] "NERAC_2016_WD"
#[1] -0.09728809
#[1] "STPAUL_2017_WD"
#[1] 0.0223498
#[1] "SZEGED_2017_OPT"
#[1] -0.06438033
cor(fdata$fvalue,fdata$dta_blup)

colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder

p1=ggplot(aes(x=fvalue,y=ft_bvs),data=fdata)+geom_point(aes(color=founder)) +
scale_color_manual(values=colorcodes[levels(fdata$variable_f),]$hex_color,labels=levels(fdata$variable_f))+


png('images/WD_0727_Factor2_eQTL_by_DTA_BLUP.png')
print(p1)
dev.off()


p2=ggplot(aes(x=fvalue,y=dta),data=fdata)+geom_point(aes(color=founder)) +
scale_color_manual(values=colorcodes[levels(fdata$variable_f),]$hex_color,labels=levels(fdata$variable_f))+

png('images/WD_0727_Factor2_eQTL_by_DTA.png')
print(p2)
dev.off()

p3=ggplot(aes(x=bv,y=ft_bvs),data=fdata)+geom_point(aes(color=founder)) +
scale_color_manual(values=colorcodes[levels(fdata$variable_f),]$hex_color,labels=levels(fdata$variable_f))+
xlab("F Genetic Values") + ylab("DTA BLUP Breeding Values") + ggtitle(sprintf('Flowering Time by Factor 2 F for %s',snp))

png('images/WD_0727_Factor2_eQTL_Fbv_by_DTA_BLUP.png')
print(p3)
dev.off()


### What SNPs are nearest to the genes loading onto this factor
pmap$end=pmap$pos+1

env1=genetable
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(pmap)
#env2$end=env2$end-1
setkey(env2,chr,pos,end)
comparison=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('chr','pos','end'),nomatch=NULL)

# 1731 have SNPs within them

find_nearest_snp=function(row){
    index=which.min(abs(row$START-pmap[pmap$chr==row$CHROM,]$pos))
    return(pmap[index,]$marker)
}

nearest_snps=sapply(seq(1,nrow(genetable)),function(x) find_nearest_snp(genetable[x,]))
genetable$SNP=nearest_snps

#Instead of doing this, I need to find the founder recombination blocks that have those SNPs
# And use those instead of the nearest SNP
all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

env1=genetable
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(all_founder_blocks)
#env2$end=env2$end-1
setkey(env2,chr,start,end)
comparison=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('chr','start','end'),nomatch=NULL)

## Calculate breeding value of DTA using SNPs near the genes in Factor 2

bvs=c()

env="ALL"
pheno="male_flowering_d6"
for(chr in 1:10){
	subtable=genetable[genetable$CHROM==chr,]
	subcomp=comparison[comparison$CHROM==chr,]
	snps=unique(subcomp$focal_snp)
	
	effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',chr,pheno,env))
	
	es=effect_sizes[effect_sizes$X_ID %in% snps,]
	

	#effect_sizes[-1]=effect_sizes[1] + effect_sizes[-1]
	
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])

	inter=intersect(exp$V1,inds)

	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
}


#########


prop_var=fread('MegaLMM/MegaLMM_WD_0727_prop_variance.txt',data.table=F)
# How much variation is enough? 
subvar=prop_var[,c('V1','Factor2')]
cutoff=0.9
test=subvar[subvar$Factor2>=cutoff,]$V1

testlist=data.frame(genes=test,stringsAsFactors=F)
fwrite(testlist,sprintf('WD_0727_Factor2_%.1f_gene_list.txt',cutoff),row.names=F,col.names=F,quote=F,sep='\t')


######### GOSeq ########
pheno_factors=c('Factor2','Factor14')

#pheno_factors=unique(pheno_df$factor)
genes=names(avg_exp)
inter2=length(intersect(names(avg_exp),annotation$Gene))
#genelength=genetable[match(genes,genetable$Gene_ID),]$LENGTH
#names(genelength)=genes

fulllist=data.frame(genes=genes,stringsAsFactors=F)
fwrite(fulllist,'WD_0727_full_gene_list.txt',row.names=F,col.names=F,quote=F,sep='\t')


ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
# Enrichment of flowering time genes?
cutoff=0.1
test=subvar[subvar$Factor2>=cutoff,]$V1
gtable=genetable[genetable$Gene_ID %in% test,]


find_nearest_snp=function(row){
    index=which.min(abs(row$START-pmap[pmap$chr==row$CHROM,]$pos))
    return(pmap[index,]$marker)
}

nearest_snps=sapply(seq(1,nrow(gtable)),function(x) find_nearest_snp(gtable[x,]))
gtable$SNP=nearest_snps

snplist=gtable[,'SNP',drop=F]
fwrite(snplist,'WD_0727_Factor2_snplist.txt',row.names=F,col.names=F,sep='\t',quote=F)

####


snplist=fread('WD_0727_Factor2_snplist.txt',data.table=F,header=F)
snplist=snplist[-1,]
snp="AX-91772415"

lambda=fread('MegaLMM/MegaLMM_WD_0727_all_Lambda_means.txt',data.table=F)
ld=fread('Factor2_rsquared.ld',data.table=F)
ld=ld[ld$SNP_A==snp,]

ld2=ld[ld$SNP_B %in% snplist,]


##### WD_0727 trans-eQTL

eqtl=fread('eqtl/results/WD_0727_trans_eQTL_scan_hits.txt',data.table=F)

pmap=c()
for(c in 1:10){
	p=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
	pmap=rbind(pmap,p)
}

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

eqtl$Gene_CHR=genetable[match(eqtl$Gene,genetable$Gene_ID),]$CHROM
eqtl$Gene_START=genetable[match(eqtl$Gene,genetable$Gene_ID),]$START
eqtl$Gene_END=genetable[match(eqtl$Gene,genetable$Gene_ID),]$END

#[1] "Zm00001d041650" "Zm00001d020311" "Zm00001d009688" "Zm00001d045677"
#[5] "Zm00001d047592"

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

eqtl$block_start=all_founder_blocks[match(eqtl$SNP,all_founder_blocks$focal_snp),]$start
eqtl$block_end=all_founder_blocks[match(eqtl$SNP,all_founder_blocks$focal_snp),]$end

fqtl$block_start=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$start
fqtl$block_end=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$end


env1=eqtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(qtl)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

f14_27=fread('eqtl/results/Factor14_trans_WD_0727_eQTL_fkeep_hits.txt',data.table=F)

goi=c("Zm00001d041650","Zm00001d009688","Zm00001d045677")
factor_groups=readRDS('MegaLMM/MegaLMM_WD_0727_factor_groups2.rds')
f2genes=factor_groups[['Factor2']]$genes
f14genes=factor_groups[['Factor14']]$genes

# nearest SNPs for genes in Factor 14

find_nearest_snp=function(row){
    index=which.min(abs(row$START-pmap[pmap$chr==row$CHROM,]$pos))
    return(pmap[index,]$marker)
}

gtable=genetable[genetable$Gene_ID %in% f14genes,]

nearest_snps=sapply(seq(1,nrow(gtable)),function(x) find_nearest_snp(gtable[x,]))
gtable$SNP=nearest_snps

snplist=gtable[,'SNP',drop=F]
snplist=rbind(snplist,data.frame(SNP=unique(f14_27$SNP)))
snp2=unique(snplist$SNP)
snp2=data.frame(SNP=snp2,stringsAsFactors=F)

fwrite(snp2,'WD_0727_Factor14_snplist.txt',row.names=F,col.names=F,sep='\t',quote=F)

ld=fread('Factor14_rsquared.ld',data.table=F)
ld=ld[ld$SNP_A==snp,]

ld2=ld[ld$SNP_B %in% snplist$V1,]


#####
# Correlation between Fvalues and phenotypes
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
rownames(exp)=exp$V1
exp=exp[,-1]
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)

metadata=metadata[metadata$experiment==time,]


geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
exp=exp[,kept_genes]

genes=names(exp)

#genes=genes[1:5]
df=fread('eqtl/results/factor_transQTL_all.txt',data.table=F)
df %>% group_by(time,Factor) %>% filter(value==max(value)) %>% arrange(time,Factor,SNP)

# A tibble: 10 × 6
# Groups:   time, Factor [3]
#   Factor     CHR        BP SNP           value time   
#   <chr>    <int>     <int> <chr>         <dbl> <chr>  
# 1 Factor16     4   5535185 AX-90856708    2.90 WD_0718
# 2 Factor16     4   5414982 AX-91597662    2.90 WD_0718
# 3 Factor16     4   5663654 AX-91849936    2.90 WD_0718
# 4 Factor14     7 171058661 AX-91065358    3.95 WD_0727
# 5 Factor14     7 170989903 AX-91743057    3.95 WD_0727
# 6 Factor14     7 170886882 PZE-107118743  3.95 WD_0727
# 7 Factor2      8 152630937 AX-91107495    3.43 WD_0727
# 8 Factor2      8 150428041 AX-91202104    3.43 WD_0727
# 9 Factor2      8 150234596 AX-91772402    3.43 WD_0727
#10 Factor2      8 150347882 AX-91772415    3.43 WD_0727

time="WD_0727"
factor="Factor2"
snp="AX-91772402"
chr="8"

time="WD_0727"
factor="Factor14"
snp="AX-91743057"
chr="7"

time="WD_0718"
factor="Factor16"
snp="AX-91597662"
chr="4"

f_all_means=fread(sprintf('MegaLMM/MegaLMM_%s_all_F_means.txt',time),data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])

inter=intersect(exp$V1,inds)

X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
colnames(X) = founders
rownames(X) = inter


founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))
fdata=data.frame(ID=inter,founder=founder,stringsAsFactors=F)


fdata$fvalue=f_all_means[match(fdata$ID,f_all_means$V1),factor]


find_qtts=function(e,p){
	pheno=phenotypes[phenotypes$Loc.Year.Treat==e,c('Genotype_code',p)]
	pheno$f=fdata[match(pheno$Genotype_code,fdata$ID),'fvalue']
	test=cor.test(pheno[,p],pheno$f,use="complete.obs")
	return(test)
}

envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:8]

all_qtts=data.frame(matrix(ncol=5,nrow=0))
names(all_qtts)=c('pheno','env','factor','r','pvalue')

for(p in phenos){
	for(e in envs){
		correlations=find_qtts(e,p)
		line=data.frame(pheno=p,env=e,factor=factor,r=unlist(correlations['estimate']),pvalue=unlist(correlations['p.value']),stringsAsFactors=F)
		all_qtts=rbind(all_qtts,line)
	}
}

ntests=48
threshold=0.05/3
all_qtts$padjust=p.adjust(all_qtts$pvalue,method='fdr')

# WD_0727 Factor 2
### total plant height was actually th emost correlated with f-values
# r of -0.17
# none of them are significant though
# Is there a way I can compare this to a null distribution rather than fdr correction?

# WD_0727 Factor 14
#### highest correlation is with male_flowering_d6 (0.15) and tkw_15 (-0.15)
## neither of them are significant after multiple test correction

# WD_0718 Factor16
#### highest correlation is with grain_yield_15 (r=0.18) and total_plant_height (r=0.16) in SZEGED_2017_OPT
# Neither are significant after fdr correction

# Raw F values and raw phenotypes are not significantly correlated
# What about phenotype genetic values?


fwrite(all_qtts,sprintf('eqtl/results/%s_%s_F_pheno_corrs.txt',time,factor),row.names=F,quote=F,sep='\t')


#### Do any factors have enrichment of FT genes? ##########
time="WD_0718"
prop_var=fread(sprintf('MegaLMM/MegaLMM_residuals_%s_prop_variance_FIXED.txt',time),data.table=F)
ft_genes=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
names(ft_genes)=c('CHR','START','END','GENE_ID')

allft=c()

factors=names(prop_var)[-1]
for(factor in factors){
	fgenes=prop_var[prop_var[,factor]>=0.1,]$V1
	tn=sum(fgenes %in% ft_genes$GENE_ID)
	ft_geneoverlap=data.frame(factor=factor,ngenes=length(fgenes),inft=tn,stringsAsFactors=F)
	allft=rbind(allft,ft_geneoverlap)
}

allgenes=prop_var$V1

null=c()
for(i in 1:nrow(allft)){
	row=allft[i,]
	ngenes=row$ngenes
	tgenes=row$inft
	alltn=c()
	for(j in 1:1000){
		draw=sample(allgenes,ngenes)
		tn=sum(draw %in% ft_genes$GENE_ID)
		alltn=c(alltn,tn)
	}
	adjust=68
	q5=quantile(alltn,1-(0.05/adjust))
	null=c(null,q5)

}
allft$null=null

allft[allft$inft>allft$null,]
#     factor ngenes inft      null
#2   Factor2   1642   83  68.53088
#4   Factor4   1126   51  49.26544
#8   Factor8   4446  158 155.79632
#10 Factor10    308   19  17.53088

f4genes=prop_var[prop_var$Factor4>0.1,c('V1','Factor4')]
f4genes=merge(f4genes,genetable,by.x='V1',by.y="Gene_ID")

rap27_1="Zm00001d010987" # Not epxressed
rap27_2="Zm00001d010988" # Not epxressed
zcn8="Zm00001d010752" # Factor 4, 8, and 14
mads69="Zm00001d042315" # Factor 17

ftf4=f4genes[f4genes$V1 %in% ft_genes$GENE_ID,]
f2=merge(ftf4,genetable,by.x='V1',by.y="Gene_ID")







ft_geneoverlap=c()
for(i in 1:length(factor_groups)){
  #inds=which(pheno_factors==f)
  #i=pheno_loc[inds]
  fgenes=factor_groups[[i]]$genes
  f=factor_groups[[i]]$factor
  ngenes=length(fgenes)
  match=intersect(fgenes,ft_genes$GENE_ID)
  if(f %in% ft_factors){
    inft=T
  }else{
    inft=F
  }
  ft_geneoverlap=rbind(ft_geneoverlap,c(f,length(match),ngenes,inft))
}
ft_geneoverlap=as.data.frame(ft_geneoverlap,stringsAsFactors=F)
names(ft_geneoverlap)=c('factor','overlap','ngenes','inft')
ft_geneoverlap$overlap=as.numeric(ft_geneoverlap$overlap)
ft_geneoverlap$ngenes=as.numeric(ft_geneoverlap$ngenes)
ft_geneoverlap$inft=factor(ft_geneoverlap$inft,levels=c('FALSE','TRUE'))
m1=lm(overlap~ngenes+inft,ft_geneoverlap)
anova(m1)



### What do Fvalues look like for Factor 14
time="WD_0727"
fvalues=fread('MegaLMM/MegaLMM_WD_0727_all_F_means.txt',data.table=F)
snp="PZE-107118743"
samples=fread(sprintf('eqtl/data/%s_samples.txt',time),data.table=F,header=F)
chr="7"

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder

#pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
X_list=lapply(X_list,function(x) x[samples$V1,])
X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
#frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
#fkeep=founders[frep2>3]

certain=apply(X,MARGIN=1,function(x) sum(x>0.75)>0)
X=X[certain,]
founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))
fvalues=fvalues[,c('V1','Factor14')]
rownames(fvalues)=fvalues$V1
fvalues=fvalues[certain,]
fvalues$founder=founder
fvalues=fvalues[order(fvalues$Factor14),]
fvalues$ind_f=factor(fvalues$V1,levels=c(fvalues$V1))
fvalues$founder_f=factor(fvalues$founder,levels=c(founders))



p1=ggplot(aes(x=ind_f,y=Factor14),data=fvalues)+geom_point(aes(color=founder_f)) +
scale_color_manual(values=colorcodes[levels(fvalues$founder_f),]$hex_color,labels=levels(fvalues$founder_f))


png('images/WD_0727_Factor14_Fvalue_by_ind.png')
print(p1)
dev.off()