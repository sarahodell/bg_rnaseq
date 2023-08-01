#!/usr/bin/env Rscript

library('data.table')

time="WD_0712"


phenotype=fread(sprintf('eqtl/noramlized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
cis=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
cis=cis[cis$time==time,]
#factor_groups=readRDS(sprintf('MegaLMM/pheno_MegaLMM_%s_factor_groups.rds',time))

pmaps=c()
all_X_list=
for(c in 1:10){
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  pmaps=rbind(pmaps,pmap)
}

adj_chr=c(5,9)

#X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
vardf=c()

h2=fread(sprintf('eqtl/data/%s_SNP_resid_h2.txt',time),data.table=F)
for(r in 1:nrow(cis)){
	g=cis[r,]$Gene
	s=cis[r,]$SNP
	chr=pmaps[pmaps$marker==s,]$chr
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_weights_results_FIXED.txt',time,chr),data.table=F)
	results=results[!is.na(results$p_value_ML),]
	pheno=phenotype[,c('ID',g)]
	varp=var(unlist(unname(phenotype[,g])))
	if(chr %in% adj_chr){
		X_list=readRDS(sprintf('phenotypes/bg%.0f_adjusted_genoprobs.rds',chr))
	}else{
		X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',chr))
	}
	X = do.call(cbind,lapply(X_list,function(x) x[,s]))
	results=results[results$Trait==g & results$X_ID==s,]
	betas=unlist(unname(results[,c(6,10:24)]))
	betas[which(is.na(betas))]=0
	cis_exp=X %*% betas
	pheno$cis_exp=cis_exp[match(pheno$ID,rownames(cis_exp)),1]
	names(pheno)=c('ID','y','cis')
	m1=lm(y~cis,data=pheno)
	r2=summary(m1)$r.squared
	betas=unlist(unname(results[,c(6,10:24)]))
	betas[-1]=betas[-1]-betas[1]
	inv_betas=2^betas
	inv_betas[-1]=inv_betas[-1]-inv_betas[1]
	inv_betas[which(is.na(inv_betas))]=0
	inv_cis_exp=X %*% inv_betas
	pheno$y_inv=2^pheno$y
	varp_inv=var(pheno$y_inv)
	pheno$cis_inv=2^pheno$cis
	m2=lm(y_inv~cis_inv,data=pheno)
	r2_inv=summary(m2)$r.squared
	#cis_prop=var_cis/varp
	vardf=rbind(vardf,c(g,chr,s,r2,varp,r2_inv,varp_inv))
}

vardf=as.data.frame(vardf,stringsAsFactors=F)
names(vardf)=c('Gene','Chr','SNP','r2','VarP','r2_inv','VarP_inv')

fwrite(vardf,sprintf('eqtl/cis/results/cis_variance_gene_exp_%s.txt',time),row.names=F,quote=F,sep='\t')


# How much variation in the whole plant phenotype a gene is clustered in is explained by cis-eQTL
pheno_df=c()
for(i in 1:nrow(cis)){
  row=cis[i,]
  g=row$Gene
  factors=which(unlist(unname(lapply(factor_groups,function(x) g %in% x$genes))))
  for(f in factors){
    if(!is.null(factor_groups[[f]]$phenotypes)){
      pheno_df=rbind(pheno_df,c(g,f))
    }
  }
}
pheno_df=as.data.frame(pheno_df,stringsAsFactors=F)
names(pheno_df)=c('Gene','Factor_index')

# Test with "BLOIS_2014_OPT-grain_yield_15"
phenotype=fread(sprintf('eqtl/results/%s_vst_counts_x_phenotypes.txt',time),data.table=F)
subdf=pheno_df[pheno_df$Factor_index==3,]

subpheno=phenotype[,c('ID',subdf$Gene,factor_groups[[3]]$phenotypes)]
names(subpheno)=c('ID','Zm00001d003592','Zm00001d004300','Zm00001d044170','y')
m3=lm(y~Zm00001d003592+Zm00001d004300+Zm00001d044170,subpheno)

subpheno$y_scaled=(subpheno$y-mean(subpheno$y,na.rm=T))/sd(subpheno$y,na.rm=T)

# R2=0.0368 so 3% of variation in yield explained by this gene
subpheno$Zm00001d003592_inv=2^subpheno$Zm00001d003592
subpheno$Zm00001d004300_inv=2^subpheno$Zm00001d004300
subpheno$Zm00001d044170_inv=2^subpheno$Zm00001d044170
m4=lm(y_scaled~Zm00001d003592_inv,subpheno)


m5=lm(y_scaled~Zm00001d004300_inv,subpheno)
# R2=0.000302, 0.03% of yield explained by this gene

m6=lm(y_scaled~Zm00001d044170_inv,subpheno)
# R2=0.003342, 0.03% of yield explained by this gene


############ QTL variance explained ############
library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('data.table')
library('dplyr')
library('stringr')


# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
#phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
#phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]
adj_chr=c("5","9")
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

qtl=fread('QTL/all_adjusted_QTL_peaks_trimmed.txt',data.table=F)
prop_vars=c()
for(i in 1:nrow(qtl)){
	row=qtl[i,]
	pheno=row$phenotype
	env=row$environment
	chr=row$CHR
	snp=row$SNP
	
	# Read in Kinship Matrix
	K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
	rownames(K)=K[,1]
	rownames(K)=gsub("-",".",rownames(K))
	K=as.matrix(K[,-1])
	colnames(K)=rownames(K)
	data=phenotypes[phenotypes$Loc.Year.Treat==env,]
	data=data[,c('Genotype_code',pheno)]
	names(data)=c('ID','y')
	#data=data.frame(ID=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
	data=data[!is.na(data$y),]
	rownames(data)=data$ID
	# Read in the haplotype group probabilities
	# Filter genotypes that are not in the K matrix
	if(chr %in% adj_chr){
		X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))

	}else{
		X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))

	}
	inds=rownames(X_list[[1]])
	inter=intersect(data$ID,inds)
	X_list=lapply(X_list,function(x) x[inter,])
	data=data[inter,]
	data$y=(data$y-mean(data$y))/sd(data$y)
	K=K[inter,inter]
	X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
	colnames(X) = founders
	rownames(X) = dimnames(X_list[[1]])[[1]]
	frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
	fkeep=founders[frep2>3]
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	betas=unlist(effect_size[,c(6:21)])
	full_betas=betas[founders]
	names(full_betas)=founders
	full_betas[!(names(full_betas) %in% fkeep)]=0
	wn=which(!is.na(betas))[1]
	intercept=names(betas)[wn]
	full_betas[names(full_betas)!=intercept]=full_betas[intercept]+full_betas[names(full_betas)!=intercept]
	bv=X %*% full_betas
	colnames(bv)=snp
	X_r=data$y-bv
	prop_var=var(bv[,snp],na.rm=T)/var(data$y,na.rm=T)
	prop_vars=c(prop_vars,prop_var)
}
qtl$prop_var=prop_vars
fwrite(qtl,'QTL/all_adjusted_QTL_peaks_trimmed.txt',row.names=F,quote=F,sep='\t')
