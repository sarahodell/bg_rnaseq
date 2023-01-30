#!/usr/bin/env Rscript

library('data.table')

time="WD_0712"


phenotype=fread(sprintf('eqtl/results/%s_vst_counts_x_phenotypes.txt',time),data.table=F)
cis=fread(sprintf('eqtl/results/%s_cis_eQTL_fkeep_hits.txt',time),data.table=F)
factor_groups=readRDS(sprintf('MegaLMM/pheno_MegaLMM_%s_factor_groups.rds',time))

pmaps=c()
all_X_list=
for(c in 1:10){
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  pmaps=rbind(pmaps,pmap)
}


#X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
vardf=c()

for(r in 1:nrow(cis)){
  g=cis[r,]$Gene
  s=cis[r,]$SNP
  chr=pmaps[pmaps$marker==s,]$chr
  results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_fkeep_results.txt',time,chr),data.table=F)
  results=results[!is.na(results$p_value_ML),]
  pheno=phenotype[,c('ID',g)]
  varp=var(unlist(unname(phenotype[,g])))
  X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',chr))
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
