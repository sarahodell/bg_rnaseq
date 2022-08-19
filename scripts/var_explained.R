#!/usr/bin/env Rscript

library('data.table')

time="WD_0712"


phenotype=fread('eqtl/normalized/WD_0712_voom_normalized_gene_counts.txt',data.table=F)
cis=fread('eqtl/results/cis_eQTL_WD_0712_all_results.txt',data.table=F)

pmaps=c()
all_X_list=
for(c in 1:10){
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  pmaps=rbind(pmaps,pmap)
}


#X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
vardf=c()

cis=cis[!is.na(cis$p_value_ML),]
for(r in 1:nrow(cis)){
  g=cis[r,]$Trait
  s=cis[r,]$X_ID
  chr=pmaps[pmaps$marker==s,]$chr

  varp=var(unlist(unname(phenotype[phenotype$V1==g,2:97])))
  X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',chr))
  X = do.call(cbind,lapply(X_list,function(x) x[,s]))
  betas=unlist(unname(cis[r,c(6,10:24)]))
  betas[which(is.na(betas))]=0

  var_cis=var(X %*% betas)

  cis_prop=var_cis/varp
  vardf=rbind(vardf,c(g,chr,s,varp,var_cis,cis_prop))
}

vardf=as.data.frame(vardf,stringsAsFactors=F)
names(vardf)=c('Gene','Chr','SNP''VarP','VarGcis','cis_prop')
