#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
factor=as.character(args[[2]])

library('data.table')
library('lme4')
library('lme4qtl')
library('emmeans')
library('ggplot2')
library('GridLMM')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


#qtl=fread(sprintf('eqtl/results/%s_pheno_residuals_trans_%s_eQTL_hits.txt',factor,time),data.table=F)
qtl=fread(sprintf('eqtl/results/%s_pheno_trans_%s_eQTL_hits.txt',factor,time),data.table=F)
#phenotype=fread(sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_all_F_means.txt',time),data.table=F)
phenotype=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_all_F_means.txt',time),data.table=F)


phenotype=phenotype[,c('V1','V1',factor)]
names(phenotype)=c('ID','ID2','y')


full_df=c()
for(c in 1:10){
  #d=fread(sprintf('eqtl/trans/results/%s_c%s_pheno_residuals_factor_trans_eQTL.txt',time,c))
  d=fread(sprintf('eqtl/trans/results/%s_c%s_pheno_factor_trans_eQTL.txt',time,c))
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  d$CHR=c
  d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
  full_df=rbind(full_df,d)
}
full_df=full_df[full_df$Trait==factor,]
full_nrow=nrow(full_df)


threshold=-log10(0.05/full_nrow)
#threshold=full_threshold
true=c()
new_value=c()
for(i in 1:nrow(qtl)){
  row=qtl[i,]
  snp=row$SNP
  chr=as.character(row$CHR)

  K=fread('../GridLMM/K_matrices/K_matrix_full.txt'),data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)
  inter=intersect(rownames(K),phenotype$ID)
  K=K[inter,inter]
  phenotype=phenotype[phenotype$ID %in% inter,]

  X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
  X_list=lapply(X_list,function(x) x[inter,])
  X = do.call(cbind,lapply(X_list,function(x) x[,snp]))

  frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
  fkeep=founders[frep2>=2]

  X_list_ordered = lapply(fkeep,function(i) X[,i,drop=F])

  null_model = GridLMM_ML(y~1 + (1|ID),phenotype,relmat=list(ID=K),ML=T,REML=F,verbose=F)

  h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
  names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
  h2_start
  V_setup=null_model$setup
  Y=as.matrix(phenotype$y)
  X_cov=null_model$lmod$X

    #X_cov is n x 0 matrix and can use full X_list_ordered
    #X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,,drop=F])

  X_list_null=NULL
    #X_list_null
  cores=1
  gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F,h2_step=1  )
  true=c(true,-log10(gwas$p_value_ML)>=threshold)
  new_value=c(new_value,-log10(gwas$p_value_ML))
}

qtl$true_sig=true
qtl$new_value=new_value
fwrite(qtl,sprintf('eqtl/results/%s_pheno_trans_%s_eQTL_hits.txt',factor,time),row.names=F,quote=F,sep='\t')
#fwrite(qtl,sprintf('eqtl/results/%s_pheno_residuals_trans_%s_eQTL_hits.txt',factor,time),row.names=F,quote=F,sep='\t')
