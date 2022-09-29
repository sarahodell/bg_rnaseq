args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
chr=as.character(args[[2]])
cores=as.numeric(args[[3]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM')
library('data.table')
library('dplyr')
#library('lme4')

phenotype=fread(sprintf('MegaLMM/MegaLMM_%s_all_F_means.txt',time),data.table=F)

factor_results=c()

for(k in 2:ncol(phenotype)){
  p=names(phenotype)[k]

  data=data.frame(ID=phenotype$V1,ID2=phenotype$V1,y=phenotype[,p],stringsAsFactors=F)
  data=data[!is.na(data$y),]
  data$y=(data$y-mean(data$y))/sd(data$y)


  K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)

  X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
  inds=rownames(X_list[[1]])
  i=intersect(inds,data$ID)

  null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

  h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
  names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
  h2_start
  V_setup=null_model$setup
  Y=as.matrix(data$y)
  X_cov=null_model$lmod$X

  X_list_ordered=lapply(X_list,function(x) x[i,])
  X_list_null=NULL

  gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
  gwas$Trait=p
  factor_results=rbind(factor_results,gwas)
}

fwrite(factor_results,sprintf('eqtl/trans/results/%s_c%s_factor_trans_eQTL.txt',time,chr),row.names=F,quote=F,sep='\t')
