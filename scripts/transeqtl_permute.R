#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
factor=as.character(args[[2]])
chr=as.character(args[[3]])
cores=as.numeric(args[[4]])
reps=as.numeric(args[[5]])

library('GridLMM')
library('data.table')
library('dplyr')
library('parallel')
library('MASS')
library('stringr')

#all_reps=c()

phenotype=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_all_F_means.txt',time),data.table=F)

#factor_results=c()
#p=names(phenotype)[factor]

data=data.frame(ID=phenotype$V1,ID2=phenotype$V1,y=phenotype[,factor],stringsAsFactors=F)
data=data[!is.na(data$y),]
data$y=(data$y-mean(data$y))/sd(data$y)
  #K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])
inter=intersect(inds,data$ID)

null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

nmarkers=dim(X_list[[1]])[2]
frep2=sapply(seq(1,nmarkers),function(i) lapply(X_list,function(j) sum(j[,i]>0.75)))
founders=names(X_list)
fkeep=apply(frep2,MARGIN=2,function(x) x>2)
markers=dimnames(X_list[[1]])[[2]]
colnames(fkeep)=markers
colnames(frep2)=markers
fgroups=unique(colSums(fkeep))

n_reps=seq(1,reps)

randomized_gwas<-function(rep){
  all_gwas=data.frame(matrix(ncol=26,nrow=0))
  names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
  for(g in fgroups){
    subm=colnames(fkeep[,colSums(fkeep)==g])
    subfkeep=fkeep[,subm]
    X_list_sub=lapply(X_list,function(x) x[inter,subm])
    len=dim(X_list_sub[[1]])[1]

    # Run GridLMM

    # randomize the order of the genotypes
    draw=sample(len,len,replace=F)
    X_list_reordered=lapply(X_list_sub,function(x) x[draw,])
    for(x in seq(1,16)){
        dimnames(X_list_reordered[[x]])[[1]]=dimnames(X_list_sub[[1]])[[1]]
    }

    pattern=apply(subfkeep,MARGIN=2,function(x) str_flatten(c(unlist(founders[x])),'-'))
    #pattern=
    fdf=data.frame(marker=subm,fpattern=pattern,stringsAsFactors=F)
    fpatterns=unique(fdf$fpattern)
    if(g==16){
      h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
      names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
      h2_start
      V_setup=null_model$setup
      Y=as.matrix(data$y)
      X_cov=null_model$lmod$X

      X_list_null=NULL
      gwas=run_GridLMM_GWAS(Y,X_cov,X_list_reordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
      gwas$Trait=factor
      names(gwas)[6:21]=founders
      gwas=gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
      all_gwas=rbind(all_gwas,gwas)
    }else{
      for(i in fpatterns){
        subm2=fdf[fdf$fpattern==i,]$marker
        subf=subfkeep[,subm2,drop=F]
        #m=marker[i]
        fk=founders[subf[,1]]
        nfk=founders[!subf[,1]]
        X_list_sub2=X_list_reordered[ - which(names(X_list_reordered) %in% nfk)]
        X_list_sub2=lapply(X_list_sub2,function(x) x[,subm2,drop=F])

        h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
        names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
        h2_start
        V_setup=null_model$setup
        Y=as.matrix(data$y)
        X_cov=null_model$lmod$X

        X_list_null=NULL

        gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub2[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
        gwas$Trait=factor
        names(gwas)[6:(6+length(fk)-1)]=fk
        gwas[,nfk]=NA
        gwas=gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
        all_gwas=rbind(all_gwas,gwas)
      }
    }
  }
  all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
  tmp=data.frame(chr=chr,replicate=rep,pval=min(all_gwas$p_value_ML))
}




#fwrite(all_gwas,sprintf('eqtl/trans/results/%s_c%s_pheno_%s_trans_results.txt',time,chr,factor),row.names=F,quote=F,sep='\t')


print(system.time({
results=mclapply(n_reps,randomized_gwas,mc.cores=cores)
}))

saveRDS(results,sprintf('eqtl/trans/permute/chr%s_%s_%s_%.0frep_min_pvalues.rds',chr,time,factor,reps))







#fwrite(factor_results,sprintf('eqtl/trans/results/%s_c%s_pheno_factor_trans_eQTL.txt',time,chr),row.names=F,quote=F,sep='\t')
#fwrite(factor_results,sprintf('eqtl/trans/results/%s_c%s_factor_trans_eQTL.txt',time,chr),row.names=F,quote=F,sep='\t')

# For each snp, figure out fkeep and separate X_list into multiple based on dimensions
# parallelize


#n_reps=seq(1,reps)

#randomized_gwas<-function(rep){
#   len=dim(X_list_full[[1]])[1]

   # Run GridLMM

   # randomize the order of the genotypes
#   draw=sample(len,len,replace=F)
#   X_list_reordered=lapply(X_list_full,function(x) x[draw,])
#   for(x in seq(1,16)){
#       dimnames(X_list_reordered[[x]])[[1]]=dimnames(X_list_full[[1]])[[1]]
#   }

#   h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
#   names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
#   h2_start
#   V_setup=null_model$setup

#   Y=as.matrix(data$y)
#   X_cov=null_model$lmod$X
#   X_list_null=NULL

#   gwas=run_GridLMM_GWAS(Y,X_cov,X_list_reordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
#   gwas=gwas[!is.na(gwas$p_value_ML),]
#   tmp=data.frame(chr=chr,replicate=rep,pval=min(gwas$p_value_ML))
#}

#print(system.time({
#results=mclapply(n_reps,randomized_gwas,mc.cores=cores)
#}))

#saveRDS(results,sprintf('test_models/chr%s_%s_x_%s_founderprobs_%.0frep_max_pvalues.rds',chr,pheno,env,reps))

#date=format(Sys.time(),'%m%d%y')


#time="WD_0720"
#chr="10"
#library('GenomicFeatures') # write a script to get a table of start and stop sites of genes from the gtf file


#options(warn=2)
# Read in Kinship Matrix
