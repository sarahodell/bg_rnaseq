#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
chr=as.character(args[[2]])
cores=as.numeric(args[[3]])
reps=as.numeric(args[[4]])

library('GridLMM')
library('data.table')
library('dplyr')
library('parallel')
library('MASS')

# Read in Kinship Matrix
K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

genetable=fread('eqtl/data/Zea_mays.B73_v4_generanges.txt',data.table=F)
genetable=genetable[genetable$TXCHROM==chr,]
genes=unique(genetable$Gene_ID)
# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts.txt',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time),data.table=F)
samples=colnames(phenotypes)[-1]
genos=metadata[match(samples,metadata$sample_name),]$genotype
testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

genes=intersect(genes,phenotypes$V1)

X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])
dhs=metadata[match(samples,metadata$sample_name),]$dh_genotype
i=intersect(dhs,inds)
K=K[i,i]


#remove(X_list)
by_gene<-function(g){
  gene=genes[g]
  snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
  if(length(snp)>1){
      snp=snp[2]
  }
  data=data.frame(sample=colnames(phenotypes)[-1],y=unlist(phenotypes[phenotypes[,1]==gene,][-1]),stringsAsFactors=F)
  data=data[!is.na(data$y),]
  data$ID=metadata[match(data$sample,metadata$sample_name),]$dh_genotype
  rownames(data)=seq(1,nrow(data))
  data=data[!is.na(data$ID),]
  data$PC1=pcs[match(pcs$sample,data$sample),]$PC1
  data$PC2=pcs[match(pcs$sample,data$sample),]$PC2
  data$PC3=pcs[match(pcs$sample,data$sample),]$PC3

  if(length(unique(data$ID))<nrow(data)){
    data=data%>%group_by(ID)%>%summarize(y=mean(y),PC1=mean(PC1),PC2=mean(PC2),PC3=mean(PC3))
  }
  data=as.data.frame(data,stringsAsFactors=F)
  rownames(data)=data$ID
  data=data[i,]
  data$ID2=data$ID
  # variance stabilize
  data$y=(data$y-mean(data$y))/sd(data$y)
  data=data[,c('ID','ID2','y','PC1','PC2','PC3')]
  null_model = GridLMM_ML(y~1+PC1+PC2+PC3+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

  h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
  names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
  h2_start
  V_setup=null_model$setup

  Y=as.matrix(data$y)
  X_cov=null_model$lmod$X
  X_list_ordered=lapply(X_list,function(x) x[i,snp,drop=F])
  X_list_null=NULL

  randomized_gwas<-function(rep){
     len=dim(X_list_ordered[[1]])[1]

     # Run GridLMM
     # randomize the order of the genotypes
     draw=sample(len,len,replace=F)
     X_list_reordered=lapply(X_list_ordered,function(x) x[draw,,drop=F])
     for(x in seq(1,16)){
         dimnames(X_list_reordered[[x]])[[1]]=dimnames(X_list_ordered[[1]])[[1]]
     }

     h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
     names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
     h2_start
     V_setup=null_model$setup

     Y=as.matrix(data$y)
     X_cov=null_model$lmod$X
     X_list_null=NULL

     #gwas=run_GridLMM_GWAS(Y,X_cov,X_list_reordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
     gwas=run_GridLMM_GWAS(Y,X_cov,X_list_reordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
     gwas$Trait=gene
     gwas=gwas[!is.na(gwas$p_value_ML),]
     gwas=gwas[!is.infinite(gwas$p_value_ML),]
     tmp=data.frame(chr=chr,gene=gene,snp=snp,replicate=rep,pval=min(gwas$p_value_ML))
  }


  n_reps=seq(1,reps)
  results=mclapply(n_reps,randomized_gwas,mc.cores=cores)

  #results$gene=gene
  #results$snp=snp
  #tmp=data.frame(chr=chr,gene=gene,snp=snp,replicate=rep,pval=min(gwas$p_value_ML))
}



n_reps=seq(1,reps)
#n_genes=seq(1,5)
n_genes=seq(1,length(genes))

print(system.time({
fullresults=mclapply(n_genes,by_gene,mc.cores=cores)
}))

resultsdf=data.frame(stringsAsFactors=F)
for(i in n_genes){
  for(j in n_reps){
    resultsdf=rbind(resultsdf,fullresults[[i]][[j]])
  }
}

fwrite(resultsdf,sprintf('eqtl/permute/%s_chr%s_ciseQTL_%.0frep_max_pvalues.txt',time,chr,reps),row.names=F,quote=F,sep='\t')
#saveRDS(fullresults,sprintf('eqtl/permute/test_models/chr%s_ciseQTL_%.0frep_max_pvalues.rds',chr,reps))



#for(g in 1:length(genes)){
#  gene=genes[g]
#  snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
#  data=data.frame(sample=colnames(phenotypes)[-1],y=unlist(phenotypes[phenotypes[,1]==gene,][-1]),stringsAsFactors=F)
#  data=data[!is.na(data$y),]
#  data$ID=metadata[match(data$sample,metadata$sample_name),]$dh_genotype
#  rownames(data)=seq(1,nrow(data))
#  data=data[!is.na(data$ID),]
#  data$PC1=pcs[match(pcs$sample,data$sample),]$PC1
#  data$PC2=pcs[match(pcs$sample,data$sample),]$PC2
#  data$PC3=pcs[match(pcs$sample,data$sample),]$PC3
#  K=K[i,i]

#  if(length(unique(data$ID))<nrow(data)){
#    data=data%>%group_by(ID)%>%summarize(y=mean(y),PC1=mean(PC1),PC2=mean(PC2),PC3=mean(PC3))
#  }
#  data=as.data.frame(data,stringsAsFactors=F)
#  rownames(data)=data$ID
#  data=data[i,]
#  data$ID2=data$ID
  # variance stabilize
#  data$y=(data$y-mean(data$y))/sd(data$y)
#  data=data[,c('ID','ID2','y','PC1','PC2','PC3')]
#  null_model = GridLMM_ML(y~1+PC1+PC2+PC3+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)
#
#  h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]3
#  names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
#  h2_start
#  V_setup=null_model$setup

#  Y=as.matrix(data$y)
#  X_cov=null_model$lmod$X

#  if(length(snp)>1){
#      X = do.call(cbind,lapply(X_list,function(x) x[,snp[2]]))
#      colnames(X) = founders
#      rownames(X) = dimnames(X_list[[1]])[[1]]
#      snp=snp[2]
#      X=X[i,]
#  }else{
#    X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
#    colnames(X) = founders
#    rownames(X) = dimnames(X_list[[1]])[[1]]
#    X=X[i,]
#  }
#    X_list_ordered=lapply(X_list,function(x) x[i,snp,drop=F])
#    frep2=apply(X,MARGIN=2,function(x) sum(x[x>0.8]))
#    fkeep=founders[frep2>2]
#    X_list_ordered=X_list_ordered[c(fkeep)]
#    X=X[,fkeep]

#    X_list_null=NULL

#    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
#    gwas$Trait=gene

#    all_res[,gene]=residuals[,gene]
#    all_gwas[g,]=new_gwas
#  }
#}#

#all_res=as.data.frame(all_res,stringsAsFactors=F)
#genenames=names(all_res)
#all_res$ID=i

#all_res=all_res[,c('ID',genenames)]
#all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)

#fwrite(all_gwas,sprintf('eqtl/cis/results/eQTL_%s_c%s_vst_results.txt',time,chr),row.names=F,quote=F,sep='\t')
