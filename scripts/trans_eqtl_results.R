#!/usr/bin/env Rscript

library('data.table')

#factors=paste0('Factor',seq(1,100))
time="WD_0720"

#files=Sys.glob(sprintf('eqtl/results/Factor*_%s_eQTL_hits.txt',time))
#all_hits=c()
#for(f in files){
#  fdf=fread(f,data.table=F)
#  all_hits=rbind(all_hits,fdf)
#}

#fwrite(all_hits,sprintf('eqtl/results/%s_all_factor_trans_eQTL_hits.txt',time),row.names=F,quote=F,sep='\t')
all_hits=fread(sprintf('eqtl/results/%s_all_factor_trans_eQTL_hits.txt',time),data.table=F)

#factors=unique(all_hits$Factor)

lambda_all_means=fread(sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means.txt',time),data.table=F)
rownames(lambda_all_means)=lambda_all_means$V1
lambda_all_means=lambda_all_means[,-1]
#lambda_all_means=lambda_all_means[,factors]

factors=names(lambda_all_means)
factor_groups=vector("list",length=ncol(lambda_all_means))
for(i in 1:length(factors)){
  factor_groups[[i]]$factor=factors[i]
  factor_groups[[i]]$genes=c(NA)
}
#v=apply(lambda_all_means,MARGIN=1,function(x) x**2)
#v2=apply(v,MARGIN=2,function(x) which(x>0.2))
#prop_var<-function(f){
#  subl=lambda_all_means[f,,drop=F]
#  gene=rownames(subl)
#  var_exp=apply(subl,MARGIN=1,function(x) x**2)
#  tot_var=sum(var_exp)
#  prop_var=var_exp/tot_var
#  fkeep=names(subl[,which(prop_var>=0.2),drop=F])
#}

prop_var=sapply(seq(1,nrow(lambda_all_means)),function(x) apply(lambda_all_means[x,,drop=F],MARGIN=1,function(i) i**2)/sum(apply(lambda_all_means[x,,drop=F],MARGIN=1,function(j) j**2)))
colnames(prop_var)=rownames(lambda_all_means)
rownames(prop_var)=names(lambda_all_means)
fkeep=sapply(seq(1,ncol(prop_var)),function(x) names(which(prop_var[,x]>=0.2)))
genes=colnames(prop_var)
for(k in 1:length(fkeep)){
  fgroups=fkeep[k]
  gene=genes[k]
  for(l in fgroups[[1]]){
    x=which(unlist(unname(lapply(factor_groups,function(x) x$factor==l))))
    if(is.na(factor_groups[[x]]$genes[1])){
      factor_groups[[x]]$genes=gene
    }else{
      factor_groups[[x]]$genes=c(factor_groups[[x]]$genes,gene)
    }
  }
}

for(i in 1:length(factor_groups)){
  print(length(factor_groups[[i]]$genes))
}
saveRDS(factor_groups,sprintf('eqtl/results/%s_factor_groupings.rds',time))
#54307147_1
#54307054_1

#for(f in 1:nrow(lambda_all_means)){
  #index=f-1
#  subl=lambda_all_means[f,,drop=F]
#  gene=rownames(subl)
#  var_exp=apply(subl,MARGIN=1,function(x) x**2)
#  tot_var=sum(var_exp)
#  prop_var=var_exp/tot_var
#  fkeep=names(subl[,which(prop_var>=0.2),drop=F])
#  for(k in fkeep){
#    x=which(unlist(unname(lapply(factor_groups,function(x) x$factor==k))))
#    if(is.na(factor_groups[[x]]$genes[1])){
#      factor_groups[[x]]$genes=gene
#    }else{
#      factor_groups[[x]]$genes=c(factor_groups[[x]]$genes,gene)
#    }
#  }
  #print(length(fkeep))
#}

#library('clusterProfiler')
