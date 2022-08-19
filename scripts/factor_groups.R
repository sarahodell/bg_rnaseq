#!/usr/bin/env Rscript

library('data.table')

time="WD_0712"

data=fread('pheno_MegaLMM_vst_residuals_WD_0712_factor_correlations.txt',data.table=F)

lambda_all_means=fread('pheno_MegaLMM_vst_residuals_WD_0712_all_Lambda_means.txt',data.table=F)
rownames(lambda_all_means)=lambda_all_means$V1
lambda_all_means=lambda_all_means[,-1]
phenotypes=rownames(lambda_all_means)[!grepl('Zm',rownames(lambda_all_means))]

count=c()
for(i in 1:nrow(data)){
  count=c(count,sum(abs(data[i,c(4:6)])>=0.75,na.rm=T))
}
strong=which(count==3)
sdata=data[strong,]

factors=names(lambda_all_means)
factor_groups=vector("list",length=nrow(lambda_all_means))
for(i in 1:length(factors)){
  factor_groups[[i]]$factor=factors[i]
  factor_groups[[i]]$genes=c(NA)
}

for(f in 1:nrow(lambda_all_means)){
  subl=lambda_all_means[f,,drop=F]
  gene=rownames(subl)
  var_exp=apply(subl,MARGIN=1,function(x) x**2)
  tot_var=sum(var_exp)
  prop_var=var_exp/tot_var
  fkeep=names(subl[,which(prop_var>=0.2),drop=F])
  for(k in fkeep){
    x=unlist(unname(lapply(factor_groups,function(x) which(x$factor==k))))
    factor_groups[[x]]$genes=c(factor_groups[[x]]$genes,gene)
  }
}

saveRDS(factor_groups,'pheno_MegaLMM_vst_residuals_WD_0712_factor_groups.rds')
