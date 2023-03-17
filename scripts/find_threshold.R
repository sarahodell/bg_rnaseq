#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
factor=as.character(args[[2]])
thresh=as.numeric(args[[3]])

library('data.table')
library('dplyr')
library('ggplot2')

options(scipen=20)

all_reps=c()
for(c in 1:10){
   r=readRDS(sprintf('eqtl/trans/permute/chr%.0f_%s_%s_1000rep_min_pvalues.rds',c,time,factor))
   df=c()
   df=sapply(seq(1,1000),function(x) rbind(df,unlist(r[[x]])))
   df=t(df)
   df=as.data.frame(df)
   names(df)=c('chr','replicate','pval')
   df=df[!is.na(df$pval),]
   all_reps=rbind(all_reps,df)
}

minp = all_reps %>% group_by(replicate) %>% summarize(pval=min(pval))
minp=as.data.frame(minp)

threshold=quantile(minp$pval,thresh,lower.tail=T)
print(threshold)
print(-log10(threshold))

method="founder_probs"

line=data.table(time=time,factor=factor,method=method,threshold=-log10(threshold))
fwrite(line,file=sprintf("eqtl/trans/threshold_%.2f_table.txt",thresh),sep=',',append=T)


#png(sprintf('GridLMM_founderprobs/permute/images/%s_x_%s_perm_1000_pval_dist.png',pheno,env))
#print(ggplot(minp,aes(x=pval)) + geom_histogram() + geom_vline(xintercept=threshold))
#dev.off()

#png(sprintf('eqtl/trans/permute/images/%s_x_%s_%.2f_perm_1000_log10pval_dist.png',pheno,env,thresh))
#print(ggplot(minp,aes(x=-log10(pval))) + geom_histogram() + geom_vline(xintercept=-log10(threshold)))
#dev.off()
