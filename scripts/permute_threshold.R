#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])

library('data.table')
library('dplyr')
library('ggplot2')

options(scipen=20)

all_reps=c()
for(c in 1:10){
   r=fread(sprintf('eqtl/permute/%s_chr%s_ciseQTL_500rep_max_pvalues.txt',time,c),data.table=F)
   all_reps=rbind(all_reps,r)
}
#test one chromosome
c=10
all_reps=fread(sprintf('eqtl/permute/%s_chr%s_ciseQTL_500rep_max_pvalues.txt',time,c),data.table=F)

thresh=0.05
minp = all_reps %>% group_by(replicate) %>% summarize(pval=min(pval))
minp=as.data.frame(minp)

threshold=quantile(minp$pval,thresh,lower.tail=T)
print(threshold)
print(-log10(threshold))
# 4.480104 for chromosome 10, probably worth doing

line=data.table(time=time,threshold=-log10(threshold))
fwrite(line,file=sprintf("ciseqtl_threshold_%.2f_table.txt",thresh),sep=',',append=T)


#png(sprintf('GridLMM_founderprobs/permute/images/%s_x_%s_perm_1000_pval_dist.png',pheno,env))
#print(ggplot(minp,aes(x=pval)) + geom_histogram() + geom_vline(xintercept=threshold))
#dev.off()

png(sprintf('eqtl/permute/%s_perm_500_log10pval_dist.png',time))
print(ggplot(minp,aes(x=-log10(pval))) + geom_histogram() + geom_vline(xintercept=-log10(threshold)))
dev.off()
