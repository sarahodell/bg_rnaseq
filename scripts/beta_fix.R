#!/usr/bin/bash
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])

library('data.table')

betas=c()
for(chr in 1:10){
  beta=fread(sprintf('eqtl/results/eQTL_%s_c%.0f_betas.txt',time,chr))
  betas=rbind(betas,beta)
}
info=betas[,c(1,2)]
betas=betas[,-c(1,2)]
beta0=unlist(betas[,1])
betas2 = as.data.frame(unlist(sapply(seq(1,nrow(betas)),function(x) as.numeric(betas[x,-1])+as.numeric(betas[x,1]))),stringsAsFactors=F)
betas2=rbind(beta0,betas2)
betas2=t(betas2)
betas=cbind(info,betas2)
names(betas)=c('Gene_ID','SNP',paste0('beta.',seq(1,16)))
fwrite(betas,sprintf('eqtl/results/eQTL_%s_all_betas.txt',time),row.names=F,quote=F,sep='\t')

all_res=c()
for(chr in 1:10){
  res=fread(sprintf('eqtl/results/eQTL_%s_c%.0f_residuals.txt',time,chr))
  if(chr==1){
    all_res=res
  }else{
    all_res=cbind(all_res,res[,2:ncol(res)])
  }
}
fwrite(all_res,sprintf('eqtl/results/cis_eQTL_%s_all_residuals.txt',time),row.names=F,quote=F,sep='\t')

#info=all_res[,c(1,2)]
#ID=all_res$ID
#genes0=unlist(genes[,1])
#genes2 = as.data.frame(unlist(sapply(seq(1,nrow(genes)),function(x) as.numeric(genes[x,-1])+as.numeric(genes[x,1]))),stringsAsFactors=F)
#betas2=rbind(beta0,betas2)
#betas2=t(betas2)
#betas=cbind(info,betas2)
