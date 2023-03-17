#!/usr/bin/bash
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])

library('data.table')

allcis=c()
times=c('WD_0712','WD_0718','WD_0720','WD_0727')
for(t in times){
  df=fread(sprintf('eqtl/results/%s_cis_eQTL_fkeep_hits.txt',t),data.table=F)
  df$time=t
  allcis=rbind(allcis,df)
}

#betas=c()
#for(chr in 1:10){
#  beta=fread(sprintf('eqtl/results/eQTL_%s_c%.0f_betas.txt',time,chr))
#  betas=rbind(betas,beta)
#}
#info=betas[,c(1,2)]
#betas=betas[,-c(1,2)]
#beta0=unlist(betas[,1])
#betas2 = as.data.frame(unlist(sapply(seq(1,nrow(betas)),function(x) as.numeric(betas[x,-1])+as.numeric(betas[x,1]))),stringsAsFactors=F)
#betas2=rbind(beta0,betas2)
#betas2=t(betas2)
#betas=cbind(info,betas2)
#names(betas)=c('Gene_ID','SNP',paste0('beta.',seq(1,16)))
#fwrite(betas,sprintf('eqtl/results/eQTL_%s_all_betas.txt',time),row.names=F,quote=F,sep='\t')

all_res=c()
for(chr in 1:10){
  res=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_fkeep_residuals.txt',time,chr),data.table=F)
  res_var=apply(res[,-1],MARGIN=2,function(x) var(x))
  k=which(is.na(res_var))
  if(length(k)!=0){
    k=k+1
    k=unname(k)
    res=res[,-c(k)]
  }
  if(chr==1){
    all_res=res
  }else{
    all_res=cbind(all_res,res[,-1])
  }
}
fwrite(all_res,sprintf('eqtl/results/cis_eQTL_%s_all_fkeep_residuals.txt',time),row.names=F,quote=F,sep='\t')

#residuals
exp=fread(sprintf('eqtl/results/cis_eQTL_%s_all_fkeep_residuals.txt',time),data.table=F)
# full gene counts
#exp=fread(sprintf('eqtl/%s_voom_normalized_gene_counts.txt',time),data.table=F)
# need to t() and add genotype codes for this to work
#meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
#snames=colnames(exp)[-1]
#gnames=meta[match(snames,meta$sample_name),]$dh_genotype
#znames=exp$V1
#colnames(exp)=c('ID',gnames)
#rownames(exp)=exp$ID
#exp=exp[,-1]
#exp=as.matrix(exp)
#exp=t(exp)
#exp=as.data.frame(exp,stringsAsFactors=F)
#exp$ID=rownames(exp)
#exp=exp[,c('ID',znames)]

phenos=c("female_flowering_d6","male_flowering_d6","total_plant_height","harvest_grain_moisture",
"grain_yield_15","tkw_15",'asi')
envs=c("BLOIS_2014_OPT","BLOIS_2017_OPT","GRANEROS_2015_OPT","NERAC_2016_WD",
"STPAUL_2017_WD","SZEGED_2017_OPT")

phenotypes=fread('../GridLMM/phenotypes_asi.csv',data.table=F)
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes2=phenotypes[phenotypes$Genotype_code %in% exp$ID,]
phenotypes2=phenotypes2[,c('Genotype_code',"Loc.Year.Treat",phenos)]
newpheno=data.frame(ID=exp$ID,stringsAsFactors=F)
for(e in envs){
  df=phenotypes2[phenotypes2$Loc.Year.Treat==e,]
  names(df)=c('ID','Loc.Year.Treat',paste0(e,'-',phenos))
  df=df[match(newpheno$ID,df$ID),]
  rownames(df)=seq(1,nrow(df))
    newpheno=cbind(newpheno,df[3:9])
}
#newpheno_scaled=as.data.frame(sapply(seq(2,ncol(newpheno)),function(x) (newpheno[,x]-mean(newpheno[,x],na.rm=T))/sd(newpheno[,x],na.rm=T)),stringsAsFactors=F)
#names(newpheno_scaled)=names(newpheno)[-1]

exp=cbind(exp,newpheno[,-1])
#residuals
fwrite(exp,sprintf('eqtl/results/%s_fkeep_residuals_x_phenotypes.txt',time),row.names=F,quote=F,sep='\t')
#gene counts
#fwrite(exp,sprintf('eqtl/results/%s_genecounts_x_phenotypes.txt',time),row.names=F,quote=F,sep='\t')

#testing
#exp=exp[,1:100]


all_res=c()
for(chr in 1:10){
  res=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_vst_results.txt',time,chr))
  if(chr==1){
    all_res=res
  }else{
    all_res=rbind(all_res,res)
  }
}
fwrite(all_res,sprintf('eqtl/results/cis_eQTL_%s_all_vst_results.txt',time),row.names=F,quote=F,sep='\t')


#info=all_res[,c(1,2)]
#ID=all_res$ID
#genes0=unlist(genes[,1])
#genes2 = as.data.frame(unlist(sapply(seq(1,nrow(genes)),function(x) as.numeric(genes[x,-1])+as.numeric(genes[x,1]))),stringsAsFactors=F)
#betas2=rbind(beta0,betas2)
#betas2=t(betas2)
#betas=cbind(info,betas2)
