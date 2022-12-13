#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
factor=as.character(args[[2]])

library('data.table')

#qtl=fread(sprintf('eqtl/results/%s_pheno_trans_%s_eQTL_hits.txt',factor,time),data.table=F)
qtl=fread(sprintf('eqtl/results/%s_pheno_residuals_trans_%s_eQTL_hits.txt',factor,time),data.table=F)

#snp="AX-9180836"
gmap=c()
for(c in 1:10){
  g=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_gmap_c%.0f.csv',c),data.table=F)
  gmap=rbind(gmap,g)
}

test_snps=c()
#threshold=0.9
for(c in 1:10){
  sub=sub=qtl[qtl$CHR==c,]
  start_ind=which.max(sub$value)
  row=sub[start_ind,]
  snp=row$SNP
  sub=sub[sub$SNP!=snp,]
  snp_cM=gmap[gmap$marker==snp,]$pos
  test_snps=c(test_snps,snp)
  while(nrow(sub)!=0){
    sub$cM=gmap[match(sub$SNP,gmap$marker),]$pos
    sub=sub[abs(snp_cM-sub$cM)>10,]
    start_ind=which.max(sub$value)
    row=sub[start_ind,]
    snp=row$SNP
    test_snps=c(test_snps,snp)
    sub=sub[sub$SNP!=snp,]
    snp_cM=gmap[gmap$marker==snp,]$pos
  }
}

kept=qtl[qtl$SNP %in% test_snps,]
kept$cM=gmap[match(kept$SNP,gmap$marker),]$pos
#fwrite(kept,sprintf('allelic/%s_%s_pheno_trans_kept_SNPS.txt',time,factor),row.names=F,quote=F,sep='\t')
fwrite(kept,sprintf('allelic/%s_%s_pheno_residuals_trans_kept_SNPS.txt',time,factor),row.names=F,quote=F,sep='\t')
