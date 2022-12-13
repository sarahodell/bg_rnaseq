#!/usr/bin/env Rscript
library('data.table')

time="WD_0712"

sigf=Sys.glob('eqtl/results/Factor[1-9]_pheno_trans_WD_0712_eQTL_hits.txt')
sigf=c(sigf,Sys.glob('eqtl/results/Factor[1-9][1-9]_pheno_trans_WD_0712_eQTL_hits.txt'))


sp1=sapply(seq(1,length(sigf)),function(x) strsplit(sigf[x],'/')[[1]][3])
sp2=sapply(seq(1,length(sigf)),function(x) strsplit(sp1[x],'_')[[1]][1])

df=data.frame(ind=seq(1,length(sigf)),factor=sp2,path=sigf,stringsAsFactors=F)
fwrite(df,'eqtl/pheno_trans_eQTL_list.txt',row.names=F,quote=F,sep=',')



sigf=Sys.glob('eqtl/results/Factor[1-9]_pheno_residuals_trans_WD_0712_eQTL_hits.txt')
sigf=c(sigf,Sys.glob('eqtl/results/Factor[1-9][1-9]_pheno_residuals_trans_WD_0712_eQTL_hits.txt'))


sp1=sapply(seq(1,length(sigf)),function(x) strsplit(sigf[x],'/')[[1]][3])
sp2=sapply(seq(1,length(sigf)),function(x) strsplit(sp1[x],'_')[[1]][1])

df=data.frame(ind=seq(1,length(sigf)),factor=sp2,path=sigf,stringsAsFactors=F)
fwrite(df,'eqtl/pheno_residuals_trans_eQTL_list.txt',row.names=F,quote=F,sep=',')
