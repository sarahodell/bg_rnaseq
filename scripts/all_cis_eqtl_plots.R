#!/usr/bin/env Rscript

library('data.table')
library('qqman')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
df=c()
for(time in times){
	for(c in 1:10){
  		d=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_weights_results.txt',time,c))
  		d$time=time
  		pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  		d$CHR=c
  		d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
  		d=d[,c('Trait','X_ID','p_value_ML','CHR','BP','time')]
  		df=rbind(df,d)
	}
}

df=df[!is.na(df$p_value_ML),]
df=df[order(df$p_value_ML),]
rownames(df)=seq(1,nrow(df))
df$p_adjusted=p.adjust(df$p_value_ML,method='fdr')
df$value=-log10(df$p_adjusted)

threshold=-log10(0.05)
sig=df[df$value>=threshold,]
sig=sig[order(sig$CHR,sig$BP),]

png('eqtl/images/all_cis_qqplot.png')
print(qqman::qq(df$p_value_ML))
dev.off()

fwrite(sig,'eqtl/results/all_cis_eQTL_weights_fdr_hits.txt',row.names=F,quote=F,sep='\t')
