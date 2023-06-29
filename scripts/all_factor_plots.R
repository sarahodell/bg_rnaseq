#!/usr/bin/env Rscript

library('data.table')
library('qqman')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
df=c()
for(time in times){
	factors=fread(sprintf('eqtl/trans/%s_factors.txt',time),header=F,data.table=F)
	for(f in factors$V1){
		for(c in 1:10){
  			d=fread(sprintf('eqtl/trans/results/%s_c%.0f_%s_trans_results.txt',time,c,f))
 			pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  			d$time=time
  			d$factor=f
  			d$CHR=c
 			d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
 			d=d[,c('Trait','X_ID','p_value_ML','time','factor','CHR','BP')]
  			df=rbind(df,d)
		}
	}
}


df=df[!is.na(df$p_value_ML),]
df=df[order(df$p_value_ML),]
rownames(df)=seq(1,nrow(df))
df$p_adjusted=p.adjust(df$p_value_ML,method='fdr')
df$value=-log10(df$p_adjusted)


threshold=-log10(0.05)
print(threshold)

sig=df[df$value>=threshold,]
sig=sig[order(sig$CHR,sig$BP),]
fwrite(sig,'eqtl/results/all_factor_trans_eqtl_fdr_hits',row.names=F,quote=F,sep='\t')

png('eqtl/images/all_factor_fdr_qqplot.png')
print(qqman::qq(df$p_value_ML))
dev.off()

