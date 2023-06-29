#!/usr/bin/env Rscript


library('qqman')
library('data.table')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")

df=c()
for(time in times){
	print(time)
	for(c in 1:10){
		print(c)
		dlist=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%.0f_weights_results.rds',time,c))
		#d=as.data.frame(do.call(rbind, dlist))
		d=rbindlist(dlist)
		d=as.data.frame(d)
		d$time=time
		rm(dlist)
		pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
		d$CHR=c
		d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
		rm(pmap)
		d=d[,c('Trait','X_ID','p_value_ML','time','CHR','BP')]
		df=rbind(df,d)
		rm(d)
	}
}



df=df[order(df$p_value_ML),]
rownames(df)=seq(1,nrow(df))
df$p_adjusted=p.adjust(df$p_value_ML,method='fdr')
df$value=-log10(df$p_adjusted)

threshold=-log10(0.05)
sig=df[df$value>=threshold,]

#png('eqtl/images/all_trans_fdr_qqplot.png')
#print(qqman::qq(df$p_value_ML))
#dev.off()

fwrite(sig,'eqtl/results/all_trans_fdr_hits.txt',row.names=F,quote=F,sep='\t')

