#!/usr/bin/env Rscript

library('qqman')
library('data.table')
library('qvalue')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")

df=c()
for(time in times){
	print(time)
	for(c in 1:10){
		print(c)
		d=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%.0f_weights_results_filtered_FIXED.txt',time,c),data.table=F)
		#dlist=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.rds',time,c))
		#d=rbindlist(dlist)
		#d=as.data.frame(d)
		d$time=time
		#rm(dlist)
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
qvalues=qvalue(df$p_value_ML,fdr.level=0.1)
summary(qvalues)
#Call:
#qvalue(p = df$p_value_ML, fdr.level = 0.1)

#pi0:	0.6890435	

#Cumulative number of significant calls:

#          <1e-04  <0.001   <0.01   <0.025    <0.05     <0.1        <1
#p-value   397041 1522578 8368944 17306335 30225060 52977748 324784652
#q-value    70587   94841  157090   228184   358672   767376 324784652
#local FDR  65877   66638   66749    66754    66757    66762     68848
write.qvalue(qvalues,'eqtl/results/trans_qvalues.txt',sep='\t')


#pdf('eqtl/trans/images/trans_qvalue_diagnostics.pdf')
#hist(qvalues)
#plot(qvalues)
#dev.off()

df$p_adjusted=qvalues$qvalues
#df$lfdr=qvalues$lfdr


#df$p_adjusted=p.adjust(df$p_value_ML,method='fdr')
df$value=-log10(df$p_adjusted)

threshold=-log10(0.1)
sig=df[df$value>=threshold,]

#png('eqtl/images/all_trans_fdr_qqplot.png')
#print(qqman::qq(df$p_value_ML))
#dev.off()

fwrite(sig,'eqtl/results/all_trans_fdr_hits_FIXED.txt',row.names=F,quote=F,sep='\t')

