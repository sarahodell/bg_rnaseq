#!/usr/bin/env Rscript

library('data.table')

fulldf=c()
for(c1 in 1:10){
	print(c1)
	for(c2 in c1:10){
		print(c2)
		ld=readRDS(sprintf('eqtl/data/founder_ld_c%.0f_c%.0f.rds',c1,c2))
		df=rbindlist(ld)
		df=as.data.frame(df)
		fulldf=rbind(fulldf,df)
	}
}
rm(df)
fulldf=as.data.frame(fulldf,stringsAsFactors=F)
fwrite(fulldf,'eqtl/data/all_founder_ld.txt',row.names=F,quote=F,sep='\t')

interld=fulldf[fulldf$CHR_A!=fulldf$CHR_B & fulldf$r2>0.5,]

fwrite(interld,'eqtl/data/founder_interchrom_ld.txt',row.names=F,quote=F,sep='\t')