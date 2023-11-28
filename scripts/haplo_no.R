#!/usr/env/bin Rscript

library('data.table')

axisdf=fread('eqtl/data/chromosome_axis.txt',data.table=F)


blocksper=c()
for(chr in 1:10){
	ibd=fread(sprintf('../ibd_segments/refinedibd/600K/bg%s_refined_ibd_blocks.txt',chr),data.table=F)
	gmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_gmap_c%s.csv',chr),data.table=F)
	nblocks=nrow(ibd)
	gsize=max(gmap$pos)
	line=data.frame(chr=chr,nblocks=nblocks,gsize=gsize,stringsAsFactors=F)
	blocksper=rbind(blocksper,line)
	
}

blocksper$per=blocksper$nblocks/blocksper$gsize

summary(blocksper$per)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  4.027   4.419   4.527   4.656   4.909   5.557 


phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
phenotypes=phenotypes[phenotypes$Loc.Year.Treat=="EXP_STPAUL_2017_WD",]

