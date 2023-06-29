#!/usr/bin/env Rscript 

#library('LDheatmap',lib='/home/sodell/R/x86_64-conda-linux-gnu/4.2')
library('data.table')
library('ggplot2')

ld=fread('eqtl/data/all_founder_ld.txt',data.table=F)

for(c1 in 1:10){
	for(c2 in c1:10){
		subld=ld[ld$CHR_A==c1 & ld$CHR_B==c2,]
		pmap1=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c1),data.table=F)
		pmap1=pmap1[order(pmap1$pos),]
		rownames(pmap1)=seq(1,nrow(pmap1))
		m1=unique(subld$SNP_A)
		pmap1=pmap1[pmap1$marker %in% m1,]
		pmap1=pmap1[order(pmap1$pos),]
		rownames(pmap1)=seq(1,nrow(pmap1))
		m1=pmap1$marker
		pmap2=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c2),data.table=F)
		pmap2=pmap2[order(pmap2$pos),]
		rownames(pmap2)=seq(1,nrow(pmap2))
		m2=unique(subld$SNP_B)
		pmap2=pmap2[pmap2$marker %in% m2,]
		pmap2=pmap2[order(pmap2$pos),]
		rownames(pmap2)=seq(1,nrow(pmap2))
		m2=pmap2$marker
		
		df=matrix(ncol=length(m2),nrow=length(m1))
		rownames(df)=m1
		colnames(df)=m2
		for(m in m1){
			sub2=subld[subld$SNP_A==m,]
			r2=sub2[match(m2,sub2$SNP_B),]$r2
			df[m,]=r2
		}
		fwrite(as.data.frame(df),sprintf('eqtl/data/founder_ld_c%.0f_c%.0f_matrix.txt',c1,c2),row.names=T,quote=F,col.names=T,sep=',')
		
		
		subld$SNP_A_f=factor(subld$SNP_A,levels=c(m1))
		subld$SNP_B_f=factor(subld$SNP_B,levels=c(m2))
		
		h1=ggplot(subld,aes(x=SNP_A_f,y=SNP_B_f,fill=r2)) + geom_tile() +
		theme_classic() + scale_fill_gradient(low="white", high="blue") +
		theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()) +
        xlab(sprintf("Chromosome %.0f",c1)) + ylab(sprintf("Chromosome %.0f",c2)) + 
        ggtitle("Founder LD")
        
		pdf(sprintf('eqtl/data/c%.0f_c%.0f_founderld_heatmap.pdf',c1,c2))
		print(h1)
		dev.off()
	}
}


allsnps=c()
for(c in 1:10){
	Xlist=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c))
	allsnps=c(allsnps,dimnames(Xlist[[1]])[[2]])

}
snps=data.frame(snps=allsnps,stringsAsFactors=F)
fwrite(snps,'eqtl/data/snplist.txt',row.names=F,quote=F,col.names=F,sep='\t')


snplist=fread('eqtl/data/snplist.txt',header=F)

### SNP LD
snpld=fread('eqtl/data/SNP_LD.ld',data.table=F)

for(c1 in 1:10){
	for(c2 in c1:10){
		subld=snpld[snpld$CHR_A==c1 & snpld$CHR_B==c2,]
		pmap1=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c1),data.table=F)
		pmap1=pmap1[order(pmap1$pos),]
		rownames(pmap1)=seq(1,nrow(pmap1))
		#Xlist1=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c1))
		#m1=dimnames(Xlist1[[1]])[[2]]
		#m1=unique(subld$SNP_A)
		#pmap1=pmap1[pmap1$marker %in% m1,]
		pmap1=pmap1[order(pmap1$pos),]
		rownames(pmap1)=seq(1,nrow(pmap1))
		#m1=pmap1$marker
		m1=intersect(pmap1$marker,subld$SNP_A)
		
		pmap2=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c2),data.table=F)
		pmap2=pmap2[order(pmap2$pos),]
		rownames(pmap2)=seq(1,nrow(pmap2))
		#Xlist2=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c2))
		#m2=dimnames(Xlist2[[1]])[[2]]
		#m2=unique(subld$SNP_B)
		#pmap2=pmap2[pmap2$marker %in% m2,]
		pmap2=pmap2[order(pmap2$pos),]
		rownames(pmap2)=seq(1,nrow(pmap2))
		#m2=pmap2$marker
		m2=intersect(pmap2$marker,subld$SNP_B)

		#subld=subld[subld$SNP_A %in% m1 & subld$SNP_B %in% m2,]

		#df=matrix(ncol=length(m2),nrow=length(m1))
		#rownames(df)=m1
		#colnames(df)=m2
		#for(m in m1){
		#	sub2=subld[subld$SNP_A==m,]
		#	r2=sub2[match(m2,sub2$SNP_B),]$R2
		#	df[m,]=r2
		#}
		#fwrite(as.data.frame(df),sprintf('eqtl/data/SNP_ld_c%.0f_c%.0f_matrix.txt',c1,c2),row.names=T,quote=F,col.names=T,sep=',')
		
		#allmarkers=union(m1,m2)
		subld$SNP_A_f=factor(subld$SNP_A,levels=c(m1))
		subld$SNP_B_f=factor(subld$SNP_B,levels=c(m2))
		
		h1=ggplot(subld,aes(x=SNP_A_f,y=SNP_B_f,fill=R2)) + geom_tile() +
		theme_classic() + scale_fill_gradient(low="white", high="blue") +
		theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()) +
        xlab(sprintf("Chromosome %.0f",c1)) + ylab(sprintf("Chromosome %.0f",c2)) + 
        ggtitle("SNP LD")
        
		png(sprintf('eqtl/data/c%.0f_c%.0f_SNP_heatmap.png',c1,c2),height=1000,width=600)
		print(h1)
		dev.off()
	}
}
