#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
c1=as.numeric(args[[1]])
c2=as.numeric(args[[2]])
cores=as.numeric(args[[3]])

library('data.table')
library('parallel')
library('MASS')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra",
           "OH43_inra","A654_inra","FV2_inra","C103_inra",
           "EP1_inra","D105_inra","W117_inra","B96","DK63",
           "F492","ND245","VA85")
 
 
adj_chr=c(5,9)
if(c1 %in% adj_chr){
	Xlist1=readRDS(sprintf('phenotypes/bg%.0f_adjusted_genoprobs.rds',c1))

}else{
	Xlist1=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c1))

}
pmap1=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c1),data.table=F)
pmap1=pmap1[order(pmap1$pos),]
rownames(pmap1)=seq(1,nrow(pmap1))
pmap1=pmap1[pmap1$marker %in% dimnames(Xlist1[[1]])[[2]],]
Xlist1=lapply(Xlist1, function(x) x[,pmap1$marker])
size1=dim(Xlist1[[1]])[2]

if(c2 %in% adj_chr){
	Xlist2=readRDS(sprintf('phenotypes/bg%.0f_adjusted_genoprobs.rds',c2))

}else{
	Xlist2=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c2))

}
pmap2=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c2),data.table=F)
pmap2=pmap2[order(pmap2$pos),]
rownames(pmap2)=seq(1,nrow(pmap2))
pmap2=pmap2[pmap2$marker %in% dimnames(Xlist2[[1]])[[2]],]
Xlist2=lapply(Xlist2, function(x) x[,pmap2$marker])
size2=dim(Xlist2[[1]])[2]
		
m_names1=dimnames(Xlist1[[1]])[[2]]
m_names2=dimnames(Xlist2[[1]])[[2]]

mcombos = expand.grid(m1 = m_names1, m2 = m_names2)
       
get_r2=function(rep){
	snp1=mcombos[rep,]$m1
	snp2=mcombos[rep,]$m2
	pos1=pmap1[pmap1$marker==snp1,]$pos
	pos2=pmap2[pmap2$marker==snp2,]$pos
	v1=c()
	for(k in 1:16){v1=c(v1,as.vector(Xlist1[[k]][,snp1]))}
	v2=c()
	for(k in 1:16){v2=c(v2,as.vector(Xlist2[[k]][,snp2]))}
	r2=cor(v1,v2)**2
	line=data.frame(CHR_A=c1,BP_A=pos1,SNP_A=snp1,CHR_B=c2,BP_B=pos2,SNP_B=snp2,r2=r2,stringsAsFactors=F)
	return(line)
}
	

n_reps=1:nrow(mcombos)

print(system.time({
results=mclapply(n_reps,get_r2,mc.cores=cores)
}))

saveRDS(results,sprintf('eqtl/data/founder_ld_c%.0f_c%.0f.rds',c1,c2))


#c1=c()
#c2=c()
#for(i in 1:10){
#	for(j in i:10){
#		c1=c(c1,i)
#		c2=c(c2,j)
#	}
#}
#combos=data.frame(c1=c1,c2=c2)
#fwrite(combos,'eqtl/data/chr_combos.txt',row.names=F,quote=F,sep=',',col.names=F)


#get_ld=function(rep){
#	full_ld=c()  
#	c1=combos[rep,]$c1
#	c2=combos[rep,]$c2       
    
#	Xlist1=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c1))
#	pmap1=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c1),data.table=F)
#	pmap1=pmap1[order(pmap1$pos),]
#	rownames(pmap1)=seq(1,nrow(pmap1))
#	pmap1=pmap1[pmap1$marker %in% dimnames(Xlist1[[1]])[[2]],]
#	Xlist1=lapply(Xlist1, function(x) x[,pmap1$marker])
#	size1=dim(Xlist1[[1]])[2]

#	Xlist2=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c2))
#	pmap2=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c2),data.table=F)
#	pmap2=pmap2[order(pmap2$pos),]
#	rownames(pmap2)=seq(1,nrow(pmap2))
#	pmap2=pmap2[pmap2$marker %in% dimnames(Xlist2[[1]])[[2]],]
#	Xlist2=lapply(Xlist2, function(x) x[,pmap2$marker])
#	size2=dim(Xlist2[[1]])[2]
		
#	m_names1=dimnames(Xlist1[[1]])[[2]]
#	m_names2=dimnames(Xlist2[[1]])[[2]]
#	for(i in 1:size1){
#		for(j in 1:size2){
#			snp1=m_names1[i]
#			pos1=pmap1[pmap1$marker==snp1,]$pos
#			snp2=m_names2[j]
#			pos2=pmap2[pmap2$marker==snp2,]$pos
#			v1=c()
#			for(k in 1:16){v1=c(v1,as.vector(Xlist1[[k]][,i]))}
#			v2=c()
#			for(k in 1:16){v2=c(v2,as.vector(Xlist2[[k]][,j]))}
#			r2=cor(v1,v2)**2
#			line=data.frame(CHR_A=c1,SNP_A=snp1,BP_A=pos1,CHR_B=c2,SNP_B=snp2,BP_B=pos2,r2=r2,stringsAsFactors=F)
#			full_ld=rbind(full_ld,line)
#		}
# 	}
#  	full_ld=as.data.frame(full_ld,stringsAsFactors=F)
#	return(full_ld)
#}



#full_ld=c()          
#
#for(c1 in 1:10){
#	Xlist1=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c1))
#	pmap1=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c1),data.table=F)
#	pmap1=pmap1[order(pmap1$pos),]
#	rownames(pmap1)=seq(1,nrow(pmap1))
#	pmap1=pmap1[pmap1$marker %in% dimnames(Xlist1[[1]])[[2]],]
#	Xlist1=lapply(Xlist1, function(x) x[,pmap1$marker])
#	size1=dim(Xlist1[[1]])[2]
#	for(c2 in c1:10){
#		Xlist2=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c2))
#		pmap2=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c2),data.table=F)
#		pmap2=pmap2[order(pmap2$pos),]
#		rownames(pmap2)=seq(1,nrow(pmap2))
#		pmap2=pmap2[pmap2$marker %in% dimnames(Xlist2[[1]])[[2]],]
#		Xlist2=lapply(Xlist2, function(x) x[,pmap2$marker])
#		size2=dim(Xlist2[[1]])[2]
#		
#		
#		m_names1=dimnames(Xlist1[[1]])[[2]]
#		m_names2=dimnames(Xlist2[[1]])[[2]]
#
#		for(i in 1:size1){
#			for(j in 1:size2){
#				snp1=m_names1[i]
#				pos1=pmap1[pmap1$marker==snp1,]$pos
#				snp2=m_names2[j]
#				pos2=pmap2[pmap2$marker==snp2,]$pos
#				v1=c()
#				for(k in 1:16){v1=c(v1,as.vector(Xlist1[[k]][,i]))}
#				v2=c()
#				for(k in 1:16){v2=c(v2,as.vector(Xlist2[[k]][,j]))}
#				
#				r2=cor(v1,v2)**2
#				line=data.frame(CHR_A=c1,SNP_A=snp1,BP_A=pos1,CHR_B=c2,SNP_B=snp2,BP_B=pos2,r2=r2,stringsAsFactors=F)
#				full_ld=rbind(full_ld,line)
#			}
#  		}
#	}
#}

#fwrite(full_ld,'founder_ld.txt',row.names=F,quote=F,sep='\t')
