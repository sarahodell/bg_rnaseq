#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')

calc_maf=function(x){
	maf=mean(x)
	if(maf>0.5){
		maf=1-maf
	}
	return(maf)
}

allmafs=c()
for(c in 1:10){
	snp=fread(sprintf('../genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%.0f_binary.csv',c),data.table=F)
	rownames(snp)=snp$ind
	psnp=fread(sprintf('../genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f.csv',c),data.table=F)
	rownames(psnp)=psnp$ind
	
	fp=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c))
	#markers=dimnames(fp[[1]])[[2]]
	#psnp=psnp[,markers]
	psnp=psnp[,-1]
	psnp=apply(psnp,MARGIN=2,function(x) ifelse(x=="A",0,1))
	inds=dimnames(fp[[1]])[[1]]
	#snp=snp[inds,markers]
	snp=snp[inds,]
	snp=snp[,-1]
	markers=intersect(colnames(snp),colnames(psnp))
	snp=snp[,markers]
	psnp=psnp[,markers]
	smaf=apply(snp,MARGIN=2,calc_maf)
	pmaf=apply(psnp,MARGIN=2,calc_maf)
	df=data.frame(chr=c,marker=markers,magicmaf=smaf,foundermaf=pmaf,stringsAsFactors=F)
	allmafs=rbind(allmafs,df)
}

fwrite(allmafs,'metadata/DH_Founder_allele_frequencies_full.txt',row.names=F,quote=F,sep='\t')

p1=ggplot(allmafs,aes(x=magicmaf)) + geom_histogram(bins=10) +
theme_classic() + xlab("MAF in MAGIC DH lines") + ylab("Frequency") +
scale_x_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
png('images/dh_afs.png')
print(p1)
dev.off()

allmafs$foundermaf_f=factor(allmafs$foundermaf,levels=c(sort(unique(allmafs$foundermaf))))


p2=ggplot(allmafs,aes(x=magicmaf)) + geom_histogram(bins=10) +
  facet_grid(foundermaf_f ~ .) + 
theme_classic() + xlab("MAF in MAGIC DH lines by Founder MAF") + ylab("Frequency") +
scale_x_continuous(limits = c(-0.01, 0.51), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))

png('images/dh_afs_by_foundermaf.png')
print(p2)
dev.off()

lmaf=allmafs[allmafs$magicmaf<0.05,]
lmaf %>% group_by(chr) %>% count

# I should refilter genotypes to drop SNPs with low MAF - does this impact founder probs at all? Probably nota
allmafs=fread('metadata/DH_Founder_allele_frequencies_full.txt',data.table=F)
kept=allmafs[allmafs$magicmaf>0.0123,]
kept=kept[kept$foundermaf>0,]
#483861 markers

get_r2=function(m1,m2){
	v1=c()
	for(k in 1:16){v1=c(v1,as.vector(fp[[1]][,k,m1]))}
	v2=c()
	for(k in 1:16){v2=c(v2,as.vector(fp[[1]][,k,m2]))}
	r2=cor(v1,v2)**2
	#line=data.frame(CHR_A=c1,BP_A=pos1,SNP_A=snp1,CHR_B=c2,BP_B=pos2,SNP_B=snp2,r2=r2,stringsAsFactors=F)
	return(r2)
}
cutoff=0.95
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

for(c in 8:9){
	fp=readRDS(sprintf('../genotypes/probabilities/geno_probs/raw/bg%.0f_genoprobs.rds',c))
	submaf=kept[kept$chr==c,]
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
	submaf$pos=pmap[match(submaf$marker,pmap$marker),]$pos
	rownames(submaf)=seq(1,nrow(submaf))
	markers=intersect(dimnames(fp[[1]])[[3]],submaf$marker)
	fp=lapply(fp,function(x) x[,,markers])
	dropped=list()
	count=1
	dropped[[count]] = list(marker=markers[1],linked=c())
	keep=c(1)
	m1=markers[1]
	size=length(markers)
	for(i in 2:size){
		m2=markers[i]
  		r2=get_r2(m1,m2)
  		if(r2<cutoff){
  			keep=c(keep,i)
  			m1=m2
  			count=count+1
  			dropped[[count]]=list(marker=m2,linked=c())
  		}else{
  			dropped[[count]]$linked=c(dropped[[count]]$linked,m2)
  		}
	}
	kept_markers=markers[keep]
	p_filtered=lapply(fp,function(x) x[,,kept_markers])
	new_fp=list()
	for(i in 1:16){
		new_fp[[i]]=p_filtered[[1]][,i,]
	}
	names(new_fp)=founders
	saveRDS(new_fp,sprintf('phenotypes/bg%.0f_updated_genoprobs.rds',c))
}


# Attempt 2


for(c in 1:10){
	fp=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c))
	fullfp=readRDS(sprintf('../genotypes/probabilities/geno_probs/raw/bg%.0f_genoprobs.rds',c))
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
	submaf=kept[kept$chr==c,]
	submaf$pos=pmap[match(submaf$marker,pmap$marker),]$pos
	submaf=submaf[order(submaf$pos),]
	rownames(submaf)=seq(1,nrow(submaf))
	
	inds=intersect(dimnames(fp[[1]])[[1]],dimnames(fullfp[[1]])[[1]])
	fp=lapply(fp,function(x) x[inds,])
	fullfp=lapply(fullfp,function(x) x[inds,,])
	markers=dimnames(fp[[1]])[[2]]
	highmaf=markers[markers %in% kept$marker]
	lowmaf=markers[!(markers %in% kept$marker)]
	if(length(lowmaf)!=0){
		dropped=readRDS(sprintf('../genotypes/probabilities/geno_probs/dropped/bg%.0f_dropped_markers_genoprobs.rds',c))
		newm=c()
		for(m in lowmaf){
			w=which(unlist(lapply(dropped,function(x) x$marker==m)))
			linked=dropped[[w]]$linked
			good=linked[linked %in% kept$marker]
			newm=c(newm,good[1])
		}
		dropped_fp=lapply(fp,function(x) x[,highmaf])
		newfp=lapply(fullfp,function(x) x[,,newm,drop=F])
		mlist=c(highmaf,newm)
		s2maf=submaf[submaf$marker %in% mlist,]
		for(k in 1:16){
			dropped_fp[[k]]=cbind(dropped_fp[[k]],newfp[[1]][,k,])
			dimnames(dropped_fp[[k]])[[2]][ncol(dropped_fp[[k]])]=newm
		}
		dropped_fp=lapply(dropped_fp,function(x) x[,s2maf$marker])
		saveRDS(dropped_fp,sprintf('phenotypes/bg%.0f_adjusted_genoprobs.rds',c))
	}
}


uniq=allmafs[allmafs$foundermaf==0.0625,]
for(c in 1:10){
	# Founder probabilities for DH lines
	fullfp=readRDS(sprintf('../genotypes/probabilities/geno_probs/raw/bg%.0f_genoprobs.rds',c))
	dimnames(fullfp[[1]])[[2]]=founders
	# Founder 600K SNPs
	psnp=fread(sprintf('../genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%.0f.csv',c),data.table=F)
	rownames(psnp)=founders
	psnp=psnp[,-1]
	psnp=apply(psnp,MARGIN=2,function(x) ifelse(x=="A",0,1))
	suniq=uniq[uniq$chr==c,]
	markers=intersect(suniq$marker,dimnames(fullfp[[1]])[[3]])
	fullfp=lapply(fullfp,function(x) x[,,markers])
	psnp=psnp[,markers]
	psums=apply(psnp,MARGIN=2,sum)
	altonly=psums[psums==1]
	markers=names(altonly)
	psnp=psnp[,markers]
	fullfp=lapply(fullfp,function(x) x[,,markers])

	#fsnp=c()
	#for(j in 1:length(psums)){
	#	p=psums[j]
	#	#print(p)
	#	sp=markers[j]
	#	if(p==1){
	#		#print(1)
	#		fs=founders[which(psnp[,sp]==1)]
	#	}else{
	#		#print(0)
	#		fs=founders[which(psnp[,sp]==0)]
	#	}
	#	fsnp=c(fsnp,fs)
	#}
	
	#apply(psnp,MARGIN=2,function(x) length(which(x==1))>1)
	fsnp=founders[unlist(unname(apply(psnp,MARGIN=2,function(x) which(x==1))))]
	names(fsnp)=markers
	snp=fread(sprintf('../genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%.0f_binary.csv',c),data.table=F)
	rownames(snp)=snp$ind
	snp=snp[,-1]
	snp=snp[,markers]
	
	inds=intersect(rownames(snp),dimnames(fullfp[[1]])[[1]])
	snp=snp[inds,]
	for(i in inds){
		print(i)
		snpi=unlist(snp[i,])
		# all the marker for which this individual has the alt allele
		snpi=names(snpi[snpi==1])
		# what is the most probable founder for this ind at those SNPs
		ifp=lapply(fullfp,function(x) x[i,,snpi])
		maxf=apply(ifp[[1]],MARGIN=2,function(x) names(which.max(x)))
		maxf=maxf[snpi]
		
		#
		matchup=fsnp[snpi]
		
		print(sum(matchup==maxf)/length(maxf)*100)
	}
	
	
}
