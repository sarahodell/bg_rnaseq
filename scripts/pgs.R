#!/usr/bin/env Rscript

library('data.table')
#library('tidyverse')
library('ggplot2')
library('parallel')
library('MASS')



  #Founders
qtl=fread('QTL/all_adjusted_QTL_SIs.txt',data.table=F)
pheno="male_flowering_d6"
env="ALL"
qtl=qtl[qtl$phenotype==pheno & qtl$environment==env,]
## DTA PGS
K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%.0f.txt',c),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
#env="EXP_STPAUL_2017_WD"

# 3 AX-90796752 instead of AX-91555799
# 8 AX-91112984

qsnps=c("AX-90834914","AX-90796752","AX-91105891","AX-91112984","AX-91145948")

all_pgs=c()
chroms=unique(qtl$CHR)
for(c in chroms){
	row1=qtl[qtl$CHR==c,]
	qsnps=unique(row1$SNP)
	fes=fread(sprintf('QTL/adjusted/Biogemma_chr%.0f_%s_x_%s_unscaled_founderprobs.txt',c,pheno,env),data.table=F)
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c))
	inds=rownames(X_list[[1]])
	inter=intersect(phenotypes$Genotype_code,inds)
	X_list=lapply(X_list,function(x) x[inter,])
	markers=colnames(X_list[[1]])
	rownames(fes)=fes$X_ID
	fes=fes[qsnps,]
	betas = fes[,founders]
	wn=apply(betas,MARGIN=1,function(x) which(!is.na(x))[1]) 
	for(x in 1:nrow(betas)){
		betas[x,-wn[x]] = betas[x,-wn[x]] + betas[x,wn[x]]
	}
	betas[is.na(betas)]=0
	pgs=sapply(qsnps,function(x) as.matrix(do.call(cbind,lapply(X_list,function(j) j[,x]))) %*% unlist(betas[x,]))

	#tmp2=as.data.frame(sapply(seq(1,nrow(betas)),function(x) tmp[x,-wn[x]] = tmp[x,-wn[x]] + tmp[x,wn[x]]))
	rownames(pgs)=inter
	ft_pgs=as.data.frame(rowSums(pgs,na.rm=T))
	if(c==3){
		all_pgs=ft_pgs
	}else{
		all_pgs=cbind(all_pgs,ft_pgs)
	}
}


allp=data.frame(inds=inter,dta=rowSums(all_pgs,na.rm=T),stringsAsFactors=F)
fwrite(allp,'QTL/Founder_ALL_DTA_polygenic_scores.txt',quote=F,sep='\t',row.names=F)
# variance 533.201

#TPH
all_pgs=c()
for(c in 1:10){
	fes=fread(sprintf('QTL/Biogemma_chr%.0f_total_plant_height_x_%s_vst_founderprobs.txt',c,env),data.table=F)
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c))
	inds=rownames(X_list[[1]])
	inter=intersect(phenotypes$ID,inds)
	X_list=lapply(X_list,function(x) x[inter,])
	markers=colnames(X_list[[1]])
	rownames(fes)=fes$X_ID
	fes=fes[markers,]
	betas = fes[,founders]
	wn=apply(betas,MARGIN=1,function(x) which(!is.na(x))[1]) 
	for(x in 1:nrow(betas)){
		betas[x,-wn[x]] = betas[x,-wn[x]] + betas[x,wn[x]]
	}
	betas[is.na(betas)]=0
	pgs=sapply(markers,function(x) as.matrix(do.call(cbind,lapply(X_list,function(j) j[,x]))) %*% unlist(betas[x,]))

	#tmp2=as.data.frame(sapply(seq(1,nrow(betas)),function(x) tmp[x,-wn[x]] = tmp[x,-wn[x]] + tmp[x,wn[x]]))
	rownames(pgs)=inter
	ft_pgs=as.data.frame(rowSums(pgs,na.rm=T))
	if(c==1){
		all_pgs=ft_pgs
	}else{
		all_pgs=cbind(all_pgs,ft_pgs)
	}
}


allp=data.frame(inds=inter,tph=rowSums(all_pgs,na.rm=T),stringsAsFactors=F)
fwrite(allp,'QTL/Founder_EXP_STPAUL_2017_WD_TPH_polygenic_scores.txt',quote=F,sep='\t',row.names=F)


## Gene eQTL
gene="Zm00001d046357"
time="WD_0712"
genechr=9
all_pgs=c()
for(c in 1:10){
	fes=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%.0f_weights_results.rds',time,c))
	w=which(unlist(lapply(fes,function(x) unique(x$Trait)==gene)))
	fes=fes[[w]]
	if(genechr==c){
		cises=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_weights_results.txt',time,c),data.table=F)
		cises=cises[cises$Trait==gene,]
		fes=rbind(cises,fes)

	}
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',c))
	inds=rownames(X_list[[1]])
	inter=intersect(phenotypes$ID,inds)
	X_list=lapply(X_list,function(x) x[inter,])
	markers=colnames(X_list[[1]])
	fes=fes[!is.na(fes$X_ID),]
	markers=intersect(markers,fes$X_ID)
	rownames(fes)=fes$X_ID
	fes=fes[markers,]
	betas = fes[,founders]
	wn=apply(betas,MARGIN=1,function(x) which(!is.na(x))[1]) 
	for(x in 1:nrow(betas)){
		betas[x,-wn[x]] = betas[x,-wn[x]] + betas[x,wn[x]]
	}
	betas[is.na(betas)]=0
	pgs=sapply(markers,function(x) as.matrix(do.call(cbind,lapply(X_list,function(j) j[,x]))) %*% unlist(betas[x,]))

	#tmp2=as.data.frame(sapply(seq(1,nrow(betas)),function(x) tmp[x,-wn[x]] = tmp[x,-wn[x]] + tmp[x,wn[x]]))
	rownames(pgs)=inter
	ft_pgs=as.data.frame(rowSums(pgs,na.rm=T))
	if(c==1){
		all_pgs=ft_pgs
	}else{
		all_pgs=cbind(all_pgs,ft_pgs)
	}
}


allp=data.frame(inds=inter,eqts=rowSums(all_pgs,na.rm=T),stringsAsFactors=F)
fwrite(allp,sprintf('QTT/Founder_%s_%s_polygenic_scores.txt',time,gene),quote=F,sep='\t',row.names=F)

dta=fread('QTL/Founder_EXP_STPAUL_2017_WD_DTA_polygenic_scores.txt',data.table=F)

allp$dta=dta[match(allp$inds,dta$ind),]$dta
#allp=allp[order(allp$)]
p1=ggplot(allp,aes(x=eqts,y=dta)) + geom_point()

png(sprintf('QTT/DTA_%s_%s_eQTS.png',time,gene))
print(p1)
dev.off()

fvalue=fread('MegaLMM/MegaLMM_WD_0727_all_F_means.txt',data.table=F)

dta$f2=fvalue[match(dta$inds,fvalue$V1),]$Factor2

p1=ggplot(dta,aes(x=f2,y=dta)) + geom_point()

allp=fread('QTL/Founder_EXP_STPAUL_2017_WD_TPH_polygenic_scores.txt',data.table=F)
allp$zscore= (allp$tph-mean(allp$tph))/sd(allp$tph)
allp$f2=fvalue[match(allp$inds,fvalue$V1),]$Factor2

p1=ggplot(allp,aes(x=f2,y=zscore)) + geom_point()


png('QTT/TPH_WD_0727_Factor2_Fvalue.png')
print(p1)
dev.off()

phenotypes$f2=fvalue[match(phenotypes$ID,fvalue$V1),]$Factor2

p1=ggplot(phenotypes,aes(x=f2,y=total_plant_height)) + geom_point()


png('QTT/TPH_pheno_WD_0727_Factor2_Fvalue.png')
print(p1)
dev.off()
