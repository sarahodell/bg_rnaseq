#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
rep=as.numeric(args[[1]])

#library('SNPRelate')
library('data.table')
library('rrBLUP')
library('lme4')

sample.id=paste0('Sim',seq(1,325))
pop_code=c(rep("DH",325))


#vcf.fn<-sprintf('merged_vcfs/MAGIC_DHSimAll_rep%.0f.vcf.gz',rep)
#snpgdsVCF2GDS(vcf.fn,sprintf("pgs/Rep%.0f.gds",rep),method="biallelic.only")
#genofile<-snpgdsOpen(sprintf("pgs/Rep%.0f.gds",rep))
#snpgdsSummary(sprintf("pgs/Rep%.0f.gds",rep))

#snpset1<-snpgdsLDpruning(genofile,ld.threshold = 0.2)
#With LD 0.2 used 18,423 markers
#snpset1.id<-unlist(snpset1)

#snprslist=read.gdsn(index.gdsn(genofile,"snp.rs.id"))[snpset1.id]
qtl=fread('QTL/all_adjusted_QTL_SIs.txt',data.table=F)
pheno="male_flowering_d6"
env="ALL"
qtl=qtl[qtl$phenotype==pheno & qtl$environment==env,]
## DTA PGS
K=fread('../GridLMM/K_matrices/K_matrix_chr10.txt',data.table=F)
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

#qsnps=c("AX-90834914","AX-90796752","AX-91105891","AX-91112984","AX-91145948")

all_pgs=c()
chroms=unique(qtl$CHR)
for(chr in chroms){
	row1=qtl[qtl$CHR==chr,]
	qsnps=unique(row1$SNP)
	fes=fread(sprintf('QTL/adjusted/Biogemma_chr%.0f_%s_x_%s_unscaled_founderprobs.txt',chr,pheno,env),data.table=F)
	rownames(fes)=fes$X_ID
	fes=fes[qsnps,]
	betas = fes[,founders]
	wn=apply(betas,MARGIN=1,function(x) which(!is.na(x))[1]) 
	for(x in 1:nrow(betas)){
		betas[x,-wn[x]] = betas[x,-wn[x]] + betas[x,wn[x]]
	}
	betas[is.na(betas)]=0
	
	X_list=readRDS(sprintf('../run_magicsim/qtl2_files/MAGIC_DHSim_rep%.0f_c%.0f_genoprobs_v3.rds',rep,chr))
	# Grab 325 of them
	# get founder probs for the markers with QTL
	X_list=X_list[[1]][,,qsnps]
	dimnames(X_list)[[2]]=founders
	new_X_list=list()
	if(length(qsnps)==1){
		for(f in 1:16){
			new_X_list[[f]]=X_list[,f,drop=F]
		}
		pgs=as.matrix(do.call(cbind,lapply(new_X_list,function(j) j[,1]))) %*% unlist(betas[x,])

	}else{
		for(f in 1:16){
			new_X_list[[f]]=X_list[,f,]
		}
		pgs=sapply(qsnps,function(x) as.matrix(do.call(cbind,lapply(new_X_list,function(j) j[,x]))) %*% unlist(betas[x,]))

	}
	
	#inds=rownames(X_list[[1]])
	#inter=intersect(phenotypes$Genotype_code,inds)
	#X_list=lapply(X_list,function(x) x[inter,])
	#markers=colnames(X_list[[1]])
	

	#tmp2=as.data.frame(sapply(seq(1,nrow(betas)),function(x) tmp[x,-wn[x]] = tmp[x,-wn[x]] + tmp[x,wn[x]]))
	rownames(pgs)=dimnames(X_list)[[1]]
	ft_pgs=as.data.frame(rowSums(pgs,na.rm=T))
	if(chr==3){
		all_pgs=ft_pgs
	}else{
		all_pgs=cbind(all_pgs,ft_pgs)
	}
}


allp=data.frame(inds=rownames(all_pgs),dta=rowSums(all_pgs,na.rm=T),stringsAsFactors=F)
fwrite(allp,sprintf('QTL/pgs_sim/Founder_ALL_DTA_rep%.0f_polygenic_scores_v3.txt',rep),quote=F,sep='\t',row.names=F)
