#!/usr/bin/env Rscript


library('data.table')
library('dplyr')

# Get cis fitted values
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

#AX-91105891
#AX-91102763
#AX-91099134
#AX-91107613
snps=c('AX-91105891','AX-91102763','AX-91099134','AX-91107613')

pheno="male_flowering_d6"
env="ALL"
chr1="8"
id1="qDTA8"
chr2="8"
id2="qDTA8"
qtl=fread('QTL/all_adjusted_QTL_peaks_trimmed.txt',data.table=F)
snp1=qtl[qtl$phenotype==pheno & qtl$environment==env & qtl$ID==id1,]$SNP
res3=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr1,pheno,env),data.table=F)
res1=res3[res3$X_ID==snp1,]

res9=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr2,pheno,env),data.table=F)
snp2=qtl[qtl$phenotype==pheno & qtl$environment==env & qtl$ID==id2,]$SNP
res2=res9[res9$X_ID==snp2,]

K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr1),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

adj_chr=c("5","9")
if(chr1 %in% adj_chr){
	X_list1=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr1))
}else{
	X_list1=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr1))
}

if(chr2 %in% adj_chr){
	X_list2=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr2))
}else{
	X_list2=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr2))
}

inter=intersect(rownames(K),dimnames(X_list1[[1]])[[1]])


X1=do.call(cbind,lapply(X_list1,function(x) x[inter,snp1]))
colnames(X1) = founders
rownames(X1) = inter
X_list_ordered1=lapply(X_list1,function(x) x[inter,snp1,drop=F])
frep1=apply(X1,MARGIN=2,function(x) round(sum(x[x>0.75])))
fkeep1=founders[frep1>3]
X_list_ordered1 = X_list_ordered1[c(fkeep1)]
betas1=unlist(res1[,fkeep1])
full_betas1=betas1[founders]
names(full_betas1)=founders
full_betas1[!(names(full_betas1) %in% fkeep1)]=0
wn1=which(!is.na(betas1))[1]
intercept1=names(betas1)[wn1]
full_betas1[names(full_betas1)!=intercept1]=full_betas1[intercept1]+full_betas1[names(full_betas1)!=intercept1]
bv1=X1 %*% full_betas1
bv1=as.data.frame(bv1,stringsAsFactors=F)
names(bv1)=snp1


X2=do.call(cbind,lapply(X_list2,function(x) x[inter,snp2]))
colnames(X2) = founders
rownames(X2) = inter
X_list_ordered2=lapply(X_list2,function(x) x[inter,snp2,drop=F])
frep2=apply(X2,MARGIN=2,function(x) round(sum(x[x>0.75])))
fkeep2=founders[frep2>3]
X_list_ordered2 = X_list_ordered2[c(fkeep2)]
betas2=unlist(res2[,fkeep2])
full_betas2=betas2[founders]
names(full_betas2)=founders
full_betas2[!(names(full_betas2) %in% fkeep2)]=0
wn2=which(!is.na(betas2))[1]
intercept2=names(betas2)[wn2]
full_betas2[names(full_betas2)!=intercept2]=full_betas2[intercept2]+full_betas2[names(full_betas2)!=intercept2]
bv2=X2 %*% full_betas2
bv2=as.data.frame(bv2,stringsAsFactors=F)
names(bv2)=snp2

cor(bv1,bv2)

cor(full_betas1,full_betas2,use="complete.obs")

# 3_2 and 9
# male flowering
#           AX-91145948
#AX-90834914   0.1977288
#  0.3081058

# 3_2 and 3_1
#            AX-91555799
#AX-90834914    0.299526

#[1] 0.6051568

# female flowering
#           AX-91146185
#AX-90834914   0.1717353

#[1] 0.3973033

#3_2 and 3_1
#            AX-91555799
#AX-90834914    0.292093

#[1] 0.5651031
snps=c('AX-91105891','AX-91102763','AX-91099134','AX-91107613')

for(i in 1:1){
	for(j in 1:nrow(new8)){
		snp1=snps[i]
		row=new8[j,]
		snp2=row$SNP
		chr=row$CHR
		pheno=row$phenotype
		env=row$environment
		pos=pmap[pmap$marker==snp,]$pos
		#results=fread(sprintf('QTL/MITE_only/results/Biogemma_chr%s_%s_x_%s_MITE_only_founderprobs.txt',chr,pheno,env),data.table=F)
		#result=results[results$X_ID==snp2,]   
	
		#betas=unlist(result[founders])
		#wn=which(!is.na(betas))[1]
		#betas[-wn]=betas[wn]+betas[-wn] 
		#betas[-1]=betas[1]+betas[-1]
		#line=c(pheno,env,snp,pos,betas)
		#new_betas=rbind(new_betas,line)
		#snp2=snps[j]
		#pheno="male_flowering_d6"
		#env="ALL"
		#chr1="8"
		#id1="qDTA8"
		#chr2="8"
		#id2="qDTA8"
		qtl=fread('QTL/all_adjusted_QTL_peaks.txt',data.table=F)#snp1=qtl[qtl$phenotype==pheno & qtl$environment==env & qtl$ID==id1,]$SNP
		res3=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno,env),data.table=F)
		res1=res3[res3$X_ID==snp1,]

		res9=fread(sprintf('QTL/MITE_only/results/Biogemma_chr%s_%s_x_%s_MITE_only_founderprobs.txt',chr,pheno,env),data.table=F)
		#snp2=qtl[qtl$phenotype==pheno & qtl$environment==env & qtl$ID==id2,]$SNP
		res2=res9[res9$X_ID==snp2,]

		phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
		phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
		phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
		phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

		data=data.frame(ID=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
		data=data[data$Loc.Year.Treat==env,]
		data=data[!is.na(data$y),]
		data$y = (data$y - mean(data$y))/sd(data$y)

		mite_prob=fread('phenotypes/mite_probabilities.txt',data.table=F)
		rownames(mite_prob)=mite_prob$ID
		mite_prob=mite_prob[data$ID,]

		has_mite=mite_prob[mite_prob$final>=0.9,]$ID

		adj_chr=c("5","9")
		if(chr %in% adj_chr){
			X_list1=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
		}else{
			X_list1=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr1))
		}

		if(chr %in% adj_chr){
			X_list2=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
		}else{
			X_list2=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr2))
		}
		inds=rownames(X_list1[[1]])
		data=data[data$ID %in% has_mite,]

		inter=intersect(data$ID,inds)

		X1=do.call(cbind,lapply(X_list1,function(x) x[inter,snp1]))
		colnames(X1) = founders
		rownames(X1) = inter
		X_list_ordered1=lapply(X_list1,function(x) x[inter,snp1,drop=F])
		frep1=apply(X1,MARGIN=2,function(x) round(sum(x[x>0.75])))
		fkeep1=founders[frep1>3]
		X_list_ordered1 = X_list_ordered1[c(fkeep1)]
		betas1=unlist(res1[,fkeep1])
		full_betas1=betas1[founders]
		names(full_betas1)=founders
		full_betas1[!(names(full_betas1) %in% fkeep1)]=0
		wn1=which(!is.na(betas1))[1]
		intercept1=names(betas1)[wn1]
		full_betas1[names(full_betas1)!=intercept1]=full_betas1[intercept1]+full_betas1[names(full_betas1)!=intercept1]
		bv1=X1 %*% full_betas1
		bv1=as.data.frame(bv1,stringsAsFactors=F)
		names(bv1)=snp1


		X2=do.call(cbind,lapply(X_list2,function(x) x[inter,snp2]))
		colnames(X2) = founders
		rownames(X2) = inter
		X_list_ordered2=lapply(X_list2,function(x) x[inter,snp2,drop=F])
		frep2=apply(X2,MARGIN=2,function(x) round(sum(x[x>0.75])))
		fkeep2=founders[frep2>3]
		X_list_ordered2 = X_list_ordered2[c(fkeep2)]
		betas2=unlist(res2[,fkeep2])
		full_betas2=betas2[founders]
		names(full_betas2)=founders
		full_betas2[!(names(full_betas2) %in% fkeep2)]=0
		wn2=which(!is.na(betas2))[1]
		intercept2=names(betas2)[wn2]
		full_betas2[names(full_betas2)!=intercept2]=full_betas2[intercept2]+full_betas2[names(full_betas2)!=intercept2]
		bv2=X2 %*% full_betas2
		bv2=as.data.frame(bv2,stringsAsFactors=F)
		names(bv2)=snp2

		print(cor(bv1,bv2,use="complete.obs"))

		print(cor(full_betas1,full_betas2,use="complete.obs"))
	}
}

full_betas=betas[founders]
names(full_betas)=founders
	full_betas[!(names(full_betas) %in% fkeep)]=0
	wn=which(!is.na(betas))[1]
	intercept=names(betas)[wn]
	full_betas[names(full_betas)!=intercept]=full_betas[intercept]+full_betas[names(full_betas)!=intercept]



factoreqtl=fread('eqtl/results/all_factor_fdr_SIs_FIXED.txt',data.table=F)
resid=fread('eqtl/results/all_residual_factor_fdr_SIs_FIXED.txt',data.table=F)
env1=factoreqtl
env1=as.data.table(env1)
env2=as.data.table(resid)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


qtl=fread('QTL/all_adjusted_QTL_all_methods.txt',data.table=F)
env1=qtl
env1=as.data.table(env1)
env2=as.data.table(factoreqtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp2=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

env1=qtl
env1=as.data.table(env1)
env2=as.data.table(resid)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp4=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


### How correlated are the fitted values for resid and whole gene
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

for(i in 1:nrow(comp)){
	row=comp[i,]
	time1=row$time
	factor1=row$factor
	time2=row$i.time
	factor2=row$i.factor
	snp1=row$SNP
	snp2=row$i.SNP
	chr=row$CHR
	f1=fread(sprintf('eqtl/trans/results/%s_residuals_c%s_%s_trans_results_FIXED.txt',time1,chr,factor1),data.table=F)
	res1=f1[f1$X_ID==snp1,]
	f2=fread(sprintf('eqtl/trans/results/%s_c%s_%s_trans_results_FIXED.txt',time2,chr,factor2),data.table=F)
	res2=f2[f2$X_ID==snp2,]
	
	adj_chr=c(5,9)
	if(chr %in% adj_chr){
		X_list1=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
	}else{
		X_list1=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	}
	inds=rownames(X_list1[[1]])
	samples=fread(sprintf('eqtl/data/%s_samples_FIXED.txt',time),data.table=F)
	inter=intersect(samples$id,inds)
	X1=do.call(cbind,lapply(X_list1,function(x) x[inter,snp1]))
	colnames(X1) = founders
	rownames(X1) = inter
	X_list_ordered1=lapply(X_list1,function(x) x[inter,snp1,drop=F])
	frep1=apply(X1,MARGIN=2,function(x) round(sum(x[x>0.75])))
	fkeep1=founders[frep1>3]
	X_list_ordered1 = X_list_ordered1[c(fkeep1)]
	betas1=unlist(res1[,fkeep1])
	full_betas1=betas1[founders]
	names(full_betas1)=founders
	full_betas1[!(names(full_betas1) %in% fkeep1)]=0
	wn1=which(!is.na(betas1))[1]
	intercept1=names(betas1)[wn1]
	full_betas1[names(full_betas1)!=intercept1]=full_betas1[intercept1]+full_betas1[names(full_betas1)!=intercept1]
	bv1=X1 %*% full_betas1
	bv1=as.data.frame(bv1,stringsAsFactors=F)
	names(bv1)=snp1

	X2=do.call(cbind,lapply(X_list1,function(x) x[inter,snp2]))
	colnames(X2) = founders
	rownames(X2) = inter
	X_list_ordered2=lapply(X_list1,function(x) x[inter,snp2,drop=F])
	frep2=apply(X2,MARGIN=2,function(x) round(sum(x[x>0.75])))
	fkeep2=founders[frep2>3]
	X_list_ordered2 = X_list_ordered2[c(fkeep2)]
	betas2=unlist(res2[,fkeep2])
	full_betas2=betas2[founders]
	names(full_betas2)=founders
	full_betas2[!(names(full_betas2) %in% fkeep2)]=0
	wn2=which(!is.na(betas2))[1]
	intercept2=names(betas2)[wn2]
	full_betas2[names(full_betas2)!=intercept2]=full_betas2[intercept2]+full_betas2[names(full_betas2)!=intercept2]
	bv2=X2 %*% full_betas2
	bv2=as.data.frame(bv2,stringsAsFactors=F)
	names(bv2)=snp2

	print(time1)
	print(factor1)
	print(time2)
	print(factor2)
	print(cor(bv1,bv2,use="complete.obs"))
	print(cor(full_betas1,full_betas2,use="complete.obs"))

}


### How correlated are the fitted values for whole factor eqtl and QTL?
for(i in 1:nrow(comp2)){
	row=comp2[i,]
	time1=row$time
	factor1=row$factor
	pheno=row$phenotype
	env=row$environment
	method=row$method
	snp1=row$SNP
	snp2=row$i.SNP
	chr=row$CHR
	f1=fread(sprintf('eqtl/trans/results/%s_c%s_%s_trans_results_FIXED.txt',time1,chr,factor1),data.table=F)
	res1=f1[f1$X_ID==snp1,]
	f2=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno,env),data.table=F)
	res2=f2[f2$X_ID==snp2,]
	
	adj_chr=c(5,9)
	if(chr %in% adj_chr){
		X_list1=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
	}else{
		X_list1=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	}
	inds=rownames(X_list1[[1]])
	samples=fread(sprintf('eqtl/data/%s_samples_FIXED.txt',time1),data.table=F)
	inter=intersect(samples$id,inds)
	X1=do.call(cbind,lapply(X_list1,function(x) x[inter,snp1]))
	colnames(X1) = founders
	rownames(X1) = inter
	X_list_ordered1=lapply(X_list1,function(x) x[inter,snp1,drop=F])
	frep1=apply(X1,MARGIN=2,function(x) round(sum(x[x>0.75])))
	fkeep1=founders[frep1>3]
	X_list_ordered1 = X_list_ordered1[c(fkeep1)]
	betas1=unlist(res1[,fkeep1])
	full_betas1=betas1[founders]
	names(full_betas1)=founders
	full_betas1[!(names(full_betas1) %in% fkeep1)]=0
	wn1=which(!is.na(betas1))[1]
	intercept1=names(betas1)[wn1]
	full_betas1[names(full_betas1)!=intercept1]=full_betas1[intercept1]+full_betas1[names(full_betas1)!=intercept1]
	bv1=X1 %*% full_betas1
	bv1=as.data.frame(bv1,stringsAsFactors=F)
	names(bv1)=snp1

	X2=do.call(cbind,lapply(X_list1,function(x) x[inter,snp2]))
	colnames(X2) = founders
	rownames(X2) = inter
	X_list_ordered2=lapply(X_list1,function(x) x[inter,snp2,drop=F])
	frep2=apply(X2,MARGIN=2,function(x) round(sum(x[x>0.75])))
	fkeep2=founders[frep2>3]
	X_list_ordered2 = X_list_ordered2[c(fkeep2)]
	betas2=unlist(res2[,fkeep2])
	full_betas2=betas2[founders]
	names(full_betas2)=founders
	full_betas2[!(names(full_betas2) %in% fkeep2)]=0
	wn2=which(!is.na(betas2))[1]
	intercept2=names(betas2)[wn2]
	full_betas2[names(full_betas2)!=intercept2]=full_betas2[intercept2]+full_betas2[names(full_betas2)!=intercept2]
	bv2=X2 %*% full_betas2
	bv2=as.data.frame(bv2,stringsAsFactors=F)
	names(bv2)=snp2

	print(time1)
	print(factor1)
	print(pheno)
	print(env)
	print(cor(bv1,bv2,use="complete.obs"))
	print(cor(full_betas1,full_betas2,use="complete.obs"))

}

### How correlated are the fitted values for residual eqtl and QTL?
for(i in 1:nrow(comp4)){
	row=comp4[i,]
	time1=row$time
	factor1=row$factor
	pheno=row$phenotype
	env=row$environment
	method=row$method
	snp1=row$SNP
	snp2=row$i.SNP
	chr=row$CHR
	f1=fread(sprintf('eqtl/trans/results/%s_residuals_c%s_%s_trans_results_FIXED.txt',time1,chr,factor1),data.table=F)
	res1=f1[f1$X_ID==snp1,]
	adj_chr=c(5,9)
	if(chr %in% adj_chr){
		X_list1=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
	}else{
		X_list1=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	}
	inds=rownames(X_list1[[1]])
	samples=fread(sprintf('eqtl/data/%s_samples_FIXED.txt',time1),data.table=F)
	inter=intersect(samples$id,inds)
	X1=do.call(cbind,lapply(X_list1,function(x) x[inter,snp1]))
	colnames(X1) = founders
	rownames(X1) = inter
	X_list_ordered1=lapply(X_list1,function(x) x[inter,snp1,drop=F])
	frep1=apply(X1,MARGIN=2,function(x) round(sum(x[x>0.75])))
	fkeep1=founders[frep1>3]
	X_list_ordered1 = X_list_ordered1[c(fkeep1)]
	betas1=unlist(res1[,fkeep1])
	full_betas1=betas1[founders]
	names(full_betas1)=founders
	full_betas1[!(names(full_betas1) %in% fkeep1)]=0
	wn1=which(!is.na(betas1))[1]
	intercept1=names(betas1)[wn1]
	full_betas1[names(full_betas1)!=intercept1]=full_betas1[intercept1]+full_betas1[names(full_betas1)!=intercept1]
	bv1=X1 %*% full_betas1
	bv1=as.data.frame(bv1,stringsAsFactors=F)
	names(bv1)=snp1

	
	if(method=="600K_SNP"){
		f2=readRDS(sprintf('../GridLMM/GridLMM_600KSNP/models/chr%s_%s_x_%s_600KSNP_ML.rds',chr,pheno,env))
		f2=f2$results
		res2=f2[f2$X_ID==snp2,]
		
		X_list=fread(sprintf('../genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
		rownames(X_list)=X_list$ind
		inter=intersect(samples$id,X_list$ind)
		X2=X_list[inter,snp2]
	
		betas2=unlist(res2[,c('beta.1','beta.2')])
		betas2[-1]=betas2[1]+betas2[-1]
		bv2=ifelse(X2==0,betas2[1],betas2[2])
		#bv2=X2 %*% full_betas2
		bv2=as.data.frame(bv2,stringsAsFactors=F)
		names(bv2)=snp2
		rownames(bv2)=inter
		fgeno=fread(sprintf('../genotypes/qtl2/Biogemma_foundergenos/Founder_genos_chr%s.csv',chr),data.table=F)
		fid=fgeno[,snp2]
		full_betas2=ifelse(fid=="A",betas2[1],betas2[2])
	}else{
		keeph=0
		for(h in 8:16){
			hapfile=sprintf('../genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h)
			if(file.exists(hapfile)){
				X_list=readRDS(hapfile)
				if(snp2 %in% dimnames(X_list[[1]])[[2]]){
					keeph=h
				}
			}
		}
		ibd=fread(sprintf('../ibd_segments/refinedibd/600K/bg%s_refined_ibd_blocks.txt',chr),data.table=F)
		pos=row$i.BP
		line=unlist(ibd[ibd$start<pos & ibd$end>pos,founders])

		f2=readRDS(sprintf('../GridLMM/GridLMM_haplotypes/models/Biogemma_chr%s_haplogrp%.0f_%s_x_%s.rds',chr,keeph,pheno,env))
		res2=f2[f2$X_ID==snp2,]
		bnames=paste0('beta.',seq(1:keeph))
		hapfile=sprintf('../genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,keeph)
		X_list=readRDS(hapfile)
		X2=do.call(cbind,lapply(X_list,function(x) x[inter,snp2]))
		colnames(X2) = bnames
		rownames(X2) = inter
		X_list_ordered2=lapply(X_list,function(x) x[inter,snp2,drop=F])
		frep2=apply(X2,MARGIN=2,function(x) round(sum(x[x>0.75])))
		fkeep2=bnames[frep2>3]
		X_list_ordered2 = X_list_ordered2[c(fkeep2)]
		betas2=unlist(res2[,fkeep2])
		full_betas2=betas2[line]
		names(full_betas2)=paste0('beta.',1:16)
		full_betas2[!(names(full_betas2) %in% fkeep2)]=0
		wn2=which(!is.na(betas2))[1]
		intercept2=names(betas2)[wn2]
		full_betas2[names(full_betas2)!=intercept2]=full_betas2[intercept2]+full_betas2[names(full_betas2)!=intercept2]
		betas2[-1]=betas2[1]+betas2[-1]
		bv2=X2 %*% betas2
		bv2=as.data.frame(bv2,stringsAsFactors=F)
		names(bv2)=snp2
	}

	print(time1)
	print(factor1)
	print(pheno)
	print(env)
	print(cor(bv1,bv2,use="complete.obs"))
	print(cor(full_betas1,full_betas2,use="complete.obs"))

}

#            AX-91102763
#AX-91105891   0.6448342
#[1] 0.938975

#            AX-91099134
#AX-91105891   0.4059834
#[1] 0.8141023

#            AX-91107613
#AX-91105891    0.768967
#[1] 0.9652284

#            AX-91099134
#AX-91102763   0.6854379
#[1] 0.9166121

#            AX-91107613
#AX-91102763   0.4587304
#[1] 0.9113466

#            AX-91107613
#AX-91099134   0.3485721
#[1] 0.7700222

pmap[pmap$marker %in% snps,]
#           marker chr       pos
#25128 AX-91099134   8 121124614
#28709 AX-91102763   8 135165296
#32144 AX-91105891   8 146537894
#34146 AX-91107613   8 153154476


bp=c()
for(i in 1:nrow(allqtl)){
	row=allqtl[i,]
	chr=row$CHR
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
	snp=row$SNP
	pos=pmap[pmap$marker==snp,]$pos
	bp=c(bp,pos)
}
allqtl$BP=bp
allqtl=allqtl[,c('time','factor','CHR','BP','SNP','left_bound_bp','alt_right_bound_bp','left_bound_snp','right_bound_snp','right_bound_bp')]
#fwrite(allqtl,'eqtl/results/all_factor_fdr_SIs_FIXED.txt',row.names=F,quote=F,sep='\t')
fwrite(allqtl,'eqtl/results/all_residual_factor_fdr_SIs_FIXED.txt',row.names=F,quote=F,sep='\t')


#allqtl=allqtl[,c('ID','phenotype','environment','CHR','BP','SNP','method','left_bound_bp','alt_right_bound_bp','left_bound_snp','right_bound_snp','right_bound_bp')]