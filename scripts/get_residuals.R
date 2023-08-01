args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
chr=as.character(args[[2]])
cores=as.numeric(args[[3]])

#date=format(Sys.time(),'%m%d%y')

library('data.table')
library('dplyr')
#library('lme4')
library('lme4qtl')
library('preprocessCore')
library('stringr')
library('parallel')
library('MASS')

#time="WD_0720"
#chr="10"
#library('GenomicFeatures') # write a script to get a table of start and stop sites of genes from the gtf file
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
all_gwas=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
all_gwas=all_gwas %>% group_by(Trait) %>% slice(which.min(p_value_ML))
all_gwas=as.data.frame(all_gwas)
all_betas=all_gwas[,c('Trait','X_ID',founders)]
#options(warn=2)
# Read in Kinship Matrix
K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
genetable=genetable[genetable$CHROM==chr,]
genes=unique(genetable$Gene_ID)
# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
metadata=fread('metadata/samples_passed_genotype_check.txt',data.table=F)


metadata=metadata[metadata$experiment==time,]

genos=phenotypes$V1
######
adj_chr=c("5","9")
if(chr %in% adj_chr){
	X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))

}else{
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	#chr 5 "AX-91671957" replaced with "AX-91671943"
}
#X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])
inter=intersect(genos,inds)

genes=intersect(genes,names(phenotypes)[-1])
rownames(phenotypes)=phenotypes$V1
phenotypes=phenotypes[,-1]
phenotypes=as.matrix(phenotypes)
phenotypes=phenotypes[inter,]
genos=rownames(phenotypes)

K=K[inter,inter]
X_list=lapply(X_list,function(x) x[inter,])


get_resids<-function(row){
	test=all_betas[row,]
	gene=test$Trait
	data=data.frame(ID=rownames(phenotypes),y=phenotypes[,gene],stringsAsFactors=F)
	data=data[!is.na(data$y),]
	rownames(data)=data$ID
	data=data[inter,]
	
	snp=test$X_ID
	X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
	colnames(X) = founders
	rownames(X) = dimnames(X_list[[1]])[[1]]
	X_list_ordered=lapply(X_list,function(x) x[,snp,drop=F])
	frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
	fkeep=founders[frep2>3]
	X_list_ordered = X_list_ordered[c(fkeep)]
	#X=X[,fkeep]
	betas=unlist(test[,fkeep])
	full_betas=betas[founders]
	names(full_betas)=founders
	full_betas[!(names(full_betas) %in% fkeep)]=0
	wn=which(!is.na(betas))[1]
	intercept=names(betas)[wn]
	full_betas[names(full_betas)!=intercept]=full_betas[intercept]+full_betas[names(full_betas)!=intercept]
	bv=X %*% full_betas
	colnames(bv)=gene
	X_r=data$y-bv
	colnames(X_r)=gene
	prop_var=var(bv[,gene],na.rm=T)/var(data$y,na.rm=T)
	line2=data.frame(gene=gene,time=time,snp=snp,prop_var=prop_var,stringsAsFactors=F)
	return(list(X_r,line2))
}

#n_reps=1:10
n_reps=1:nrow(all_betas)

print(system.time({
results=mclapply(n_reps,get_resids,mc.cores=cores)
}))

all_res=do.call(cbind,lapply(results,function(x) x[[1]]))
all_res=as.data.frame(all_res,stringsAsFactors=F)
rownames(all_res)=inter
fwrite(all_res,sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_residuals_FIXED.txt',time,chr),row.names=T,quote=F,sep='\t')

all_props=do.call(rbind,lapply(results,function(x) x[[2]]))
all_props=as.data.frame(all_props,stringsAsFactors=F)
fwrite(all_props,sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_prop_var_FIXED.txt',time,chr),row.names=F,quote=F,sep='\t')






