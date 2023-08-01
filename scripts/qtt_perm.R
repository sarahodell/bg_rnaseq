args=commandArgs(trailingOnly=T)
row=as.numeric(args[[1]])
cores=as.numeric(args[[2]])


library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('data.table')
library('dplyr')
library('lme4qtl')
library('preprocessCore',lib='/home/sodell/R/x86_64-conda-linux-gnu/4.2')
library('stringr')
library('parallel')
library('MASS')

#comp=fread('QTT/QTL_cis_eQTL_overlap.txt',data.table=F)
comp=fread('QTT/cis_eQTL_STPAUL_QTL_overlaps.txt',data.table=F)
test=comp[row,]
#test
gene=test$Trait
time=test$time
chr=test$CHR
env=test$environment
pheno=test$phenotype
qsnp=test$SNP
xid=test$X_ID
qtlid=test$ID

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
exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
#metadata=fread('metadata/BG_completed_sample_list_FIXED.txt',data.table=F)
metadata=fread('metadata/samples_passed_genotype_check.txt',data.table=F)

pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates_2.txt',time),data.table=F)

metadata=metadata[metadata$experiment==time,]
#metadata=metadata[metadata$read==1,]

genos=exp$V1
######
adj_chr=c(5,9)
if(chr %in% adj_chr){
	X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))

}else{
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
}
inds=rownames(X_list[[1]])
inter=intersect(genos,inds)
X_list=lapply(X_list,function(x) x[inter,])

genes=intersect(genes,colnames(exp)[-1])
rownames(exp)=exp$V1
exp=exp[,-1]
exp=as.matrix(exp)
exp=exp[inter,]
genos=rownames(exp)

K=K[inter,inter]

plate=metadata[match(genos,metadata$dh_genotype),]$plate
df=data.frame(ID=genos,plate=plate,stringsAsFactors=F)
df$plate=as.factor(df$plate)

# permute

testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

allweights=fread(sprintf('eqtl/normalized/%s_voom_weights_2.txt',time),data.table=F)
allweights=allweights[,c('V1',inter)]

phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

eqtl_perm=function(draw){
	#snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
	gweights=unlist(allweights[allweights$V1==gene,draw])
	gweights=gweights[draw]
	data=data.frame(ID=rownames(exp),y=exp[,gene],stringsAsFactors=F)
	data=data[!is.na(data$y),]

    data$PC1=pcs[match(data$ID,pcs$V1),]$PC1
    data$PC2=pcs[match(data$ID,pcs$V1),]$PC2
    data$PC3=pcs[match(data$ID,pcs$V1),]$PC3

    rownames(data)=data$ID
    data=data[draw,]
    data$ID2=data$ID
    data=data[,c('ID','ID2','y','PC1','PC2','PC3')]
    newK=K
	rownames(newK)=draw
	colnames(newK)=draw
    null_model = GridLMM_ML(y~1+PC1+PC2+PC3+(1|ID),data=data,weights=gweights,relmat=list(ID=newK),ML=T,REML=F,verbose=F)
	
    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    Y=as.matrix(data$y)
    X_cov=null_model$lmod$X
    X_list_null=NULL
	X_list_reordered=lapply(X_list,function(x) x[draw,])
	for(x in seq(1,16)){
   		dimnames(X_list_reordered[[x]])[[1]]=dimnames(X_list[[1]])[[1]]
	}
    X = do.call(cbind,lapply(X_list_reordered,function(x) x[,xid]))
    colnames(X) = founders
    rownames(X) = dimnames(X_list_reordered[[1]])[[1]]
    X_list_reordered=lapply(X_list_reordered,function(x) x[,xid,drop=F])
    frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
    fkeep=founders[frep2>3]
   	X_list_reordered = X_list_reordered[c(fkeep)]
	X=X[,fkeep]
	gwas=run_GridLMM_GWAS(Y,X_cov,X_list_reordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=1,verbose=F)
    gwas$Trait=gene
	betas=unlist(gwas[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	names(betas)=fkeep
	return(betas)
}



qtl_perm=function(draw){
	data=data.frame(ID=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
	data=data[!is.na(data$y),]
	#data$y = (data$y - mean(data$y))/sd(data$y)
	rownames(data)=data$ID
	data=data[draw,]
	newK=K
	rownames(newK)=draw
	colnames(newK)=draw
	null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=newK),ML=T,REML=F,verbose=F)

	h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
	names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
	h2_start
	V_setup=null_model$setup

	Y=as.matrix(data$y)
	X_cov=null_model$lmod$X
	X_list_null=NULL
	
	X_list_reordered=lapply(X_list,function(x) x[draw,])
	for(x in seq(1,16)){
   		dimnames(X_list_reordered[[x]])[[1]]=dimnames(X_list[[1]])[[1]]
	}
   	X = do.call(cbind,lapply(X_list_reordered,function(x) x[,qsnp]))
    colnames(X) = founders
    rownames(X) = dimnames(X_list_reordered[[1]])[[1]]
    X_list_reordered=lapply(X_list_reordered,function(x) x[,qsnp,drop=F])
    frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
    fkeep=founders[frep2>3]
   	X_list_reordered = X_list_reordered[c(fkeep)]
	X=X[,fkeep]
  	gwas=run_GridLMM_GWAS(Y,X_cov,X_list_reordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=1,verbose=F)
    #names(gwas)[6:21]=founders
    
    betas=unlist(gwas[,6:21])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	names(betas)=fkeep
	return(betas)	
}


perm_cor=function(rep){
	draw=sample(inter,length(inter))
	eqtl_betas=eqtl_perm(draw)
	qtl_betas=qtl_perm(draw)
	ctest=cor.test(eqtl_betas,qtl_betas,use="complete.obs")
	r=ctest$estimate
	p=ctest$p.value
	line=data.frame(rep=rep,r=r,pvalue=p,stringsAsFactors=F)
	return(line)
}

n_reps=1:1000

print(system.time({
results=mclapply(n_reps,perm_cor,mc.cores=cores)
}))
d=rbindlist(results)
d=as.data.frame(d)
fwrite(d,sprintf('QTT/permute/%s_%s_%s_%s_random_permutation.txt',gene,time,pheno,env))