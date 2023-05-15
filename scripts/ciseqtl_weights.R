args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
chr=as.character(args[[2]])
cores=as.numeric(args[[3]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('data.table')
library('dplyr')
#library('lme4')
library('lme4qtl')
library('preprocessCore')

#time="WD_0720"
#chr="10"
#library('GenomicFeatures') # write a script to get a table of start and stop sites of genes from the gtf file


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
phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time),data.table=F)

geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
phenotypes=phenotypes[,c('V1',kept_genes)]


metadata=metadata[metadata$experiment==time,]

genos=phenotypes$V1
######

X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
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

testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

allweights=fread(sprintf('eqtl/normalized/%s_voom_weights.txt',time),data.table=F)
allweights=allweights[,c('V1',inter)]


all_gwas=data.frame(matrix(ncol=29,nrow=0))
names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',"B73_inra","PC1","PC2","PC3",founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')

#all_res=data.frame(matrix(ncol=length(genes),nrow=length(i)))
#names(all_res)=genes
for(g in 1:length(genes)){
  gene=genes[g]
  gweights=unlist(allweights[allweights$V1==gene,inter])

  snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
  if(length(snp)!=0){
    data=data.frame(ID=rownames(phenotypes),y=phenotypes[,gene],stringsAsFactors=F)
    data=data[!is.na(data$y),]

    data$PC1=pcs[match(data$ID,pcs$V1),]$PC1
    data$PC2=pcs[match(data$ID,pcs$V1),]$PC2
    data$PC3=pcs[match(data$ID,pcs$V1),]$PC3

    rownames(data)=data$ID
    data=data[inter,]
    data$ID2=data$ID
    data=data[,c('ID','ID2','y','PC1','PC2','PC3')]
    null_model = GridLMM_ML(y~1+PC1+PC2+PC3+(1|ID),data=data,weights=gweights,relmat=list(ID=K),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    Y=as.matrix(data$y)
    X_cov=null_model$lmod$X
    if(length(snp)>1){
        X = do.call(cbind,lapply(X_list,function(x) x[,snp[1]]))
        colnames(X) = founders
        rownames(X) = dimnames(X_list[[1]])[[1]]
        snp=snp[1]
        #X=X[i,]
    }else{
    	X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
    	colnames(X) = founders
    	rownames(X) = dimnames(X_list[[1]])[[1]]
    }
    X_list_ordered=lapply(X_list,function(x) x[,snp,drop=F])
    frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
    fkeep=founders[frep2>3]
    X_list_ordered = X_list_ordered[c(fkeep)]

    X=X[,fkeep]

    X_list_null=NULL

    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
    gwas$Trait=gene

    if(!("B73_inra" %in% fkeep)){
      end=10+length(fkeep)-2
      betas=gwas[,c(6,10:end)]
      betas=unlist(unname(betas))
      betas[-1]=betas[1]+betas[-1]
      names(betas)=fkeep
      new_gwas=gwas[,1:5]
      ncol1=data.frame(matrix(ncol=1,nrow=1))
      names(ncol1)='B73_inra'
      new_gwas=cbind(new_gwas,ncol1)
      new_gwas=cbind(new_gwas,gwas[,7:9])
      new_gwas=cbind(new_gwas,gwas[,c(6,10:end)])
      names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra',"PC1","PC2","PC3",fkeep)
      fdrop=founders[!(founders %in% fkeep)]
      fdrop=fdrop[fdrop!="B73_inra"]
      nacol=data.frame(matrix(ncol=length(fdrop),nrow=1))
      names(nacol)=fdrop
      new_gwas=cbind(new_gwas,nacol)
      new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
      new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',"B73_inra","PC1","PC2","PC3",founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
    }else{
      end=10+length(fkeep)-2
      betas=gwas[,c(6,10:end)]
      betas=unlist(unname(betas))
      betas[-1]=betas[1]+betas[-1]
      names(betas)=fkeep
      new_gwas=gwas[,1:5]
      new_gwas=cbind(new_gwas,betas[1])
      new_gwas=cbind(new_gwas,gwas[,7:9])
      new_gwas=cbind(new_gwas,gwas[,10:end])
      names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra',"PC1","PC2","PC3",fkeep[-1])
      fdrop=founders[!(founders %in% fkeep)]
      nacol=data.frame(matrix(ncol=length(fdrop),nrow=1))
      names(nacol)=fdrop
      new_gwas=cbind(new_gwas,nacol)
      new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
      new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',"B73_inra","PC1","PC2","PC3",founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
    }
    all_gwas=rbind(all_gwas,new_gwas)
  }
}

all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
fwrite(all_gwas,sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results.txt',time,chr),row.names=F,quote=F,sep='\t')


