args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
chr=as.character(args[[2]])
cores=as.numeric(args[[3]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM')
library('data.table')
library('dplyr')
library('lme4')
library('lme4qtl')

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


# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
#phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts.txt',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time),data.table=F)
samples=pcs$sample
genos=metadata[match(samples,metadata$sample_name),]$dh_genotype
pcs$ID=genos
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
pc_list=c("PC1","PC2","PC3")

heritability=c()

for(k in 1:length(pc_list)){
  pc=pc_list[k]

  data=data.frame(ID=pcs$ID,ID2=pcs$ID,y=pcs[,pc],stringsAsFactors=F)
  data=data[!is.na(data$y),]
  #data$y=(data$y-mean(data$y))/sd(data$y)
  #K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
  K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)

  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)

  X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
  inds=rownames(X_list[[1]])
  i=intersect(inds,data$ID)
  K=K[i,i]
  if(length(unique(data$ID))<nrow(data)){
    data=data%>%group_by(ID)%>%summarize(y=mean(y))
  }
  data=as.data.frame(data,stringsAsFactors=F)
  rownames(data)=data$ID
  data=data[i,]
  null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

  h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
  names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
  h2_start
  V_setup=null_model$setup
  Y=as.matrix(data$y)
  X_cov=null_model$lmod$X

  X_list_ordered=lapply(X_list,function(x) x[i,])
  X_list_null=NULL

  gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
  gwas$Trait=pc
  h2=gwas$setup$h2_start
  hinfo=data.frame(method="Founder_probs",phenotype=pc,environment=time,chr=chr,h2=h2,stringsAsFactors=F)
  fwrite(hinfo,'eqtl/results/PC_heritabilities.txt',quote=F,sep='\t',row.names=F,append=T)

  fwrite(gwas,sprintf('eqtl/results/%s_c%s_%s_eQTL.txt',time,chr,pc),row.names=F,quote=F,sep='\t')
}


#fwrite(factor_results,sprintf('eqtl/trans/results/%s_c%s_factor_trans_eQTL.txt',time,chr),row.names=F,quote=F,sep='\t')
