#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
factors=as.numeric(args[[2]])
cores=as.numeric(args[[3]])
library('data.table')
library('GridLMM')
library('ggplot2')

time="WD_0712"
factor=6
cores=2

qtl=fread(sprintf('eqtl/results/Factor%.0f_pheno_trans_%s_eQTL_hits.txt',factor,time),data.table=F)
#snp="AX-9180836"
gmap=c()
for(c in 1:10){
  g=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_gmap_c%.0f.csv',c),data.table=F)
  gmap=rbind(gmap,g)
}

test_snps=c()
#threshold=0.9
for(c in 1:10){
  sub=sub=qtl[qtl$CHR==c,]
  start_ind=which.max(sub$value)
  row=sub[start_ind,]
  snp=row$SNP
  sub=sub[sub$SNP!=snp,]
  snp_cM=gmap[gmap$marker==snp,]$pos
  test_snps=c(test_snps,snp)
  while(nrow(sub)!=0){
    sub$cM=gmap[match(sub$SNP,gmap$marker),]$pos
    sub=sub[abs(snp_cM-sub$cM)>10,]
    start_ind=which.max(sub$value)
    row=sub[start_ind,]
    snp=row$SNP
    test_snps=c(test_snps,snp)
    sub=sub[sub$SNP!=snp,]
    snp_cM=gmap[gmap$marker==snp,]$pos
  }
}

kept=qtl[qtl$SNP %in% test_snps,]
kept$cM=gmap[match(kept$SNP,gmap$marker),]$pos
fwrite(kept,sprintf('allelic/%s_Factor%.0f_kept_SNPS.txt',time,factor),row.names=F,quote=F,sep='\t')
#qtl=fread(sprintf('eqtl/results/%s_cis_eQTL_hits.txt',time),data.table=F)
phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time),data.table=F)
genos=metadata[match(pcs$V1,metadata$sample_name),]$dh_genotype
#pcs$ID=genos
#pcs=pcs[,c('ID','PC1','PC2','PC3')]

samples=phenotypes$ID

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

genes=names(phenotypes)[-1]

all_effect_sizes=c()

for(q in 1:nrow(kept)){
  row=qtl[q,]
  chr=as.character(row$CHR)
  snp=row$SNP

  K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)

  X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
  inds=rownames(X_list[[1]])
  dhs=metadata[match(samples,metadata$sample_name),]$dh_genotype

  i=intersect(rownames(K),phenotypes$ID)
  K=K[i,i]

  #pcs=pcs[pcs$ID %in%i,]
  phenotypes=phenotypes[phenotypes$ID %in% i,]
  #all_gwas=data.frame(matrix(ncol=26,nrow=length(genes)))
  #names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')

  all_res=data.frame(matrix(ncol=16,nrow=length(genes)))
  names(all_res)=founders
  rownames(all_res)=genes
  for(g in 1:length(genes)){
    gene=genes[g]
    data=data.frame(ID=phenotypes$ID,y=phenotypes[,gene],stringsAsFactors=F)
    data=data[!is.na(data$y),]
    data=as.data.frame(data,stringsAsFactors=F)
    rownames(data)=data$ID
    data=data[i,]
    data$ID2=data$ID
    # variance stabilize
    data$y=(data$y-mean(data$y))/sd(data$y)
    data=data[,c('ID','ID2','y')]
      #rownames(data)=data$ID
    null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    Y=as.matrix(data$y)
    X_cov=null_model$lmod$X

    #X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
    #colnames(X) = founders
    #rownames(X) = dimnames(X_list[[1]])[[1]]
    #X=X[i,]
    X_list_ordered=lapply(X_list,function(x) x[i,snp,drop=F])

    X_list_null=NULL

    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
    gwas$Trait=gene
    betas=gwas[,6:21]
    betas=unlist(unname(betas))
    betas[-1]=betas[1]+betas[-1]
    names(betas)=founders
    all_res[gene,]=betas
    #all_gwas[g,]=gwas
  }

  all_res=as.data.frame(all_res,stringsAsFactors=F)
  genenames=rownames(all_res)
  all_res$SNP=snp
  all_res$factor=row$Trait
  all_effect_sizes=rbind(all_effect_sizes,all_res)
  #all_res$ID=i
  #all_res=all_res[,c('ID',genenames)]


}
all_effect_sizes=as.data.frame(all_effect_sizes,stringsAsFactors=F)
fwrite(all_effect_sizes,sprintf('eqtl/results/%s_factor%.0f_trans_peak_effect_sizes.txt',time,factor),row.names=F,quote=F,sep='\t')

snp=kept$SNP[1]
sub=all_effect_sizes[all_effect_sizes$SNP==snp,c(founders)]

pca=prcomp(t(sub), center = TRUE)
pcs=as.data.frame(pca$x,stringsAsFactors=F)
pcs$founder=rownames(pcs)

var_explained <- pca$sdev^2/sum(pca$sdev^2)

p1=ggplot(aes(x=PC1,y=PC2),data=pcs) + geom_point() + ggtitle(sprintf('PCA Expression Effect Sizes for SNP %s',snp))
png(sprintf('eqtl/trans/images/%s_Factors%.0f_Gene_PCA_SNP_%s.png',time,factor,snp))
print(p1)
dev.off()

#all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
#fwrite(all_res,sprintf('eqtl/results/single_trans_eQTL_%s_%s_betas.txt',time,factor),row.names=F,quote=F,sep='\t')
#fwrite(all_gwas,sprintf('eqtl/results/single_trans_eQTL_%s_%s_results.txt',time,factor),row.names=F,quote=F,sep='\t')
