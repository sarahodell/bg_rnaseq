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


options(warn=2)
# Read in Kinship Matrix
K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

genetable=fread('eqtl/data/Zea_mays.B73_v4_generanges.txt',data.table=F)
genetable=genetable[genetable$TXCHROM==chr,]
genes=unique(genetable$Gene_ID)
# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread(sprintf('star/vst_%s_gene_counts_results2.csv',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
samples=colnames(phenotypes)[-1]
genos=metadata[match(samples,metadata$sample_name),]$genotype
testsnps=readRDS(sprintf('eqtl/gene_focal_snps_c%s.rds',chr))
founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)


genes=intersect(genes,phenotypes$V1)
#genes=genes[1:10]
#ID      BGA_ID  CODE_HYBRIDE
#linenames=fread('metadata/Lines_Names_BALANCE.txt',data.table=F)
#dh_genotype=linenames[match(metadata$genotype,linenames$CODE_HYBRIDE),]$BGA_ID
#missing=unique(metadata$genotype[is.na(dh_genotype)]) # 12 genotypes that are not in the file
# "9RGTCONEXXION" "9P9486"        "9A277-H3"      "9LG30315"
# "9DK440"        "9P9400"        "9DKC4117"      "9LG30500"
# "9PR38N86"      "9LG30444"      "9DKC4490"      "9ESGALLERY"
#metadata$dh_genotype=dh_genotype
all_gwas=c()
all_res=c()
for(g in 1:length(genes)){
  gene=genes[g]
  snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
  if(length(snp)!=0){
    data=data.frame(sample=colnames(phenotypes)[-1],y=unlist(phenotypes[phenotypes[,1]==gene,][-1]),stringsAsFactors=F)
    data=data[!is.na(data$y),]
    data$ID=metadata[match(data$sample,metadata$sample_name),]$dh_genotype
    rownames(data)=seq(1,nrow(data))
    data=data[!is.na(data$ID),]
    X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
    inds=rownames(X_list[[1]])

    i=intersect(data$ID,inds)

    founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
    K=K[i,i]
    if(length(unique(data$ID))<nrow(data)){
      data=data%>%group_by(ID)%>%summarize(y=mean(y))
    }
    data=as.data.frame(data,stringsAsFactors=F)
    rownames(data)=data$ID
    data=data[i,]
    data$ID2=data$ID
    data=data[,c('ID','ID2','y')]
      #rownames(data)=data$ID
    null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    Y=as.matrix(data$y)
    X_cov=null_model$lmod$X
    X_list_ordered=lapply(X_list,function(x) x[i,snp,drop=F])

    X_list_null=NULL

    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
    gwas$Trait=gene

    if(length(snp)>1){
        betas=unlist(unname(gwas[2,6:21]))
        X = do.call(cbind,lapply(X_list,function(x) x[,snp[2]]))
        colnames(X) = founders
        rownames(X) = dimnames(X_list[[1]])[[1]]
        snp=snp[2]
        X=X[i,]
    }else{
      betas=unlist(unname(gwas[6:21]))
      X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
      colnames(X) = founders
      rownames(X) = dimnames(X_list[[1]])[[1]]
      X=X[i,]
    }
    #frep2=apply(X,MARGIN=2,function(x) sum(x[x>0.8]))
    #fkeep=founders[frep2>2]
    #X=X[,fkeep]
    # drop one of them to test
    # to get residuals, get difference of fitted values

    betas[-1]=betas[1]+betas[-1]
    X_r=X %*% betas
    residuals=data$y - X_r
    residuals=as.data.frame(residuals,stringsAsFactors=F)
    names(residuals)=gene
    #saveRDS(gwas,sprintf('results/eqtl_gridlmm_chr%s_%s_x_%s_results.rds',chr,time,env))


    #m4 = relmatLmer(y ~ 0 + X + (1|ID),data=data,relmat = list(ID=K))
    #residuals=as.data.frame(summary(m4)$residuals,stringsAsFactors=F)
    #names(residuals)=gene
      #se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
      #names(se4)=c('value','se','tvalue')
      #rownames(se4)=founders
      #se4$founder=rownames(se4)

      #se4$variable_f=factor(se4$founder,levels=se4$founder)
    if(is.null(all_res)){
      all_res=residuals
      all_gwas=gwas
      #all_betas=c(gene,snp,betas)
    }else{
      all_res=cbind(all_res,residuals)
      all_gwas=rbind(all_gwas,gwas)
      #all_betas=rbind(all_betas,c(gene,snp,betas))
    }
  }
}

all_res=as.data.frame(all_res,stringsAsFactors=F)
genenames=names(all_res)
all_res$ID=i
#names(all_res)=c(genes,'ID')
all_res=all_res[,c('ID',genenames)]


all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
#all_betas=as.data.frame(all_betas,stringsAsFactors=F)
#names(all_betas)=c('Gene_ID','SNP',paste0('beta.',seq(1,16)))

fwrite(all_res,sprintf('eqtl/results/eQTL_%s_c%s_residuals.txt',time,chr),row.names=F,quote=F,sep='\t')
fwrite(all_gwas,sprintf('eqtl/results/eQTL_%s_c%s_results.txt',time,chr),row.names=F,quote=F,sep='\t')


    # Run GridLMM
    #null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

    #h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    #names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    #h2_start
    #V_setup=null_model$setup
    #Y=as.matrix(data$y)
    #X_cov=null_model$lmod$X
    #X_list_ordered=lapply(X_list,function(x) x[i,snp,drop=F])

    #X_list_null=NULL

    #gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
    # save residuals to a file
    #saveRDS(gwas,sprintf('results/eqtl_gridlmm_chr%s_%s_x_%s_results.rds',chr,time,env))
    #saveRDS(gwas_adjusted,sprintf('models/Biogemma_chr%s_%s_x_%s_founderprobs_adjusted.rds',chr,pheno,env))
    #hinfo=data.frame(method="Founder_probs",phenotype=pheno,environment=env,chr=chr,h2=h2_start,hap=NA,stringsAsFactors=F)
    #fwrite(hinfo,'../heritabilities.txt',quote=F,sep='\t',row.names=F,append=T)

    #rbind all results together for all genes
    #keep track of number of tests


    #betas[-1]=betas[1]+betas[-1]
    #X_r=X %*% betas

    #resid_model = GridLMM_ML(y~0 + X_r+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)
    #h2_start=resid_model$results[,grepl('.ML',colnames(resid_model$results),fixed=T),drop=FALSE]
    #names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    #h2_start
    #V_setup=resid_model$setup
    #Y=as.matrix(data$y)
    #X_cov=resid_model$lmod$X
    #X_list_ordered=lapply(X_list,function(x) x[i,snp,drop=F])
    #X_list_null=NULL

    #gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)

    #if(sum(is.na(betas))<3){
      # should just use X , not X_r from GridLMM


#fwrite(all_res,sprintf('eqtl/results/cis_eQTL_%s_c%s_residuals.txt',time,chr),row.names=F,quote=F,sep='\t')
#gtf <- makeTxDbFromGFF(file="/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.chr.gtf",format="gtf")
#allgenes=keys(gtf,keytype="TXNAME")
#gtable=select(gtf,keys=allgenes,columns=c("TXCHROM","TXNAME","TXSTART","TXEND"),keytype="TXNAME")
#gtable$Gene_ID=sapply(seq(1,nrow(gtable)),function(x) strsplit(gtable$TXNAME[x],"_T00")[[1]][[1]])
#fwrite(gtable,'metadata/Zea_mays.B73_v4_generanges.txt',row.names=F,quote=F,sep='\t')
#genes=phenotypes$V1 # list of unique genes

# Make a table of founder prob sites and their ranges
#chroms=c("1","2","3","4","5","6","7","8","9","10")
#for(chr in chroms){
#  X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
#  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
#  snps=colnames(X_list[[1]])
#  founder_blocks=c()
#  for(snp in snps){
#    m=pmap[pmap$marker==snp,]$pos
#    fdropped=readRDS(sprintf('../genotypes/probabilities/geno_probs/dropped/bg%s_dropped_markers_genoprobs.rds',chr))
#    dropped_markers=fdropped[[which(unlist(lapply(fdropped,function(x) x$marker==snp)))]]$linked
#    if(!is.null(dropped_markers)){
#      sub=pmap[pmap$marker %in% dropped_markers,]
#      rownames(sub)=seq(1,dim(sub)[1])
#      right=which.max(sub$pos)
#      left_bound=m
#      right_bound=sub[right,]$pos
#    }else{
#      left_bound=m
#      right_bound=m+1
#    }
#    if(snp==snps[1]){
#      founder_blocks=c(chr,snp,left_bound,right_bound)
#    }else{
#      founder_blocks=rbind(founder_blocks,c(chr,snp,left_bound,right_bound))
#    }
#  }
#  founder_blocks=as.data.frame(founder_blocks,stringsAsFactors=F)
#  names(founder_blocks)=c("chr","focal_snp","start","end")
#  founder_blocks$chr=as.numeric(founder_blocks$chr)
#  founder_blocks$start=as.numeric(founder_blocks$start)
#  founder_blocks$end=as.numeric(founder_blocks$end)
#  fwrite(founder_blocks,sprintf('metadata/founder_recomb_blocks_c%s.txt',chr),row.names=F,quote=F,sep='\t')
#}

# grab founder recombination block that encompasses the gene - overlap so may include more than one...
#for(chr in chroms){
#  genetable=fread('eqtl/data/Zea_mays.B73_v4_generanges.txt',data.table=F)#
#  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
#  genetable=genetable[genetable$TXCHROM==chr,]
#  genes=unique(genetable$Gene_ID)
#  testsnps=vector("list",length=length(genes))
#  for(i in 1:length(genes)){
#    gene=genes[i]
#    subtable=as.data.table(genetable[genetable$Gene_ID==gene,])
#    env2=as.data.table(founder_blocks)
#    setkey(subtable,TXSTART,TXEND)
#    comparison=foverlaps(env2,subtable,by.y=c('TXSTART','TXEND'),by.x=c('start','end'),nomatch=NA)
#    comparison=comparison[!is.na(comparison$Gene_ID),]
#    snps=unique(comparison$focal_snp)
#    testsnps[[i]]=list(gene=gene,focal_snps=snps)
#  }
#  saveRDS(testsnps,sprintf('eqtl/gene_focal_snps_c%s.rds',chr))
#}

  # find location of gene in annotation?

  # for each site in fp, test association between gene expression and site



#write output to file
#what is my significance threshold?

#phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
#phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
#phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]
