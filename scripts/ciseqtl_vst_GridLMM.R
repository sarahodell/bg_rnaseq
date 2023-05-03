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
#data = phenotypes[,c('V1'),drop=F]
#names(data)=c('ID')

genos=phenotypes$V1
######

X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])
#dhs=metadata[match(samples,metadata$sample_name),]$dh_genotype
i=intersect(genos,inds)

genes=intersect(genes,names(phenotypes)[-1])
rownames(phenotypes)=phenotypes$V1
phenotypes=phenotypes[,-1]
phenotypes=as.matrix(phenotypes)
phenotypes=phenotypes[i,]
genos=rownames(phenotypes)

K=K[i,i]

plate=metadata[match(genos,metadata$dh_genotype),]$plate
df=data.frame(ID=genos,plate=plate,stringsAsFactors=F)
df$plate=as.factor(df$plate)
#Separate individuals by plate and quantile normalize separately
#Ynorm=c()
#plates=unique(df$plate)
#for(p in plates){
#  pinds=df[df$plate==p,]$ID
#  subY=phenotypes[pinds,]
#  subYnorm=normalize.quantiles(as.matrix(subY))
#  rownames(subYnorm)=rownames(subY)
#  colnames(subYnorm)=colnames(subY)
#  Ynorm=rbind(Ynorm,subYnorm)
#}

#Ynorm=as.matrix(Ynorm)
#phenotypes=Ynorm
#genos=metadata[match(samples,metadata$sample_name),]$genotype
testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")



# Quantile normalization of genes
#pnorm=normalize.quantiles(phenotypes)
#rownames(pnorm)=rownames(phenotypes)
#colnames(pnorm)=colnames(phenotypes)
#phenotypes=pnorm
#genes=genes[1:20]
#ID      BGA_ID  CODE_HYBRIDE
#linenames=fread('metadata/Lines_Names_BALANCE.txt',data.table=F)
#dh_genotype=linenames[match(metadata$genotype,linenames$CODE_HYBRIDE),]$BGA_ID
#missing=unique(metadata$genotype[is.na(dh_genotype)]) # 12 genotypes that are not in the file
# "9RGTCONEXXION" "9P9486"        "9A277-H3"      "9LG30315"
# "9DK440"        "9P9400"        "9DKC4117"      "9LG30500"
# "9PR38N86"      "9LG30444"      "9DKC4490"      "9ESGALLERY"
#metadata$dh_genotype=dh_genotype


nob73=c()

all_gwas=data.frame(matrix(ncol=29,nrow=length(genes)))
names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',"B73_inra","PC1","PC2","PC3",founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')

all_res=data.frame(matrix(ncol=length(genes),nrow=length(i)))
names(all_res)=genes
for(g in 1:length(genes)){
  gene=genes[g]
  snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
  if(length(snp)!=0){
    data=data.frame(ID=rownames(phenotypes),y=phenotypes[,gene],stringsAsFactors=F)
    data=data[!is.na(data$y),]
    #data$ID=metadata[match(data$sample,metadata$sample_name),]$dh_genotype
    #rownames(data)=seq(1,nrow(data))
    #data=data[!is.na(data$ID),]
    data$PC1=pcs[match(data$ID,pcs$V1),]$PC1
    data$PC2=pcs[match(data$ID,pcs$V1),]$PC2
    data$PC3=pcs[match(data$ID,pcs$V1),]$PC3
    #K=K[i,i]
    #if(length(unique(data$ID))<nrow(data)){
    #  data=data%>%group_by(ID)%>%summarize(y=mean(y),PC1=mean(PC1),PC2=mean(PC2),PC3=mean(PC3))
    #}
    #data=as.data.frame(data,stringsAsFactors=F)
    rownames(data)=data$ID
    data=data[i,]
    data$ID2=data$ID

    #plate=metadata[match(data$ID,metadata$dh_genotype),]$plate
    #data$plate=as.factor(plate)
    # variance stabilize
    data$y=(data$y-mean(data$y))/sd(data$y)
    data=data[,c('ID','ID2','y','PC1','PC2','PC3')]
      #rownames(data)=data$ID
    null_model = GridLMM_ML(y~1+PC1+PC2+PC3+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    Y=as.matrix(data$y)
    X_cov=null_model$lmod$X
    if(length(snp)>1){
        #betas=unlist(unname(gwas[2,6:21]))
        X = do.call(cbind,lapply(X_list,function(x) x[,snp[1]]))
        colnames(X) = founders
        rownames(X) = dimnames(X_list[[1]])[[1]]
        snp=snp[1]
        X=X[i,]
    }else{
      #betas=unlist(unname(gwas[6:21]))
      X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
      colnames(X) = founders
      rownames(X) = dimnames(X_list[[1]])[[1]]
      X=X[i,]
    }
    X_list_ordered=lapply(X_list,function(x) x[i,snp,drop=F])
    frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
    fkeep=founders[frep2>2]
    X_list_ordered = X_list_ordered[c(fkeep)]
    #frep2=apply(X,MARGIN=2,function(x) sum(x[x>0.8]))
    #fkeep=founders[frep2>2]
    #X_list_ordered=X_list_ordered[c(fkeep)]
    X=X[,fkeep]

    X_list_null=NULL

    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
    gwas$Trait=gene
    #betas=gwas[,c(6,10:24)]
    #betas=unlist(unname(betas))
    #betas[-1]=betas[1]+betas[-1]
    if(!("B73_inra" %in% fkeep)){
      nob73=c(nob73,gene)
      end=10+length(fkeep)-2
      betas=gwas[,c(6,10:end)]
      betas=unlist(unname(betas))
      betas[-1]=betas[1]+betas[-1]
      #names(gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',fkeep,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
      names(betas)=fkeep
      new_gwas=gwas[,1:5]
      #new_gwas=cbind(new_gwas,betas[1,"B73_inra"])
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
      #new_X=cbind(X,data[,c('PC1','PC2','PC3')])
      #new_X=new_X[,c(fkeep[1],'PC1','PC2','PC3',fkeep[-1])]
    }else{
      end=10+length(fkeep)-2
      betas=gwas[,c(6,10:end)]
      betas=unlist(unname(betas))
      betas[-1]=betas[1]+betas[-1]
      #names(gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',fkeep,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
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
      #new_X=cbind(X,data[,c('PC1','PC2','PC3')])
      #new_X=new_X[,c('B73_inra','PC1','PC2','PC3',fkeep[-1])]
    }

    # drop one of them to test
    # to get residuals, get difference of fitted values
    #betas=unlist(unname(betas))
    #betas[-1]=betas[1]+betas[-1]
    #all_betas=gwas[,c(6:end)]
    #all_betas=unlist(unname(all_betas))
    #all_betas[-1]=all_betas[1]+all_betas[-1]

    X_r=X %*% betas
    #pc_betas=gwas[,7:9]
    #pc_betas=unlist(unname(pc_betas))
    #pc_betas[-1]=pc_betas[1]+pc_betas[-1]
    #pc_X = as.matrix(data[,c('PC1','PC2','PC3')]) %*% pc_betas
    residuals=data$y - X_r
    #residuals=residuals - pc_X
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

    all_res[,gene]=residuals[,gene]
    all_gwas[g,]=new_gwas
      #all_betas=rbind(all_betas,c(gene,snp,betas))
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

fwrite(all_res,sprintf('eqtl/cis/results/eQTL_%s_c%s_vst_residuals.txt',time,chr),row.names=F,quote=F,sep='\t')
fwrite(all_gwas,sprintf('eqtl/cis/results/eQTL_%s_c%s_vst_results.txt',time,chr),row.names=F,quote=F,sep='\t')


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
# founder_blocks=as.data.frame(founder_blocks,stringsAsFactors=F)
#  names(founder_blocks)=c("chr","focal_snp","start","end")
#  founder_blocks$chr=as.numeric(founder_blocks$chr)
#  founder_blocks$start=as.numeric(founder_blocks$start)
#  founder_blocks$end=as.numeric(founder_blocks$end)
#  fwrite(founder_blocks,sprintf('metadata/founder_recomb_blocks_c%s.txt',chr),row.names=F,quote=F,sep='\t')
#}

# grab founder recombination block that encompasses the gene - overlap so may include more than one...
#for(chr in chroms){
#  genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)#
#  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
#  genetable=genetable[genetable$CHROM==chr,]
#  genes=unique(genetable$Gene_ID)
#  testsnps=vector("list",length=length(genes))
#  for(i in 1:length(genes)){
#    gene=genes[i]
#    subtable=as.data.table(genetable[genetable$Gene_ID==gene,])
#    env2=as.data.table(founder_blocks)
#    setkey(env2,start,end)
    #setkey(subtable,TXSTART,TXEND)
#    comparison=foverlaps(subtable,env2,by.x=c('START','END'),by.y=c('start','end'),nomatch=NA)
#    comparison=comparison[!is.na(comparison$Gene_ID),]
#    snps=unique(comparison$focal_snp)
#    snps=snps[!is.na(snps)]
#    if(length(snps)==0){
#      s_dist=which.min(abs(subtable$START-founder_blocks$end))
#      e_dist=which.min(abs(founder_blocks$start-subtable$END))
#      choose=which.min(c(subtable$START-founder_blocks$end[s_dist],founder_blocks$start[e_dist]-subtable$END))
#      snps=c(founder_blocks$focal_snp[s_dist],founder_blocks$focal_snp[e_dist])
#    }
#  testsnps[[i]]=list(gene=gene,focal_snps=snps)

# }
#  saveRDS(testsnps,sprintf('eqtl/gene_focal_snps_c%s.rds',chr))
#}

  # find location of gene in annotation?

  # for each site in fp, test association between gene expression and site



#write output to file
#what is my significance threshold?

#phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
#phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
#phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]
