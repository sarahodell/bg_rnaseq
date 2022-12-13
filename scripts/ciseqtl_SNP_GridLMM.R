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

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
genetable=genetable[genetable$CHROM==chr,]
genes=unique(genetable$Gene_ID)
# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts.txt',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
samples=colnames(phenotypes)[-1]

genes=intersect(genes,phenotypes$V1)

X=fread(sprintf('../genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F,stringsAsFactors=F)
rownames(X)=X$ind
X=X[,-1]
inds=rownames(X)

dhs=metadata[match(samples,metadata$sample_name),]$dh_genotype
i=intersect(dhs,inds)
K=K[i,i]
X=X[i,]

pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)

all_gwas=c()

for(g in 1:length(genes)){
  gene=genes[g]
  gene_start=genetable[genetable$Gene_ID==gene,]$START
  gene_end=genetable[genetable$Gene_ID==gene,]$END
  inside=which(abs(gene_start-pmap$pos)<1000|abs(pmap$pos-gene_end)<1000|(pmap$pos>gene_start & pmap$pos<gene_end))
  snp=pmap$marker[inside]
  test=X[,names(X)[2],drop=F]
  test=test[i,,drop=F]
  snp=snp[snp %in% colnames(X)]

  #snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
  if(length(snp)!=0){
    subX=X[,snp,drop=F]
    data=data.frame(sample=colnames(phenotypes)[-1],y=unlist(phenotypes[phenotypes[,1]==gene,][-1]),stringsAsFactors=F)
    data=data[!is.na(data$y),]
    data$ID=metadata[match(data$sample,metadata$sample_name),]$dh_genotype
    rownames(data)=seq(1,nrow(data))
    data=data[!is.na(data$ID),]

    if(length(unique(data$ID))<nrow(data)){
      data=data%>%group_by(ID)%>%summarize(y=mean(y))
    }
    data=as.data.frame(data,stringsAsFactors=F)
    rownames(data)=data$ID
    data=data[i,]
    data$ID2=data$ID

    i=intersect(rownames(X),data$ID)
    subX=subX[i,,drop=F]

    subX=as.matrix(subX)
    if(length(snp)==1){
      subX=cbind(subX,test)
    }
    #subX=
    data$y=(data$y-mean(data$y))/sd(data$y)
    data=data[,c('ID','ID2','y')]

    gwas = GridLMM_GWAS(
                            formula = y~1+(1|ID),
                            test_formula = ~1,
                            reduced_formula = ~0,
                            data = data,
                            weights = NULL,
                            X = subX,
                            X_ID = 'ID',
                            h2_start = NULL,
                            h2_step = 0.01,
                            max_steps = 100,
                            relmat = list(ID=K),
                            centerX = FALSE,
                            scaleX = FALSE,
                            fillNAX = FALSE,
                            method = 'REML',
                            mc.cores = cores,
                            verbose = FALSE
    )

    gwas_res=gwas$results
    if(length(snp)==1){
      gwas_res=gwas_res[gwas_res$X_ID!="AX-91157322",]
    }
    gwas_res$Trait=gene
    all_gwas=rbind(all_gwas,gwas_res)

  }
}


all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','beta.1','beta.2','REML_logLik','ID.RMEL','F,1','n_steps','Df_X','p_value_REML')

fwrite(all_gwas,sprintf('eqtl/cis/results/eQTL_%s_c%s_SNP_results.txt',time,chr),row.names=F,quote=F,sep='\t')
