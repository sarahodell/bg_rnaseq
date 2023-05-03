#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
chr=as.character(args[[2]])
cores=as.numeric(args[[3]])


library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('data.table')
library('dplyr')
library('preprocessCore')
library('stringr')
library('parallel')
library('MASS')

# Read in Kinship Matrix
K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
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
metadata=metadata[metadata$read==1,]
#data = phenotypes[,c('V1'),drop=F]
#names(data)=c('ID')

genos=phenotypes$V1
######

X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])
#dhs=metadata[match(samples,metadata$sample_name),]$dh_genotype
inter=intersect(genos,inds)

genes=intersect(genes,names(phenotypes)[-1])
rownames(phenotypes)=phenotypes$V1
phenotypes=phenotypes[,-1]
phenotypes=as.matrix(phenotypes)
phenotypes=phenotypes[inter,]
genos=rownames(phenotypes)

K=K[inter,inter]


founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

n_reps=seq(1,length(genes))
#n_reps=seq(1,5)

transeqtl_gwas=function(rep){
	all_gwas=data.frame(matrix(ncol=29,nrow=0))
	names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
	gene=genes[rep]
	gene_chrom=genetable[genetable$Gene_ID==gene,]$CHROM
	data=data.frame(ID=rownames(phenotypes),y=phenotypes[,gene],stringsAsFactors=F)
	data=data[!is.na(data$y),]
	data$PC1=pcs[match(data$ID,pcs$V1),]$PC1
	data$PC2=pcs[match(data$ID,pcs$V1),]$PC2
	data$PC3=pcs[match(data$ID,pcs$V1),]$PC3
	rownames(data)=data$ID
	data=data[inter,]
	data$ID2=data$ID
	data=data[,c('ID','ID2','y','PC1','PC2','PC3')]
	null_model = GridLMM_ML(y~1+PC1+PC2+PC3+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)
	h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
	names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
	h2_start
	V_setup=null_model$setup
	Y=as.matrix(data$y)
	X_cov=null_model$lmod$X
	X_list_null=NULL
	nmarkers=dim(X_list[[1]])[2]
	frep2=sapply(seq(1,nmarkers),function(i) lapply(X_list,function(j) sum(j[,i]>0.75)))
	founders=names(X_list)
	fkeep=apply(frep2,MARGIN=2,function(x) x>2)
	markers=dimnames(X_list[[1]])[[2]]
	colnames(fkeep)=markers
	colnames(frep2)=markers
	fgroups=unique(colSums(fkeep))
	
	for(g in fgroups){
		subm=colnames(fkeep[,colSums(fkeep)==g])
		if(chr==gene_chrom){
  			snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
  			if(sum(snp %in% subm)!=0){
  				subm=subm[subm!=snp]
  			}
  		}
  		subfkeep=fkeep[,subm]
  		X_list_sub=lapply(X_list,function(x) x[inter,subm])
  		if(g==16){
      		gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
      		gwas$Trait=gene
      		names(gwas)[c(6,10:24)]=founders
      		names(gwas)[7:9]=c('PC1','PC2','PC3')
      		gwas=gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
      		all_gwas=rbind(all_gwas,gwas)
      	}else{
      		pattern=apply(subfkeep,MARGIN=2,function(x) str_flatten(c(unlist(founders[x])),'-'))
    		fdf=data.frame(marker=subm,fpattern=pattern,stringsAsFactors=F)
    		fpatterns=unique(fdf$fpattern)
      		for(i in fpatterns){
        		subm2=fdf[fdf$fpattern==i,]$marker
        		subf=subfkeep[,subm2,drop=F]
        		fk=founders[subf[,1]]
        		nfk=founders[!subf[,1]]
        		X_list_sub2=X_list_sub[ - which(names(X_list_sub) %in% nfk)]
        		X_list_sub2=lapply(X_list_sub2,function(x) x[,subm2,drop=F])

        		gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub2[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
        		gwas$Trait=gene
        		if(!("B73_inra" %in% fk)){
          			end=10+length(fk)-2
          			new_gwas=gwas[,1:5]
          			ncol1=data.frame(matrix(ncol=1,nrow=nrow(gwas)))
          			names(ncol1)='B73_inra'
          			new_gwas=cbind(new_gwas,ncol1)
         	 		new_gwas=cbind(new_gwas,gwas[,7:9])
          			new_gwas=cbind(new_gwas,gwas[,c(6,10:end)])
          			names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',fk)
          			fdrop=nfk[nfk!="B73_inra"]
          			nacol=data.frame(matrix(ncol=length(fdrop),nrow=nrow(gwas)))
          			names(nacol)=fdrop
          			new_gwas=cbind(new_gwas,nacol)
          			new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
          			new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
        		}else{
          			end=10+length(fk)-2
          			new_gwas=gwas[,1:5]
          			new_gwas=cbind(new_gwas,gwas[,6])
          			new_gwas=cbind(new_gwas,gwas[,7:9])
          			new_gwas=cbind(new_gwas,gwas[,10:end])
          			names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',fk[-1])
          			nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
          			names(nacol)=nfk
          			new_gwas=cbind(new_gwas,nacol)
          			new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
          			new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
          		}
        		all_gwas=rbind(all_gwas,new_gwas)
        	}
        }
    }
	all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
	tmp=all_gwas
	return(tmp)
}



#fwrite(all_gwas,sprintf('eqtl/trans/results/%s_c%s_pheno_%s_trans_results.txt',time,chr,factor),row.names=F,quote=F,sep='\t')


print(system.time({
results=mclapply(n_reps,transeqtl_gwas,mc.cores=cores)
}))

saveRDS(results,sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_fkeep_results.rds',chr,time))


