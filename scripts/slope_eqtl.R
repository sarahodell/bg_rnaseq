#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])
cores=as.numeric(args[[2]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('data.table')
library('dplyr')
#library('lme4')
#library('lme4qtl')
library('preprocessCore')
library('stringr')
library('parallel')
library('MASS')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
genes=unique(genetable$Gene_ID)

xs=c(18,20,27)
phenotypes=fread('eqtl/expression_slopes.txt',data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)

#data = phenotypes[,c('V1'),drop=F]
#names(data)=c('ID')

genos=phenotypes$ind
######
X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])

inter=intersect(genos,inds)
K=K[inter,inter]

genes=unique(phenotypes$gene)



n_reps=seq(1,length(genes))
#n_reps=seq(1,5)

transeqtl_gwas=function(rep){
	gene=genes[rep]
	gene_chrom=genetable[genetable$Gene_ID==gene,]$CHROM

	data=phenotypes[phenotypes$gene==gene,]
	rownames(data)=data$ind
	data=data[inter,]
	
  	#data=data.frame(ID=rownames(phenotypes),y=phenotypes[,gene],stringsAsFactors=F)
  	data=data[!is.na(data$slope),]


  	data$ID2=data$ind
  	names(data)=c('ID','gene','WD_0718','WD_0720','WD_0727','y','ID2')
  	data=data[,c('ID','ID2','y')]
  	# sd-scale slopes?
  	#data$y=(data$y-mean(data$y))/sd(data$y)
  
  	null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

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

	all_gwas=data.frame(matrix(ncol=26,nrow=0))
  	names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')

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
      		names(gwas)[c(6:21)]=founders
      		gwas=gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
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
          			end=6+length(fk)-1
          			new_gwas=gwas[,1:5]
          			ncol1=data.frame(matrix(ncol=1,nrow=nrow(gwas)))
          			names(ncol1)='B73_inra'
          			new_gwas=cbind(new_gwas,ncol1)
          			new_gwas=cbind(new_gwas,gwas[,6:end])
          			names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra',fk)
          			fdrop=nfk[nfk!="B73_inra"]
          			nacol=data.frame(matrix(ncol=length(fdrop),nrow=nrow(gwas)))
          			names(nacol)=fdrop
          			new_gwas=cbind(new_gwas,nacol)
          			new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
          			new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
        		}else{
          			end=6+length(fk)-1
          			new_gwas=gwas[,1:5]
          			new_gwas=cbind(new_gwas,gwas[,6:end])
          			names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',fk)
          			nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
          			names(nacol)=nfk
          			new_gwas=cbind(new_gwas,nacol)
          			new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
          			new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
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

saveRDS(results,sprintf('eqtl/trans/results/slope_eQTL_c%s_fkeep_results.rds',chr))






























































































































