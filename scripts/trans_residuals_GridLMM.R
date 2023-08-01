#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
factor=as.character(args[[2]])
chr=as.character(args[[3]])
cores=as.numeric(args[[4]])

library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('data.table')
library('dplyr')
library('parallel')
library('MASS')
library('stringr')

#all_reps=c()

#drop=list('Factor19'="EB.09S.H.00429",'Factor38'="EB.09S.H.00424",
#'Factor48'="EB.09S.H.00494",'Factor49'="EB.09S.H.00402")

phenotype=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_all_F_means_FIXED.txt',time),data.table=F)

if(time=="WD_0712" & factor=="Factor5"){
	drop_ind="EB.09S.H.00417"
	phenotype[phenotype$V1==drop_ind,factor]=-1.4175
}

if(time=="WD_0712" & factor=="Factor23"){
	drop_ind="EB.09S.H.00005"
	phenotype[phenotype$V1==drop_ind,factor]=1.9960830
}

#factor_results=c()
#p=names(phenotype)[factor]
#metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
metadata=fread('metadata/samples_passed_genotype_check.txt',data.table=F)
metadata=metadata[metadata$experiment==time,]
plate=metadata[match(phenotype$V1,metadata$dh_genotype),]$plate

data=data.frame(ID=phenotype$V1,ID2=phenotype$V1,plate=plate,y=phenotype[,factor],stringsAsFactors=F)
data=data[!is.na(data$y),]
data$plate=as.factor(data$plate)

#data$y=(data$y-mean(data$y))/sd(data$y)
  #K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

adj_chr=c("5","9")
if(chr %in% adj_chr){
	X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))

}else{
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
}
#if(factor %in% names(drop)){
#  drop_ind=drop[factor]
#  X_list=lapply(X_list,function(x) x[rownames(x)!=drop_ind,])
#}

inds=rownames(X_list[[1]])
inter=intersect(inds,data$ID)
rownames(data)=data$ID
data=data[inter,]
K=K[inter,inter]
X_list=lapply(X_list,function(x) x[inter,])

null_model = GridLMM_ML(y~1+plate+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

nmarkers=dim(X_list[[1]])[2]
frep2=sapply(seq(1,nmarkers),function(i) lapply(X_list,function(j) sum(j[,i]>0.75)))
founders=names(X_list)
fkeep=apply(frep2,MARGIN=2,function(x) x>3)
markers=dimnames(X_list[[1]])[[2]]
colnames(fkeep)=markers
colnames(frep2)=markers
fgroups=unique(colSums(fkeep))
plates=paste0('plate',unique(data$plate))

dimname=c('Trait','X_ID','s2','ML_logLik','ID.ML',plates,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
all_gwas=data.frame(matrix(ncol=length(dimname),nrow=0))
names(all_gwas)=dimname

for(g in fgroups){
	subm=colnames(fkeep[,colSums(fkeep)==g,drop=F])
	subfkeep=fkeep[,subm,drop=F]
	X_list_sub=lapply(X_list,function(x) x[,subm,drop=F])
	if(g==16){
		h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
		names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
		h2_start
		V_setup=null_model$setup
		Y=as.matrix(data$y)
		X_cov=null_model$lmod$X
		X_list_null=NULL
		gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
		gwas$Trait=factor
		endp=6+length(plates)-1
		end=endp+15
		names(gwas)[6:endp]=plates
		names(gwas)[(endp+1):end]=founders[-1]
		gwas[,plates[-1]]=gwas[,plates[1]] + gwas[,plates[-1]]
		gwas$B73_inra=apply(gwas[,plates],MARGIN=1,mean)
		gwas=gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',plates,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
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
			h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
			names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
			h2_start
			V_setup=null_model$setup
			Y=as.matrix(data$y)
			X_cov=null_model$lmod$X
			X_list_null=NULL
			gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub2[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
			gwas$Trait=factor
			if(!("B73_inra" %in% fk)){
				endp=6+length(plates)-1
				end=endp+(length(fk)-1)
				names(gwas)[6:endp]=plates
				names(gwas)[(endp+1):end]=fk[-1]
				gwas[,plates[-1]]=gwas[,plates[1]] + gwas[,plates[-1]]
				gwas[,fk[1]]=apply(gwas[,plates],MARGIN=1,mean)
				end=endp+length(fk)-1
				new_gwas=gwas
				#names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',plates,fk[-1])
				nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
				names(nacol)=nfk
				new_gwas=cbind(new_gwas,nacol)
				#new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
				new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',plates,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
			}else{
				endp=6+length(plates)-1
				end=endp+(length(fk)-1)
				names(gwas)[6:endp]=plates
				names(gwas)[(endp+1):end]=fk[-1]
				gwas[,plates[-1]]=gwas[,plates[1]] + gwas[,plates[-1]]
				gwas$B73_inra=apply(gwas[,plates],MARGIN=1,mean)
				end=endp+length(fk)-1
				new_gwas=gwas
				#names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',plates,fk[-1])
				nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
				names(nacol)=nfk
				new_gwas=cbind(new_gwas,nacol)
				#new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
				new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',plates,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
        	}
      		all_gwas=rbind(all_gwas,new_gwas)
    	}
    }
}
all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
fwrite(all_gwas,sprintf('eqtl/trans/results/%s_residuals_c%s_%s_trans_results_FIXED.txt',time,chr,factor),row.names=F,quote=F,sep='\t')







#fwrite(factor_results,sprintf('eqtl/trans/results/%s_c%s_pheno_factor_trans_eQTL.txt',time,chr),row.names=F,quote=F,sep='\t')
#fwrite(factor_results,sprintf('eqtl/trans/results/%s_c%s_factor_trans_eQTL.txt',time,chr),row.names=F,quote=F,sep='\t')

# For each snp, figure out fkeep and separate X_list into multiple based on dimensions
# parallelize


#n_reps=seq(1,reps)

#randomized_gwas<-function(rep){
#   len=dim(X_list_full[[1]])[1]

   # Run GridLMM

   # randomize the order of the genotypes
#   draw=sample(len,len,replace=F)
#   X_list_reordered=lapply(X_list_full,function(x) x[draw,])
#   for(x in seq(1,16)){
#       dimnames(X_list_reordered[[x]])[[1]]=dimnames(X_list_full[[1]])[[1]]
#   }

#   h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
#   names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
#   h2_start
#   V_setup=null_model$setup

#   Y=as.matrix(data$y)
#   X_cov=null_model$lmod$X
#   X_list_null=NULL

#   gwas=run_GridLMM_GWAS(Y,X_cov,X_list_reordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
#   gwas=gwas[!is.na(gwas$p_value_ML),]
#   tmp=data.frame(chr=chr,replicate=rep,pval=min(gwas$p_value_ML))
#}

#print(system.time({
#results=mclapply(n_reps,randomized_gwas,mc.cores=cores)
#}))

#saveRDS(results,sprintf('test_models/chr%s_%s_x_%s_founderprobs_%.0frep_max_pvalues.rds',chr,pheno,env,reps))

#date=format(Sys.time(),'%m%d%y')


#time="WD_0720"
#chr="10"
#library('GenomicFeatures') # write a script to get a table of start and stop sites of genes from the gtf file


#options(warn=2)
# Read in Kinship Matrix
