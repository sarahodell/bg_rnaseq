#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
chr=as.character(args[[3]])
cores=as.numeric(args[[4]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('data.table')
library('dplyr')
library('stringr')

# Read in Kinship Matrix
K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

data=data.frame(ID=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
data=data[data$Loc.Year.Treat==env,]
data=data[!is.na(data$y),]
data$y = (data$y - mean(data$y))/sd(data$y)

# Read in the haplotype group probabilities
# Filter genotypes that are not in the K matrix
X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])

inter=intersect(data$ID,inds)
X_list=lapply(X_list,function(x) x[inter,])
#founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


#Make B73 the first in the list so that it is the one that is dropped
#names(X_list)=founders
#X_list=X_list[new_founders]
K=K[inter,inter]
rownames(data)=data$ID
data=data[inter,]
# Run GridLMM
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
fkeep=apply(frep2,MARGIN=2,function(x) x>3)
markers=dimnames(X_list[[1]])[[2]]
colnames(fkeep)=markers
colnames(frep2)=markers
fgroups=unique(colSums(fkeep))
	
all_gwas=data.frame(matrix(ncol=26,nrow=0))
names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')

for(g in fgroups){
	subm=colnames(fkeep[,colSums(fkeep)==g,drop=F])
	subfkeep=fkeep[,subm,drop=F]
  	X_list_sub=lapply(X_list,function(x) x[,subm,drop=F])
  	if(g==16){
  		gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
      	names(gwas)[6:21]=founders
      	all_gwas=rbind(all_gwas,gwas)
    }else{
    	pattern=apply(subfkeep,MARGIN=2,function(x) str_flatten(c(unlist(founders[x])),'-'))
    	fdf=data.frame(marker=subm,fpattern=pattern,stringsAsFactors=F)
    	fpatterns=unique(fdf$fpattern)
    	for(i in fpatterns){
      		subm2=fdf[fdf$fpattern==i,,drop=F]$marker
        	subf=subfkeep[,subm2,drop=F]
        	fk=founders[subf[,1]]
        	nfk=founders[!subf[,1]]
        	X_list_sub2=X_list_sub[ - which(names(X_list_sub) %in% nfk)]
        	X_list_sub2=lapply(X_list_sub2,function(x) x[,subm2,drop=F])
        	gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub2[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
        	end=6+length(fk)-1
          	new_gwas=gwas[,1:5]
          	new_gwas=cbind(new_gwas,gwas[,6:end])
          	names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',fk)
          	nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
          	names(nacol)=nfk
          	new_gwas=cbind(new_gwas,nacol)
          	new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
          	new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
          	all_gwas=rbind(all_gwas,new_gwas)
        }
    }
}

#X_list_ordered=lapply(X_list,function(x) x[i,])



#gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)

fwrite(all_gwas,sprintf('QTL/Biogemma_chr%s_%s_x_%s_vst_founderprobs.txt',chr,pheno,env),row.names=F,quote=F,sep='\t')

# Convert all very high and very low probabilities to 1 and 0, respectively
#X_list_full = lapply(X_list_ordered,function(x) sapply(seq(1,dim(x)[2]), function(i) ifelse(x[,i]>=0.95,1,ifelse(x[,i]<=0.05,0,x[,i]))))
#dimnames(X_list_full[[1]])[[2]]=dimnames(X_list_ordered[[1]])[[2]]

#gwas_adjusted=gwas
#sums=lapply(X_list_full,function(x) colSums(x))
#for(i in 1:16){
#    s=sums[[i]]
#    t=dim(X_list_full[[i]])[1]-2
#    l=2
#    grab=which(s>t,s)
#    grab=c(grab,which(s<l,s))
#    grab=sort(grab)
#    beta=sprintf('beta.%.0f',seq(1,16))
#    gwas_adjusted[grab,beta]=0
#    gwas_adjusted[grap,'p_value_ML']=0.99
#    print(grab)
#}

#saveRDS(gwas_adjusted,sprintf('models/Biogemma_chr%s_%s_x_%s_founderprobs_adjusted.rds',chr,pheno,env))
#hinfo=data.frame(method="Founder_probs",phenotype=pheno,environment=env,chr=chr,h2=h2_start,hap=NA,stringsAsFactors=F)
#fwrite(hinfo,'../heritabilities.txt',quote=F,sep='\t',row.names=F,append=T)
