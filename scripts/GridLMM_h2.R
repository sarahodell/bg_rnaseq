#!/usr/bin/env Rscript


library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('data.table')
library('dplyr')
library('stringr')
library('ggplot2')


founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

# Read in Kinship Matrix
K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
#phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:9]

chr=10
X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])

h2_df=c()

for(env in envs){
	for(pheno in phenos){
		data=data.frame(ID=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
		data=data[data$Loc.Year.Treat==env,]
		data=data[!is.na(data$y),]
		
		inter=intersect(data$ID,inds)
		#new_X_list=lapply(X_list,function(x) x[inter,])
		newK=K[inter,inter]
		rownames(data)=data$ID
		data=data[inter,]
		# Run GridLMM
		null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=newK),ML=T,REML=F,verbose=F)
		h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
		names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
		#h2_start
		
		line=data.frame(phenotype=pheno,environment=env,h2=h2_start,stringsAsFactors=F)
		h2_df=rbind(h2_df,line)
	}
}

h2_df=as.data.frame(h2_df,stringsAsFactors=F)
names(h2_df)=c('phenotype','environment','h2')

p1=ggplot(data=h2_df,aes(x=phenotype,y=h2,group=environment)) + geom_point(aes(color=environment)) + geom_line(aes(color=environment)) +
 xlab("Phenotype") +
ylab(expression("h"^2))

png('QTL/images/h2_by_env.png',width=800,height=600)
print(p1)
dev.off()

