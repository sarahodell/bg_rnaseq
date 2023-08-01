#!/usr/bin/env Rscript

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

hits=fread(sprintf('QTL/adjusted/%s_%s_QTL_scan_0.100000_hits.txt',pheno,env),data.table=F)
hits=hits[hits$CHR==chr,]
peaks=fread(sprintf('QTL/adjusted/%s_%s_QTL_scan_0.10_peaks.txt',pheno,env),data.table=F)
thresh=0.10
threshtable=fread(sprintf('QTL/adjusted_threshold_%.2f_table.txt',thresh),data.table=F)
threshold=threshtable[threshtable$method=="founder_probs" & threshtable$environment==env & threshtable$phenotype==pheno,]$threshold

peaks=peaks[peaks$CHR==chr,]
loc=which.max(peaks$value)

covar=peaks$SNP[loc]

X=do.call(cbind,lapply(X_list,function(x) x[,covar]))
colnames(X)=paste0(colnames(X),'-cov')
founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))
data$founder=founder
data$founder=factor(data$founder,levels=c(paste0(founders,'-cov')))
X=X[,colSums(X)>3]
K=K[inter,inter]
rownames(data)=data$ID
data=data[inter,]
# Run GridLMM
X_list_new=lapply(X_list,function(x) x[, - which(dimnames(x)[[2]]==covar)])
null_model = GridLMM_ML(y~X+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)
#null_model = GridLMM_ML(y~founder+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)
h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
h2_start
V_setup=null_model$setup
Y=as.matrix(data$y)
X_cov=null_model$lmod$X
X_list_null=NULL
nmarkers=dim(X_list_new[[1]])[2]
frep2=sapply(seq(1,nmarkers),function(i) lapply(X_list_new,function(j) sum(j[,i]>0.75)))
founders=names(X_list_new)
fkeep=apply(frep2,MARGIN=2,function(x) x>3)
markers=dimnames(X_list_new[[1]])[[2]]
colnames(fkeep)=markers
colnames(frep2)=markers
fgroups=unique(colSums(fkeep))
cov.names=paste0(colnames(X),'-cov')
all_gwas=data.frame(matrix(ncol=41,nrow=0))
cov.names=colnames(X_cov)
dimname=c('Trait','X_ID','s2','ML_logLik','ID.ML',cov.names,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
all_gwas=data.frame(matrix(ncol=length(dimname),nrow=0))
names(all_gwas)=dimname
for(g in fgroups){
	subm=colnames(fkeep[,colSums(fkeep)==g,drop=F])
	subfkeep=fkeep[,subm,drop=F]
	X_list_sub=lapply(X_list_new,function(x) x[,subm,drop=F])
	if(g==16){
		gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
		endp=6+length(cov.names)-1
		end=endp+15
		names(gwas)[6:endp]=cov.names
		names(gwas)[(endp+1):end]=founders[-1]
		gwas[,cov.names[-1]]=gwas[,cov.names[1]] + gwas[,cov.names[-1]]
		gwas$B73_inra=apply(gwas[,cov.names],MARGIN=1,mean)
		gwas=gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',cov.names,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
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
			if(!("B73_inra" %in% fk)){
				endp=6+length(cov.names)-1
				end=endp+(length(fk)-1)
				names(gwas)[6:endp]=cov.names
				names(gwas)[(endp+1):end]=fk[-1]
				gwas[,cov.names[-1]]=gwas[,cov.names[1]] + gwas[,cov.names[-1]]
				gwas[,fk[1]]=apply(gwas[,cov.names],MARGIN=1,mean)
				end=endp+length(fk)-1
				new_gwas=gwas
				nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
				names(nacol)=nfk
				new_gwas=cbind(new_gwas,nacol)
				new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',cov.names,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
			}else{
				endp=6+length(cov.names)-1
				end=endp+(length(fk)-1)
				names(gwas)[6:endp]=cov.names
				names(gwas)[(endp+1):end]=fk[-1]
				gwas[,cov.names[-1]]=gwas[,cov.names[1]] + gwas[,cov.names[-1]]
				gwas$B73_inra=apply(gwas[,cov.names],MARGIN=1,mean)
				end=endp+length(fk)-1
				new_gwas=gwas
				nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
				names(nacol)=nfk
				new_gwas=cbind(new_gwas,nacol)
				new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',cov.names,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
			}
			all_gwas=rbind(all_gwas,new_gwas)
		}
	}
}
sig=all_gwas[-log10(all_gwas$p_value_ML)>=threshold,]
print(dim(sig))

pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
sig$CHR=chr
sig$BP=pmap[match(sig$X_ID,pmap$marker),]$pos

#X_list_ordered=lapply(X_list,function(x) x[i,])



#gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)

fwrite(all_gwas,sprintf('QTL/Biogemma_chr%s_%s_x_%s_vst_founderprobs.txt',chr,pheno,env),row.names=F,quote=F,sep='\t')
