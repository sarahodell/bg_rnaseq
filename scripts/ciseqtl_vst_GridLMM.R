args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
chr=as.character(args[[2]])
cores=as.numeric(args[[3]])

#date=format(Sys.time(),'%m%d%y')

library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('data.table')
library('dplyr')
#library('lme4')
library('lme4qtl')
library('preprocessCore',lib='/home/sodell/R/x86_64-conda-linux-gnu/4.2')
library('stringr')


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
phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list_FIXED.txt',data.table=F)
pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates_2.txt',time),data.table=F)

#geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
#kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
#phenotypes=phenotypes[,c('V1',kept_genes)]

metadata=metadata[metadata$experiment==time,]
metadata=metadata[metadata$read==1,]

genos=phenotypes$V1
######

X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])
inter=intersect(genos,inds)
X_list=lapply(X_list,function(x) x[inter,])

genes=intersect(genes,names(phenotypes)[-1])
rownames(phenotypes)=phenotypes$V1
phenotypes=phenotypes[,-1]
phenotypes=as.matrix(phenotypes)
phenotypes=phenotypes[inter,]
genos=rownames(phenotypes)

K=K[inter,inter]

plate=metadata[match(genos,metadata$dh_genotype),]$plate.x
df=data.frame(ID=genos,plate=plate,stringsAsFactors=F)
df$plate=as.factor(df$plate)
#Separate individuals by plate and quantile normalize separately
Ynorm=c()
plates=unique(df$plate)
for(p in plates){
  pinds=df[df$plate==p,]$ID
  subY=phenotypes[pinds,]
  subYnorm=normalize.quantiles(as.matrix(subY))
  rownames(subYnorm)=rownames(subY)
  colnames(subYnorm)=colnames(subY)
  Ynorm=rbind(Ynorm,subYnorm)
}

Ynorm=as.matrix(Ynorm)
phenotypes=Ynorm
#genos=metadata[match(samples,metadata$sample_name),]$genotype
testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

allweights=fread(sprintf('eqtl/normalized/%s_voom_weights_2.txt',time),data.table=F)
allweights=allweights[,c('V1',inter)]

dimname=c('Trait','X_ID','s2','ML_logLik','ID.ML',plates,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
all_gwas=data.frame(matrix(ncol=length(dimname),nrow=0))
names(all_gwas)=dimname

for(g in 1:length(genes)){
	gene=genes[g]
	gweights=unlist(allweights[allweights$V1==gene,inter])
	snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
	if(length(snp)!=0){
		data=data.frame(ID=rownames(phenotypes),y=phenotypes[,gene],stringsAsFactors=F)
		data=data[!is.na(data$y),]
    	#data$PC1=pcs[match(data$ID,pcs$V1),]$PC1
    	#data$PC2=pcs[match(data$ID,pcs$V1),]$PC2
    	#data$PC3=pcs[match(data$ID,pcs$V1),]$PC3
    	plate=metadata[match(data$ID,metadata$dh_genotype),]$plate.x
		data$plate=as.factor(plate)
    	rownames(data)=data$ID
    	data=data[inter,]
    	data$ID2=data$ID
    	data=data[,c('ID','ID2','y','plate')]
    	null_model = GridLMM_ML(y~1+plate+(1|ID),data=data,weights=gweights,relmat=list(ID=K),ML=T,REML=F,verbose=F)
    	h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    	names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    	h2_start
    	V_setup=null_model$setup
    	Y=as.matrix(data$y)
    	X_cov=null_model$lmod$X
    	X_list_null=NULL
    	# More than one snp in cis
    	if(length(snp)>1){
    		frep2=sapply(snp,function(i) lapply(X_list,function(j) sum(j[,i]>0.75)))
    		fkeep=apply(frep2,MARGIN=2,function(x) x>3)
    		colnames(fkeep)=snp
    		colnames(frep2)=snp
    		fgroups=unique(colSums(fkeep))
    		plates=as.character(paste0('plate',unique(data$plate)))
    		# Only one number of low rep founders
    		for(g in fgroups){
    			subm=colnames(fkeep[,colSums(fkeep)==g,drop=F])
    			subfkeep=fkeep[,subm,drop=F]
    			X_list_sub=lapply(X_list,function(x) x[inter,subm,drop=F])
    			if(g==16){
    				gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=1,verbose=F)
    				gwas$Trait=gene
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
    					gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub2[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
    					gwas$Trait=gene
    					endp=6+length(plates)-1
    					end=endp+(length(fk)-1)
    					names(gwas)[6:endp]=plates
    					names(gwas)[(endp+1):end]=fk[-1]
    					gwas[,plates[-1]]=gwas[,plates[1]] + gwas[,plates[-1]]
    					gwas[,fk[1]]=apply(gwas[,plates],MARGIN=1,mean)
    					end=endp+length(fk)-1
    					new_gwas=gwas
    					nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
    					names(nacol)=nfk
    					new_gwas=cbind(new_gwas,nacol)
    					new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',plates,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
    					all_gwas=rbind(all_gwas,new_gwas)	
    				}
      			}
    		}
    	# Only one snp in cis
    	}else{
    		X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
    		colnames(X) = founders
    		rownames(X) = dimnames(X_list[[1]])[[1]]
    		X_list_ordered=lapply(X_list,function(x) x[,snp,drop=F])
    		frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
    		fkeep=founders[frep2>3]
    		nfk=founders[!(founders %in% fkeep)]
    		X_list_ordered = X_list_ordered[c(fkeep)]
    		X=X[,fkeep]
    		gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
    		gwas$Trait=gene
    		plates=as.character(paste0('plate',unique(data$plate)))
    		if(!("B73_inra" %in% fkeep)){
    			endp=6+length(plates)-1
    			end=endp+(length(fkeep)-1)
    			names(gwas)[6:endp]=plates
    			names(gwas)[(endp+1):end]=fkeep[-1]
    			gwas[,plates[-1]]=gwas[,plates[1]] + gwas[,plates[-1]]
    			gwas[,fkeep[1]]=apply(gwas[,plates],MARGIN=1,mean)
    			end=endp+length(fkeep)-1
    			new_gwas=gwas
    			nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
    			names(nacol)=nfk
    			new_gwas=cbind(new_gwas,nacol)
    			new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',plates,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
    		}else{
    			endp=6+length(plates)-1
    			end=endp+(length(fkeep)-1)
    			names(gwas)[6:endp]=plates
    			names(gwas)[(endp+1):end]=fkeep[-1]
    			gwas[,plates[-1]]=gwas[,plates[1]] + gwas[,plates[-1]]
    			gwas$B73_inra=apply(gwas[,plates],MARGIN=1,mean)
    			end=endp+length(fkeep)-1
    			new_gwas=gwas
    			nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
    			names(nacol)=nfk
    			new_gwas=cbind(new_gwas,nacol)
    			new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',plates,founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
    		}
    		all_gwas=rbind(all_gwas,new_gwas)
    	}
    }
}

all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
fwrite(all_gwas,sprintf('eqtl/cis/results/eQTL_%s_c%s_quantnorm_results_FIXED.txt',time,chr),row.names=F,quote=F,sep='\t')
