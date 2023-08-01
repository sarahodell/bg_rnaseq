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
library('preprocessCore')
library('stringr')
library('parallel')
library('MASS')

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
phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
metadata=fread('metadata/samples_passed_genotype_check.txt',data.table=F)

#metadata=fread('metadata/BG_completed_sample_list_FIXED.txt',data.table=F)
pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates_2.txt',time),data.table=F)

#geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
#kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
#phenotypes=phenotypes[,c('V1',kept_genes)]


metadata=metadata[metadata$experiment==time,]

genos=phenotypes$V1
######
adj_chr=c("5","9")
if(chr %in% adj_chr){
	X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))

}else{
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	#chr 5 "AX-91671957" replaced with "AX-91671943"
}
#X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])
inter=intersect(genos,inds)

genes=intersect(genes,names(phenotypes)[-1])
rownames(phenotypes)=phenotypes$V1
phenotypes=phenotypes[,-1]
phenotypes=as.matrix(phenotypes)
phenotypes=phenotypes[inter,]
genos=rownames(phenotypes)

K=K[inter,inter]
X_list=lapply(X_list,function(x) x[inter,])

testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

allweights=fread(sprintf('eqtl/normalized/%s_voom_weights_2.txt',time),data.table=F)
allweights=allweights[,c('V1',inter)]


all_gwas=data.frame(matrix(ncol=29,nrow=0))
names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',"B73_inra","PC1","PC2","PC3",founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')

all_res=data.frame(matrix(ncol=length(genes),nrow=length(inter)))
names(all_res)=genes
rownames(all_res)=inter

ciseqtl_gwas<-function(rep){
	gene=genes[rep]
	gweights=unlist(allweights[allweights$V1==gene,inter])
	snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
	if(length(snp)!=0){
		data=data.frame(ID=rownames(phenotypes),y=phenotypes[,gene],stringsAsFactors=F)
		data=data[!is.na(data$y),]
		data$PC1=pcs[match(data$ID,pcs$V1),]$PC1
		data$PC2=pcs[match(data$ID,pcs$V1),]$PC2
		data$PC3=pcs[match(data$ID,pcs$V1),]$PC3
		rownames(data)=data$ID
		data=data[inter,]
		data$ID2=data$ID
		data=data[,c('ID','ID2','y','PC1','PC2','PC3')]
		null_model = GridLMM_ML(y~1+PC1+PC2+PC3+(1|ID),data=data,weights=gweights,relmat=list(ID=K),ML=T,REML=F,verbose=F)
		h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
		names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
		V_setup=null_model$setup
		Y=as.matrix(data$y)
		X_cov=null_model$lmod$X
		X_list_null=NULL
		all_gwas=c()
		if(length(snp)>1){
			frep2=sapply(snp,function(i) lapply(X_list,function(j) sum(j[,i]>0.75)))
			fkeep=apply(frep2,MARGIN=2,function(x) x>3)
			colnames(fkeep)=snp
			colnames(frep2)=snp
			fgroups=unique(colSums(fkeep))
			subgwas=c()
			for(f in fgroups){
				subm=colnames(fkeep[,colSums(fkeep)==f,drop=F])
				subfkeep=fkeep[,subm,drop=F]
				X_list_sub=lapply(X_list,function(x) x[inter,subm,drop=F])
				if(f==16){
					gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
					gwas$Trait=gene
					names(gwas)[c(6,10:24)]=founders
					names(gwas)[7:9]=c('PC1','PC2','PC3')
					new_gwas=gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
					all_gwas=rbind(all_gwas,new_gwas)
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
						end=10+length(fk)-2
						names(gwas)[c(6,10:end)]=fk
						names(gwas)[7:9]=c('PC1','PC2','PC3')
						new_gwas=gwas[,1:end]
						nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
						names(nacol)=nfk
						new_gwas=cbind(new_gwas,nacol)
						new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
						new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
						all_gwas=rbind(all_gwas,new_gwas)
					}
				}
				subgwas=rbind(subgwas,new_gwas)
				highest=which.min(subgwas$p_value_ML)
				betas=unlist(subgwas[highest,founders])
				fk=founders[fkeep[,highest]]
				betas=betas[fk]
				wn=which(!is.na(betas))[1]
				betas[-wn]=betas[-wn]+betas[wn]
				X = do.call(cbind,lapply(X_list,function(x) x[,snp[highest]]))
				colnames(X) = founders
				rownames(X) = dimnames(X_list[[1]])[[1]]
				X=X[,fk]	
			}
			high2=which.min(all_gwas$p_value_ML)[1]
			snp=snp[highest]
		}else{
			X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
			colnames(X) = founders
			rownames(X) = dimnames(X_list[[1]])[[1]]
			X_list_ordered=lapply(X_list,function(x) x[,snp,drop=F])
			frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
			fkeep=founders[frep2>3]
			X_list_ordered = X_list_ordered[c(fkeep)]
			X=X[,fkeep]
			gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
			gwas$Trait=gene
			if(!("B73_inra" %in% fkeep)){
				end=10+length(fkeep)-2
				new_gwas=gwas[,1:5]
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
				betas=unlist(new_gwas[,fkeep])
				wn=which(!is.na(betas))[1]
				betas[-wn]=betas[-wn]+betas[wn]
			}else{
				end=10+length(fkeep)-2
				new_gwas=gwas[,1:5]
				new_gwas=cbind(new_gwas,gwas[6])
				new_gwas=cbind(new_gwas,gwas[,7:9])
				new_gwas=cbind(new_gwas,gwas[,10:end])
				names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra',"PC1","PC2","PC3",fkeep[-1])
				fdrop=founders[!(founders %in% fkeep)]
				nacol=data.frame(matrix(ncol=length(fdrop),nrow=1))
				names(nacol)=fdrop
				new_gwas=cbind(new_gwas,nacol)
				new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
				new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',"B73_inra","PC1","PC2","PC3",founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
				betas=unlist(new_gwas[,fkeep])
				wn=which(!is.na(betas))[1]
				betas[-wn]=betas[wn]+betas[-wn]
			}
			all_gwas=rbind(all_gwas,new_gwas)
		}
		bv=X %*% betas
		colnames(bv)=gene
		X_r=data$y-bv
		prop_var=var(bv[,gene],na.rm=T)/var(data$y,na.rm=T)
		line=data.frame(gene=gene,time=time,snp=snp,prop_var=prop_var,stringsAsFactors=F)
		colnames(X_r)=gene
		#all_res[inter,gene]=X_r
	}
	return(list(all_gwas,X_r,line))
}

n_reps=1:length(genes)

print(system.time({
results=mclapply(n_reps,ciseqtl_gwas,mc.cores=cores)
}))

all_res=do.call(cbind,lapply(results,function(x) x[[2]]))
all_res=as.data.frame(all_res,stringsAsFactors=F)
fwrite(all_res,sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_residuals_FIXED.txt',time,chr),row.names=F,quote=F,sep='\t')


all_gwas=do.call(rbind,lapply(results,function(x) x[[1]]))
all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
fwrite(all_gwas,sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),row.names=F,quote=F,sep='\t')

all_props=do.call(rbind,lapply(results,function(x) x[[3]]))
all_props=as.data.frame(all_props,stringsAsFactors=F)
fwrite(all_props,sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_prop_var_FIXED.txt',time,chr),row.names=F,quote=F,sep='\t')

#saveRDS(results,sprintf('eqtl/cis/results/cis_eQTL_%s_c%s_weights_results_FIXED.rds',time,chr))


#genes=genes[1:10]
#for(g in 1:length(genes)){
#	
#}





