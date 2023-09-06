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
phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
metadata=fread('metadata/samples_passed_genotype_check.txt',data.table=F)

pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates_2.txt',time),data.table=F)
metadata=metadata[metadata$experiment==time,]
genos=phenotypes$V1
######
adj_chr=c("5","9")
if(chr %in% adj_chr){
	X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))

}else{
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
}
inds=rownames(X_list[[1]])
inter=intersect(genos,inds)

genes=intersect(genes,names(phenotypes)[-1])
rownames(phenotypes)=phenotypes$V1
phenotypes=phenotypes[,-1]
phenotypes=as.matrix(phenotypes)
phenotypes=phenotypes[inter,]
genos=rownames(phenotypes)

K=K[inter,inter]

#cutoff=0.2
#ld=fread('eqtl/data/all_founder_ld.txt',data.table=F)
X_list=lapply(X_list,function(x) x[inter,])


founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
#testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
allweights=fread(sprintf('eqtl/normalized/%s_voom_weights_2.txt',time),data.table=F)
allweights=allweights[,c('V1',inter)]
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
n_reps=seq(1,length(genes))

print(length(n_reps))

#n_reps=seq(1,5)

transeqtl_gwas=function(rep){
	gene=genes[rep]
	all_gwas=data.frame(matrix(ncol=29,nrow=0))
	names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
	gene_chrom=genetable[genetable$Gene_ID==gene,]$CHROM
	testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',gene_chrom))
	snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
	#if(length(snp)>1){
	#	snp=snp[1]
	#}
	#subld=ld[(ld$SNP_A==snp & ld$CHR_B==chr & ld$r2<=cutoff &ld$SNP_B!=snp),]
	#if(nrow(subld)==0){
	#	subld=ld[(ld$SNP_B==snp & ld$CHR_A==chr & ld$r2<=cutoff &ld$SNP_A!=snp),]
	#	tsnps=unique(subld$SNP_A)
	#	tsnps=intersect(tsnps,dimnames(X_list[[1]])[[2]])
	#}else{
	#	tsnps=unique(subld$SNP_B)
	#	tsnps=intersect(tsnps,dimnames(X_list[[1]])[[2]])
	#}
	#X_list_ld=lapply(X_list,function(x) x[,tsnps])
	gweights=unlist(allweights[allweights$V1==gene,inter])
	data=data.frame(ID=rownames(phenotypes),y=phenotypes[,gene],stringsAsFactors=F)
	data=data[!is.na(data$y),]
	data$PC1=pcs[match(data$ID,pcs$V1),]$PC1
	data$PC2=pcs[match(data$ID,pcs$V1),]$PC2
	data$PC3=pcs[match(data$ID,pcs$V1),]$PC3
	rownames(data)=data$ID
	data=data[inter,]
	data$ID2=data$ID
	data=data[,c('ID','ID2','y','PC1','PC2','PC3')]
	null_model = GridLMM_ML(y~1+PC1+PC2+PC3+(1|ID),data,weights=gweights,relmat=list(ID=K),ML=T,REML=F,verbose=F)
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
	for(g in fgroups){
		subm=colnames(fkeep[,colSums(fkeep)==g,drop=F])
		if(chr==gene_chrom){
  			#snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
  			if(sum(snp %in% subm)!=0){
  				subm=subm[!(subm %in% snp)]
  			}
  		}
  		if(length(subm)!=0){
  			subfkeep=fkeep[,subm,drop=F]
  			X_list_sub=lapply(X_list,function(x) x[inter,subm,drop=F])
  			if(g==16){
      			gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=1,verbose=F)
      			gwas$Trait=gene
      			names(gwas)[6]="B73_inra"
      			names(gwas)[7:9]=c("PC1","PC2","PC3")
        		end=9+(length(founders)-1)
				names(gwas)[10:end]=founders[-1]
				new_gwas=gwas	
				new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',"B73_inra","PC1","PC2","PC3",founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
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
        			gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub2[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=1,verbose=F)
        			gwas$Trait=gene
        			names(gwas)[7:9]=c("PC1","PC2","PC3")
					end=9+(length(fk)-1)
					names(gwas)[6]=fk[1]
					names(gwas)[10:end]=fk[-1]
					new_gwas=gwas					
					nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
					names(nacol)=nfk
					new_gwas=cbind(new_gwas,nacol)
					new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',"B73_inra","PC1","PC2","PC3",founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
        			all_gwas=rbind(all_gwas,new_gwas)
        		}
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


#for(g in 1:length(n_reps)){
#	b[[n_reps[g]]]=results[[g]]
#}
saveRDS(results,sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_FIXED.rds',time,chr))

#saveRDS(b,'eqtl/trans/results/trans_eQTL_WD_0718_c6_weights_results_full.rds')

