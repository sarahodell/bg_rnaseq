#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
factor=as.character(args[[2]])

library('data.table')
#library('lme4')
library('lme4qtl')
library('emmeans')
library('ggplot2')
library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')

library('dplyr')
library('preprocessCore')
library('stringr')
library('parallel')
library('MASS')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

time="WD_0727"
factor="Factor14"

transqtl=fread('eqtl/results/WD_0727_trans_eQTL_scan_hits.txt',data.table=F)
factorqtl=fread('eqtl/results/Factor14_trans_WD_0727_eQTL_fkeep_hits.txt',data.table=F)
prop_var=fread('MegaLMM/MegaLMM_WD_0727_prop_variance.txt',data.table=F)

summary(prop_var[,'Factor14'])
gois=c("Zm00001d041650","Zm00001d009688","Zm00001d045677")

#Correlation between -log10(pvalue) of transeqtl and lambda**2 values

#prop_var[prop_var$V1 %in% gois,'Factor14']
#[1] 0.7775218 0.5596646 0.7550769

# Zm00001d021130 has the highest lambda**2

# What about just lambda values?

lambda_all_means=fread('MegaLMM/MegaLMM_WD_0727_all_Lambda_means.txt',data.table=F)


#exp=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_all_F_means.txt',time),data.table=F)

gois=unique(transqtl$Gene)
plot_list=list()
count=1
for(g in gois){
	sub=exp[,c('V1',g)]
	names(sub)=c('ind','y')
	sub=sub[order(sub$y),]
	rownames(sub)=seq(1,nrow(sub))
	sub$ind_f=factor(sub$ind,levels=c(unique(sub$ind)))
	p1=ggplot(aes(x=ind_f,y=y),data=sub) + geom_point()
	plot_list[[count]]=p1
	count=count+1
}
pdf('eqtl/images/trans_eqtl_by_ind.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])	
}
dev.off()


time="WD_0727"

transeqtl=fread('eqtl/results/WD_0727_trans_eQTL_scan_hits.txt',data.table=F)
indices=c(4,5,9,13,14)
transeqtl=transeqtl[indices,]
factorqtl=fread('eqtl/results/Factor14_trans_WD_0727_eQTL_fkeep_hits.txt',data.table=F)

#transeqtl=fread('eqtl/results/%s_cis_eQTL_fkeep_hits.txt',time),data.table=F)
norm=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
rownames(norm)=norm$V1
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder

plot_list=list()
count=1
for(i in 1:nrow(transeqtl)){
  gene=transeqtl[i,]$Gene
  if(gene %in% names(norm)){
  	chr=transeqtl[i,]$CHR
    gene_chr=genetable[genetable$Gene_ID==gene,]$CHROM
    X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',chr))
    testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%.0f.rds',chr))
    inter=intersect(norm$V1,dimnames(X_list[[1]])[[1]])
    norm=norm[inter,]
	snp=transeqtl[i,]$SNP
    #snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
    X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	
    colnames(X) = founders
    rownames(X) = inter
    
    frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
    fkeep=founders[frep2>2]
    founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))

    #print("True")
    ex=data.frame(ID=norm$V1,exp=norm[,gene],founder=founder,stringsAsFactors=F)
    ex=ex[ex$ID!="",]
    ex=ex[order(ex$exp),]
    rownames(ex)=seq(1,nrow(ex))
    ex$ID_f=factor(ex$ID,levels=c(unique(ex$ID)))
    ex$founder_f=factor(ex$founder,levels=c(unique(ex$founder)))
    ex=ex[ex$founder %in% fkeep,]
    #subex=subset(ex, ID %in% drop_ind)
    p=ggplot(aes(x=ID_f,y=exp),data=ex) + geom_point() +
    geom_point(aes(color=founder_f)) +
    scale_color_manual(values=colorcodes[levels(ex$founder_f),]$hex_color,labels=levels(ex$founder_f))+
     xlab('Sample') +
     ylab('Expression (log2CPM)') + geom_hline(yintercept=1)
     
     p2=ggplot(aes(x=founder_f,y=exp),data=ex) + geom_boxplot(aes(color=founder_f)) + 
     scale_color_manual(values=colorcodes[levels(ex$founder_f),]$hex_color,labels=levels(ex$founder_f))+
     xlab('Sample') +
     ylab('Expression (log2CPM)') + geom_hline(yintercept=1)
  
     #exmean=ex %>% group_by(founder) %>% summarize(avg_exp=mean(exp))
    plot_list[[count]]=p2
    count=count+1
  }
}

pdf('images/founder_eqtl_by_ind.pdf')
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()


## pvalues of genes in Factor 14 tested against AX-91508580
snp="AX-91508580"
fgenes=prop_var[prop_var$Factor14>=0.1,]$V1
results=readRDS('eqtl/trans/results/trans_eQTL_1_cWD_0727_fkeep_results.rds')
tgenes=which(unlist(lapply(results,function(x) unique(x$Trait))) %in% fgenes)
results=results[tgenes]
df=do.call(rbind,lapply(results,function(x) x[x$X_ID==snp,]))
df$log10p=-log10(df$p_value_ML)

df$lambda2=prop_var[match(df$Trait,prop_var$V1),]$Factor14
df$lambda=lambda_all_means[match(df$Trait,lambda_all_means$V1),]$Factor14

cor(df$log10p,df$lambda2)
#[1] 0.3322874

cor(df$log10p,abs(df$lambda))
#[1] 0.4344795

gois=c("Zm00001d041650","Zm00001d009688","Zm00001d045677")

p2=ggplot(df,aes(x=abs(lambda),y=log10p)) + geom_point() +
geom_point(data = df %>% filter(Trait %in% gois), color = "red")

png('eqtl/images/lamba_by_log10p.png')
print(p2)
dev.off()
#lambda_all_means[lambda_all_means$V1 %in% gois,'Factor14']
#[1] -0.3952658 -0.2256223 -0.4550680
#quantile(lambda_all_means[,'Factor14'],0.05)
#        5% 
#-0.1320597 

#qtl=fread(sprintf('eqtl/results/%s_pheno_residuals_trans_%s_eQTL_hits.txt',factor,time),data.table=F)
#qtl=fread(sprintf('eqtl/results/%s_trans_%s_eQTL_hits.txt',factor,time),data.table=F)
#phenotype=fread(sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_all_F_means.txt',time),data.table=F)
phenotype=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_all_F_means.txt',time),data.table=F)


phenotype=phenotype[,c('V1','V1',factor)]
names(phenotype)=c('ID','ID2','y')


full_df=c()
for(c in 1:10){
  #d=fread(sprintf('eqtl/trans/results/%s_c%s_pheno_residuals_factor_trans_eQTL.txt',time,c))
  d=fread(sprintf('eqtl/trans/results/%s_c%s_pheno_factor_trans_eQTL.txt',time,c))
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  d$CHR=c
  d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
  full_df=rbind(full_df,d)
}
full_df=full_df[full_df$Trait==factor,]
full_nrow=nrow(full_df)


threshold=-log10(0.05/full_nrow)
#threshold=full_threshold
true=c()
new_value=c()
for(i in 1:nrow(qtl)){
  row=qtl[i,]
  snp=row$SNP
  chr=as.character(row$CHR)

  K=fread('../GridLMM/K_matrices/K_matrix_full.txt'),data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)
  inter=intersect(rownames(K),phenotype$ID)
  K=K[inter,inter]
  phenotype=phenotype[phenotype$ID %in% inter,]

  X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
  X_list=lapply(X_list,function(x) x[inter,])
  X = do.call(cbind,lapply(X_list,function(x) x[,snp]))

  frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
  fkeep=founders[frep2>=2]

  X_list_ordered = lapply(fkeep,function(i) X[,i,drop=F])

  null_model = GridLMM_ML(y~1 + (1|ID),phenotype,relmat=list(ID=K),ML=T,REML=F,verbose=F)

  h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
  names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
  h2_start
  V_setup=null_model$setup
  Y=as.matrix(phenotype$y)
  X_cov=null_model$lmod$X

    #X_cov is n x 0 matrix and can use full X_list_ordered
    #X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,,drop=F])

  X_list_null=NULL
    #X_list_null
  cores=1
  gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F,h2_step=1  )
  true=c(true,-log10(gwas$p_value_ML)>=threshold)
  new_value=c(new_value,-log10(gwas$p_value_ML))
}

qtl$true_sig=true
qtl$new_value=new_value
fwrite(qtl,sprintf('eqtl/results/%s_pheno_trans_%s_eQTL_hits.txt',factor,time),row.names=F,quote=F,sep='\t')
#fwrite(qtl,sprintf('eqtl/results/%s_pheno_residuals_trans_%s_eQTL_hits.txt',factor,time),row.names=F,quote=F,sep='\t')


###### check that these trans eQTL are real
#### 1) Do they still show up if frep cutoff is 3, instead of 2?
time="WD_0727"

transqtl=fread(sprintf('eqtl/results/%s_trans_eQTL_scan_hits.txt',time),data.table=F)
qgenes=unique(transqtl$Gene)


genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
genes=unique(genetable$Gene_ID)

phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time),data.table=F)

geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
phenotypes=phenotypes[,c('V1',kept_genes)]


genes=intersect(genes,names(phenotypes)[-1])
rownames(phenotypes)=phenotypes$V1
phenotypes=phenotypes[,-1]
phenotypes=as.matrix(phenotypes)
phenotypes=phenotypes[inter,]
genos=rownames(phenotypes)

metadata=metadata[metadata$experiment==time,]
metadata=metadata[metadata$read==1,]

founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


n_m=4716
n_genes=31879
adjust=n_genes*n_m
threshold=-log10(0.05/adjust)
print(threshold)

full_gwas=data.frame(matrix(ncol=29,nrow=0))
names(full_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')

for(r in 1:nrow(transqtl)){
	row=transqtl[r,]
	chr=as.character(row$CHR)
	snp=row$SNP
	value=row$value
	gene=row$Gene
	print(gene)
	print(snp)
	K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
	rownames(K)=K[,1]
	rownames(K)=gsub("-",".",rownames(K))
	K=as.matrix(K[,-1])
	colnames(K)=rownames(K)
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])
	inter=intersect(genos,inds)
	K=K[inter,inter]
	X_list=lapply(X_list,function(x) x[inter,])
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

    X_list_ordered=lapply(X_list,function(x) x[,snp,drop=F])
    X = do.call(cbind,lapply(X_list_ordered,function(x) x[,snp]))
    colnames(X) = founders
    rownames(X) = dimnames(X_list[[1]])[[1]]

    frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
    fkeep=founders[frep2>3]
    X_list_ordered = X_list_ordered[c(fkeep)]
    X=X[,fkeep]
    X_list_null=NULL
    gwas2=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
    gwas2$Trait=gene
	print("one SNP")
	print(-log10(gwas2$p_value_ML))
	all_gwas=data.frame(matrix(ncol=29,nrow=0))
	names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','PC1','PC2','PC3',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
	nmarkers=dim(X_list[[1]])[2]
	frep2=sapply(seq(1,nmarkers),function(k) lapply(X_list,function(j) sum(j[,k]>0.75)))
	founders=names(X_list)
	fkeep=apply(frep2,MARGIN=2,function(x) x>3)
	markers=dimnames(X_list[[1]])[[2]]
	colnames(fkeep)=markers
	colnames(frep2)=markers
	fgroups=unique(colSums(fkeep))
	for(g in fgroups){
		subm=colnames(fkeep[,colSums(fkeep)==g])
		if(chr==gene_chrom){
			cissnp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
			if(sum(cissnp %in% subm)!=0){
				subm=subm[subm!=cissnp]
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
	print("all SNPs")
	test=all_gwas[all_gwas$X_ID==snp,]
	print(-log10(test$p_value_ML))
	if(-log10(test$p_value_ML)>=threshold){
		print(TRUE)
	}else{
		print(FALSE)
	}
	full_gwas=rbind(full_gwas,test)
}






