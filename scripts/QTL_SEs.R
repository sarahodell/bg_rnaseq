#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
#factor=as.character(args[[2]])

library('data.table')
library('lme4')
library('lmerTest')
library('pbkrtest')
#library('PLS205')
library('lme4qtl',lib='/home/sodell/R/x86_64-conda-linux-gnu/4.2')
library('emmeans',lib='/home/sodell/R/x86_64-conda-linux-gnu/4.2')
library('multcomp',lib='/home/sodell/R/x86_64-conda-linux-gnu/4.2')
library('multcompView',lib='/home/sodell/R/x86_64-conda-linux-gnu/4.2')
library('ggplot2')
library('GridLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('preprocessCore')



colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
#has_mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)

#qtl=fread(sprintf('allelic/%s_%s_pheno_trans_kept_SNPS.txt',time,factor),data.table=F)
#qtl=fread(sprintf('eqtl/results/%s_cis_eQTL_fkeep_hits.txt',time),data.table=F)
qtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits.txt',data.table=F)
qtl=qtl[qtl$time==time,]
#qtl=fread(sprintf('eqtl/results/WD_0727_trans_eQTL_scan_hits.txt',time),data.table=F)
phenotype=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)

#phenotype=fread(sprintf('eqtl/normalized/WD_0712_voom_normalized_gene_counts_formatted.txt',time),data.table=F)

K = fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table = F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

inter=intersect(rownames(K),phenotype$V1)
K=K[inter,inter]
phenotype=phenotype[phenotype$V1 %in% inter,]

pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
genos=phenotype$V1


genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)


geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
phenotype=phenotype[,c('V1',kept_genes)]

metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
metadata=metadata[metadata$experiment==time,]

######
genos=phenotype$V1
rownames(phenotype)=phenotype$V1
phenotype=phenotype[,-1]
phenotype=as.matrix(phenotype)

allweights=fread(sprintf('eqtl/normalized/%s_voom_weights.txt',time),data.table=F)
allweights=allweights[,c('V1',inter)]

#Ynorm=c()
#plates=unique(data$plate)
#for(p in plates){
#  pinds=data[data$plate==p,]$ID
#  subY=phenotype[pinds,]
#  subYnorm=normalize.quantiles(as.matrix(subY))
#  rownames(subYnorm)=rownames(subY)
#  colnames(subYnorm)=colnames(subY)
#  Ynorm=rbind(Ynorm,subYnorm)
#}

#Ynorm=as.matrix(Ynorm)
#phenotype=Ynorm
#genos=metadata[match(pcs$sample,metadata$sample_name),]$dh_genotype
#pcs$ID=genos



######### Try just using lme4 without kinship matrix and input that into emmeans and cld
#qtl=fread(sprintf('eqtl/results/%s_cis_eQTL_hits.txt',time),data.table=F)

#phenotype=fread('eqtl/normalized/WD_0712_voom_normalized_gene_counts.txt',data.table=F)
#metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
#cs=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time),data.table=F)
#samples=names(phenotype)[-1]
#genos=metadata[match(samples,metadata$sample_name),]$dh_genotype
#phenotype=as.data.frame(t(phenotype),stringsAsFactors=F)
#names(phenotype)=c(phenotype[1,])
#phenotype=phenotype[-1,]
#phenotype$ID=genos
ses=list()
ses_plot=list()
count=1
for(i in 1:nrow(qtl)){
	chr=as.character(qtl[i,]$CHR)
	snp=qtl[i,]$X_ID
	pos=qtl[i,]$BP
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results.txt',time,chr),data.table=F)
	founder_probs = readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	founder_probs=lapply(founder_probs,function(x) x[inter,])
	genetable2=genetable[genetable$CHROM==chr,]
	genes=unique(genetable2$Gene_ID)
	genes=intersect(genes,names(phenotype)[-1])
	inter=intersect(rownames(founder_probs[[1]]),genos)
	X = do.call(cbind,lapply(founder_probs,function(x) x[,snp]))
	y=qtl[i,]$Trait
	gweights=unlist(allweights[allweights$V1==y,inter])
	subpheno = phenotype[,y,drop=F]
	subpheno=as.data.frame(subpheno,stringsAsFactors=F)
	subpheno$ID=rownames(subpheno)
	names(subpheno)=c('y','ID')
	subpheno=subpheno[!is.na(subpheno$y),]
	subpheno=subpheno[inter,]
	frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
	fkeep=founders[frep2>3]
	X=X[,fkeep]
	subpheno=cbind(subpheno,X)
	data=subpheno
	data$PC1=pcs[match(data$ID,pcs$sample),]$PC1
	data$PC2=pcs[match(data$ID,pcs$sample),]$PC2
	data$PC3=pcs[match(data$ID,pcs$sample),]$PC3
	m1=lm(y~PC1+PC2+PC3,data)
	resids=summary(m1)$residuals
	subpheno$resids=resids
	f=as.formula(paste('resids',paste(c(0,fkeep,"(1|ID)"),collapse=" + "),sep=" ~ "))
	m3=relmatLmer(f,data=subpheno,relmat = list(ID = K),weights=gweights)
	for(f in fkeep){
		atlist=list()
		for(f2 in fkeep){
			if(f==f2){
				atlist[[f2]]=1
			}else{
				atlist[[f2]]=0
			}
		}
    	refgrid=ref_grid(m3,at=atlist)
    	em=emmeans(refgrid,specs=f,lmer.df="kenward-roger")
    	if(f==fkeep[1]){
    		allem=em
    	}else{
    		allem=rbind(allem,em)
    	}
	}
	cld=cld(allem,Letters=letters)
	cld=as.data.frame(cld,stringsAsFactors=F)
	cld=cld[order(cld$emmean),]
	rownames(cld)=1:nrow(cld)
	variable_f=sapply(seq(1,length(fkeep)),function(x) fkeep[which(cld[x,fkeep]==1)])
	cld$variable_f=factor(variable_f,levels=variable_f)
	fgroups=cld$.group
	qtl$gene_snp=paste0(qtl$Trait,'_',qtl$X_ID)
	gs=qtl[i,]$gene_snp
	print(gs)
	results$gene_snp=paste0(results$Trait,'_',results$X_ID)
	row=results[results$gene_snp==gs,]
	betas=unlist(row[,fkeep])
	betas[-1]=betas[-1]+betas[1]
	rownames(cld)=cld$variable_f
	em=cld[fkeep,]
	print(cor(em$emmean,unname(betas)))
	ses[[count]]=list(gene=y,values=cld,tukey_res=fgroups,time=time,fkeep=fkeep,colsum=colSums(X))
	p1=ggplot(cld,aes(x=variable_f,y=emmean,color=variable_f)) + geom_point() +
 	 geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),data=cld) +
  	geom_text(aes(label=.group),vjust=-5,color="black",size=6)+
  	ylim(min(cld$lower.CL)-0.25,max(cld$upper.CL)+0.25) +
  	scale_color_manual(values=colorcodes[levels(cld$variable_f),]$hex_color,labels=levels(cld$variable_f))+
  	theme(axis.text.x=element_text(size=10,angle=45)) +
  	xlab("Founder") + ylab("Expression (log2CPM)") +
  	labs(title=sprintf("%s cis-eQTL Founder Effect Sizes (lme4qtl)",y))
  	ses_plot[[count]]=p1
  	count=count+1
}
saveRDS(ses,sprintf('eqtl/results/%s_ciseQTL_founder_lme4qtl_es.rds',time))
pdf(sprintf('eqtl/images/%s_ciseQTL_founder_lme4qtl_es.pdf',time))
for(i in 1:length(ses_plot)){
  print(ses_plot[[i]])
}
dev.off()


# Exp by Ind

ciseqtl=fread(sprintf('eqtl/results/%s_cis_eQTL_weights_fdr_hits.txt',time),data.table=F)
norm=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)


plot_list=list()
count=1
for(i in 1:nrow(ciseqtl)){
	gene=ciseqtl[i,]$Gene
	ex=data.frame(ID=norm$V1,exp=unlist(norm[,gene]),stringsAsFactors=F)
	ex=ex[ex$ID %in% inter,]
	ex=ex[order(ex$exp),]
	rownames(ex)=seq(1,nrow(ex))
	ex$ID_f=factor(ex$ID,levels=c(unique(ex$ID)))
	p=ggplot(aes(x=ID_f,y=exp),data=ex) + geom_point() +
		xlab('Sample') +
		ylab('Expression (log2CPM)') + geom_hline(yintercept=1)
	plot_list[[count]]=p
	count=count+1
}

pdf(sprintf('eqtl/images/%s_eqtl_by_ind.pdf',time))
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()


##### Trans ######

colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
#has_mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)


qtl=fread(sprintf('eqtl/results/%s_trans_eQTL_weights_hits.txt',time),data.table=F)

# for WD_0712
#indices=c(1,2,6,8,10)

# for WD_0718
#indices=c(1,2,6)

# for WD_0720
#indices=c(1)

#for WD_0727
indices=c(1,3,8,14)





qtl=qtl[indices,]
rownames(qtl)=seq(1:nrow(qtl))
#qtl=fread(sprintf('eqtl/results/WD_0727_trans_eQTL_scan_hits.txt',time),data.table=F)
phenotype=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)

#phenotype=fread(sprintf('eqtl/normalized/WD_0712_voom_normalized_gene_counts_formatted.txt',time),data.table=F)

K = fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table = F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

inter=intersect(rownames(K),phenotype$V1)
K=K[inter,inter]
phenotype=phenotype[phenotype$V1 %in% inter,]

pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time),data.table=F)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
genos=phenotype$V1

allweights=fread(sprintf('eqtl/normalized/%s_voom_weights.txt',time),data.table=F)
allweights=allweights[,c('V1',inter)]

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)


geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
phenotype=phenotype[,c('V1',kept_genes)]

metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
metadata=metadata[metadata$experiment==time,]

######
genos=phenotype$V1
rownames(phenotype)=phenotype$V1
phenotype=phenotype[,-1]
phenotype=as.matrix(phenotype)

ses=list()
ses_plot=list()
count=1

chroms=unique(qtl$CHR)
for(c in chroms){
	subqtl=qtl[qtl$CHR==c,]
	subqtl$gene_snp=paste0(subqtl$Gene,'_',subqtl$SNP)

	qgenes=unique(subqtl$Gene)
	chr=as.character(c)
	results=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results.rds',time,chr))
  	tgenes=which(unlist(lapply(results,function(x) unique(x$Trait))) %in% qgenes)
  	results=results[tgenes]
  	founder_probs = readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	for(i in 1:nrow(subqtl)){
		gene=subqtl[i,]$Gene
		tgenes=which(unlist(lapply(results,function(x) unique(x$Trait)))==gene)
		res=results[[tgenes]]
  		snp=subqtl[i,]$SNP
  		pos=subqtl[i,]$BP
  		genetable2=genetable[genetable$CHROM==chr,]
  		genes=unique(genetable2$Gene_ID)
 		genes=intersect(genes,colnames(phenotype))
 		inter=intersect(rownames(founder_probs[[1]]),genos)
 		X = do.call(cbind,lapply(founder_probs,function(x) x[inter,snp]))
		y=subqtl[i,]$Gene
		gweights=unlist(allweights[allweights$V1==y,inter])
  		subpheno = phenotype[,y,drop=F]
 		subpheno=as.data.frame(subpheno,stringsAsFactors=F)
  		subpheno$ID=rownames(subpheno)
  		names(subpheno)=c('y','ID')
  		subpheno=subpheno[!is.na(subpheno$y),]
  		subpheno=subpheno[inter,]
  		frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
  		fkeep=founders[frep2>3]
  		X=X[,fkeep]
  		subpheno=cbind(subpheno,X)

  		data=subpheno
  		data$PC1=pcs[match(data$ID,pcs$sample),]$PC1
  		data$PC2=pcs[match(data$ID,pcs$sample),]$PC2
  		data$PC3=pcs[match(data$ID,pcs$sample),]$PC3
  		
  		m1=lm(y~PC1+PC2+PC3,data)
  		resids=summary(m1)$residuals
  		subpheno$resids=resids
  		subpheno$resid_inter=subpheno$resids + summary(m1)$coefficients[1,1]
   	 	f=as.formula(paste('resid_inter',paste(c(1,fkeep[-1],"(1|ID)"),collapse=" + "),sep=" ~ "))
  
  		m3=relmatLmer(f,data=subpheno,relmat = list(ID = K),weights=gweights)

  		for(f in fkeep){
  			if(f==fkeep[1]){
  				atlist=vector("list",length(fkeep))
  				atlist=lapply(atlist,function(x) x=0)
  				atlist[[1]]=1
  				names(atlist)=c('(Intercept)',fkeep[-1])
  			 	refgrid=ref_grid(m3,at=atlist)
  			 	em=emmeans(refgrid,specs=fkeep[2],lmer.df="kenward-roger")
				allem=em
  			}else{
  				atlist=vector("list",length(fkeep))
  				atlist=lapply(atlist,function(x) x=0)
  				names(atlist)=c('(Intercept)',fkeep[-1])
  				atlist[[1]]=1
    			atlist[[f]]=1

    			refgrid=ref_grid(m3,at=atlist)
    			em=emmeans(refgrid,specs=f,lmer.df="kenward-roger")
      			allem=rbind(allem,em)
  			}	
  		}
  		
  		cld=cld(allem,Letters=letters)
  		cld=as.data.frame(cld,stringsAsFactors=F)
  		iloc=which(cld[,fkeep[2]]=="0")
  		cld[,fkeep[1]]=""
		cld[iloc,fkeep[1]]="1"
		cld[iloc,fkeep[2]]=""

  		cld=cld[order(cld$emmean),]
		cld=cld[,c(fkeep,"emmean","SE","df","lower.CL","upper.CL",".group")]
  		rownames(cld)=1:nrow(cld)
  		variable_f=sapply(seq(1,length(fkeep)),function(x) fkeep[which(cld[x,fkeep]=="1")])
  		cld$variable_f=factor(variable_f,levels=variable_f)
  		fgroups=cld$.group

  		gs=subqtl[i,]$gene_snp
  		print(gs)
  		
  		res$gene_snp=paste0(res$Trait,'_',res$X_ID)
  		row=res[res$gene_snp==gs,]
  		betas=unlist(row[,fkeep])
  		betas[-1]=betas[-1]+betas[1]
  		rownames(cld)=cld$variable_f
  		em=cld[fkeep,]
  		print(cor(em$emmean,unname(betas)))

  		ses[[count]]=list(gene=y,values=em,tukey_res=fgroups,time=time,fkeep=fkeep,colsum=colSums(X))
     	
     	p1=ggplot(cld,aes(x=variable_f,y=emmean,color=variable_f)) + geom_point() +
  		geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),data=cld) +
  		geom_text(aes(label=.group),vjust=-5,color="black",size=6)+
 		ylim(min(cld$lower.CL)-0.25,max(cld$upper.CL)+0.25) +
  		scale_color_manual(values=colorcodes[levels(cld$variable_f),]$hex_color,labels=levels(cld$variable_f))+
  		theme(axis.text.x=element_text(size=10,angle=45)) +
  		xlab("Founder") + ylab("Expression (log2CPM)") +
  		labs(title=sprintf("%s trans-eQTL Founder Effect Sizes (lme4qtl)",y))
  		
  		#png('test.png')
  		#print(p1)
  		#dev.off()
  		
  		
  		ses_plot[[count]]=p1
  		count=count+1
	}
}

saveRDS(ses,sprintf('eqtl/results/%s_transeQTL_founder_lme4qtl_es.rds',time))

pdf(vfsprintf('eqtl/images/%s_transeQTL_founder_lme4qtl_es.pdf',time))
for(j in 1:length(ses_plot)){
  print(ses_plot[[j]])
}
dev.off()



####### Factor trans-eQTL
time="WD_0712"

colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
#has_mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)

qtl=fread('eqtl/results/all_residual_factor_fdr_peaks_FIXED.txt',data.table=F)#indices=c(1,5)




qtl=qtl[qtl$time==time,]
rownames(qtl)=seq(1:nrow(qtl))
#qtl=fread(sprintf('eqtl/results/WD_0727_trans_eQTL_scan_hits.txt',time),data.table=F)
phenotype=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_all_F_means_FIXED.txt',time),data.table=F)

#phenotype=fread(sprintf('eqtl/normalized/WD_0712_voom_normalized_gene_counts_formatted.txt',time),data.table=F)

K = fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table = F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

inter=intersect(rownames(K),phenotype$V1)
K=K[inter,inter]
phenotype=phenotype[phenotype$V1 %in% inter,]

#pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time),data.table=F)
metadata=fread('metadata/samples_passed_genotype_check.txt',data.table=F)
genos=phenotype$V1

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

metadata=metadata[metadata$experiment==time,]

######
genos=phenotype$V1
rownames(phenotype)=phenotype$V1
phenotype=phenotype[,-1]
phenotype=as.matrix(phenotype)

ses=list()
ses_plot=list()
count=1

factors=unique(qtl$Trait)

adj_chr=c("5","9")

chroms=unique(qtl$CHR)
for(c in chroms){
	subqtl=qtl[qtl$CHR==c,]
	subqtl$gene_snp=paste0(subqtl$Trait,'_',subqtl$X_ID)
	chr=as.character(c)
	if(chr %in% adj_chr){
		founder_probs = readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
	}else{
		founder_probs = readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	}
	for(i in 1:nrow(subqtl)){
		factor=subqtl[i,]$Trait
		res=fread(sprintf('eqtl/trans/results/%s_residuals_c%s_%s_trans_results_FIXED.txt',time,chr,factor),data.table=F)
  		snp=subqtl[i,]$X_ID
  		#pos=subqtl[i,]$BP
 		#genes=intersect(genes,colnames(phenotype))
 		inter=intersect(rownames(founder_probs[[1]]),genos)
 		X = do.call(cbind,lapply(founder_probs,function(x) x[inter,snp]))
		y=subqtl[i,]$Trait
		
  		subpheno = phenotype[,y,drop=F]
 		subpheno=as.data.frame(subpheno,stringsAsFactors=F)
  		subpheno$ID=rownames(subpheno)
  		names(subpheno)=c('y','ID')
  		subpheno=subpheno[!is.na(subpheno$y),]
  		subpheno=subpheno[inter,]
  		frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
  		fkeep=founders[frep2>3]
  		X=X[,fkeep]
  		subpheno=cbind(subpheno,X)

  		data=subpheno
  		data$plate=as.factor(metadata[match(data$ID,metadata$dh_genotype),]$plate)
  		m1=lm(y~plate,data)
  		resids=summary(m1)$residuals
  		subpheno$resids=resids
   	 	f=as.formula(paste('resids',paste(c(0,fkeep,"(1|ID)"),collapse=" + "),sep=" ~ "))
  
  		m3=relmatLmer(f,data=subpheno,relmat = list(ID = K))

  		for(f in fkeep){
    		atlist=list()
    		for(f2 in fkeep){
      			if(f==f2){
        			atlist[[f2]]=1
      			}else{
        			atlist[[f2]]=0
     	 		}
    		}
    		refgrid=ref_grid(m3,at=atlist)
    		em=emmeans(refgrid,specs=f,lmer.df="kenward-roger")
    		if(f==fkeep[1]){
      			allem=em
    		}else{
      			allem=rbind(allem,em)
    		}
  		}
  		cld=cld(allem,Letters=letters)
  		cld=as.data.frame(cld,stringsAsFactors=F)
  		cld=cld[order(cld$emmean),]
  		rownames(cld)=1:nrow(cld)
  		variable_f=sapply(seq(1,length(fkeep)),function(x) fkeep[which(cld[x,fkeep]==1)])
  		cld$variable_f=factor(variable_f,levels=variable_f)
  		fgroups=cld$.group

  		gs=subqtl[i,]$gene_snp
  		print(gs)
  		
  		res$gene_snp=paste0(res$Trait,'_',res$X_ID)
  		row=res[res$gene_snp==gs,]
  		betas=unlist(row[,fkeep])
  		betas[-1]=betas[-1]+betas[1]
  		rownames(cld)=cld$variable_f
  		em=cld[fkeep,]
  		print(cor(em$emmean,unname(betas)))

  		ses[[count]]=list(gene=y,values=em,tukey_res=fgroups,time=time,fkeep=fkeep,colsum=colSums(X))
     	
     	p1=ggplot(cld,aes(x=variable_f,y=emmean,color=variable_f)) + geom_point() +
  		geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),data=cld) +
  		geom_text(aes(label=.group),vjust=-5,color="black",size=6)+
 		ylim(min(cld$lower.CL)-0.25,max(cld$upper.CL)+0.25) +
  		scale_color_manual(values=colorcodes[levels(cld$variable_f),]$hex_color,labels=levels(cld$variable_f))+
  		theme(axis.text.x=element_text(size=10,angle=45)) +
  		xlab("Founder") + ylab("Expression (log2CPM)") +
  		labs(title=sprintf("%s %s c%s trans-eQTL Founder Effect Sizes (lme4qtl)",time, y,c))
  		
  		#png('test.png')
  		#print(p1)
  		#dev.off()
  		
  		
  		ses_plot[[count]]=p1
  		count=count+1
	}
}

saveRDS(ses,sprintf('eqtl/results/%s_residuals_transeQTL_founder_lme4qtl_es.rds',time))

pdf(sprintf('eqtl/images/%s_residuals_factor_transeQTL_founder_lme4qtl_es.pdf',time))
for(j in 1:length(ses_plot)){
  print(ses_plot[[j]])
}
dev.off()


### test
#plate=metadata[match(subpheno$ID,metadata$dh_genotype),]$plate
#subpheno$plate=as.factor(plate)

#m1=lm(y~PC1+PC2+PC3,data)
#resids=summary(m1)$residuals
#subpheno$resids=resids


#founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))
#subpheno$founder=factor(founder,levels=fkeep)
#m4 = relmatLmer(y ~ (1|plate) + (1|founder:plate)+ (1|ID) + 0 + founder ,
#data=subpheno,relmat = list(ID=K),REML=T,control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",
#   check.nobs.vs.nRE="ignore"))


#em=emmeans(m4,specs='founder',lmer.df="kenward-roger")
#cld=cld(em)
#cld=cld[order(cld$emmean),]
#rownames(cld)=1:nrow(cld)
#variable_f=sapply(seq(1,length(fkeep)),function(x) fkeep[which(cld[x,fkeep]==1)])
#cld$variable_f=factor(variable_f,levels=variable_f)

#p1=ggplot(cld,aes(x=variable_f,y=emmean,color=variable_f)) + geom_point() +
#  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),data=cld) +
#   geom_text(aes(label=.group),vjust=-5,color="black",size=6)+
#   ylim(min(cld$lower.CL)-0.25,max(cld$upper.CL)+0.25) +
#   scale_color_manual(values=colorcodes[levels(cld$variable_f),]$hex_color,labels=levels(cld$variable_f))+
#   theme(axis.text.x=element_text(size=10,angle=45)) +
   #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#   xlab("Founder") + ylab("Expression (log2CPM)") +
#   labs(title=sprintf("%s cis-eQTL Founder Effect Sizes (lme4qtl)",y))
   #dev.off()


#png('test.png')
#print(p1)
#dev.off()

#results=fread('eqtl/cis/results/eQTL_WD_0712_c1_fkeep_results.txt',data.table=F)
#qtl$gene_snp=paste0(qtl$Gene,'_',qtl$SNP)
#results$gene_snp=paste0(results$Trait,'_',results$X_ID)
#gs=qtl[1,]$gene_snp
#row=results[results$gene_snp==gs,]
#betas=row[,fkeep]
#em=as.data.frame(em,stringsAsFactors=F
#cor(em$emmean,unlist(unname(betas)))
#   (is_crossed(founder~plate,subpheno))
#   [1] FALSE
#               plate
#   founder      7 8
#     B73_inra   5 3
#     CO255_inra 4 0
#     FV252_inra 1 2
#     A654_inra  5 0
#     FV2_inra   4 2
#     C103_inra  3 5
#     EP1_inra   8 1
#     D105_inra  1 6
#     W117_inra  2 4
#     DK63       3 2
#     F492       5 1
#     ND245      4 2
#     VA85       2 3


#m1 = lm(y ~ 0+ X,data=subpheno)
#m2=lmer(y ~ 0 + X + (1|ID),data=subpheno,control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",
#   check.nobs.vs.nRE="ignore"))


##### multcompl code ####
#get_letters <- function( n, Letters=c(letters, LETTERS), separator="." ){

#  n.complete <- floor(n / length(Letters))        # number of complete sets of Letters
#  n.partial <- n %% length(Letters)               # number of additional Letters
#  lett <- character()
#  separ=""
#  if( n.complete > 0 ){
#    for( i in 1:n.complete ){
#      lett <- c(lett, paste(separ, Letters, sep="") )
#      separ <- paste( separ, separator, sep="" )
#    }
#  }
#  if(n.partial > 0 )
#    lett <- c(lett, paste(separ, Letters[1:n.partial], sep="") )
#  return(lett)
#}

#sweepLetters <- function(mat, start.col=1, Letters=c(letters, LETTERS), separator="."){

#  stopifnot( all(start.col %in% 1:ncol(mat)) )
#  locked <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))          # 1 indicates that another letter dependes on this entry
#  cols <- 1:ncol(mat)
#  cols <- cols[c( start.col, cols[-start.col] )]
#  if( any(is.na(cols) ) )
#    cols <- cols[-which(is.na(cols))]

#  for( i in cols){
#    tmp <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))
#    tmp[which(mat[,i]),] <- mat[which(mat[,i]),]                        # get items of those rows which are TRUE in col "i"
#    one <- which(tmp[,i]==1)

#    if( all(apply(tmp[,-i,drop=FALSE], 1, function(x) return( any(x==1) ))) ){     # there is at least one row "l" where mat[l,i] is the only item which is TRUE i.e. no item can be removed in this column
#      next
#    }
#    for( j in one ){                                                    # over all 1's
#      if( locked[j,i] == 1 ){                                           # item is locked
#        next
#      }
#      chck <- 0
#      lck <- list()
#      for( k in one ){
#        if( j==k ){
#          next
#        }
#        else{                                                           # pair j-k
#          rows <- tmp[c(j,k),]
#          dbl <- rows[1,] & rows[2,]
#          hit <- which(dbl)
#          hit <- hit[-which(hit==i)]
#          dbl <- rows[1,-i,drop=FALSE] & rows[2,-i,drop=FALSE]
#          if( any(dbl) ){
#            chck <- chck + 1
#            lck[[chck]] <- list(c(j,hit[length(hit)]), c(k,hit[length(hit)]))      # record items which have to be locked, use last column if multiple hits
#          }
#        }
#      }
#      if( (chck == (length(one)-1)) && chck != 0 ){                     # item is redundant
#        for( k in 1:length(lck) ){                                      # lock items
#          locked[ lck[[k]][[1]][1], lck[[k]][[1]][2] ] <- 1
#          locked[ lck[[k]][[2]][1], lck[[k]][[2]][2] ] <- 1
#        }
#        mat[j,i] <- FALSE                                               # delete redundant entry
#      }
#    }
#    if(all(mat[,i]==FALSE)){                                           # delete column where each entry is FALSE and restart
#      mat <- mat[,-i,drop=FALSE]
#      colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
#      return(sweepLetters(mat, Letters=Letters, separator=separator))
#    }
#  }
#  onlyF <- apply(mat, 2, function(x) return(all(!x)))
#  if( any(onlyF) ){                                                     # There are columns with just FALSE entries
#    mat <- mat[,-which(onlyF),drop=FALSE]
#    colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
#  }
#  return( mat )
#}

#insert_absorb <- function( x, Letters=c(letters, LETTERS), separator=".", decreasing = FALSE,
#                           comps = NULL, lvl_order){

#  obj_x <- deparse(substitute(x))
#  if (is.null(comps)) {
#      namx <- names(x)
#      namx <- gsub(" ", "", names(x))
#      if(length(namx) != length(x))
#          stop("Names required for ", obj_x)
#      split_names <- strsplit(namx, "-")
#      stopifnot( sapply(split_names, length) == 2 )
#      comps <- t(as.matrix(as.data.frame(split_names)))
#  }
#  rownames(comps) <- names(x)
#  lvls <- lvl_order
#  n <- length(lvls)
#  lmat <- array(TRUE, dim=c(n,1), dimnames=list(lvls, NULL) )

#  if( sum(x) == 0 ){                                                        # no differences
#    ltrs <- rep(get_letters(1, Letters=Letters, separator=separator), length(lvls) )
#    names(ltrs) <- lvls
#    colnames(lmat) <- ltrs[1]
#    msl <- ltrs
#    ret <- list(Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat)
#    class(ret) <- "multcompLetters"
#    return(ret)
#  }
#  else{
#    signifs <- comps[x,,drop=FALSE]

#    absorb <- function(m){
#      for(j in 1:(ncol(m)-1)){
#        for(k in (j+1):ncol(m)){
#          if( all(m[which(m[,k]),k] & m[which(m[,k]),j]) ){                 # column k fully contained in column j
#            m <- m[,-k, drop=FALSE]
#            return(absorb(m))
#          }
#          else if( all(m[which(m[,j]),k] & m[which(m[,j]),j]) ){            # column j fully contained in column k
#            m <- m[,-j, drop=FALSE]
#            return(absorb(m))
#          }
#        }
#      }
#      return(m)
#    }
#    for( i in 1:nrow(signifs) ){                                            # insert
#      tmpcomp <- signifs[i,]
#      wassert <- which(lmat[tmpcomp[1],] & lmat[tmpcomp[2],])               # which columns wrongly assert nonsignificance
#      if(any(wassert)){
#        tmpcols <- lmat[,wassert,drop=FALSE]
#        tmpcols[tmpcomp[2],] <- FALSE
#        lmat[tmpcomp[1],wassert] <- FALSE
#        lmat <- cbind(lmat, tmpcols)
#        colnames(lmat) <- get_letters( ncol(lmat), Letters=letters,
#                                       separator=separator)
#        if(ncol(lmat) > 1){                                                 # absorb columns if possible
#          lmat <- absorb(lmat)
#          colnames(lmat) <- get_letters( ncol(lmat),  Letters=letters,
#                                         separator=separator )
#        }
#      }
#    }
#  }
#  lmat <- lmat[,order(apply(lmat, 2, sum))]
#  lmat <- sweepLetters(lmat)                                                                  # 1st
#  lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x))))))]                # reorder columns
#  colnames(lmat) <- get_letters( ncol(lmat),  Letters=letters,
#                                 separator=separator)
#  lmat <- lmat[,order(apply(lmat, 2, sum))]                                                   # 2nd sweep
#  lmat <- sweepLetters(lmat)
#  lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x))))))]                # reorder columns
#  colnames(lmat) <- get_letters( ncol(lmat),  Letters=letters,
#                                 separator=separator)
#  ltrs <- apply(lmat,1,function(x) return(paste(names(x)[which(x)], sep="", collapse="") ) )
#  msl <- matrix(ncol=ncol(lmat), nrow=nrow(lmat))                                             # prepare monospaced letters
#  for( i in 1:nrow(lmat) ){
#    msl[i,which(lmat[i,])] <- colnames(lmat)[which(lmat[i,])]
#    absent <- which(!lmat[i,])
#    if( length(absent) < 2 ){
#      if( length(absent) == 0 )
#        next
#      else{
#        msl[i,absent] <- paste( rep(" ", nchar(colnames(lmat)[absent])), collapse="" )
#      }
#    }
#    else{
#      msl[i,absent] <- unlist( lapply( sapply( nchar(colnames(lmat)[absent]),
#                                               function(x) return(rep( " ",x)) ),
#                                       paste, collapse="") )
#    }
#  }
#  msl <- apply(msl, 1, paste, collapse="")
#  names(msl) <- rownames(lmat)
#  ret <- list( Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat,
#               aLetters = letters, aseparator = separator )
#  class(ret) <- "multcompLetters"
#  return(ret)
#}

##################


#phenotype=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_all_F_means.txt',time),data.table=F)
#phenotype=phenotype[,c('V1','V1',factor)]
#names(phenotype)=c('ID','ID2','y')

#ses=list()
#ses_plot=list()
#count=1
#for(i in 1:nrow(qtl)){
#  chr=as.character(qtl[i,]$CHR)
#  snp=qtl[i,]$SNP
#  pos=qtl[i,]$BP
#  results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_results2.txt',time,chr),data.table=F)
#  founder_probs = readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
#  K = fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
#  rownames(K) = K[,1]
#  K = as.matrix(K[,-1])
#  K=K[inter,inter]
#  X = do.call(cbind,lapply(founder_probs,function(x) x[inter,snp]))
  #colnames(X) = founders
  #rownames(X) = dimnames(founder_probs[[1]])[[1]]
#  y=qtl[i,]$Gene
#  subpheno = phenotype[,c('ID',y)]
#  names(subpheno)=c('ID','y')
#  subpheno=subpheno[!is.na(subpheno$y),]
#  rownames(subpheno)=subpheno$ID
#  subpheno=subpheno[inter,]
  #phenotype$y=phenotype$y-mean(phenotype$y)
#  frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
#  fkeep=founders[frep2>2]
#  X=X[,fkeep]
#  m4 = relmatLmer(y ~ 0 + X + (1|ID),data=subpheno,relmat = list(ID=K),REML=T)

#  se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)

#  rownames(se4)=fkeep
#  names(se4)=c('value','se','tvalue')
#  se4$founder=rownames(se4)
#  se4=se4[order(se4$value),]
#  se4$variable_f=factor(se4$founder,levels=se4$founder)
#  glm_betas=unlist(unname(results[results$Trait==y,fkeep]))
#  glm_betas[-1]=glm_betas[1]+glm_betas[-1]
#  corr=cor(unlist(glm_betas),se4$value)

#  t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE){
#    if( equal.variance==FALSE )
#    {
#        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
#        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
#    } else
#    {
        # pooled standard deviation, scaled by the sample sizes
#        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
#        df <- n1+n2-2
#    }
#    t <- (m1-m2-m0)/se
#    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
#    names(dat) <- c("Difference", "StdError", "t", "pvalue")
#    return(dat)
#  }
#  N <- nrow(subpheno) # total sample size
#  k <- length(fkeep) # number of treatments
#  whichf=unlist(unname(apply(X,1,function(x) which.max(x))))
#  subpheno$maxf=fkeep[whichf]
#  n=table(subpheno$maxf)

#  threshold=0.05/((k*(k-1))/2)
#  fkeep=se4$founder
#  tukey_res=c()
#  for(a in 1:(k-1)){
#    for(b in ((a+1):k)){
#      f1=fkeep[a]
#      f2=fkeep[b]
#      m1=se4[se4$founder==f1,]$value
#      s1=se4[se4$founder==f1,]$se
#      m2=se4[se4$founder==f2,]$value
#      s2=se4[se4$founder==f2,]$se
#      n1=unname(n[f1])
#      n2=unname(n[f2])
#      result=t.test2(m1,m2,s1,s2,n1,n2)
#      nam=names(result)
#      result=c(f1,f2,result)
#      names(result)=c('founder1','founder2',nam)
#      tukey_res=rbind(tukey_res,result)
#    }
#  }
#  tukey_res=as.data.frame(tukey_res,stringsAsFactors=F)
#  tukey_res$Difference=as.numeric(tukey_res$Difference)
#  tukey_res$StdError=as.numeric(tukey_res$StdError)
#  tukey_res$t=as.numeric(tukey_res$t)
#  tukey_res$pvalue=as.numeric(tukey_res$pvalue)
#  tukey_res$bonf_sig=tukey_res$pvalue<=threshold
#  rownames(tukey_res)=seq(1,nrow(tukey_res))
#  tukey_res$fdr_sig=p.adjust(tukey_res$pvalue,method="fdr")<=0.05

#  sigres=tukey_res[tukey_res$fdr_sig==T,]
#  all_comps=c()
#  for(a in 1:(k-1)){
#    for(b in (a+1):k){
#      all_comps=c(all_comps,paste0(fkeep[a],'-',fkeep[b]))
#    }
#  }
  #lvl_order=all_comps
#  lvl_order=se4$founder
#  sigcomps=sapply(seq(1,nrow(sigres)),function(x) paste0(sigres[x,]$founder1,'-',sigres[x,]$founder2))
#  x=all_comps %in% sigcomps
#  names(x)=all_comps

#  fgroups=insert_absorb(x=x,Letters=letters,lvl_order=fkeep)$Letters
#  se4$letter=unlist(unname(fgroups))
#  ses[[count]]=list(gene=y,values=se4,tukey_res=fgroups,time=time,cor=corr,fkeep=fkeep,colsum=colSums(X))
    #se4[2,]$value=se4[-1,]$value + se4[1,]$value
    #se4=se4[order(se4$value),]
    #GRIDLMM
  #png(sprintf('GridLMM/effect_sizes/founder_ES/%s_BLUP_founder_effect_sizes_lme4qtl.png',name),width=1000,height=800)
#  p1=ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) + geom_point() +
#  geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
#  geom_text(aes(label=letter),vjust=-5,color="black",size=10)+
#  scale_color_manual(values=colorcodes[levels(se4$variable_f),]$hex_color,labels=levels(se4$variable_f))+
#  theme(axis.text.x=element_text(size=10)) +
#  xlab("Founder") + ylab("Effect Size (log2CPM)") +
#  labs(title=sprintf("%s cis-eQTL Founder Effect Sizes",y))
  #dev.off()


#  png('test.png')
#  print(p1)
#  dev.off()

#  ses_plot[[count]]=p1
#  count=count+1
#}
#saveRDS(ses,sprintf('eqtl/results/%s_ciseQTL_founder_es.rds',time))
#pdf(sprintf('eqtl/images/%s_ciseQTL_founder_es.pdf',time))
#for(i in 1:length(ses_plot)){
#  print(ses_plot[[i]])
#}
#dev.off()
