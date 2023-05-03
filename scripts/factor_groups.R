#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')

time="WD_0727"

#data=fread(sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_factor_correlations.txt',time),data.table=F)
data=fread(sprintf('MegaLMM/MegaLMM_%s_factor_correlations.txt',time),data.table=F)

lambda_all_means=fread(sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means.txt',time),data.table=F)
#lambda_all_means=fread(sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_all_Lambda_means.txt',time),data.table=F)
rownames(lambda_all_means)=lambda_all_means$V1
lambda_all_means=lambda_all_means[,-1]

count=c()
for(i in 1:nrow(data)){
  count=c(count,sum(abs(data[i,c(4:6)])>=0.9,na.rm=T))
}
strong=which(count==3)
sdata=data[strong,]


lambda=lambda_all_means[,'Factor2',drop=F]
lambda$lambda2=lambda$Factor2**2
sublambda=lambda[lambda$lambda2>=0.1,]
lgenes=rownames(lambda[lambda$lambda2>=0.1,])
lgenes=data.frame(genes=lgenes,stringsAsFactors=F)
fwrite(lgenes,'WD_0727_Factor2_lambda2_genelist.txt',row.names=F,quote=F,sep='\t',col.names=F)

true_nft=length(intersect(ft_genelist$V4,lgenes$genes))

genes=rownames(lambda)
ngenes=nrow(lgenes)
allgenes=nrow(lambda)
nfts=c()
for(i in 1:10000){
	draw=sample(seq(1,allgenes),ngenes)
	nft=intersect(ft_genelist$V4,genes[draw])
	nfts=c(nfts,length(nft))
}

#phenotypes=rownames(lambda_all_means)#[grepl('_',rownames(lambda_all_means))]

# Get a prop_var matrix from lambdas

prop_matrix=as.data.frame(t(apply(lambda_all_means,MARGIN=1,function(x) x^2/sum(x^2))),stringsAsFactors=F)
fwrite(prop_matrix,sprintf('MegaLMM/MegaLMM_%s_prop_variance.txt',time),row.names=T,quote=F,sep='\t')

factors=names(lambda_all_means)
factor_groups=vector("list",length=length(factors))
for(i in 1:length(factors)){
  factor_groups[[i]]$factor=factors[i]
  factor_groups[[i]]$genes=c()
  #factor_groups[[i]]$phenotypes=c()
}



for(f in 1:nrow(lambda_all_means)){
  subl=lambda_all_means[f,,drop=F]
  gene=rownames(subl)
  var_exp=apply(subl,MARGIN=1,function(x) x**2)
  tot_var=sum(var_exp)
  prop_var=var_exp/tot_var
  fkeep=names(subl[,which(var_exp>=0.1),drop=F])
  for(k in fkeep){
    x=which(unlist(unname(lapply(factor_groups,function(x) x$factor==k))))
    #if(gene %in% phenotypes){
    #  factor_groups[[x]]$phenotypes=c(factor_groups[[x]]$phenotypes,gene)
    #}else{
    factor_groups[[x]]$genes=c(factor_groups[[x]]$genes,gene)
    #}
  }
}
names(factor_groups)=factors
#saveRDS(factor_groups,sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_factor_groups.rds',time))
saveRDS(factor_groups,sprintf('MegaLMM/MegaLMM_%s_factor_groups2.rds',time))

factor_groups=readRDS(sprintf('MegaLMM/pheno_MegaLMM_%s_factor_groups.rds',time))

pheno_factors=c()
pheno=c()
prop_vars=c()
ngenes=c()
for(f in 1:length(factor_groups)){
  if(length(factor_groups[[f]]$phenotype)>0){
    print(f)
    factor=factor_groups[[f]]$factor
    for(i in 1:length(factor_groups[[f]]$phenotype)){
      p=factor_groups[[f]]$phenotype[i]
      subl=lambda_all_means[p,,drop=F]
      var_exp=apply(subl,MARGIN=1,function(x) x**2)
      tot_var=sum(var_exp)
      prop_var=var_exp/tot_var
      pvar=prop_var[which(rownames(prop_var)==factor)]
      pheno_factors=c(pheno_factors,factor_groups[[f]]$factor)
      pheno=c(pheno,p)
      prop_vars=c(prop_vars,pvar)
    }
  }
}

pheno_df=data.frame(factor=pheno_factors,phenotype=pheno,prop_var=prop_vars,stringsAsFactors=F)
fwrite(pheno_df,sprintf('MegaLMM/pheno_MegaLMM_%s_sig_factors.txt',time),row.names=F,quote=F,sep='\t')
#fwrite(pheno_df,sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_sig_factors.txt',time),row.names=F,quote=F,sep='\t')

####### Without phenotypes ##############
#########################################

#data=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_factor_correlations.txt',time),data.table=F)
data=fread(sprintf('MegaLMM/MegaLMM_%s_factor_correlations.txt',time),data.table=F)

lambda_all_means=fread(sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means.txt',time),data.table=F)
#lambda_all_means=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_all_Lambda_means.txt',time),data.table=F)
rownames(lambda_all_means)=lambda_all_means$V1
lambda_all_means=lambda_all_means[,-1]

count=c()
for(i in 1:nrow(data)){
  count=c(count,sum(abs(data[i,c(4:6)])>=0.75,na.rm=T))
}
strong=which(count==3)
sdata=data[strong,]

phenotypes=rownames(lambda_all_means)[grepl('_',rownames(lambda_all_means))]


factors=names(lambda_all_means)
factor_groups=vector("list",length=length(factors))
for(i in 1:length(factors)){
  factor_groups[[i]]$factor=factors[i]
  factor_groups[[i]]$genes=c()
  factor_groups[[i]]$phenotypes=c()
}

for(f in 1:nrow(lambda_all_means)){
  subl=lambda_all_means[f,,drop=F]
  gene=rownames(subl)
  var_exp=apply(subl,MARGIN=1,function(x) x**2)
  tot_var=sum(var_exp)
  prop_var=var_exp/tot_var
  fkeep=names(subl[,which(prop_var>=0.1),drop=F])
  for(k in fkeep){
    x=which(unlist(unname(lapply(factor_groups,function(x) x$factor==k))))
    if(gene %in% phenotypes){
      factor_groups[[x]]$phenotypes=c(factor_groups[[x]]$phenotypes,gene)
    }else{
      factor_groups[[x]]$genes=c(factor_groups[[x]]$genes,gene)
    }
  }
}

#saveRDS(factor_groups,sprintf('MegaLMM_residuals_%s_factor_groups.rds',time))
saveRDS(factor_groups,sprintf('MegaLMM/MegaLMM_%s_factor_groups.rds',time))

pheno_df=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_sig_factors.txt',time),data.table=F)
nFactors=length(unique(pheno_df$factor))
nTraits=length(unique(pheno_df$phenotype))
factorTraitvar=matrix(0, nrow = nFactors, ncol = nTraits)
rownames(factorTraitvar)=unique(pheno_df$factor)
colnames(factorTraitvar)=unique(pheno_df$phenotype)
for(i in 1:nrow(pheno_df)){
  row=pheno_df[i,]
  f=row$factor
  p=row$phenotype
  factorTraitvar[f,p]=round(row$prop_var,2)
}
#pdf('images/factor_trait_propvar.pdf')
#par(mar = c(10, 8.5, 4, 4));
#labeledHeatmap(Matrix = factorTraitvar,
#xLabels = colnames(factorTraitvar),
#yLabels = rownames(factorTraitvar),
#ySymbols = rownames(factorTraitvar),
#colorLabels = FALSE,
#colors = blueWhiteRed(50),
#textMatrix = factorTraitvar,
#setStdMargins = FALSE,
#cex.text = 0.5,
#zlim = c(-1,1),
#main = paste("Factor-trait membership"))
#dev.off()


# Dummy data
#x <- LETTERS[1:20]
#y <- paste0("var", seq(1,20))
#data <- expand.grid(X=x, Y=y)
#data$Z <- runif(400, 0, 5)
nFactors=length(unique(pheno_df$factor))
factors=unique(pheno_df$factor)
nTraits=length(unique(pheno_df$phenotype))
phenos=unique(pheno_df$phenotype)

phenosf=c("BLOIS_2014_OPT-asi",
"BLOIS_2017_OPT-asi",
"GRANEROS_2015_OPT-asi",
"NERAC_2016_WD-asi",
"STPAUL_2017_WD-asi",
"SZEGED_2017_OPT-asi",
"BLOIS_2014_OPT-female_flowering_d6",
"BLOIS_2017_OPT-female_flowering_d6",
"GRANEROS_2015_OPT-female_flowering_d6",
"NERAC_2016_WD-female_flowering_d6",
"STPAUL_2017_WD-female_flowering_d6",
"SZEGED_2017_OPT-female_flowering_d6",
"BLOIS_2014_OPT-male_flowering_d6",
"BLOIS_2017_OPT-male_flowering_d6",
"GRANEROS_2015_OPT-male_flowering_d6",
"STPAUL_2017_WD-male_flowering_d6",
"SZEGED_2017_OPT-male_flowering_d6",
"BLOIS_2014_OPT-grain_yield_15",
"BLOIS_2017_OPT-grain_yield_15",
"GRANEROS_2015_OPT-grain_yield_15",
"NERAC_2016_WD-grain_yield_15",
"STPAUL_2017_WD-grain_yield_15",
"SZEGED_2017_OPT-grain_yield_15",
"BLOIS_2014_OPT-harvest_grain_moisture",
"BLOIS_2017_OPT-harvest_grain_moisture",
"GRANEROS_2015_OPT-harvest_grain_moisture",
"NERAC_2016_WD-harvest_grain_moisture",
"STPAUL_2017_WD-harvest_grain_moisture",
"SZEGED_2017_OPT-harvest_grain_moisture",
"BLOIS_2014_OPT-tkw_15",
"BLOIS_2017_OPT-tkw_15",
"GRANEROS_2015_OPT-tkw_15",
"NERAC_2016_WD-tkw_15",
"STPAUL_2017_WD-tkw_15",
"SZEGED_2017_OPT-tkw_15",
"BLOIS_2014_OPT-total_plant_height",
"BLOIS_2017_OPT-total_plant_height",
"GRANEROS_2015_OPT-total_plant_height",
"NERAC_2016_WD-total_plant_height",
"STPAUL_2017_WD-total_plant_height",
"SZEGED_2017_OPT-total_plant_height")

new_phenodf=c()
for(f in factors){
  for(p in phenos){
    #val=pheno_df[pheno_df$factor==f & pheno_df$phenotype==p,]
    val=prop_matrix[p,f]
    line=c(f,p,val)
    new_phenodf=rbind(new_phenodf,line)
  }
}
new_phenodf=as.data.frame(new_phenodf,stringsAsFactors=F)
names(new_phenodf)=c('factor','phenotype','prop_var')
new_phenodf$factor_f=factor(new_phenodf$factor,levels=factors)
new_phenodf$phenos_f=factor(new_phenodf$phenotype,levels=phenosf)
new_phenodf$prop_var=as.numeric(new_phenodf$prop_var)
# Heatmap
p1=ggplot(new_phenodf, aes(factor_f, phenos_f, fill= prop_var)) +
scale_fill_gradient(low="white", high="darkblue") +
  geom_tile(color="black") +
  theme_classic()+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(text = element_text(size=10)) + xlab("Factor") + ylab("Environment-Phenotype") +
  ggtitle("Proportion Variance")

pdf(sprintf('images/%s_factor_trait_propvar.pdf',time))
print(p1)
dev.off()

cis=fread(sprintf('eqtl/results/%s_cis_eQTL_fkeep_hits.txt',time),data.table=F)
pheno_factors=unique(pheno_df$factor)

pheno_loc=c()
for(f in pheno_factors){
  indices=which(unlist(unname(lapply(factor_groups,function(x) x$factor==f))))
  pheno_loc=c(pheno_loc,indices)
}

cis_genes=c()
cis_factors=c()

for(factor in pheno_factors){
  test=factor_groups[[factor]]$genes
  overlap=cis$Gene %in% test
  if(sum(overlap)>0){
    print(factor)
    cis_gene=cis[which(overlap),]$Gene
    for(c in cis_gene){
      cis_genes=c(cis_genes,c)
      cis_factors=c(cis_factors,factor)
    }
  }
}

ft_genes=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
names(ft_genes)=c('CHR','START','END','GENE_ID')

ft_df=pheno_df[pheno_df$phenotype %in% phenotypes[grep('flowering',phenotypes)],]
ft_factors=unique(ft_df$factor)

ft_geneoverlap=c()
for(i in 1:length(factor_groups)){
  #inds=which(pheno_factors==f)
  #i=pheno_loc[inds]
  fgenes=factor_groups[[i]]$genes
  f=factor_groups[[i]]$factor
  ngenes=length(fgenes)
  match=intersect(fgenes,ft_genes$GENE_ID)
  if(f %in% ft_factors){
    inft=T
  }else{
    inft=F
  }
  ft_geneoverlap=rbind(ft_geneoverlap,c(f,length(match),ngenes,inft))
}
ft_geneoverlap=as.data.frame(ft_geneoverlap,stringsAsFactors=F)
names(ft_geneoverlap)=c('factor','overlap','ngenes','inft')
ft_geneoverlap$overlap=as.numeric(ft_geneoverlap$overlap)
ft_geneoverlap$ngenes=as.numeric(ft_geneoverlap$ngenes)
ft_geneoverlap$inft=factor(ft_geneoverlap$inft,levels=c('FALSE','TRUE'))
m1=lm(overlap~ngenes+inft,ft_geneoverlap)
anova(m1)
#Analysis of Variance Table

#Response: overlap
#          Df Sum Sq Mean Sq   F value  Pr(>F)
#ngenes     1  56213   56213 1099.1462 < 2e-16 ***
#inft       1    212     212    4.1411 0.04483 *
#Residuals 89   4552      51
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(m1)

#Call:
#lm(formula = overlap ~ ngenes + inft, data = ft_geneoverlap)

#Residuals:
#     Min       1Q   Median       3Q      Max
#-23.0984  -2.9110   0.2802   1.8528  22.2714

#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept) -1.4600888  0.9650060  -1.513   0.1338
#ngenes       0.0248140  0.0008632  28.745   <2e-16 ***
#inftTRUE     5.7096393  2.8057780   2.035   0.0448 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 7.151 on 89 degrees of freedom
#Multiple R-squared:  0.9254,	Adjusted R-squared:  0.9237
#F-statistic: 551.6 on 2 and 89 DF,  p-value: < 2.2e-16

#> summary(m1)$coef
#               Estimate   Std. Error   t value     Pr(>|t|)
#(Intercept) -1.46008876 0.9650060065 -1.513036 1.338136e-01
#ngenes       0.02481401 0.0008632376 28.745286 8.057256e-47
#inftTRUE     5.70963933 2.8057779527  2.034958 4.483028e-02


cis_overlap=data.frame(cis_gene=cis_genes,factor=cis_factors,stringsAsFactors=F)
fwrite(cis_overlap,sprintf('eqtl/results/%s_cis_pheno_factor_overlap.txt',time),row.names=F,quote=F,sep='\t')

cis_eqtl_overlap=fread('eqtl/results/WD_0712_eQTL_QTL_overlap.txt',data.table=F)
qtt=intersect(cis_eqtl_overlap$Gene,cis_overlap$cis_gene)



##### No. of genes per pheno_factor
factor_groups=readRDS(sprintf('MegaLMM/pheno_MegaLMM_%s_factor_groups.rds',time))

pheno_df=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_sig_factors.txt',time),data.table=F)
pheno_factors=unique(pheno_df$factor)
ngenes=sapply(pheno_factors,function(x) length(factor_groups[[x]]$genes))

genes_df=data.frame(factor=names(ngenes),ngenes=ngenes,stringsAsFactors=F)
genes_df$factor_f=factor(genes_df$factor,levels=c(genes_df$factor))

p1=ggplot(genes_df,aes(x=factor_f,y=ngenes)) + geom_bar(stat='identity') +
geom_text(aes(label=ngenes), vjust=-0.3, size=3)+ xlab("Factor") +
ylab("Number of Genes") + ggtitle(sprintf('Number of Genes in Phenotype Factors %s',time)) +
theme(axis.text.x=element_text(angle=-45,hjust=-0.25))

png(sprintf('images/factor_gene_counts_%s.png',time))
print(p1)
dev.off()
