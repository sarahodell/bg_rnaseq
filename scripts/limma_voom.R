#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])

library('data.table')
library('edgeR')
library('dplyr')
library('ggplot2')
# Make gene expression count matrix of all samples
#timep="WD_0712"
time="WD_0727"
mat=fread(sprintf('star/%s_updated_gene_counts.txt',time),data.table=F)
rownames(mat)=mat$Gene_ID
mat=mat[5:nrow(mat),]
mat=mat[,-1]
genes=rownames(mat)
mat=t(mat)
meta=fread('metadata/samples_passed_genotype_check.txt',data.table=F)

#meta=fread('metadata/BG_completed_sample_list_FIXED.txt',data.table=F)
#meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
meta=meta[meta$experiment==time,]
snames=rownames(mat)
gnames=meta[match(snames,meta$sample_name),]$dh_genotype
mat=as.data.frame(mat,stringsAsFactors=F)
mat$gnames=gnames
mat=mat[mat$gnames!="",]
plate=meta[match(mat$gnames,meta$dh_genotype),]$plate

gnames=mat$gnames

#mat=t(mat)
#mat2=mat %>% group_by(gnames) %>% summarise(across(everything(), list(mean)))
mat1 = mat %>% select(-gnames) %>% as.matrix
rownames(mat1)=gnames
mat=t(mat1)
#rownames(mat1) = mat2$gnames
coldata=data.frame(genotype=gnames,plate=plate,stringsAsFactors=F)
coldata$plate=as.factor(coldata$plate)
#dds <- DESeqDataSetFromMatrix(countData = mat,
#                              colData = coldata,
#                              design = ~ 1)

#sample_table=fread('samples.tsv',data.table=F)
#samples=sample_table[sample_table$include==TRUE,]$sample

#convert sequencing sample names to sample names from experiment
#conversion=fread('raw_reads/metadata/BG_sequencing_sample_conversion_table.txt',data.table=F)

# Preprocessing

#counts=read.delim('batch_1_raw_readcount_matrix.txt',row.names=1)
#mat=t(mat)
d0=DGEList(mat)
d0=calcNormFactors(d0)

summary(d0$samples$lib.size/1e6)

# WD_0712
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8677  1.4455  2.2734  3.4449  3.6523 17.5589 

#WD_0712
cutoff=5 #WD_0712
nind=15

#WD_0718
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9011  2.2289  2.7520  3.1207  3.7236 10.1522

#updated WD_0718
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9159  2.3899  2.9855  3.5358  3.9623 19.0705 
#cutoff=5# maybe higher cutoff for WD_0718?
#nind=15

#WD_0720
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.8933  1.6238  2.0577  2.7422  2.9605 13.3744 

#cutoff=5 # WD_0720
#nind=15

#WD_0725
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.8923  1.5518  1.9795  3.2689  3.8135 41.6707  

#cutoff=5 # WD_075
#nind=15

#WD_0727
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8357  1.3540  1.7517  2.7997  2.9597 12.0689
cutoff=5
nind=15

maxes=apply(cpm(d0),1,max)
# fewer than 2 samples with cpm higher than cutoff
drop=which(apply(cpm(d0),1,function(x) sum(x>=cutoff)<nind))
d=d0[-drop,]
dim(d)

snames=colnames(mat)



# Voom transformation and calculation of variance weights
#design_matrix=fread('vgt1_founder_test_WD_0712.txt',data.table=F)
#Add design matrix based off of founder identity at vgt1
#founder=design_matrix$founder
genotype=coldata$genotype
mm=model.matrix(~1+plate,coldata)
png(sprintf('images/voom_%s_2.png',time))
y=voom(d,mm,plot=T)
dev.off()
saveRDS(y,sprintf('eqtl/normalized/%s_voom_results_2.rds',time))
norm=as.data.frame(y$E)
fwrite(norm,sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_2.txt',time),row.names=T,quote=F,sep='\t')
weights=as.data.frame(y$weights)
rownames(weights)=rownames(norm)
names(weights)=names(norm)
fwrite(weights,sprintf('eqtl/normalized/%s_voom_weights_2.txt',time),row.names=T,quote=F,sep='\t')


y=readRDS(sprintf('eqtl/normalized/%s_voom_results_2.rds',time))

ym=as.matrix(y$E)
ym=t(ym)

pca <- prcomp(ym, center = TRUE)
 pcs=as.data.frame(pca$x,stringsAsFactors=F)
var_explained <- pca$sdev^2/sum(pca$sdev^2)

pc_mm=pcs[,1:3]
pc_mm$sample=rownames(pc_mm)
pc_mm$plate=coldata[match(pc_mm$sample,coldata$genotype),]$plate
fwrite(pc_mm,sprintf('eqtl/normalized/%s_PCA_covariates_2.txt',time),row.names=T,quote=F,sep='\t')


p1=ggplot(aes(x=PC1,y=PC2),data=pc_mm) + geom_point(aes(color=plate)) +
 xlab(sprintf('PC1 %.2f %%',var_explained[1]*100)) +
 ylab(sprintf('PC2 %.2f %%',var_explained[2]*100))

png(sprintf('images/PC1_PC2_%s_2.png',time))
print(p1)
dev.off()

norm=as.data.frame(t(norm),stringsAsFactors=F)
fwrite(norm,sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),row.names=T,quote=F,sep='\t')

#### Drop outliers and format



#WD_0712 Outliers
#pc_mm[pc_mm$PC1<= -75,]
#                      PC1        PC2        PC3         sample plate
#EB.10H.H.00062 -109.72253  25.547987  24.466926 EB.10H.H.00062     7
#EB.09S.H.00041 -113.79115 119.614185 -45.485948 EB.09S.H.00041     7
#EB.09S.H.00081 -117.98430  24.018115  -3.460107 EB.09S.H.00081     7
#EB.09S.H.00133  -90.25365  83.598869 -43.484926 EB.09S.H.00133     7
#EB.09S.H.00528  -82.10763  17.714362 -10.307004 EB.09S.H.00528     7
#EB.09S.H.00107 -100.09700 -58.018899 -70.415520 EB.09S.H.00107     7
#EB.09S.H.00074  -95.30411  11.119079  54.135581 EB.09S.H.00074     7
#EB.09S.H.00086 -118.16719  32.183642  -2.490881 EB.09S.H.00086     7
#EB.09S.H.00116  -93.67894  -5.223644  54.689449 EB.09S.H.00116     7
#EB.09S.H.00460  -89.83154 -23.936969 -57.914722 EB.09S.H.00460     7
#drop_ind=c('EB.09S.H.00041')
drop_ind=c('EB.10H.H.00025')
norm=norm[,!names(norm)%in% drop_ind]

norm=as.data.frame(t(norm),stringsAsFactors=F)
fwrite(norm,'eqtl/normalized/WD_0712_voom_normalized_gene_counts_formatted_FIXED.txt',row.names=T,quote=F,sep='\t')

#WD_0718 Outliers
#pc_mm[pc_mm$PC1<= -50,]
#EB.09S.H.00390  -68.86051 -34.485184  -8.983030 EB.09S.H.00390
#EB.09S.H.00377 -146.65385 -32.092430 -21.437155 EB.09S.H.00377
#EB.09S.H.00353  -55.95032  17.391344  21.902194 EB.09S.H.00353
#EB.09S.H.00185  -52.83590   2.284391  -2.432465 EB.09S.H.00185
#EB.10H.H.00054  -96.17086   1.761571 -22.847481 EB.10H.H.00054
#EB.09S.H.00433  -59.54090  -5.009809 -41.097642 EB.09S.H.00433
#EB.09S.H.00529  -68.28829 -61.468211  -7.014567 EB.09S.H.00529
#EB.09S.H.00264  -85.86447 -30.142220 -11.475979 EB.09S.H.00264

#norm=norm[,!names(norm)%in% drop_ind]
#fwrite(norm,'eqtl/WD_0718_voom_normalized_gene_counts.txt',row.names=F,quote=F,sep='\t')

drop_ind=c('EB.09S.H.00377')

#drop_ind=c('EB.10H.H.00043')
norm=norm[,!names(norm)%in% drop_ind]

norm=as.data.frame(t(norm),stringsAsFactors=F)
fwrite(norm,'eqtl/normalized/WD_0718_voom_normalized_gene_counts_formatted_FIXED.txt',row.names=T,quote=F,sep='\t')

#WD_0720 Outliers

#pc_mm[pc_mm$PC1>= 100,]
#                    PC1       PC2          PC3         sample plate
#EB.09S.H.00373 132.6266 -20.23686  -0.08339753 EB.09S.H.00373    16
#EB.09S.H.00519 132.3031 -22.01188  -7.92784132 EB.09S.H.00519    14
#EB.09S.H.00428 140.9779 -18.18513 -21.16529231 EB.09S.H.00428    14
#EB.09S.H.00093 123.0088  49.48852  -7.67474676 EB.09S.H.00093    14



#norm=norm[,!names(norm)%in% drop_ind]
norm=as.data.frame(t(norm),stringsAsFactors=F)
fwrite(norm,'eqtl/normalized/WD_0720_voom_normalized_gene_counts_formatted_FIXED.txt',row.names=T,quote=F,sep='\t')


#norm=norm[,!names(norm)%in% drop_ind]
#fwrite(norm,'eqtl/WD_0720_voom_normalized_gene_counts.txt',row.names=F,quote=F,sep='\t')
norm=as.data.frame(t(norm),stringsAsFactors=F)
fwrite(norm,'eqtl/normalized/WD_0727_voom_normalized_gene_counts_formatted_FIXED.txt',row.names=T,quote=F,sep='\t')


#mm=model.matrix(~1 + genotype)
norm=fread(sprintf('eqtl/%s_voom_normalized_gene_counts.txt',time),data.table=F)
resid=fread(sprintf('eqtl/results/eQTL_%s_c9_residuals.txt',time),data.table=F)
i=intersect(norm$V1,colnames(resid)[-1])
norm=norm[norm$V1 %in% i,]
rownames(norm)=norm$V1
norm=norm[i,]
norm=norm[,-1]

rownames(resid)=resid$ID
resid=resid[,-1]
resid=resid[,i]

var_resid=apply(resid,MARGIN=2,function(x) var(x))
var_norm=apply(norm,MARGIN=1,function(x) var(x))

data=data.frame(pos=seq(1,length(var_resid)),var_resid=var_resid,var_norm=var_norm,stringsAsFactors=F)
data$gene=rownames(data)

p1=ggplot(aes(x=pos,y=var_norm),data=data) + geom_point() + ggtitle(sprintf('%s Normalized Variance by gene',time))

p2=ggplot(aes(x=pos,y=var_resid),data=data) + geom_point() + ggtitle(sprintf('%s Residual Variance by gene',time))

png(sprintf('images/%s_norm_var_c9.png',time))
print(p1)
dev.off()

png(sprintf('images/%s_resid_var_c9.png',time))
print(p2)
dev.off()

# Fitting linear models in limma
#fit=lmFit(y,mm)
#saveRDS(fit,sprintf('%s_model.rds',time))
#coefficients=as.data.frame(coef(fit))
#names(coefficients)=c("A632_usa","A654_inra","B73_inra","B96",
#"C103_inra","CO255_inra","D105_inra","DK63","EP1_inra",
#"F492","FV2_inra","FV252_inra","ND245","OH43_inra","VA85","W117_inra")
#fwrite(coefficients,sprintf('%s_normalized_coefficients.txt',time),row.names=T,quote=F,sep='\t')

# Checking what gene counts looked like before filtering for the cis-eQTL we did find
#
#
#
#
######
# Original filtering
mat=fread(sprintf('star/%s_updated_gene_counts.txt',time),data.table=F)
rownames(mat)=mat$Gene_ID
mat=mat[,-1]
meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
snames=colnames(mat)
gnames=meta[match(snames,meta$sample_name),]$genotype
mat=as.matrix(mat)
coldata=data.frame(sample=snames,genotype=gnames)

d0=DGEList(mat)
d0=calcNormFactors(d0)

cutoff=1
drop=which(apply(cpm(d0),1,max)<cutoff)
d=d0[-drop,]
dim(d)

snames=colnames(mat)

genotype=coldata$genotype
mm=model.matrix(~1,coldata)
y=voom(d,mm,plot=T)
norm=as.data.frame(y$E)

snames=names(norm)
gnames=meta[match(snames,meta$sample_name),]$dh_genotype
norm=as.data.frame(t(norm),stringsAsFactors=F)
norm$ID=gnames
#####

maxes=apply(cpm(d0),1,max)

ciseqtl=fread('eqtl/results/all_cis_eQTL_fkeep_hits.txt',data.table=F)
ciseqtl=ciseqtl[ciseqtl$time==time,]

drop_ind=c('EB.09S.H.00041','EB.09S.H.00133','EB.10H.H.00087')
norm=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)


plot_list=list()
count=1
for(i in 1:nrow(ciseqtl)){
  gene=ciseqtl[i,]$Gene
  if(gene %in% names(norm)){

    #print("True")
    ex=data.frame(ID=norm$V1,exp=unlist(norm[,gene]),stringsAsFactors=F)
    ex=ex[ex$ID %in% inter,]
    ex=ex[order(ex$exp),]
    rownames(ex)=seq(1,nrow(ex))
    ex$ID_f=factor(ex$ID,levels=c(unique(ex$ID)))
    #subex=subset(ex, ID %in% drop_ind)
    X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
    colnames(X) = founders
    rownames(X) = dimnames(X_list[[1]])[[1]]
    X=X[inter,]
    founders=apply(X,MARGIN=1,function(x) names(x[which.max(x)]))
    ex$founder=founders[match(ex$ID,names(founders))]
    p=ggplot(aes(x=ID_f,y=exp),data=ex) + 
    geom_point(aes(color=founder)) +
     xlab('Sample') +
     ylab('Expression (log2CPM)') + geom_hline(yintercept=1)
    plot_list[[count]]=p
    count=count+1
  }
}

pdf('images/founder_eqtl_by_ind.pdf')
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()

# Testing about presence and absence of MITE
#has=design_matrix$mite
#has=factor(has,levels=c("MITE","NO_MITE"))
#mm2=model.matrix(~0+has)
#y2=voom(d,mm2,plot=T)

# Fitting linear models in limma
#fit2=lmFit(y2,mm2)
#saveRDS(fit2,sprintf('eqtl/%s_limma_voom_model.rds',time))
#coefficients2=as.data.frame(coef(fit2))
#names(coefficients2)=c("MITE","NO_MITE")
#fwrite(coefficients2,sprintf('%s_normalized_coefficients.txt',time),row.names=T,quote=F,sep='\t')


#contr=makeContrasts(hasNO_MITE-hasMITE,levels=colnames(coef(fit2)))

#tmp=contrasts.fit(fit2,contr)
#tmp=eBayes(tmp)

#top.table=topTable(tmp,sort.by="P",n="Inf")
#length(which(top.table$adj.P.Val < 0.05))

#top.table$Gene <- rownames(top.table)
#top.table <- top.table[,c("Gene", names(top.table)[1:6])]
#write.table(top.table, file = "NO_MITE_v_MITE_DEGs.txt", row.names = F, sep = "\t", quote = F)


# Correlation of plate with PC
time="WD_0727"
pcs=fread(sprintf('eqtl/normalized/%s_PCA_covariates_2.txt',time),data.table=F)
meta=fread('metadata/samples_passed_genotype_check.txt',data.table=F)

#meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
meta=meta[meta$experiment==time,]
#meta=meta[meta$read==1,]
meta=meta[match(pcs$V1,meta$dh_genotype),]

cor(meta$plate,pcs$PC1)
cor(meta$plate,pcs$PC2)
cor(meta$plate,pcs$PC3)

# WD_0712
#[1] 0.8159085
#[1] -0.08469473
#[1] 0.2859349

#WD0718
#[1] -0.03237676
#[1] -0.5439166
#[1] -0.04923342

#WD0720
#[1] 0.2180963
#[1] 0.1039951
#[1] -0.1608267

#WD_0727
#[1] 0.02979704
#[1] 0.04257913
#[1] -0.3597453