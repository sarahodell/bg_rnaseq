#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])

library('data.table')
library('edgeR')

# Make gene expression count matrix of all samples
#timep="WD_0712"
mat=fread(sprintf('star/%s_gene_counts.txt',time),data.table=F)
rownames(mat)=mat$Gene_ID
mat=mat[,-1]
meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
snames=colnames(mat)
gnames=meta[match(snames,meta$sample_name),]$genotype
mat=as.matrix(mat)
coldata=data.frame(sample=snames,genotype=gnames)

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

cutoff=1
drop=which(apply(cpm(d0),1,max)<cutoff)
d=d0[-drop,]
dim(d)

snames=colnames(mat)


# Voom transformation and calculation of variance weights
#design_matrix=fread('vgt1_founder_test_WD_0712.txt',data.table=F)
#Add design matrix based off of founder identity at vgt1
#founder=design_matrix$founder
genotype=coldata$genotype
mm=model.matrix(~1,coldata)
y=voom(d,mm,plot=T)
saveRDS(y,sprintf('eqtl/%s_voom_results.rds',time))
norm=as.data.frame(y$E)
fwrite(norm,sprintf('eqtl/%s_voom_normalized_gene_counts.txt',time),row.names=T,quote=F,sep='\t')
weights=as.data.frame(y$weights)
rownames(weights)=rownames(norm)
names(weights)=names(norm)
fwrite(weights,sprintf('eqtl/%s_voom_weights.txt',time),row.names=T,quote=F,sep='\t')

ym=as.matrix(y)
ym=t(ym)

pca <- prcomp(ym, center = TRUE)
pcs=as.data.frame(pca$x,stringsAsFactors=F)
var_explained <- pca$sdev^2/sum(pca$sdev^2)

pc_mm=pcs[,1:3]
pc_mm$sample=rownames(pc_mm)
fwrite(pc_mm,sprintf('eqtl/%s_PCA_covariates.txt',time),row.names=T,quote=F,sep='\t')
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
