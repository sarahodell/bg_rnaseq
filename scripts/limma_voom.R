#!/usr/bin/env Rscript

library('data.table')
library('edgeR')

# Make gene expression count matrix of all samples

sample_table=fread('samples.tsv',data.table=F)
samples=sample_table[sample_table$include==TRUE,]$sample

#convert sequencing sample names to sample names from experiment
conversion=fread('raw_reads/metadata/BG_sequencing_sample_conversion_table.txt',data.table=F)

# Preprocessing

counts=read.delim('batch_1_raw_readcount_matrix.txt',row.names=1)
d0=DGEList(counts)
d0=calcNormFactors(d0)

cutoff=1
drop=which(apply(cpm(d0),1,max)<cutoff)
d=d0[-drop,]
dim(d)

snames=colnames(counts)


# Voom transformation and calculation of variance weights

#Add design matrix based off of founder identity at vgt1
group=c()
mm=model.matrix(~0+group)
y=voom(d,mm,plot=T)

# Fitting linear models in limma
fit=lmFit(y,mm)

coefficients=as.data.frame(coef(fit))
fwrite(coefficients,'batch_1_normalized_coefficients.txt',row.names=T,col.names=F,quote=F,sep='\t')
