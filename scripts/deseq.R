#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
timep=as.character(args[[1]])
#chr=as.character(args[[2]])
#cores=as.numeric(args[[3]])

library('data.table')
library('DESeq2')

#timep='WD_0720'
# WD_0718, WD_0720, WD_0727, WD_0712
mat=fread(sprintf('star/%s_gene_counts.txt',timep),data.table=F)
rownames(mat)=mat$Gene_ID
mat=mat[,-1]
meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
snames=colnames(mat)
gnames=meta[match(snames,meta$sample_name),]$genotype
mat=as.matrix(mat)
coldata=data.frame(sample=snames,genotype=gnames)

dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = coldata,
                              design = ~ 1)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#dds2 = DESeq(dds)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

vsd <- vst(dds, blind=TRUE)

cts = fread(sprintf('star/vst_%s_gene_counts_results2.csv',timep),data.table=F)
rownames(cts)=cts$V1
cts=cts[,-1]
cts=t(cts)
pca <- prcomp(cts, center = TRUE)
pcs=as.data.frame(pca$x,stringsAsFactors=F)
snames=rownames(cts)
gnames=meta[match(snames,meta$sample_name),]$genotype
#mat=as.matrix(mat)
coldata=data.frame(sample=snames,genotype=gnames,stringsAsFactors=F)
coldata$pc1=pcs[match(coldata$sample,rownames(pcs)),]$PC1
coldata$pc2=pcs[match(coldata$sample,rownames(pcs)),]$PC2

cts=t(cts)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ pc1 + pc2)

#normalized_counts <- counts(vsd, normalized=TRUE)
#write.csv(as.data.frame(normalized_counts),
#          file=sprintf("star/vst_%s_gene_counts_results2.csv",timep))
write.csv(as.data.frame(assay(vsd)),
          file=sprintf("star/vst_%s_gene_counts_results2.csv",timep))

#ntd <- normTransform(dds)
library("vsn")
#p1=meanSdPlot(assay(ntd))
#png(sprintf('images/%s_normTransform_meanSdPlot.png',timep))
#print(p1)
#dev.off()

p2=meanSdPlot(assay(vsd))
png(sprintf('images/%s_varStabilized_meanSdPlot2.png',timep))
print(p2)
dev.off()
