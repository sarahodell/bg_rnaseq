#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')

trans=fread('eqtl/results/all_trans_fdr_SIs_FIXED.txt',data.table=F)
dim(trans)
length(unique(trans$gene))

trans %>% group_by(time) %>% summarize(length(unique(gene)))
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

trans=merge(trans,genetable,by.x='gene',by.y='Gene_ID')

dim(trans[trans$CHR!=trans$CHROM,])
trans$gene_snp=paste0(trans$gene,'-',trans$SNP)

ntimes=trans %>% group_by(gene_snp) %>% count()

ngenes2=trans %>% group_by(SNP,time) %>% count()

p=ggplot(ngenes2,aes(x=n)) + geom_histogram() + xlab("Number of Genes per Variant") +
ylab("Frequency") + ggtitle('Distal-eQTL Hotpsots')

png('paper_figures/distal_eQTL_hotspots_hist.png')
print(p)
dev.off()