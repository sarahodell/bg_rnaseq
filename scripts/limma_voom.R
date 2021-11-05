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
design_matrix=fread('vgt1_founder_test_WD_0712.txt',data.table=F)
#Add design matrix based off of founder identity at vgt1
founder=design_matrix$founder

mm=model.matrix(~0+founder)
y=voom(d,mm,plot=T)

# Fitting linear models in limma
fit=lmFit(y,mm)
saveRDS(fit,'vgt1_founder_model.rds')
coefficients=as.data.frame(coef(fit))
names(coefficients)=c("A632_usa","A654_inra","B73_inra","B96",
"C103_inra","CO255_inra","D105_inra","DK63","EP1_inra",
"F492","FV2_inra","FV252_inra","ND245","OH43_inra","VA85","W117_inra")
fwrite(coefficients,'batch_1_normalized_coefficients.txt',row.names=T,quote=F,sep='\t')


# Testing about presence and absence of MITE
has=design_matrix$mite
has=factor(has,levels=c("MITE","NO_MITE"))
mm2=model.matrix(~0+has)
y2=voom(d,mm2,plot=T)

# Fitting linear models in limma
fit2=lmFit(y2,mm2)
saveRDS(fit2,'vgt1_MITE_model.rds')
coefficients2=as.data.frame(coef(fit2))
names(coefficients2)=c("MITE","NO_MITE")
fwrite(coefficients2,'batch_1_normalized_coefficients_MITE.txt',row.names=T,quote=F,sep='\t')


contr=makeContrasts(hasNO_MITE-hasMITE,levels=colnames(coef(fit2)))

tmp=contrasts.fit(fit2,contr)
tmp=eBayes(tmp)

top.table=topTable(tmp,sort.by="P",n="Inf")
length(which(top.table$adj.P.Val < 0.05))

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "NO_MITE_v_MITE_DEGs.txt", row.names = F, sep = "\t", quote = F)
