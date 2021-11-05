#!/usr/bin/env Rscript

library('data.table')
library('tibble')

reads=fread('limma_results/batch_1_raw_readcount_matrix.txt',data.table=F)
rownames(reads)=reads$V1
reads=reads[,2:dim(reads)[2]]
ind=rownames_to_column(ind,"marker")
zcn8="Zm00001d010752"
rap27="Zm00001d010987"
mads69="Zm00001d042315"
goi=c(zcn8,rap27,mads69)

reads=reads[goi,]

zcn8_rap27=cor(unlist(unname(reads[goi[1],])),unlist(unname(reads[goi[2],])))
zcn8_mads69=cor(unlist(unname(reads[goi[1],])),unlist(unname(reads[goi[3],])))
mads69_rap27=cor(unlist(unname(reads[goi[3],])),unlist(unname(reads[goi[2],])))
