#!/usr/bin/env Rscript

library('data.table')
library('tibble')

#reads=fread('limma_results/batch_1_raw_readcount_matrix.txt',data.table=F)
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
  reads=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts.txt',time),data.table=F)
  rownames(reads)=reads$V1
  reads=reads[,-1]
  samples=colnames(reads)
  metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
  genos=metadata[match(samples,metadata$sample_name),]$genotype
  colnames(reads)=genos

  zcn8="Zm00001d010752"
  rap27="Zm00001d010987"
  mads69="Zm00001d042315"
  goi=c(zcn8,rap27,mads69)

  reads=reads[goi,]
  print(time)
  zcn8_rap27=cor(unlist(unname(reads[goi[1],])),unlist(unname(reads[goi[2],])))
  print(zcn8_rap27)
  zcn8_mads69=cor(unlist(unname(reads[goi[1],])),unlist(unname(reads[goi[3],])))
  print(zcn8_mads69)
  mads69_rap27=cor(unlist(unname(reads[goi[3],])),unlist(unname(reads[goi[2],])))
  print(mads69_rap27)
}
