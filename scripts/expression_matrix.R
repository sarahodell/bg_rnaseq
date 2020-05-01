#!/usr/bin/env Rscript

library('data.table')

# Make gene expression count matrix of all samples

sample_table=fread('samples.tsv',data.table=F)
samples=sample_table[sample_table$include==TRUE,]$sample

#convert sequencing sample names to sample names from experiment
conversion=fread('raw_reads/metadata/BG_sequencing_sample_conversion_table.txt',data.table=F)

read_counts=c()
# need to figure out strandedness of data
gene_names=c()
wanted_samples=c()
others=c("Blanco_lab","")
for(s in samples){
  if(!(conversion[conversion$seq_name==sprintf('%s_R1_001',s),]$experiment %in% others)){
    reads=fread(sprintf('%s_pass2/ReadsPerGene.out.tab',s),data.table=F,header=F)
    readcol=reads[5:dim(reads)[1],'V2']
    read_counts=cbind(read_counts,readcol)
    gene_names=reads$V1[5:dim(reads)[1]]
    wanted_samples=c(wanted_samples,s)
  }
}

full_samples=sprintf('%s_R1_001',wanted_samples)
#full_samples=conversion[conversion$experiment != "Blanco_lab"]$seq_name
new_names=conversion[match(full_samples,conversion$seq_name),]$sample

read_counts=as.data.frame(read_counts)
names(read_counts)=new_names
rownames(read_counts)=gene_names

fwrite(read_counts,'batch_1_raw_readcount_matrix.txt',row.names=T,col.names=T,sep='\t',quote=F)
