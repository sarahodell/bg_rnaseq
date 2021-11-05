#!/usr/bin/env Rscript

library('data.table')
library('stringr')
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


test_samples=fread('raw_reads/metadata/expression_testing_samples_WD_0712.txt',data.table=F)
test_samples$genotype_2 = gsub('-','.',test_samples$genotype_2)
mite_founder=fread('raw_reads/metadata/vgt1_max_prob_founder.txt',data.table=F)
#mite_founder$ID_2=test_samples[match(mite_founder$ID,test_samples$genotype_2),]$sample

test_samples$founder=mite_founder[match(test_samples$genotype_2,mite_founder$ID),]$founder
test_samples=unique(test_samples[,c('sample','founder')])

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
mite=c("MITE","NO_MITE","MITE","MITE","NO_MITE","MITE","MITE","MITE","MITE","MITE","MITE","NO_MITE","MITE","MITE","MITE","NO_MITE")

design_matrix=test_samples[,c('sample','founder')]
design_matrix$mite=mite[match(design_matrix$founder,founders)]
design_matrix=design_matrix[!(is.na(design_matrix$founder)),]
design_matrix=design_matrix[design_matrix$sample %in% names(read_counts),]
rownames(design_matrix)=seq(1,dim(design_matrix)[1])

read_counts=read_counts[,c(design_matrix$sample)]
fwrite(read_counts,'batch_1_raw_readcount_matrix.txt',row.names=T,col.names=T,sep='\t',quote=F)


fwrite(design_matrix,'vgt1_founder_test.txt',row.names=F,quote=F,sep='\t')
