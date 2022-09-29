#!/usr/bin/env Rscript

library('data.table')
library('stringr')
library('dplyr')
library('edgeR')
# Make gene expression count matrix of all samples

time="WD_0712"
zcn8="Zm00001d010752"
rap27="Zm00001d010987"
mads69="Zm00001d042315"
goi=c(zcn8,rap27,mads69)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
mite=c("NO_MITE","MITE","MITE","MITE","NO_MITE","MITE","MITE","MITE","MITE","MITE","MITE","NO_MITE","MITE","MITE","MITE","NO_MITE")
mite_late=c("FV252_inra","C103_inra","F492","A632_usa")

reads=fread(sprintf('star/%s_gene_counts.txt',time),data.table=F)
rownames(reads)=reads$Gene_ID
reads=reads[,-1]
genes=rownames(reads)

samples=colnames(reads)
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
genos=metadata[match(samples,metadata$sample_name),]$genotype
#drop duplicates because they'll be dropped anyway
drop=which(duplicated(genos))
if(length(drop)!=0){
  genos=genos[-drop]
  reads=reads[,-drop]
}
reads=t(reads)
reads=as.data.frame(reads,stringsAsFactors=F)

#reads$ID=genos
#reads=t(reads)

#meanreads=reads %>% group_by(ID) %>% mean()
rownames(reads)=genos
#reads=reads[,goi]


rownames(reads)=genos
reads$ID=rownames(reads)
reads$ID2=metadata[match(reads$ID,metadata$genotype),]$dh_genotype


mite_founder=fread('metadata/vgt1_max_prob_founder.txt',data.table=F)
#mite_founder$ID_2=test_samples[match(mite_founder$ID,test_samples$genotype_2),]$sample

reads$founder=mite_founder[match(reads$ID2,mite_founder$ID),]$founder
reads=reads[!is.na(reads$founder),]
#reads=unique(reads[,c('sample','founder')])


design_matrix=reads[,c('ID2','founder')]
design_matrix$mite=mite[match(design_matrix$founder,founders)]

design_matrix=design_matrix[!(is.na(design_matrix$founder)),]
design_matrix=design_matrix[design_matrix$ID %in% reads$ID2,]
rownames(design_matrix)=seq(1,dim(design_matrix)[1])

#read=reads[,c(design_matrix$ID2)]
fwrite(reads,sprintf('limma_results/%s_full_readcount_matrix.txt',time),row.names=F,col.names=T,sep='\t',quote=F)
fwrite(design_matrix,sprintf('limma_results/%s_vgt1_founder_test.txt',time),row.names=F,quote=F,sep='\t')


design_matrix=fread(sprintf('limma_results/%s_vgt1_founder_test.txt',time),data.table=F)
reads=fread(sprintf('limma_results/%s_full_readcount_matrix.txt',time),data.table=F)
rownames(reads)=reads$ID2
d=matrix(as.numeric(unlist(reads[,genes])),nrow=nrow(reads))
rownames(d)=reads$ID2
colnames(d)=genes
rownames(design_matrix)=design_matrix$ID2
d=t(d)

has=design_matrix$mite
has=factor(has,levels=c("MITE","NO_MITE"))
mm2=model.matrix(~0+has)
#d=t(d)
y2=voom(d,mm2,plot=T)
fit=lmFit(y2,mm2)
saveRDS(fit,sprintf('limma_results/%s_MITE_full_model.rds',time))
contr <- makeContrasts(hasMITE-hasNO_MITE, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
fwrite(top.table,sprintf('limma_results/%s_top_MITE_DEGs.txt',time),row.names=T,quote=F,sep='\t')


f=design_matrix$founder
f=factor(f,levels=founders)

mm1=model.matrix(~0+f)
#d=t(d)
y1=voom(d,mm1,plot=T)
fit=lmFit(y1,mm1)
saveRDS(fit,sprintf('limma_results/%s_founder_model.rds',time))

#MITE_early vs MITE_late
design_matrix2=design_matrix[design_matrix$mite=="MITE",]
design_matrix2$mite=ifelse(design_matrix2$founder %in% mite_late,"MITE_LATE","MITE_EARLY")

has=design_matrix$mite
has=factor(has,levels=c("MITE_EARLY","MITE_LATE"))
mm3=model.matrix(~0+has)
#d=t(d)
y3=voom(d,mm3,plot=T)
fit=lmFit(y3,mm3)
saveRDS(fit,sprintf('limma_results/%s_MITE_late_full_model.rds',time))
contr <- makeContrasts(hasMITE_EARLY-hasMITE_LATE, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
fwrite(top.table,sprintf('limma_results/%s_top_MITE_late_DEGs.txt',time),row.names=T,quote=F,sep='\t')


sig=top.table[top.table$adj.P.Val<=0.05,]
fwrite(sig,sprintf('limma_results/%s_top_MITE_late_sig_DEGS.txt',time),row.names=T,quote=F,sep='\t')
#sample_table=fread('samples.tsv',data.table=F)
#samples=sample_table[sample_table$include==TRUE,]$sample

#convert sequencing sample names to sample names from experiment
#conversion=fread('raw_reads/metadata/BG_sequencing_sample_conversion_table.txt',data.table=F)

#read_counts=c()
# need to figure out strandedness of data
#gene_names=c()
#wanted_samples=c()
#others=c("Blanco_lab","")
#for(s in samples){
#  if(!(conversion[conversion$seq_name==sprintf('%s_R1_001',s),]$experiment %in% others)){
#    reads=fread(sprintf('%s_pass2/ReadsPerGene.out.tab',s),data.table=F,header=F)
#    readcol=reads[5:dim(reads)[1],'V2']
#    read_counts=cbind(read_counts,readcol)
#    gene_names=reads$V1[5:dim(reads)[1]]
#    wanted_samples=c(wanted_samples,s)
#  }
#}

#full_samples=sprintf('%s_R1_001',wanted_samples)
#full_samples=conversion[conversion$experiment != "Blanco_lab"]$seq_name
#new_names=conversion[match(full_samples,conversion$seq_name),]$sample

#read_counts=as.data.frame(read_counts)
#names(read_counts)=new_names
#rownames(read_counts)=gene_names


#test_samples=fread('raw_reads/metadata/expression_testing_samples_WD_0712.txt',data.table=F)
#test_samples$genotype_2 = gsub('-','.',test_samples$genotype_2)



# Fitting linear models in limma
#fit=lmFit(y,mm)
#saveRDS(fit,sprintf('%s_model.rds',time))
#coefficients=as.data.frame(coef(fit))
#names(coefficients)=c("A632_usa","A654_inra","B73_inra","B96",
#"C103_inra","CO255_inra","D105_inra","DK63","EP1_inra",
#"F492","FV2_inra","FV252_inra","ND245","OH43_inra","VA85","W117_inra")
#fwrite(coefficients,sprintf('%s_normalized_coefficients.txt',time),row.names=T,quote=F,sep='\t')




# Testing about presence and absence of MITE


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
