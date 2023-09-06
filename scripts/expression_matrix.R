#!/usr/bin/env Rscript

library('data.table')
library('stringr')
library('dplyr')
library('edgeR')
# Make gene expression count matrix of all samples

time="WD_0712"

zcn8="Zm00001d010752"
rap27_1="Zm00001d010987"
rap27_2="Zm00001d010988"
mads69="Zm00001d042315"
goi=c(zcn8,rap27_1,rap27_2,mads69)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
mite=c("NO_MITE","MITE","MITE","MITE","NO_MITE","MITE","MITE","MITE","MITE","MITE","MITE","NO_MITE","MITE","MITE","MITE","NO_MITE")
mite_late=c("FV252_inra","C103_inra","F492","A632_usa")

mat=fread(sprintf('star/%s_updated_gene_counts.txt',time),data.table=F)
rownames(mat)=mat$Gene_ID
mat=mat[5:nrow(mat),]
mat=mat[,-1]
genes=rownames(mat)
mat=t(mat)
meta=fread('metadata/samples_passed_genotype_check.txt',data.table=F)

#meta=fread('metadata/BG_completed_sample_list_FIXED.txt',data.table=F)
#meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
meta=meta[meta$experiment==time,]
snames=rownames(mat)
gnames=meta[match(snames,meta$sample_name),]$dh_genotype
mat=as.data.frame(mat,stringsAsFactors=F)
mat$gnames=gnames
mat=mat[mat$gnames!="",]
plate=meta[match(mat$gnames,meta$dh_genotype),]$plate

gnames=mat$gnames

#mat=t(mat)
#mat2=mat %>% group_by(gnames) %>% summarise(across(everything(), list(mean)))
mat1 = mat %>% select(-gnames) %>% as.matrix
rownames(mat1)=gnames
mat=t(mat1)
#rownames(mat1) = mat2$gnames
coldata=data.frame(genotype=gnames,plate=plate,stringsAsFactors=F)
coldata$plate=as.factor(coldata$plate)
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


summary(d0$samples$lib.size/1e6)

cutoff=5
nind=15

maxes=apply(cpm(d0),1,max)
# fewer than 2 samples with cpm higher than cutoff
drop=which(apply(cpm(d0),1,function(x) sum(x>=cutoff)<nind))
d=d0[-drop,]
dim(d)

snames=colnames(mat)

mat=as.data.frame(mat,stringsAsFactors=F)

#reads$ID=genos
#reads=t(reads)

#meanreads=reads %>% group_by(ID) %>% mean()
#rownames(mat)=genos
#reads=reads[,goi]



#mat$ID=rownames(mat)
#mat$ID2=metadata[match(mat$ID,metadata$genotype),]$dh_genotype
#mat$plate=metadata[match(samples,metadata$sample_name),]$plate
#mat$plate=as.factor(mat$plate)
mite_prob=fread('phenotypes/mite_probabilities.txt',data.table=F)
rownames(mite_prob)=mite_prob$ID
mite_prob=mite_prob[mite_prob$final>0.9 | mite_prob$final<0.1,]
mite_prob$final=round(mite_prob$final)
#has_mite=mite_prob[mite_prob$final>=0.9,]$ID
#mite_founder=fread('metadata/vgt1_max_prob_founder.txt',data.table=F)

#mite_founder$mite=mite[match(mite_founder$founder,founders)]


#mite_prob=fread('../GridLMM/mite_probabilities.txt',data.table=F)
#mite_founder$final=mite_prob[match(mite_founder$ID,mite_prob$ID),]$final
#rownames(mite_founder)=mite_founder$ID
#mite_founder$ID_2=test_samples[match(mite_founder$ID,test_samples$genotype_2),]$sample

coldata$founder=mite_prob[match(coldata$genotype,mite_prob$ID),]$founder_AX91102985
coldata$mite=mite_prob[match(coldata$genotype,mite_prob$ID),]$final

coldata=coldata[!is.na(coldata$founder),]
d=d[,coldata$genotype]
#coldata=coldata[coldata$mite>=0.9 | coldata$mite<0.1,]
#reads=unique(reads[,c('sample','founder')])


#design_matrix=mat[,c('ID','plate','founder')]

#coldata$mite=mite[match(coldata$founder,founders)]


#read=reads[,c(design_matrix$ID2)]
#fwrite(mat,sprintf('limma_results/%s_full_readcount_matrix.txt',time),row.names=F,col.names=T,sep='\t',quote=F)
fwrite(coldata,sprintf('limma_results/%s_vgt1_founder_test.txt',time),row.names=F,quote=F,sep='\t')


design_matrix=fread(sprintf('limma_results/%s_vgt1_founder_test.txt',time),data.table=F)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
qtl=fread('QTL/all_adjusted_QTL_support_intervals.txt',data.table=F)

design_matrix$mite=factor(design_matrix$mite,levels=c("0","1"))
mm2=model.matrix(~0+plate+mite,design_matrix)
#d=t(d)
y2=voom(d,mm2,plot=T)
fit=lmFit(y2,mm2)
saveRDS(fit,sprintf('limma_results/%s_MITE_full_model.rds',time))
contr <- makeContrasts(mite1-mite0, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
top.table=top.table[top.table$`adj.P.Val`<=0.05,]

top.table$Gene_ID=rownames(top.table)
top.table=merge(top.table,genetable,by="Gene_ID")

env1=top.table
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


sum(top.table$Gene_ID %in% ft_genelist$V4)
fwrite(top.table,sprintf('limma_results/%s_top_MITE_DEGs.txt',time),row.names=T,quote=F,sep='\t')


present=founders[founders %in% design_matrix$founder]
#f=design_matrix$founder
#f=factor(f,levels=founders)
design_matrix$founder_f=factor(design_matrix$founder,levels=present)
mm1=model.matrix(~0+plate+founder_f,design_matrix)
#d=t(d)
y1=voom(d,mm1,plot=T)
fit=lmFit(y1,mm1)
saveRDS(fit,sprintf('limma_results/%s_founder_model.rds',time))

#MITE_early vs MITE_late
design_matrix2=design_matrix[design_matrix$mite=="1",]
design_matrix2$mite=ifelse(design_matrix2$founder %in% mite_late,"MITE_LATE","MITE_EARLY")

#has=design_matrix$mite
design_matrix2$mite=factor(design_matrix2$mite,levels=c("MITE_EARLY","MITE_LATE"))
d2=d[,design_matrix2$genotype]
mm3=model.matrix(~0+plate+mite,design_matrix2)
#d=t(d)
y3=voom(d2,mm3,plot=T)
fit=lmFit(y3,mm3)
saveRDS(fit,sprintf('limma_results/%s_MITE_late_full_model.rds',time))
contr <- makeContrasts(miteMITE_EARLY-miteMITE_LATE, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
top.table=top.table[top.table$`adj.P.Val`<=0.05,]



top.table$Gene_ID=rownames(top.table)
top.table=merge(top.table,genetable,by="Gene_ID")

sum(top.table$Gene_ID %in% ft_genelist$V4)

fwrite(top.table,sprintf('limma_results/%s_top_MITE_late_DEGs.txt',time),row.names=T,quote=F,sep='\t')

env1=top.table
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


time="WD_0720"

top.table=fread(sprintf('limma_results/%s_top_MITE_late_DEGs.txt',time),data.table=F)

allq=fread('QTL/MITE_only/all_MITE_only_QTL_peaks.txt',data.table=F)
env1=top.table
env1=as.data.table(env1)
env2=as.data.table(allq)
setkey(env2,CHR,leftmost,alt_rightmost)
comp=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)

allq=fread('QTL/MITE_only/all_MITE_only_QTL_support_intervals.txt',data.table=F)
env1=top.table
env1=as.data.table(env1)
env2=as.data.table(allq)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


# Shared across timepoitns?

alldegs=c()
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	deg=fread(sprintf('limma_results/%s_top_MITE_DEGs.txt',time),data.table=F)
	deg$time=time
	alldegs=rbind(alldegs,deg)
}
degs=unique(alldegs$Gene_ID)
degs=data.frame(degs=degs,stringsAsFactors=F)
fwrite(degs,sprintf('limma_results/%s_top_MITE_DEGs.txt',time),row.names=F,quote=F,sep='\t')

env1=alldegs
env1=as.data.table(env1)
env2=as.data.table(allqtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

unique(comp$ID)
#[1] "qDTA8"   "qDTS8"   "qTPH8"   "qTKW7_2" "qHGM7"   "qDTA7"   "qTKW7_1"
#[8] "qASI10" 

#46 of the 61 MITE_noMITE DEGs overlap with qtl SIs
# 4 of these are not on CHR 8 - "Zm00001d020692" "Zm00001d001183" "Zm00001d019102" "Zm00001d025407"
# Overlap with qTKW7_1, qTKW7_2, qHGM7, qDTA7, qASI10


alldegs=c()
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	deg=fread(sprintf('limma_results/%s_top_MITE_late_DEGs.txt',time),data.table=F)
	deg$time=time
	alldegs=rbind(alldegs,deg)
}
degs=unique(alldegs$Gene_ID)
degs=data.frame(degs=degs,stringsAsFactors=F)
fwrite(degs,sprintf('limma_results/%s_top_MITE_late_DEGs.txt',time),row.names=F,quote=F,sep='\t')
# 30 DEGs for Early-Late MITE

factoreqtl=fread('eqtl/results/all_residual_factor_fdr_SIs_FIXED.txt',data.table=F)
env1=alldegs
env1=as.data.table(env1)
env2=as.data.table(factoreqtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


yas=c("Zm00001d027874","Zm00001d031092","Zm00001d033602",
"Zm00001d006835","Zm00001d013676","Zm00001d018255","Zm00001d022109")

allqtl=fread('QTL/MITE_only/all_MITE_only_QTL_support_intervals.txt',data.table=F)
env1=alldegs
env1=as.data.table(env1)
env2=as.data.table(allqtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

qtl=fread('QTL/all_adjusted_QTL_all_methods.txt',data.table=F)
env1=alldegs
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


unique(comp$ID)
#"qDTA8" "qTPH8" "qTPH3" "qDTS8"

# Zm00001d042673 overlaps with qTPH3


#sig=top.table[top.table$adj.P.Val<=0.05,]
#fwrite(sig,sprintf('limma_results/%s_top_MITE_late_sig_DEGS.txt',time),row.names=T,quote=F,sep='\t')
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
