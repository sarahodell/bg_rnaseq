#!/usr/bin/env Rscript

library('data.table')
library('dplyr')


time="WD_0720"
X_list=readRDS('../genotypes/probabilities/geno_probs/bg10_filtered_genotype_probs.rds')
inds=rownames(X_list[[1]])


times=c("WD_0712","WD_0718","WD_0720","WD_0727")

for(time in times){
	phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	genos=phenotypes$V1

	inter=intersect(genos,inds)
	print(length(inter))
	metadata=fread('metadata/samples_passed_genotype_check.txt',data.table=F)
	metadata=metadata[metadata$experiment==time,]
	#metadata=metadata[metadata$read==1,]
	metadata=metadata[metadata$dh_genotype %in% inter,]
	print(dim(metadata))
	allsamples=metadata$dh_genotype
	allpaths=metadata$bam_file
	data=data.frame(path=allpaths,sample=allsamples,stringsAsFactors=F)
	fwrite(data,sprintf('final_bams/%s_path_samples_FIXED.txt',time),row.names=F,quote=F,sep=',',col.names=F)
}

#WD_0712 83
#WD_0718 143
#WD_0720 220
#WD_0727 194

#tsamples=data.frame(samples=inter,stringsAsFactors=F)
#fwrite(tsamples,sprintf('eqtl/data/%s_samples.txt',time),row.names=F,col.names=F,sep=',',quote=F)

metadata=fread('metadata/samples_passed_genotype_check.txt',data.table=F)
metadata=metadata[metadata$experiment==time,]
#metadata=metadata[metadata$read==1,]
metadata=metadata[metadata$dh_genotype %in% inter,]
print(dim(metadata))
#metadata
dirs=unique(metadata$batch)

#dirs=c('18048-85-lane12-16','batch_2','18048-85-04-11','18048-85-01-03','batch_1')
all_files=c()
#d=dirs[2] #batch_2 only
  #sample=unique(samples[samples$batch==d,]$sample_name)
  
#allpaths=c()
#allsamples=c()
#for(d in dirs){
#	submeta=metadata[metadata$batch==d,]
#	filenames=unique(submeta$sample_name)
#	paths=paste0('final_bams/',d,'/',filenames,'.Aligned.sortedByCoord.MKDup.Processed.out.bam')
# 	samplenames=submeta[match(filenames,submeta$sample_name),]$dh_genotype
#	allpaths=c(allpaths,paths)
#	allsamples=c(allsamples,samplenames)
#}
allsamples=metadata$dh_genotype
allpaths=metadata$bam_file
data=data.frame(path=allpaths,sample=allsamples,stringsAsFactors=F)
fwrite(data,sprintf('final_bams/%s_path_samples_FIXED.txt',time),row.names=F,quote=F,sep=',',col.names=F)


#### Get Gene Bed Files
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

allhits=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
#allhits=fread('eqtl/results/all_eQTL_hits.txt',data.table=F)
allhits=allhits[allhits$time==time,]
allhits=allhits[allhits$class!='factor',]

genebed=data.frame(Gene_ID=unique(allhits$Trait),stringsAsFactors=F)
#WD_0720
for(time in times){
	genebed=data.frame(Gene_ID=c("Zm00001d011152","Zm00001d010752","Zm00001d010987","Zm00001d010988","Zm00001d042315"))
	#genebed=rbind(genebed,addons)
	if(time!="WD_0712"){
		genebed=rbind(genebed,c('Zm00001d011189'))
		genebed=rbind(genebed,c('Zm00001d010946'))
	}
	genebed$CHR=genetable[match(genebed$Gene_ID,genetable$Gene_ID),]$CHROM
	genebed$START=genetable[match(genebed$Gene_ID,genetable$Gene_ID),]$START
	genebed$END=genetable[match(genebed$Gene_ID,genetable$Gene_ID),]$END
	genebed=genebed[,c('CHR','START','END','Gene_ID')]
	fwrite(genebed,sprintf('eqtl/results/%s_genes_FIXED.bed',time),row.names=F,quote=F,col.names=F,sep='\t')

}

#if(time=="WD_0727"){
#	genebed=rbind(genebed,c('Zm00001d011513'))
#	genebed=rbind(genebed,c('Zm00001d022126'))
#	genebed=rbind(genebed,c('Zm00001d022141'))
#}



# For each of our samples, what founder prob to they have at this gene?
time="WD_0718"

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

genebed=fread(sprintf('eqtl/results/%s_genes_FIXED.bed',time),data.table=F)
samples=fread(sprintf('eqtl/data/%s_samples_FIXED.txt',time),data.table=F)
for(i in 1:nrow(genebed)){
	chr=genebed[i,]$V1
	start=genebed[i,]$V2
	end=genebed[i,]$V3
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',chr))
	X_list=lapply(X_list,function(x) x[samples$id,])
	#X_list_ordered=lapply(X_list,function(x) x[,snp,drop=F])
    
	pmap=pmap[pmap$marker %in% colnames(X_list[[1]]),]
	pmap$end=pmap$pos+1
	lower=pmap[pmap$pos<=start,]
	leftindex=which.min(start-lower$pos)
	leftmarker=lower[leftindex,]
	
	#upper=pmap[pmap$pos>start,]
	#rightindex=which.min(upper$pos-start)
	#upper=upper[rightindex,]
	
	X = do.call(cbind,lapply(X_list,function(x) x[,leftmarker$marker]))
	frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
    fkeep=founders[frep2>3]
    #X=X[,fkeep]
    certain=apply(X,MARGIN=1,function(x) sum(x>0.75)>0)
    X=X[certain,]
	founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))
	gene=genebed[i,]$V4
	newdf=data.frame(gene=gene,chr=chr,start=start,end=end,marker=leftmarker$marker,sample=rownames(X),founder=founder,stringsAsFactors=F)
	fwrite(newdf,sprintf('eqtl/results/%s_%s_founder_ident_FIXED.txt',time,gene),row.names=F,quote=F,sep='\t')
}


times=c("WD_0712","WD_0718","WD_0720","WD_0727")
fullbed=c()
for(time in times){
	bed=fread(sprintf('eqtl/results/%s_genes_FIXED.bed',time),data.table=F)
	bed$time=time
	fullbed=rbind(fullbed,bed)
}

fullbed=as.data.frame(fullbed,stringsAsFactors=F)
fwrite(fullbed,'eqtl/results/all_eQTL_genes_FIXED.bed',row.names=F,quote=F,col.names=F,sep='\t')

genes=unique(fullbed$V4)
#genes=c(genes,)
gtf=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46.gtf',data.table=F)
subgtf=gtf[grepl(paste(genes, collapse='|'), gtf$V9),]

fwrite(subgtf,'eqtl/results/all_eQTL_genes_FIXED.gtf',row.names=F,quote=F,sep='\t',col.names=F)


# For each gene, group by founder at location
time="WD_0727"

genebed=fread(sprintf('eqtl/results/%s_genes_FIXED.bed',time),data.table=F)
data=fread(sprintf('final_bams/%s_path_samples_FIXED.txt',time),data.table=F,header=F)
for(i in 1:nrow(genebed)){
	row=genebed[i,]
	query=paste0(row$V1,':',row$V2,'-',row$V3)
	gene=row$V4
	newdf=fread(sprintf('eqtl/results/%s_%s_founder_ident_FIXED.txt',time,gene),data.table=F)
	for(f in unique(newdf$founder)){
		subsamp=unique(newdf[newdf$founder==f,]$sample)
		#bamlist=data.frame(path=data[data$V2 %in% subsamp,]$V1)
		bamlist=data.frame(path=sprintf('/home/sodell/projects/biogemma/expression/final_bams/alignment_check/%s_%s.sorted.bam',time,subsamp))

		#outfile="final_bams/alignment_check/${time}_${sample}.bam"
		#bamlist$path=sprintf(,bamlist$path)
		fwrite(bamlist,sprintf('final_bams/alignment_check/%s_%s_%s.bam.list',time,gene,f),row.names=F,col.names=F,quote=F,sep=',')
	}
}
