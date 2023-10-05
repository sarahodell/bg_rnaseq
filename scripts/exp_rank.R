#!/usr/bin/env Rscript

library('data.table')

time1="WD_0712"


log_inverse=function(x){
  	return(2^x)
}
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
#exp=exp[,c('V1',top5k)]
rownames(exp)=exp$V1
exp=exp[,-1]

unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
avg_exp = apply(unlog,2,mean)
names(avg_exp)=names(exp)
avg_exp[avg_exp<1]=0
	# re-log the input data
avg_logexp=log2(avg_exp)
avg_logexp[is.infinite(avg_logexp)]=0

expdf=data.frame(Gene_ID=names(avg_logexp),avg_logexp=avg_logexp,stringsAsFactors=F)
expdf=expdf[order(expdf$avg_logexp),]
expdf$rank=seq(nrow(expdf),1)

expdf=merge(expdf,genetable,by="Gene_ID")


inds=rownames(exp)
genes=names(exp)
find_rank=function(gene){
	sub=exp[order(exp[,gene],decreasing=TRUE),gene,drop=F]
	r=match(inds,rownames(sub))
	return(r)
}
ranks=as.data.frame(unlist(sapply(genes,function(x) find_rank(x))),stringsAsFactors=F)
rownames(ranks)=inds
colnames(ranks)=genes

fwrite(ranks,sprintf('eqtl/results/%s_all_exp_ranks.txt',time1),row.names=T,quote=F,sep='\t')