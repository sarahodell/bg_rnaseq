#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])
time1=as.character(args[[2]])

library('data.table')
library('ggplot2')
library('dplyr')

# Disregulation smile plot

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10


#time1="WD_0712"

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

#expdf=expdf[order(expdf$avg_logexp),]
#expdf$rank=seq(nrow(expdf),1)
#expdf=merge(expdf,genetable,by="Gene_ID")


#allrares=c()

#rcount=fread(sprintf('eqtl/results/%s_%s_5kb_rare_counts_lower_exp.txt',time1,chr),data.table=F)
rcount=fread(sprintf('eqtl/results/%s_%s_5kb_rare_counts.txt',time1,chr),data.table=F)
rcount$time=time1
rcount$chr=chr
adj_chr=c(5,9)
testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
if(chr %in% adj_chr){
	X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
}else{
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
}
results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time1,chr),data.table=F)
bv=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_bv_FIXED.txt',time1,chr),data.table=F)
rownames(bv)=bv$V1
max_fs=c()
allbetas=c()
beta_ranks=c()
add_ranks=c()
for(i in 1:nrow(rcount)){
	row1=rcount[i,]
	gene=row1$Gene_ID
	id=row1$ID
	add_rank=which(bv[order(bv[,gene],decreasing=TRUE),c('V1',gene)]$V1==id)
	snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps[1]
	fprobs=unlist(lapply(X_list,function(x) x[id,snp]))
	if(max(fprobs)>0.75){
		max_f=names(which.max(fprobs))
		res=results[results$X_ID==snp & results$Trait==gene,]
		betas=unlist(res[,founders])
		wn=which(!is.na(betas))[1]
		betas[-wn]=betas[-wn]+betas[wn]
		betas=sort(betas,decreasing=TRUE)
		beta_rank=match(max_f,names(betas))
		beta1=betas[max_f]
	}else{
		max_f=NA
		beta1=NA
		beta_rank=NA
		#add_rank=NA
	}
	max_fs=c(max_fs,max_f)
	allbetas=c(allbetas,beta1)
	beta_ranks=c(beta_ranks,beta_rank)
	add_ranks=c(add_ranks,add_rank)
}
#max_fs=unlist(sapply(seq(1,nrow(rcount)),function(x) get_maxf(x)))
rcount$max_f=max_fs
rcount$beta=allbetas
rcount$beta_rank=beta_ranks
rcount$add_rank=add_ranks

#allrares=rbind(allrares,rcount)
#allrares=as.data.frame(allrares,stringsAsFactors=F)

fwrite(rcount,sprintf('eqtl/results/rare_counts_%s_c%s_max_f.txt',time1,chr),row.names=F,quote=F,sep='\t')

#fwrite(rcount,sprintf('eqtl/results/rare_counts_%s_c%s_max_f_lower_exp.txt',time1,chr),row.names=F,quote=F,sep='\t')

