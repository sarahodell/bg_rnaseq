#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])
time1=as.character(args[[2]])
cores=as.numeric(args[[3]])

library('data.table')
library('dplyr')
library('stringr')
library('ggplot2')
library('parallel')
library('MASS')


hapsnps=fread(sprintf('datasets/hapmap/chr%s_hapmap_biogemma_rare_alleles_ref.txt',chr),data.table=F)
names(hapsnps)=c("SNP","CHR","POS","MAJ","MIN")
hapsnps=hapsnps[,c("SNP","CHR","POS","MAJ","MIN")]

founder=fread(sprintf('datasets/hapmap/chr%s_founder_rare_alleles.txt',chr),data.table=F)
names(founder)=c("SNP","CHR","POS","REF","ALT","A632_usa","A654_inra","B73_inra","B96","C103_inra","CO255_inra","D105_inra","DK63","EP1_inra","F492", "FV252_inra", "FV2_inra","MBS847","ND245","OH43_inra","VA85","W117_inra")

founders=c("MBS847","B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
founder=founder[,c("SNP","CHR","POS","REF","ALT",founders)]

founder[,founders]=apply(founder[,founders],MARGIN=2,function(x) ifelse(x=="0/0",0,ifelse(x=="1/1",1,ifelse(x=="2/2",2,NA))))
# remove multiallelic sites?
founder$ALT2=NA
mult=which(grepl(',',founder$ALT))
alts=founder$ALT[mult]
alt1=sapply(seq(1,length(alts)),function(x) strsplit(alts[x],',')[[1]][[1]])
alt2=sapply(seq(1,length(alts)),function(x) strsplit(alts[x],',')[[1]][[2]])
founder[mult,]$ALT=alt1
founder[mult,]$ALT2=alt2

founder=founder[,c("SNP","CHR","POS","REF","ALT","ALT2",founders)]

overlap=merge(hapsnps,founder,by.x='POS',by.y='POS')

# remove reference minor alleles 
b73alts=unique(overlap[!is.na(overlap$B73_inra) & overlap$B73_inra!=0,]$B73_inra)
overlap=overlap[overlap$B73_inra!=1,]

# remove multi-allelic sites
overlap=overlap[is.na(overlap$ALT2),]

min_founder=readRDS(sprintf('datasets/hapmap/minor_allele_founders_c%s.rds',chr))

tmpdf=fread(sprintf('eqtl/results/eQTL_%s_freq_chi_data.txt',time1),data.table=F)
tmpdf$gene_time_snp=paste0(tmpdf$Trait,'-',tmpdf$time,'-',tmpdf$X_ID)

#top5k=unique(tmpdf[tmpdf$rank<=5000,]$Trait)
#top5k=unique(tmpdf[tmpdf$rank>5000,]$Trait)

#ranks=fread(sprintf('eqtl/results/%s_all_exp_ranks.txt',time1),data.table=F)
#ranks=fread(,sprintf('eqtl/results/%s_top5k_exp_ranks.txt',time1),data.table=F)
#rownames(ranks)=ranks$V1

upstream=fread('eqtl/data/upstream_gene_list.txt',data.table=F)
upstream=upstream[upstream$CHROM==chr,]
chr5k=intersect(tmpdf$Trait,upstream$Gene_ID)

#ranks=ranks[,chr5k]

#bimbam=fread(sprintf('datasets/chr%s_biogemma_rare_allele_probs.txt',chr),data.table=F)

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
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


samples=fread(sprintf('eqtl/data/%s_samples_FIXED.txt',time1),data.table=F)
inds=samples$id
rmelt=as.data.frame(expand.grid(chr5k,inds),stringsAsFactors=F)
names(rmelt)=c("Gene_ID","ID")
cutoff=0.75

n_reps=1:nrow(rmelt)


#n_reps=100:200
get_rare_counts=function(rep){
	row1=rmelt[rep,]
	gene=as.character(row1$Gene_ID)
	id=as.character(row1$ID)
	#r=ranks[id,gene]
	rstart=upstream[upstream$Gene_ID==gene,]$UP_START
	rend=upstream[upstream$Gene_ID==gene,]$UP_END
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
		if(max_f %in% names(betas)){
			beta_rank=match(max_f,names(betas))
			beta1=betas[max_f]
			t1=overlap[overlap$POS>=rstart & overlap$POS<=rend,]
			t1=t1[!is.na(t1$POS),]
			t1=t1[,c('SNP.x','ALT','REF',max_f)]
			if(nrow(t1)!=0){
				rcount=sum(t1[,max_f],na.rm=T)
			}else{
				rcount=0
			}
		}else{
			beta1=Inf
			beta_rank=NA
			rcount=0
		}	
	}else{
		max_f=NA
		beta1=NA
		beta_rank=NA
		rcount=NA
		#add_rank=NA
	}

	line=data.frame(time=time1,chr=chr,Gene_ID=gene,ID=id,snp=snp,add_rank=add_rank,rare_count=rcount,max_f=max_f,beta=beta1,beta_rank=beta_rank,stringsAsFactors=F)
	rownames(line)=1
	return(line)
}

print(system.time({
dresults=mclapply(n_reps,get_rare_counts,mc.cores=cores)
}))
d=rbindlist(dresults)
#saveRDS(results,sprintf('eqtl/trans/permute/chr%s_%s_%s_%.0frep_min_pvalues.rds',chr,time,factor,reps))
#all_props=do.call(rbind,lapply(results,function(x) x))
d=as.data.frame(d,stringsAsFactors=F)

#fwrite(d,sprintf('eqtl/results/%s_%s_5kb_rare_counts_v3.txt',time1,chr),row.names=T,quote=F,sep='\t')
fwrite(d,sprintf('eqtl/results/%s_%s_5kb_rare_counts_all_v3.txt',time1,chr),row.names=T,quote=F,sep='\t')

#fwrite(d,sprintf('eqtl/results/%s_%s_5kb_rare_counts_lower_exp.txt',time1,chr),row.names=F,quote=F,sep='\t')
