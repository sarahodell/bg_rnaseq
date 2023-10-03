#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')
#chr="10"

#bimbam=fread(sprintf('datasets/chr%s_biogemma_rare_allele_probs.txt',chr),data.table=F)

#inds=rownames(ranks)
#rmelt=as.data.frame(expand.grid(chr5k,inds),stringsAsFactors=F)
#names(rmelt)=c("Gene_ID","ID")
cutoff=0.75

#n_reps=1:nrow(rmelt)

get_rare_counts=function(id){
	#row1=rmelt[rep,]
	#gene=as.character(row1$Gene_ID)
	#id=as.character(row1$ID)
	#r=ranks[id,gene]
	
	#rend=genetable[genetable$Gene_ID==gene,]$START
	#rstart=rend-5000
	#t=bimbam[bimbam$pos>=rstart & bimbam$pos<=rend,]
	subt=t[,c('marker','alt1','ref',id)]
	t[,id]=t[,id]/2
	rcount=sum(t[,id]>=cutoff)
	max_f=names(which.max(unlist(lapply(X_list,function(x) x[id,snp]))))
	line=data.frame(time=time1,factor=factor1,ID=id,rare_count=rcount,founder=max_f,stringsAsFactors=F)
	return(line)
}

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

factoreqtl=fread('eqtl/results/all_residual_factor_fdr_SIs_FIXED.txt',data.table=F)

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}
factoreqtl$block_start=all_founder_blocks[match(factoreqtl$SNP,all_founder_blocks$focal_snp),]$start
factoreqtl$block_end=all_founder_blocks[match(factoreqtl$SNP,all_founder_blocks$focal_snp),]$end


factordf=c()

for(i in 1:nrow(factoreqtl)){
	row1=factoreqtl[i,]
	time1=row1$time
	factor1=row1$factor
	chr=row1$CHR
	
	adj_chr=c("5","5")

	if(chr %in% adj_chr){
		X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))

	}else{
		X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	#chr 5 "AX-91671957" replaced with "AX-91671943"
	}
	snp=row1$SNP
	results=fread(sprintf('eqtl/trans/results/%s_residuals_c%s_%s_trans_results_FIXED.txt',time1,chr,factor1),data.table=F)
	results=results[results$X_ID==snp,]
	betas=results[,founders]
	bimbam=fread(sprintf('datasets/chr%s_biogemma_rare_allele_probs.txt',chr),data.table=F)
	Fvalue=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_all_F_means_FIXED.txt',time1),data.table=F)
	Fvalue=Fvalue[,c('V1',factor1)]
	Fvalue$rank=match(Fvalue[,factor1],sort(Fvalue[,factor1],decreasing=TRUE))
	rstart=row1$block_start
	rend=row1$block_end
	#rstart=row1$left_bound_bp
	#rend=row1$alt_right_bound_bp
	t=bimbam[bimbam$pos>=rstart & bimbam$pos<=rend,]
	inds=Fvalue$V1
	counts=as.data.frame(t(sapply(inds,function(x) get_rare_counts(x))),stringsAsFactors=F)
	counts$rare_count=as.numeric(counts$rare_count)
	Fvalue=merge(Fvalue,counts,by.x='V1',by.y='ID')
	names(Fvalue)=c('ID','Fvalue','F_rank','time','factor','rare_count','founder')
	Fvalue$chr=chr
	Fvalue$f_beta=unlist(betas[match(Fvalue$founder,names(betas))])
	betas=sort(unlist(betas),decreasing=TRUE)
	Fvalue$beta_rank=match(Fvalue$founder,names(betas))
	factordf=rbind(factordf,Fvalue)
}
factordf=as.data.frame(factordf,stringsAsFactors=F)
factordf$time_factor=paste0(factordf$time,'-',factordf$factor)
fwrite(factordf,'eqtl/results/residual_factor_rare_count.txt',row.names=F,quote=F,sep='\t')

factordf=factordf[factordf$founder!="B73_inra",]

founder_rank=factordf %>% group_by(founder,time_factor) %>% summarize(avg_rank=mean(beta_rank),avg_rare=mean(rare_count))
#
m2=lm(avg_rare~time_factor+poly(avg_rank,2),founder_rank)
m1=lm(rare_count~time_factor+poly(F_rank,2),factordf)

p1=ggplot(factordf,aes(x=as.factor(F_rank),y=rare_count))+ geom_point()
png('eqtl/images/factor_disregulation.png')
print(p1)
dev.off()


p1=ggplot(founder_rank,aes(x=avg_rank,y=avg_rare))+ geom_point()
png('eqtl/images/founder_factor_disregulation.png')
print(p1)
dev.off()
# What founder allele do each of these individuals have at these positions


