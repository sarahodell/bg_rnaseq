#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('ggplot2')
library('lme4')
library('lmerTest')

# Z scores of expression
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10

totalrares=c()
for(time1 in times){
	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
	allrares=allrares[allrares$max_f!="B73_inra",]
	totalrares=rbind(totalrares,allrares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

ft_df=fread('eqtl/data/FT_genelist.txt',data.table=F)
ftgenes=ft_df$Gene_ID
ft_rares=totalrares[totalrares$Gene_ID %in% ftgenes,]
ft_rares=ft_rares[!is.na(ft_rares$beta_rank),]
print(length(unique(ft_rares$Gene_ID)))
# 267 genes
### Beta Z-scores
totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]
totf=totalrares %>% group_by(gene_time_founder) %>% summarize(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))

beta_z=totf%>% group_by(gene_time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))

beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
# Extreme bins are >=3/<=-3 SD
bbreaks=c(min(beta_z$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_z$beta_z)+0.1)
beta_z$bin=cut(beta_z$beta_z, breaks=bbreaks, label=F)
beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)

alldf=fread('eqtl/data/founder_frequency.txt',data.table=F)
alldf$snp_founder=paste0(alldf$snp,'-',alldf$variable)

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

beta_merge=merge(beta_z,eqtl,by='gene_time')
beta_merge$snp_founder=paste0(beta_merge$X_ID,'-',beta_merge$max_f)
# 769 gene time, 233 genes
beta_merge$freq=alldf[match(beta_merge$snp_founder,alldf$snp_founder),]$value



fwrite(beta_merge,'eqtl/results/local_eqtl_effect_size_frequency_beta_z.txt',row.names=F,quote=F,sep='\t')

### For all of them, what is the gene's correlation with flowering time?
fts=c("male_flowering_d6","female_flowering_d6")
cutoff=0.75
localcomp=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
localcomp=localcomp[localcomp$phenotype %in% fts,]
localcomp=localcomp[abs(localcomp$r)>=cutoff,]

tmp=merge(localcomp,beta_merge,by='gene_time')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
pheno="male_flowering_d6"
env="EXP_STPAUL_2017_WD"
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
rownames(phenotypes)=phenotypes$Genotype_code
gts=unique(beta_merge$gene_time)
beta_merge$pvalue=1
beta_merge$cor=0

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
# correlation of local additive genetic effect and phenotype (should I do BLUPs?) I think correlation will be higher in 
# environment where tissue was sampled
for(time1 in times){
	for(chr in 1:10){
		exp=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_bv_FIXED.txt',time1,chr),data.table=F)
		sub0=beta_merge[beta_merge$time.x==time1 & beta_merge$chr==chr,]
		gts=unique(sub0$Gene_ID)
		exp=exp[,c('V1',gts)]
		rownames(exp)=exp$V1
		exp=exp[,-1]
		inter=intersect(rownames(phenotypes),rownames(exp))
		subpheno=phenotypes[inter,]
		exp=exp[inter,]
		for(i in gts){
			test=cor.test(exp[,i],subpheno$male_flowering_d6,use="complete.obs")
			r=test$estimate
			p=test$p.value
			beta_merge[beta_merge$time.x==time1 & beta_merge$Gene_ID==i,]$pvalue=p
			beta_merge[beta_merge$time.x==time1 & beta_merge$Gene_ID==i,]$cor=r
		}
		#cors=c(cors,r)
		#pvals=c(pvals,p)
	}

	#exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
	
}

fwrite(beta_merge,'eqtl/results/local_eqtl_effect_size_frequency_beta_z.txt',row.names=F,quote=F,sep='\t')

# Which alleles are associated with later flowering? (Higher exp for positive r, lower exp for negative r)

high=beta_merge[abs(beta_merge$cor)>0.1,]
# "Zm00001d034045" "Zm00001d048474" "Zm00001d051135"


# If correlation is negative, look for lowest expressed (negative beta)
# if correlation is positive, look for highest expressed (positive beta)

over_under_dir=c()
over_under_div=c()

over=0.0625 + 0.01
under=0.0625 - 0.01
gts=unique(high$gene_time)
for(t in gts){
	sub=high[high$gene_time==t,]
	cor=unique(sub$cor)
	sub0=sub %>% group_by(max_f) %>% summarize(freq=mean(freq),beta_z=mean(beta_z),bin=unique(bin),rare_count=mean(rare_count))
	sub0=as.data.frame(sub0,stringsAsFactors=F)
	sub0$abs_betaz=abs(sub0$beta_z)

	# Number of extreme(z-score > 1 or <-1 ) that are over or under represented

	if(cor<0){ # Negative
		test1=sub0[sub0$beta_z<= -1,]
	}else{ #Positive
		test1=sub0[sub0$beta_z>= 1,]
	}
	test2=sub0[sub0$abs_betaz>=1,]
	ov_dir=nrow(test1[test1$freq>=over,])
	un_dir=nrow(test1[test1$freq<=under,])
	tot_dir=nrow(test1)
	
	ov_div=nrow(test2[test2$freq>=over,])
	un_div=nrow(test2[test2$freq<=under,])
	tot_div=nrow(test2)
	
	line_dir=data.frame(gene_time=t,over=ov_dir,under=un_dir,total=tot_dir)
	over_under_dir=rbind(over_under_dir,line_dir)
	
	line_div=data.frame(gene_time=t,over=ov_div,under=un_div,total=tot_div)
	over_under_div=rbind(over_under_div,line_div)
}

over_under_dir=as.data.frame(over_under_dir,stringsAsFactors=F)
over_under_dir$perc=over_under_dir$over/over_under_dir$total
fwrite(over_under_dir,'eqtl/results/directional_FT_freq_counts.txt',row.names=F,quote=F,sep='\t')

over_under_div=as.data.frame(over_under_div,stringsAsFactors=F)
over_under_div$perc=over_under_div$over/over_under_div$total

fwrite(over_under_div,'eqtl/results/diversifying_FT_freq_counts.txt',row.names=F,quote=F,sep='\t')

#Directional - how often is the extreme founder allele associated with late FT over represented?


# Diversifying
##### Which alleles for these genes are the extreme one? Are they at higher frequency

