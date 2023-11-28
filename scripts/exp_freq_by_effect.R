#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('ggplot2')
library('lme4')
library('lmerTest')

# Z scores of expression
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10


# Just looking at the top 5000 highest expressed genes
totalrares=c()
for(time1 in times){
	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
	allrares=allrares[allrares$max_f!="B73_inra",]
	totalrares=rbind(totalrares,allrares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)


#totalrares=fread('eqtl/results/all_rare_counts_max_f_all_exp.txt',data.table=F)
#totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))

#ft_df=fread('eqtl/data/FT_genelist.txt',data.table=F)
#ftgenes=ft_df$Gene_ID
#ft_rares=totalrares[totalrares$Gene_ID %in% ftgenes,]
#ft_rares=ft_rares[!is.na(ft_rares$beta_rank),]
#print(length(unique(ft_rares$Gene_ID)))
# 267 genes
### Beta Z-scores
totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]
totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))

check=totf %>% group_by(gene_time_founder) %>% reframe(nrares=length(unique(rare_count)),nobs=length(gene_time_founder))

beta_z=totf%>% group_by(gene_time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T),abs_beta_z=abs((beta-mean(beta,na.rm=T))/sd(beta,na.rm=T)))

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



fwrite(beta_merge,'eqtl/results/local_eqtl_effect_size_frequency_beta_z_top5k.txt',row.names=F,quote=F,sep='\t')

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

fwrite(beta_merge,'eqtl/results/local_eqtl_effect_size_frequency_beta_z_all_genes.txt',row.names=F,quote=F,sep='\t')

# Which alleles are associated with later flowering? (Higher exp for positive r, lower exp for negative r)
summary(abs(beta_merge$cor))


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
fwrite(over_under_dir,'eqtl/results/directional_FT_freq_counts_all_genes.txt',row.names=F,quote=F,sep='\t')

over_under_div=as.data.frame(over_under_div,stringsAsFactors=F)
over_under_div$perc=over_under_div$over/over_under_div$total

fwrite(over_under_div,'eqtl/results/diversifying_FT_freq_counts_all_genes.txt',row.names=F,quote=F,sep='\t')

#Directional - how often is the extreme founder allele associated with late FT over represented?


# Diversifying
##### Which alleles for these genes are the extreme one? Are they at higher frequency

over_under_div=fread('eqtl/results/diversifying_FT_freq_counts_all_genes.txt',data.table=F)

# For each gene_time, rank by zscore and make same "dysregulation plot", but with avg frequency on the y-axis

totalrares=c()
for(time1 in times){
	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
	allrares=allrares[allrares$max_f!="B73_inra",]
	totalrares=rbind(totalrares,allrares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
totalrares$gene_time_id=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$ID)

# first I need to break up by gene_time_founder
totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

beta_merge=fread('eqtl/results/local_eqtl_effect_size_frequency_beta_z_all_genes.txt',data.table=F)
beta_merge=beta_merge[!is.na(beta_merge$beta_z),]
# Extreme bins are >=3/<=-3 SD
bbreaks=c(min(beta_merge$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_merge$beta_z)+0.1)
beta_merge$bin=cut(beta_merge$beta_z, breaks=bbreaks, label=F)
#beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
#beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

avg_z=beta_merge %>% group_by(bin) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)


p3=ggplot(aes(x=bin,y=avg_freq),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/all_gene_allele_frequency_bin_smooth_smile.png')
print(p3)
dev.off()


#### What about just FT genes & genes with high correlation with FT
high=beta_merge[abs(beta_merge$cor)>0.2,]
ft_corgenes=unique(high$Gene_ID)
ft_df=fread('eqtl/data/FT_genelist.txt',data.table=F)

# 101 genes that are in the gene list and correlated with FT

bbreaks=c(min(high$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_merge$beta_z)+0.1)
#beta_merge$bin=cut(beta_merge$beta_z, breaks=bbreaks, label=F)
#beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
#beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

avg_z=high %>% group_by(bin) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)


p3=ggplot(aes(x=bin,y=avg_freq),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/FT_gene_allele_frequency_bin_smooth_smile.png')
print(p3)
dev.off()

### Highly confident of involvement in FT

ft_corgenes=intersect(ft_df$Gene_ID,ft_corgenes)
# 101 genes
vhigh=high[high$Gene_ID %in% ft_corgenes,]

avg_z=vhigh %>% group_by(bin) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)


p3=ggplot(aes(x=bin,y=avg_freq),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/conf_FT_gene_allele_frequency_bin_smooth_smile.png')
print(p3)
dev.off()

### Frequency plots of candidate genes

cand_beta=beta_merge[beta_merge$Gene_ID %in% cand$Trait & beta_merge$time.x=="WD_0712",]

alldf=fread('eqtl/data/founder_frequency.txt',data.table=F)
alldf$snp_founder=paste0(alldf$snp,'-',alldf$variable)

cand=fread('QTT/sig_candidate_genes.txt',data.table=F)

theme_set(theme_classic())
#theme_update(text=element_text(family="Times"))
theme_update(plot.caption = element_text(hjust = 0))
theme_update(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
theme_update(plot.title = element_text(size=20),axis.title=element_text(size=20))
theme_update(panel.background=element_blank())
theme_update(plot.caption=element_text(size=20))

plot_list=list()
count=1
cands=cand$Trait
for(cg in cands){
	sub1=cand_beta[cand_beta$Gene_ID==cg,]
	bv_r=unique(sub1$cor)
	p1=ggplot(data=sub1,aes(x=beta,y=freq)) + geom_point(aes(color=max_f)) +
	geom_smooth(method='lm', formula= y~x) +
	#geom_hline(yintercept=0.0625,color='black') +
	xlab("Founder Effect Size")+ ylab("Founder Allele Frequency") +
	ggtitle(sprintf("%s, bv r=%.2f",cg,bv_r))
	
	plot_list[[count]]=p1
	count=count+1
}

pdf('QTT/images/candidate_gene_z_by_freq.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()
# For this list of candidate genes

# Does a linear model fit better than a quadratic?
# For linear model, is slope in the same direction as correlation?
#ft_corgenes=intersect(ft_df$Gene_ID,ft_corgenes)
# 101 genes
#vhigh=high[high$Gene_ID %in% ft_corgenes,]
ft_df=fread('eqtl/data/FT_genelist.txt',data.table=F)
high=beta_merge[abs(beta_merge$cor)>0.2,]
ft_corgenes=unique(high$Gene_ID)

# Is this more than we'd expect - for a random set of genes,
# Is there a different ratio of quadratic to linear? (More quadratic )

#ft_corgenes=intersect(ft_df$Gene_ID,ft_corgenes)
# 101 genes
#vhigh=high[high$Gene_ID %in% ft_corgenes,]



#vhigh=vhigh[vhigh$time.x=="WD_0712",]
vhigh=beta_merge[beta_merge$Gene_ID %in% ft_corgenes,]
# 1730 genes
# 6592 gene_times

comp_df=c()
#vhigh=vhigh[vhigh$time.x=="WD_0712",]
#vhigh=beta_merge[beta_merge$Gene_ID %in% ft_df$Gene_ID,]
# 593 genes
gene_times=unique(vhigh$gene_time)
# 2,131
for(gt in gene_times){
	sub1=vhigh[vhigh$gene_time==gt,]
	gene=unique(sub1$Gene_ID)
	time1=unique(sub1$time.x)
	
	m0=lm(freq~1,sub1)
	m1=lm(freq~beta_z,sub1)
	# Linear model
	m1sum=summary(m1)
	m1_pval= m1sum$coefficients['beta_z',4]
	m1_slope= m1sum$coefficients['beta_z',1]
	# Is slope in same direction as correlation?
	bv_r=unique(sub1$cor)
	if(bv_r>0 & m1_slope>0){
		consistent=TRUE
	}else{
		consistent=FALSE
	}
	# Quadratic model
	m2=lm(freq~poly(beta_z,2),sub1)
	m2sum=summary(m2)
	m2_pval= m2sum$coefficients[3,4]
	m2_curve= m2sum$coefficients[3,1]
	
	#Model comparison
	modelcomp=anova(m0,m1,m2)
	modelcomp_res1=modelcomp$`Pr(>F)`[2]
	modelcomp_res2=modelcomp$`Pr(>F)`[3]
	if(modelcomp_res2<=0.05){
		result="quadratic"
	}else if(modelcomp_res1<=0.05){
		result="linear"
	}else{
		result="neither"
	}
	
	
	line=data.frame(gene_time=gt,Gene_ID=gene,time=time1,m1_pvalue=m1_pval,m1_slope=m1_slope,bv_cor=bv_r,consistent=consistent,m2_pvalue=m2_pval,m2_curve=m2_curve,modelcomp_m1_pvalue=modelcomp_res1,modelcomp_m2_pvalue=modelcomp_res2,better=result)
	comp_df=rbind(comp_df,line)
}

table(comp_df$better)



truedf=data.frame(better=c('linear','p_quadratic','n_quadratic','neither'),n=c(349,149,279,5815),stringsAsFactors=F)

#   linear   neither quadratic 
#      104      1887       140 
comp_df=as.data.frame(comp_df,stringsAsFactors=F)

fwrite(comp_df,'eqtl/results/div_dir_model_comparison.txt',row.names=F,quote=F,sep='\t')
# Linear is significant

linear=comp_df[comp_df$better=="linear",]

# Quadratic is significant

quad=comp_df[comp_df$better=="quadratic",]

dim(quad[quad$m2_curve>0,])
# 60 out of 140 are positively curved (diversifying)

bonf=0.05/nrow(comp_df)
dim(comp_df[comp_df$m1_pvalue<=bonf,])
dim(comp_df[comp_df$m2_pvalue<=bonf,])

plot_list=list()
count=1
cands=quad$gene_time
for(cg in cands){
	sub1=vhigh[vhigh$gene_time==cg,]
	bv_r=unique(sub1$cor)
	p1=ggplot(data=sub1,aes(x=beta_z,y=freq)) + geom_point(aes(color=max_f)) +
	geom_smooth(method='lm', formula = y ~ x + I(x^2)) +
	#geom_hline(yintercept=0.0625,color='black') +
	xlab("Founder Effect Size")+ ylab("Founder Allele Frequency") +
	ggtitle(sprintf("%s, bv r=%.2f",cg,bv_r))
	
	plot_list[[count]]=p1
	count=count+1
}


#png('QTT/test.png')
#print(p1)
#dev.off()

pdf('QTT/images/quadratic_candidate_gene_beta_by_freq.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()

# Is this more than we'd expect - for a random set of genes,
# Is there a different ratio of quadratic to linear? (More quadratic )

ft_corgenes=intersect(ft_df$Gene_ID,ft_corgenes)
# 101 genes
vhigh=high[high$Gene_ID %in% ft_corgenes,]


comp_df=c()
#vhigh=vhigh[vhigh$time.x=="WD_0712",]
vhigh=beta_merge[beta_merge$Gene_ID %in% ft_df$Gene_ID,]
# 593 genes
gene_times=unique(vhigh$gene_time)

rand=beta_merge[!(beta_merge$gene_time %in% gene_times),]

n=length(gene_times)
set.seed(100)
n_reps=1:100

rep_comp_df=c()
for(rep in n_reps){
	draw=sample(rand$gene_time,n)
	for(gt in draw){
		sub1=rand[rand$gene_time==gt,]
		gene=unique(sub1$Gene_ID)
		time1=unique(sub1$time.x)
		m0=lm(freq~1,sub1)
		m1=lm(freq~beta_z,sub1)
		# Linear model
		m1sum=summary(m1)
		m1_pval= m1sum$coefficients['beta_z',4]
		m1_slope= m1sum$coefficients['beta_z',1]
		# Is slope in same direction as correlation?
		bv_r=unique(sub1$cor)
		if(bv_r>0 & m1_slope>0){
			consistent=TRUE
		}else{
			consistent=FALSE
		}
		# Quadratic model
		m2=lm(freq~poly(beta_z,2),sub1)
		m2sum=summary(m2)
		m2_pval= m2sum$coefficients[3,4]
		m2_curve= m2sum$coefficients[3,1]
		#Model comparison
		modelcomp=anova(m0,m1,m2)
		modelcomp_res1=modelcomp$`Pr(>F)`[2]
		modelcomp_res2=modelcomp$`Pr(>F)`[3]
		if(modelcomp_res2<=0.05){
			result="quadratic"
		}else if(modelcomp_res1<=0.05){
			result="linear"
		}else{
			result="neither"
		}
		line=data.frame(rep=rep,gene_time=gt,Gene_ID=gene,time=time1,m1_pvalue=m1_pval,m1_slope=m1_slope,bv_cor=bv_r,consistent=consistent,m2_pvalue=m2_pval,m2_curve=m2_curve,modelcomp_m1_pvalue=modelcomp_res1,modelcomp_m2_pvalue=modelcomp_res2,better=result)
		rep_comp_df=rbind(rep_comp_df,line)
	}
}
rep_comp_df=as.data.frame(rep_comp_df,stringsAsFactors=F)
# 2,131


fwrite(rep_comp_df,'eqtl/results/div_dir_model_comparison_perms.txt',row.names=F,quote=F,sep='\t')

truedf=data.frame(better=c('linear','p_quadratic','n_quadratic','neither'),true_n=c(349,149,279,5815),stringsAsFactors=F)


rep_comp_df=c()
for(rep in 1:100){
	p=fread(sprintf('eqtl/results/div_dir_model_comparison_perm%s.txt',rep),data.table=F)
	rep_comp_df=rbind(rep_comp_df,p)
}
rep_comp_df=as.data.frame(rep_comp_df,stringsAsFactors=F)

rep_sum=rep_comp_df %>% group_by(rep) %>% count(better,sort=TRUE)

# 
chi_res=c()

for(rep in 1:100){
	t0=rep_sum[rep_sum$rep==rep,]
	t1=truedf
	t1$rep_n=t0[match(t1$better,t0$better),]$n
	#t1=t1[t1$better!='neither',]
	test=chisq.test(t1$true_n,t1$rep_n,simulate.p.value=TRUE)
	line=data.frame(rep=rep,pvalue=test$p.value,stat=test$parameter)
	chi_res=rbind(chi_res,line)
}
chi_res=as.data.frame(chi_res,stringsAsFactors=F)

# None of them are significant
# Get average # of each time from reps
rep_sum %>% group_by(better) %>% reframe(avg_n=mean(n))
chi_df=data.frame(better=c('linear','quadratic_p','quadratic_n','neither'),
true_count=c(104,60,80,1887),
perm_count=c())

