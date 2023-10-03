#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')
library('emmeans')
library('lme4')
library('lmerTest')
library('PLS205')
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

alldf=fread('eqtl/data/founder_frequency.txt',data.table=F)
localcomp=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
# Separate out FT genes and genes overlapping FT QTL

ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
ftgenes=ft_genelist$V4
fts=c("male_flowering_d6","female_flowering_d6")
ftgenes=c(ftgenes,unique(localcomp[localcomp$phenotype %in% fts,]$Trait))
ftgenes=unique(ftgenes)
# 2437 genes
nftgenes=unique(eqtl[!(eqtl$Trait %in% ftgenes),]$Trait)
# 16675 genes

# For local-eQTL related to FT, do more extreme founder effect sizes tend to be at higher frequency

fteqtl=eqtl[eqtl$Trait %in% ftgenes,]
fteqtl=merge(fteqtl,alldf,by.x='X_ID',by.y='snp')
fsnps=unique(fteqtl$X_ID)
# 705 markers
# For local-eQTL, do more extreme founder effect sizes tend to be at lower frequency
nft=eqtl[!(eqtl$Trait %in% ftgenes),]
# eQTL are not linked with FT genes
nft=nft[!(nft$X_ID %in% fsnps),] #12945

nft=merge(nft,alldf,by.x='X_ID',by.y='snp')

nsnps = nft %>% group_by(X_ID) %>% summarize(n=length(unique(Trait)))

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
alltmp=c()
for(time1 in times){
	tmpdf=fread(sprintf('eqtl/results/eQTL_%s_freq_chi_data.txt',time1),data.table=F)
	alltmp=rbind(alltmp,tmpdf)
}
#time1="WD_0712"

tmpdf=alltmp
tmpdf=tmpdf[,1:42]

# Founder with highest deviation from avg effect sizes
max_beta=function(row1){
	avg=row1$avg_logexp
	betas=unlist(row1[,founders])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	betas=betas[!is.na(betas)]
	diff=betas-mean(betas)
	w=which.max(abs(diff))
	return(diff[w])
}

mbetas=sapply(seq(1,nrow(tmpdf)),function(x) max_beta(tmpdf[x,]))
tmpdf$max_betas=unlist(mbetas)
tmpdf$max_founder=names(mbetas)

tmpdf$gene_time_snp=paste0(tmpdf$Trait,'-',tmpdf$time,'-',tmpdf$X_ID)
eqtl$gene_time_snp=paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
tmpdf$eqtl_sig=tmpdf$gene_time_snp %in% eqtl$gene_time_snp

total_tests=73248
bonf=-log10(0.05 / total_tests)
tmpdf$chi_sig=-log10(tmpdf$p_chi)>=bonf


fwrite(tmpdf,'eqtl/results/eQTL_all_freq_chi_data.txt',row.names=F,quote=F,sep='\t')


max_beta2=funtion(x){

}


# How many eQTL are in seg distorted regions?
# 19403 eQTL in X2 regions, 1451 markers

# For these, are more the extreme founders at lower frequency?
# For each snp, what is the
# for each gene, rank 
chidf=tmpdf[tmpdf$eqtl_sig==TRUE & tmpdf$chi_sig==TRUE,]

# in rank, higher values associated with larger B diff
allranks=c()
csnps=unique(chidf$X_ID)
for(snp in csnps){
	sub=chidf[chidf$X_ID==snp,]
	mgenes=unique(sub$gene_time_snp)
	tmp=alldf[alldf$snp==snp,]
	rownames(tmp)=tmp$variable
	bcor=c()
	for(m in mgenes){
		sub2=sub[sub$gene_time_snp==m,]
		avg=sub2$avg_logexp
		betas=unlist(sub2[,founders])
		wn=which(!is.na(betas))[1]
		betas[-wn]=betas[-wn]+betas[wn]
		betas=betas[!is.na(betas)]
		
		#tmp=tmp[names(betas),]
		#tmp$beta=betas
		diff=abs(betas-mean(betas))
		rowm=paste0(m,'_rank')
		#tmp[,m]=diff[match(rownames(tmp),names(betas))]
		tmp[,rowm]=match(rownames(tmp),names(sort(abs(diff))))
		#tmp$diff=diff
		#diff=avg-betas
		#w=which.max(abs(diff))
	}
	#apply(tmp[,5:ncol(tmp),drop=F],MARGIN=1,function(x) mean(x,na.rm=T))
	tmp$avg_rank=unlist(apply(tmp[,5:ncol(tmp),drop=F],MARGIN=1,function(x) mean(x,na.rm=T)))
	#tmp$sd_rank
	tmp=tmp[,c('snp','variable','value','chr','avg_rank')]
	tmp$freq_rank=match(tmp$variable,tmp[order(tmp$value,decreasing=TRUE),]$variable)
	allranks=rbind(allranks,tmp)
}

allranks=as.data.frame(allranks,stringsAsFactors=F)
fwrite(allranks,'eqtl/results/all_chisq_rank_exp.txt',row.names=F,quote=F,sep='\t')

allranks$freq_f=factor(allranks$freq_rank,levels=c(seq(1,16)))

p2=ggplot(allranks,aes(x=freq_f,y=avg_rank)) + geom_boxplot()

png('eqtl/images/freq_by_exp_rank.png')
print(p2)
dev.off()

snp="AX-91558396"

# for each site, how often is the most extreme founder allele the one at the lowest frequency?


# 47010 eQTL in non-X2 regions

# what is the correlation between eQTL sig and max B diff
cor(abs(tmpdf$max_betas),-log10(tmpdf$p_value_ML))
#0.5133101
p1=ggplot(tmpdf,aes(x=abs(tmpdf$max_betas),y=-log10(tmpdf$p_value_ML))) + geom_point()

png('eqtl/images/pval_by_bdiff.png')
print(p1)
dev.off()
# X2 pvalue
cor(abs(tmpdf$max_betas),-log10(tmpdf$p_chi))
# not correlated at all -0.03819631

p1=ggplot(tmpdf,aes(x=abs(max_betas),y=-log10(p_chi))) + geom_point()

png('eqtl/images/x2_pval_by_bdiff.png')
print(p1)
dev.off()

# Grab marker - what is the correlation between founder freq and effect on expression
snp="AX-91558396"
sub=tmpdf[tmpdf$X_ID==snp,]


############################


alldf=c()
adj_chr=c(5,9)
for(chr in 1:10){
	#pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
	if(chr %in% adj_chr){
		X_list=readRDS(sprintf('phenotypes/bg%.0f_adjusted_genoprobs.rds',chr))
	}else{
		X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',chr))
	}
	snps=dimnames(X_list[[1]])[[2]]
	freq=as.data.frame(t(unlist(sapply(snps,function(x) colSums(do.call(cbind,lapply(X_list,function(y) y[,x])))/nrow(X_list[[1]])))))
	freq$snp=rownames(freq)
	freqm=reshape2::melt(freq,'snp')
	freqm$chr=chr
	#frep=sapply(snps,function(x) apply(do.call(cbind,lapply(X_list,function(y) y[,x])),MARGIN=2,function(z) round(sum(z>0.75))))
	#fvar=apply(frep,MARGIN=2,function(x) var(x[x>3]))
	#n_est=apply(frep,MARGIN=2,function(x) length(x[x>3]))
	#df=data.frame(chr=chr,snp=snps,fvar=fvar,n_est=n_est,stringsAsFactors=F)
	alldf=rbind(alldf,freqm)
}

fwrite(alldf,'eqtl/data/founder_frequency.txt',row.names=F,quote=F,sep='\t')

minhaps=fread('../genotypes/probabilities/haplotype_probs/RefinedIBD_600K/min_haps.txt',data.table=F)
hdf=c()
baselist=minhaps$V2
for(chr in 1:10){
	for(h in baselist[chr]:16){
		hapfile=sprintf('../genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%s_probs.rds',chr,h)
		if(file.exists(hapfile)){
			X_list=readRDS(hapfile)
			snps=dimnames(X_list[[1]])[[2]]
			freq=as.data.frame(t(unlist(sapply(snps,function(x) colSums(do.call(cbind,lapply(X_list,function(y) y[,x])))/nrow(X_list[[1]])))))
			colnames(freq)=paste0('haplotype',seq(1:h))
			freq$snp=rownames(freq)
			freqm=reshape2::melt(freq,'snp')
			freqm$nhaps=h
			freqm$chr=chr
		}
	}
	#pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
	
	#frep=sapply(snps,function(x) apply(do.call(cbind,lapply(X_list,function(y) y[,x])),MARGIN=2,function(z) round(sum(z>0.75))))
	#fvar=apply(frep,MARGIN=2,function(x) var(x[x>3]))
	#n_est=apply(frep,MARGIN=2,function(x) length(x[x>3]))
	#df=data.frame(chr=chr,snp=snps,fvar=fvar,n_est=n_est,stringsAsFactors=F)
	hdf=rbind(hdf,freqm)
}

fwrite(hdf,'eqtl/data/haplotype_frequency.txt',row.names=F,quote=F,sep='\t')

alldf=fread('eqtl/data/haplotype_frequency.txt',data.table=F)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

qtl=fread('QTL/all_adjusted_QTL_SIs.txt',data.table=F)
qtl$pheno_env_ID=paste0(qtl$phenotype,'-',qtl$environment,'-',qtl$ID)
#qtlmerge=merge(qtl,alldf,by.x="SNP",by.y="snp")


#alldf=alldf[alldf$snp %in% qtl$SNP,]
names(alldf)=c('snp',"haplotype","frequency","nhaps","chr")
# Get founder effect sizes for each qtl SNP

plot_list=list()
count=1
#qsnps=qtl$SNP
newdf=c()
ft=c("male_flowering_d6","female_flowering_d6")
for(i in 1:nrow(qtl)){
	row=qtl[i,]
	pheno=row$phenotype
	env=row$environment
	chr=row$CHR
	snp=row$SNP
	id=row$ID
	pei=row$pheno_env_ID
	ibd=fread(sprintf('../ibd_segments/refinedibd/600K/bg%s_refined_ibd_blocks.txt',chr),data.table=F)
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
	spos=pmap[pmap$marker==snp,]$pos
	subhdf=alldf[alldf$chr==chr,]
	subhdf$pos=pmap[match(subhdf$snp,pmap$marker),]$pos
	subhdf$dist=subdf$pos-spos
	subhdf=subhdf[subhdf$dist>0,]
	hsnp=subhdf[which.min(subhdf$dist),]$snp
	hpos=subhdf[which.min(subhdf$dist),]$pos
	fgroups=ibd[ibd$start<=hpos & ibd$end>hpos,founders]
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,founders])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	subdf=alldf[alldf$snp==hsnp,]
	subdf$effect_size=effect_size[founders]
	subdf$hapgrp=unlist(fgroups)
	
	newdf=data.frame(founder=founders,effect_size=effect_size,hapgrp=unlist(fgroups),stringsAsFactors=F)
	
	if(pheno %in% ft | id=="qHGM3_2"){
		subdf$ft=TRUE
	}else{
		subdf$ft=FALSE
	}
	subdf$pei=pei
	subdf$phenotype=pheno
	subdf$environment=env
	subdf$ID=id
	 p1=ggplot(subdf,aes(x=frequency,y=effect_size)) + geom_point(aes(color=founder)) +
	geom_smooth(method='lm', formula= y~x) + 
	xlab("Founder Allele Frequency") + ylab("QTL Effect Size") +
	ggtitle(sprintf("%s %s",pei,snp))
	plot_list[[count]]=p1
	count=count+1
	newdf=rbind(newdf,subdf)
}

pdf('QTL/images/effect_size_by_haplotype_freq.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()



## For QTL

alldf=fread('eqtl/data/founder_frequency.txt',data.table=F)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

qtl=fread('QTL/all_adjusted_QTL_peaks_trimmed.txt',data.table=F)
qtl$pheno_env_ID=paste0(qtl$phenotype,'-',qtl$environment,'-',qtl$ID)
#qtlmerge=merge(qtl,alldf,by.x="SNP",by.y="snp")


alldf=alldf[alldf$snp %in% qtl$SNP,]
names(alldf)=c('snp',"founder","frequency","chr")
# Get founder effect sizes for each qtl SNP

plot_list=list()
count=1
#qsnps=qtl$SNP
newdf=c()
ft=c("male_flowering_d6","female_flowering_d6")
for(i in 1:nrow(qtl)){
	row=qtl[i,]
	pheno=row$phenotype
	env=row$environment
	chr=row$CHR
	snp=row$SNP
	id=row$ID
	pei=row$pheno_env_ID
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,founders])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	subdf=alldf[alldf$snp==snp,]
	subdf$effect_size=effect_size[subdf$founder]
	if(pheno %in% ft | id=="qHGM3_2"){
		subdf$ft=TRUE
	}else{
		subdf$ft=FALSE
	}
	subdf$pei=pei
	subdf$phenotype=pheno
	subdf$environment=env
	subdf$ID=id
	 p1=ggplot(subdf,aes(x=frequency,y=effect_size)) + geom_point(aes(color=founder)) +
	geom_smooth(method='lm', formula= y~x) + 
	xlab("Founder Allele Frequency") + ylab("QTL Effect Size") +
	ggtitle(sprintf("%s %s",pei,snp))
	plot_list[[count]]=p1
	count=count+1
	newdf=rbind(newdf,subdf)
}

pdf('QTL/images/effect_size_by_freq.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()



### For FT, grab largest effect size
ft=newdf[newdf$ft==TRUE,]
topft=ft %>% group_by(pei) %>% slice_max(effect_size, n = 1)
topft=as.data.frame(topft,stringsAsFactors=F)

nft=newdf[newdf$ft==FALSE,]
topnft=nft %>% group_by(pei) %>% slice_max(abs(effect_size), n = 1)
topnft=as.data.frame(topnft,stringsAsFactors=F)
topnft$effect_size=abs(topnft$effect_size)

testdf=rbind(topft,topnft)
testdf=as.data.frame(testdf,stringsAsFactors=F)

testdf2=testdf %>% group_by(ID_snp) %>% mutate(effect_size=mean(effect_size)) %>% distinct(frequency,.keep_all=T)
newdf$ID_snp=paste0(newdf$ID,'-',newdf$snp)



mod=lmer(frequency ~ (1|ID) + factor(ft),data=testdf2)
anova(mod,ddf="K")

em=emmeans(mod,specs="ft")
contrasts=contrast(em,'trt.vs.ctrl',ref=c('FALSE'))

# contrast     estimate      SE   df t.ratio p.value
# TRUE - FALSE   0.0129 0.00987 20.2   1.307  0.2059
 
 
# Evidence of stabilizing selection 
# one-sided t-test
# Is frequency of non-FT alleles of large effect size less than the frequeny of FT alleles of large effect sizes?
nft_freq=testdf2[testdf2$ft==FALSE,]$frequency
ft_freq=testdf2[testdf2$ft==TRUE,]$frequency

test=t.test(nft_freq,ft_freq,alternative="less")

#	Welch Two Sample t-test
#
#data:  nft_freq and ft_freq
#t = -1.7877, df = 38.07, p-value = 0.04089
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#          -Inf -0.0007223553
#sample estimates:
# mean of x  mean of y 
#0.04913212 0.06181125
testdf2=as.data.frame(testdf2,stringsAsFactors=F)
p1=ggplot(data=testdf2,aes(y=frequency,x=ft)) + geom_boxplot() + geom_jitter() + xlab("QTL Phenotype") + ylab("Frequency of Max Beta Founder Allele ") +
scale_x_discrete(labels=c("Non-FT","FT"))

png('QTL/images/freq_by_FT.png')
print(p1)
dev.off()

newdf$ID_founder=interaction(newdf$ID,newdf$founder)
newdf$env_founder=interaction(newdf$envirnment,newdf$founder)
newdf$environment=as.factor(newdf$envirnment)
# create multiple linear model
# read dataset
ids=unique(newdf$ID)

all_trends=c()
testdf=data.frame(ids=ids,stringsAsFactors=F)
pvalues=c()
for(id in ids){
	tmpdf=newdf[newdf$ID==id,]
	tmpdf$environment=as.factor(tmpdf$envirnment)
	tmpdf$env_founder=interaction(tmpdf$envirnment,tmpdf$founder)
	if(length(unique(tmpdf$environment))==1){
		lm_fit <- lm(effect_size ~ I(frequency), data=tmpdf)
		results=anova(lm_fit)
	}else{
		 lm_fit=lmer(effect_size~I(frequency)+(1|founder),data=tmpdf)
		 results=anova(lm_fit,ddf='K',type='I')
		#envmean = tmpdf %>% group_by(founder) %>% mutate(effect_size=mean(effect_size)) %>% distinct(founder,.keep_all=T)
		#lm_fit <- lm(effect_size ~ I(frequency) , data=envmean)
	}

	#lm_fit <- lmer(effect_size ~ (1|environment) + I(frequency) , data=tmpdf)
	#results=anova(lm_fit,ddf="K",type="I")
	
	
	pvalue=results['I(frequency)',]$`Pr(>F)`
	predicted_data <- data.frame(frequency = seq(0,0.2,by=0.025)) 
	trend2_means = emmeans(lm_fit,spec = 'frequency',at=list(frequency = predicted_data$frequency))
	trend2_CIs = summary(trend2_means,infer = c(T,F),level = 0.95,ddf='K')
	trend2_CIs_table = as.data.frame(trend2_CIs)
	trend2_CIs_table$ID=id
	trend2_CIs_table$pvalue=pvalue
	all_trends=rbind(all_trends,trend2_CIs_table)
	pvalues=c(pvalues,pvalue)
}
testdf$pvalue=pvalues
testdf$adjusted_p=p.adjust(testdf$pvalue,method='fdr')

testdf[testdf$adjusted_p<=0.05,]

ftdf=newdf[newdf$phenotype %in% c("male_flowering_d6","female_flowering_d6"),]
ftdf$environment_ID=interaction(ftdf$environment,ftdf$ID)
ftdf$ID_founder=interaction(ftdf$ID,ftdf$founder)
ftdf$env_founder=interaction(ftdf$environment,ftdf$founder)
ftdf$environment=as.factor(ftdf$environment)

ftmod=lmer(effect_size~ ID + (1|ID_founder) + I(frequency),data=ftdf)

ftmod2=lmer(effect_size~ (1|ID) + (1|environment_ID) + (1|ID_founder) + I(frequency),data=ftdf)
# What about absolute effect size?
all_trends=c()
testdf=data.frame(ids=ids,stringsAsFactors=F)
pvalues=c()
for(id in ids){
	tmpdf=newdf[newdf$ID==id,]
	tmpdf$environment=as.factor(tmpdf$envirnment)
	tmpdf$env_founder=interaction(tmpdf$envirnment,tmpdf$founder)
	if(length(unique(tmpdf$environment))){
		lm_fit <- lm(abs(effect_size) ~ I(frequency) , data=tmpdf)
	}else{
		lm_fit <- lm(abs(effect_size) ~ environment + I(frequency) , data=tmpdf)
	}

	#lm_fit <- lmer(effect_size ~ (1|environment) + I(frequency) , data=tmpdf)
	#results=anova(lm_fit,ddf="K",type="I")
	
	results=anova(lm_fit)
	pvalue=results['I(frequency)',]$`Pr(>F)`
	predicted_data <- data.frame(frequency = seq(0,0.2,by=0.025)) 
	trend2_means = emmeans(lm_fit,spec = 'frequency',at=list(frequency = predicted_data$frequency))
	trend2_CIs = summary(trend2_means,infer = c(T,F),level = 0.95,ddf='K')
	trend2_CIs_table = as.data.frame(trend2_CIs)
	trend2_CIs_table$ID=id
	trend2_CIs_table$pvalue=pvalue
	all_trends=rbind(all_trends,trend2_CIs_table)
	pvalues=c(pvalues,pvalue)
}
testdf$pvalue=pvalues
testdf$adjusted_p=p.adjust(testdf$pvalue,method='fdr')

testdf[testdf$adjusted_p<=0.05,]
#      ids       pvalue   adjusted_p  adjusted_p2
#1 qDTS3_2 3.260865e-06 3.586951e-05 1.482211e-07
#5 qDTA3_2 2.370248e-09 5.214546e-08 1.077386e-10

# save predictions of the model in the new data frame 
# together with variable you want to plot against
predicted_df <- data.frame(mpg_pred = predict(lm_fit, df), hp=df$hp)

# this is the predicted line of multiple linear regression
ggplot(data = df, aes(x = mpg, y = hp)) + 
  geom_point(color='blue') +
  geom_line(color='red',data = predicted_df, aes(x=mpg_pred, y=hp))