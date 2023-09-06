#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')
library('emmeans')
library('lme4')
library('lmerTest')
library('PLS205')

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

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