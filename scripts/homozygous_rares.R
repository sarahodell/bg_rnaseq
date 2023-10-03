#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('stringr')
library('ggplot2')
#library('glmnet')
#library('lme4')
library('lmerTest')
library('lme4qtl')

chr=10
K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)
inds=rownames(K)

chroms=1:10

allcounts=data.frame(inds=inds,stringsAsFactors=F)
for(chr in chroms){
	min_founder=readRDS(sprintf('datasets/hapmap/minor_allele_founders_c%s.rds',chr))
	bimbam=fread(sprintf('datasets/chr%s_biogemma_rare_allele_probs.txt',chr),data.table=F)
	testrare=which(unlist(lapply(min_founder,function(x) "MBS847" %in% x)))
	sites=names(min_founder)[testrare]
	bimbam=bimbam[bimbam$marker %in% sites,]
	#names(bimbam)=c('marker','alt1','ref',rownames(K),'pos')
	bimbam$chr=chr
	bimbam=bimbam[,c('marker','chr','pos','alt1','ref',rownames(K))]
	#For each individual, how many sites are they homozygous for rare alleles?
	cutoff=0.75
	rcount=apply(bimbam[,inds],MARGIN=2,function(x) sum(x>cutoff))
	fwrite(bimbam,sprintf('datasets/chr%s_tester_rare_allele_probs.txt',chr),row.names=F,quote=F,sep='\t')
	coln=paste0('chr',chr)	
	allcounts[,coln]=unlist(rcount)
}


allcounts=as.data.frame(allcounts,stringsAsFactors=F)
allcounts$total=rowSums(allcounts[,2:11])
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

phenotypes=phenotypes[phenotypes$Loc.Year.Treat=="ALL",]
allcounts$grain_yield_15=phenotypes[match(allcounts$inds,phenotypes$Genotype_code),]$grain_yield_15
allcounts$tkw_15=phenotypes[match(allcounts$inds,phenotypes$Genotype_code),]$tkw_15


fwrite(allcounts,'datasets/homozygous_rare_count.txt',row.names=F,quote=F,sep='\t') 


## How many rare alleles
bimbam=fread('datasets/all_tester_rare_allele_probs.txt',data.table=F)
upstream=fread('eqtl/data/upstream_gene_list.txt',data.table=F)
bimbam$pos1=bimbam$pos+1

env1=as.data.table(bimbam)
env2=as.data.table(upstream)
setkey(env2,CHROM,UP_START,UP_END)
comp=foverlaps(env1,env2,by.x=c('chr','pos','pos1'),by.y=c('CHROM','UP_START','UP_END'),nomatch=NULL)
sites=unique(comp$marker)
bimbam=bimbam[bimbam$marker %in% sites,]

cutoff=0.75
rcount=apply(bimbam[,inds],MARGIN=2,function(x) sum(x>cutoff))
allcounts$UPSTREAM_5kb_count=rcount
fwrite(allcounts,'datasets/homozygous_rare_count.txt',row.names=F,quote=F,sep='\t') 


#x=allcounts$UPSTREAM_5kb_count
#y=allcounts$grain_yield_15

allcounts=fread('datasets/homozygous_rare_count.txt',data.table=F)
rownames(allcounts)=allcounts$inds
model=relmatLmer(grain_yield_15 ~ UPSTREAM_5kb_count + (1|inds), allcounts, relmat = list(inds = K))

lme4qtl::VarProp(model)
lme4::VarCorr(model)
anova(model,ddf='k')

m0 <- relmatLmer(grain_yield_15 ~ (1|inds), allcounts, relmat = list(inds = K))
#m1 <- relmatLmer(trait1 ~ AGE + (1|ID), dat40, relmat = list(ID = kin2))
anova(m0, model)

#refitting model(s) with ML (instead of REML)

#Data: allcounts
#Models:
#m0: grain_yield_15 ~ (1 | inds)
#model: grain_yield_15 ~ UPSTREAM_5kb_count + (1 | inds)
#      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#m0       3 2027.2 2038.5 -1010.6   2021.2                         
#model    4 2017.3 2032.4 -1004.6   2009.3 11.907  1  0.0005594 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#model2=relmatLmer(grain_yield_15 ~ total + (1|inds), allcounts, relmat = list(inds = K))

#model3=lm(grain_yield_15 ~ UPSTREAM_5kb_count,data=allcounts)

#model4=lm(grain_yield_15 ~ total,data=allcounts)

#
#ridgefit=glmnet(x,y,alpha=0)

#plot(ridgefit)


p1=ggplot(aes(x=UPSTREAM_5kb_count,y=grain_yield_15),data=allcounts) + geom_point()+
xlab("# Homozygous Rares Alleles 5kb upstream") + ylab("Grain Yield")
png('eqtl/images/yield_by_homo_5kb_upstream.png')
print(p1)
dev.off()




#### I number of rare alleles possessed by founder associated with founder freqeuncy?
chr=10

rare_blocks=c()
for(chr in 1:10){
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
	sites=bimbam[bimbam$chr==chr,]$pos
	ounder=founder[founder$POS %in% sites,]
	# remove reference minor alleles 
	founder=founder[founder$B73_inra!=1,]
	# remove multi-allelic sites
	founder=founder[is.na(founder$ALT2),]
	# within windows, how many rare alleles do each founder have?
	founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
	founders2=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
	block_df=c()
	for(i in 1:nrow(founder_blocks)){
		row1=founder_blocks[i,]
		start=row1$start
		end=row1$end
		subf=founder[founder$POS>=start & founder$POS<end,]
		df=data.frame(founder=founders2,rfcount=colSums(subf[,founders2],na.rm=T),stringsAsFactors=F)
		df=cbind(df,row1)
		block_df=rbind(block_df,df)
	}
	rare_blocks=rbind(rare_blocks,block_df)
}

rare_blocks=as.data.frame(rare_blocks,stringsAsFactors=F)
rare_blocks=rare_blocks[rare_blocks$founder!="B73_inra",]

fwrite(rare_blocks,'datasets/founder_blocks_homo_rares.txt',row.names=F,quote=F,sep='\t')

rare_blocks=fread('datasets/founder_blocks_homo_rares.txt',data.table=F)
rare_blocks$founder_snp=paste0(rare_blocks$founder,'-',rare_blocks$focal_snp)
# frequency data
alldf=fread('eqtl/data/founder_frequency.txt',data.table=F)
alldf$founder_snp=paste0(alldf$variable,'-',alldf$snp)
alldf=alldf[alldf$variable!="B73_inra",]
alldf=alldf[alldf$snp %in% unique(rare_blocks$focal_snp),]

rare_blocks$frequency=alldf[match(rare_blocks$founder_snp,alldf$founder_snp),]$value
# chisq data
chidf=c()
for(chr in 1:10){
	chi=fread(sprintf('../selection/founder_probs/bg%s_founder_chisq_results.txt',chr),data.table=F)
	chi$chr=chr
	chidf=rbind(chidf,chi)
}

rare_blocks$p_chi=chidf[match(rare_blocks$focal_snp,chidf$marker),]$p_chi

ntests=75456
bonf=0.05/ntests
sig_blocks=rare_blocks[rare_blocks$p_chi<=bonf,]
#### What is the relationship between frequency and # of rare alleles and frequency of founder

#model=lmer(frequency ~ rfcount + (1|focal_snp) + (1|founder),sig_blocks)

# For each site, is the relationship between rfcount and frequency significantly less than 0?
snps=unique(sig_blocks$focal_snp)

pvals=c()
rs=c()
for(snp in snps){
	sblock=sig_blocks[sig_blocks$focal_snp==snp,]
	sblock=sblock[complete.cases(sblock),]
	test=cor.test(sblock$frequency,sblock$rfcount,alternative="less")
	p=test$p.value
	r=test$estimate
	pvals=c(pvals,p)
	rs=c(rs,r)
}

res=data.frame(snp=snps,r=rs,pvalue=pvals,stringsAsFactors=F)
res=res[complete.cases(res),]
snps=unique(res$snp)



sig_blocks=sig_blocks[sig_blocks$focal_snp %in% snps,]
# Do the founders under 1/16 frequency have more rare alleles than those over 1/16 frequency?

even=1/16

# does it make sense to mean standardize rfcounts first?

sig_blocks$under=sig_blocks$frequency<(even-0.01)
sig_blocks$over=sig_blocks$frequency>(even+0.01)

sig_blocks=sig_blocks[complete.cases(sig_blocks),]
under_rares=sig_blocks[sig_blocks$under==TRUE,]$rfcount
over_rares=sig_blocks[sig_blocks$over==TRUE,]$rfcount

t.test(under_rares,over_rares,alternative="greater")
#	Welch Two Sample t-test
#
#data:  under_rares and over_rares
#t = 3.0577, df = 17200, p-value = 0.001117
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
# 5.290255      Inf
#sample estimates:
#mean of x mean of y 
# 55.41118  43.96114 
#p1=ggplot(sig_blocks,aes(x=under,y=rfcount)) + geom_boxplot() + geom_jitter() + xlab('Under/Over Founder Frequency') + 
#ylab("# of Rare Alleles in Founder") + ggtitle("Founder Homozygous Rare Allele Count per LD block")
#png('datasets/frequency_hmzg_rare_count.png')
#print(p1)
#dev.off()

blocks2=sig_blocks %>% group_by(focal_snp) %>% summarize(rf_std=(rfcount-mean(rfcount))/sd(rfcount),frequency=frequency,over=over,under=under)
blocks2=as.data.frame(blocks2,stringsAsFactors=F)

blocks2$under=blocks2$frequency<even
blocks2=blocks2[complete.cases(blocks2),]
under_rares=blocks2[blocks2$under==TRUE,]$rf_std
over_rares=blocks2[blocks2$over==TRUE,]$rf_std

t.test(under_rares,over_rares,alternative="greater")
#		Welch Two Sample t-test
#
#data:  under_rares and over_rares
#t = 3.9135, df = 17779, p-value = 4.565e-05
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
# 0.03066414        Inf
#sample estimates:
#  mean of x   mean of y 
# 0.02021114 -0.03268740 

pblocks2=blocks2[blocks2$under==TRUE | blocks2$over==TRUE,]
p1=ggplot(pblocks2,aes(x=under,y=rf_std)) + geom_boxplot() + xlab('Under/Over Founder Frequency') + 
ylab("# of Rare Alleles in Founder (Standardized)") + ggtitle("Founder Homozygous Rare Allele Count per LD block")
png('datasets/frequency_hmzg_rare_count_std.png')
print(p1)
dev.off()


# There are signfiic
model=lmer(frequency~rf_std + (1|focal_snp),pblocks2)
anova(model,ddf='K')
#Type III Analysis of Variance Table with Kenward-Roger's method
#          Sum Sq   Mean Sq NumDF DenDF F value  Pr(>F)  
#rf_std 0.0044357 0.0044357     1 23183  4.7918 0.02861 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

model2=lmer(frequency~rfcount + (1|focal_snp),sig_blocks)
anova(model2,ddf='K')
#Type III Analysis of Variance Table with Kenward-Roger's method
#          Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#rfcount 0.003208 0.003208     1 5697.7  3.4567 0.06305 .
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

####

# What about sites that are 5kb upstream of a gene

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

allbimbams=c()
for(chr in chroms){
	bimbam=fread(sprintf('datasets/chr%s_tester_rare_allele_probs.txt',chr),data.table=F)
	#bimbam$chr=chr
	#names(bimbam)=c('marker','chr','pos','alt1','ref',rownames(K))
	allbimbams=rbind(allbimbams,bimbam)
}
allbimbams=as.data.frame(allbimbams,stringsAsFactors=F)
fwrite(allbimbams,'datasets/all_tester_rare_allele_probs.txt',row.names=F,quote=F,sep='\t')

# Which side is upstream?
genetable$UPSTREAM=genetable$START-5000
gtf=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46.gtf',data.table=F)
gtf=gtf[gtf$V3=="gene",]
# Grab gene ID
gtf$ID=unlist(sapply(seq(1,nrow(gtf)),function(x) strsplit(gtf$V9[x],'\"')[[1]][2]))

genetable$STRAND=gtf[match(genetable$Gene_ID,gtf$ID),]$V7

genetable$UPSTREAM=sapply(seq(1,nrow(genetable)),function(x) ifelse(genetable$STRAND[x]=='+',genetable$START[x]-5000,genetable$END[x]+5000))

names(genetable)=c("Gene_ID","CHROM","START","END","STRAND","UPSTREAM_5kb")
fwrite(genetable,'eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',row.names=F,quote=F,sep='\t')

upstream=genetable[,c("Gene_ID","CHROM")]
upstream$UP_START=sapply(seq(1,nrow(genetable)),function(x) ifelse(genetable$STRAND[x]=='+',genetable$UPSTREAM_5kb[x],genetable$END[x]))
upstream$UP_END=sapply(seq(1,nrow(genetable)),function(x) ifelse(genetable$STRAND[x]=='+',genetable$START[x],genetable$UPSTREAM_5kb[x]))
fwrite(upstream,'eqtl/data/upstream_gene_list.txt',row.names=F,quote=F,sep='\t')