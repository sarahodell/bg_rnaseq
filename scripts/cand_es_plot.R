#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')


rap27_1="Zm00001d010987"
rap27_2="Zm00001d010988"
zcn8="Zm00001d010752"
mads69="Zm00001d042315"
cand=fread('QTT/sig_candidate_genes.txt',data.table=F)
cand$pei=paste0(cand$phenotype,'-',cand$environment,'-',cand$ID)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
plot_list=list()
count=1


# Make data table of gene and effect size

time="WD_0712"

plotdf=c()
for(i in 1:nrow(cand)){
	row=cand[i,]
	pheno=row$phenotype
	env=row$environment
	chr=row$CHR
	gene=row$Trait
	qsnp=row$SNP
	id=row$ID
	pei=row$pheno_env_ID
	gts=row$gene_time_SNP
	time=row$time
	gene=row$Trait
	esnp=row$X_ID
	
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	results=results[results$X_ID==esnp & results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	
	df=data.frame(founder=founders,qtl_es=effect_size,eqtl_es=betas,stringsAsFactors=F)
	df$gene=gene
	df$pei=pei
	df$pvalue=p
	df$r=r
	df$time=time
	plotdf=rbind(df,plotdf)
	#p1=ggplot(df,aes(x=eqtl_es,y=qtl_es)) + geom_point(aes(color=founder)) +
	#xlab(sprintf("local-eQTL %s Expression Effect (log2cpm)",gene)) + ylab(" Male Flowering Effect (gdd)") +
	#ggtitle(sprintf("%s and %s,",gts,pei),subtitle=sprintf('r=%.2f',r)) +
	#theme(plot.title = element_text(size = 9))
	#plot_list[[count]]=p1
	#count=count+1
	
}

##### rap27_1 ####3
eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)


pheno="male_flowering_d6"
env="EXP_STPAUL_2017_WD"
qsnp="AX-91771656"
id="qDTA8"
chr=8

gene=rap27_1
# WD_0720 only
# local-eQTL in T20
testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
esnp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps


effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr,pheno,env),data.table=F)
effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
effect_size=unlist(effect_size[,c(6:21)])
wn=which(!is.na(effect_size))[1]
effect_size[-wn]=effect_size[-wn]+effect_size[wn]

time="WD_0720"
results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
results=results[results$X_ID==esnp & results$Trait==gene,]
betas=unlist(results[,c(6,10:24)])
wn=which(!is.na(betas))[1]
betas[-wn]=betas[-wn]+betas[wn]

test=cor.test(effect_size,betas,use="complete.obs")
r=test$estimate
p=test$p.value
	
df=data.frame(founder=founders,qtl_es=effect_size,eqtl_es=betas,stringsAsFactors=F)
df$gene=gene
df$pei=pei
df$pvalue=p
df$r=r
df$time=time

plotdf=rbind(plotdf,df)

##### rap27_2 ####3

pheno="male_flowering_d6"
env="EXP_STPAUL_2017_WD"
qsnp="AX-91771656"
id="qDTA8"
chr=8

gene=rap27_2

testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
esnp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps


effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr,pheno,env),data.table=F)
effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
effect_size=unlist(effect_size[,c(6:21)])
wn=which(!is.na(effect_size))[1]
effect_size[-wn]=effect_size[-wn]+effect_size[wn]

# which timpeoint is it expressed in?
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	print(time)
	print(gene %in% results$Trait)
}
# None of them
results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
results=results[results$X_ID==esnp & results$Trait==gene,]
betas=unlist(results[,c(6,10:24)])
wn=which(!is.na(betas))[1]
betas[-wn]=betas[-wn]+betas[wn]

test=cor.test(effect_size,betas,use="complete.obs")
r=test$estimate
p=test$p.value
	
df=data.frame(founder=founders,qtl_es=effect_size,eqtl_es=betas,stringsAsFactors=F)
df$gene=gene
df$pei=pei
df$pvalue=p
df$r=r
df$time=time

plotdf=rbind(plotdf,df)

###### zcn8 #####

pheno="male_flowering_d6"
env="EXP_STPAUL_2017_WD"
qsnp="AX-91771656"
id="qDTA8"
chr=8

gene=zcn8

testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
esnp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps


effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr,pheno,env),data.table=F)
effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
effect_size=unlist(effect_size[,c(6:21)])
wn=which(!is.na(effect_size))[1]
effect_size[-wn]=effect_size[-wn]+effect_size[wn]

# which timpeoint is it expressed in?
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	print(time)
	print(gene %in% results$Trait)
}
# All of them, Which has the highest r?
# local-eQTL? in T12, T20, & T27 - 

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	results=results[results$X_ID==esnp & results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	print(r)
}
#WD_0712
#0.2311321 
# WD_0718 
#-0.1762191 
#  WD_0720
#0.1556671 
#  WD_0727
#-0.002351313 
time="WD_0712"
results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
results=results[results$X_ID==esnp & results$Trait==gene,]
betas=unlist(results[,c(6,10:24)])
wn=which(!is.na(betas))[1]
betas[-wn]=betas[-wn]+betas[wn]
test=cor.test(effect_size,betas,use="complete.obs")
r=test$estimate
p=test$p.value

df=data.frame(founder=founders,qtl_es=effect_size,eqtl_es=betas,stringsAsFactors=F)
df$gene=gene
df$pei=pei
df$pvalue=p
df$r=r
df$time=time

plotdf=rbind(plotdf,df)

####### zmmads69 ######
pheno="male_flowering_d6"
env="EXP_STPAUL_2017_WD"
qsnp="PZE-103093413"
id="qDTA3_2"
chr=3

gene=mads69

testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
esnp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps


effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr,pheno,env),data.table=F)
effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
effect_size=unlist(effect_size[,c(6:21)])
wn=which(!is.na(effect_size))[1]
effect_size[-wn]=effect_size[-wn]+effect_size[wn]

# which timpeoint is it expressed in?
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	print(time)
	print(gene %in% results$Trait)
}
# All of them, Which has the highest r?
# local eQTL in all timepoints - increases in effect size over time
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	results=results[results$X_ID==esnp & results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	print(r)
}
#WD_0712
#0.4254668 
# WD_0718 
#-0.2588646 
#  WD_0720
#-0.5000195
#  WD_0727
#-0.02005968 
time="WD_0712" # positive in T12
results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
results=results[results$X_ID==esnp & results$Trait==gene,]
betas=unlist(results[,c(6,10:24)])
wn=which(!is.na(betas))[1]
betas[-wn]=betas[-wn]+betas[wn]
test=cor.test(effect_size,betas,use="complete.obs")
r=test$estimate
p=test$p.value

df=data.frame(founder=founders,qtl_es=effect_size,eqtl_es=betas,stringsAsFactors=F)
df$gene=gene
df$pei=pei
df$pvalue=p
df$r=r
df$time=time
plotdf=rbind(plotdf,df)


time="WD_0720" #Negative in T20?
results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
results=results[results$X_ID==esnp & results$Trait==gene,]
betas=unlist(results[,c(6,10:24)])
wn=which(!is.na(betas))[1]
betas[-wn]=betas[-wn]+betas[wn]
test=cor.test(effect_size,betas,use="complete.obs")
r=test$estimate
p=test$p.value

df=data.frame(founder=founders,qtl_es=effect_size,eqtl_es=betas,stringsAsFactors=F)
df$gene=gene
df$pei=pei
df$pvalue=p
df$r=r
df$time=time

plotdf=rbind(plotdf,df)

fwrite(plotdf,'paper_figures/data/FT_gene_effect_sizes.txt',row.names=F,quote=F,sep='\t')
#pdf('QTT/images/candidate_gene_effect_size_correlations.pdf')
#for(i in 1:length(plot_list)){#
#	print(plot_list[[i]])
#}
#dev.off()