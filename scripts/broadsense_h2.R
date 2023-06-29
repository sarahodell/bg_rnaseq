#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")

allgenes=c()
for(time in times){
	norm=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
	allgenes=union(allgenes,names(norm)[-1])
}
#19554

calc_cor=function(col){
	d1=norm1[,col,drop=F]
	names(d1)=c('exp1')
	d1$PC1=pc1[match(rownames(d1),pc1$V1),]$PC1
    d1$PC2=pc1[match(rownames(d1),pc1$V1),]$PC2
    d1$PC3=pc1[match(rownames(d1),pc1$V1),]$PC3
	m1=lm(exp1 ~ PC1 + PC2 + PC3,d1)
	resid1=summary(m1)$residuals
	
	d2=norm2[,col,drop=F]
	names(d2)=c('exp2')
	d2$PC1=pc1[match(rownames(d2),pc1$V1),]$PC1
    d2$PC2=pc1[match(rownames(d2),pc1$V1),]$PC2
    d2$PC3=pc1[match(rownames(d2),pc1$V1),]$PC3
	m2=lm(exp2 ~ PC1 + PC2 + PC3,d2)
	resid2=summary(m2)$residuals
	r2=cor(resid1,resid2,use="complete.obs")**2
	return(r2)
}

mincutoff=1e-3

data=as.data.frame(matrix(nrow=length(allgenes),ncol=11))
names(data)=c('gene','WD_0712-WD_0718','WD_0712-WD_0720','WD_0712-WD_0727','WD_0718-WD_0720','WD_0718-WD_0727','WD_0720-WD_0727','WD_0712_snp_h2','WD_0718_snp_h2',"WD_0720_snp_h2","WD_0727_snp_h2")
data$gene=allgenes

for(time in times){
	column=paste0(time,'_snp_h2')
	#snp=fread(sprintf('eqtl/data/%s_SNP_resid_h2.txt',time),data.table=F)

	snp=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
	data[,column]=snp[match(data$gene,snp$gene),]$h2
}

data$mean_snp_h2=apply(data[,c('WD_0712_snp_h2','WD_0718_snp_h2',"WD_0720_snp_h2","WD_0727_snp_h2")],MARGIN=1,mean)
data$sd_snp_h2=apply(data[,c('WD_0712_snp_h2','WD_0718_snp_h2',"WD_0720_snp_h2","WD_0727_snp_h2")],MARGIN=1,sd)

data=fread('eqtl/data/h2_correlations.txt',data.table=F)

for(time in times){
	column=paste0(time,'_resid_snp_h2')
	snp=fread(sprintf('eqtl/data/%s_SNP_resid_h2.txt',time),data.table=F)

	#snp=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
	data[,column]=snp[match(data$gene,snp$gene),]$h2
}

for(i in 1:3){
	time1=times[i]
	print(time1)
	norm1=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time1),data.table=F)
	pc1=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time1),data.table=F)
	rownames(norm1)=norm1$V1
	for(j in (i+1):4){
		time2=times[j]
		print(time2)
		norm2=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time2),data.table=F)
		rownames(norm2)=norm2$V1
		pc2=fread(sprintf('eqtl/normalized/%s_PCA_covariates.txt',time2),data.table=F)
		sharedi=intersect(rownames(norm1),rownames(norm2))
		print(length(sharedi))
		norm1=norm1[sharedi,]
		norm2=norm2[sharedi,]
		if(length(sharedi)!=0){
			sharedg=intersect(names(norm1)[-1],names(norm2)[-1])
			#print(length(sharedg))
			norm1=norm1[,sharedg]
			norm2=norm2[,sharedg]
			r2=sapply(seq(1,length(sharedg)),function(x) calc_cor(x))
			subdata=data.frame(gene=sharedg,r2=r2,stringsAsFactors=F)
			column=sprintf('%s-%s',time1,time2)
			print(column)
			data[,column]=subdata[match(data$gene,subdata$gene),]$r2
		}	
	}
}

data$mean_bs_h2=apply(data[,c('WD_0712-WD_0718','WD_0718-WD_0720','WD_0720-WD_0727')],MARGIN=1,mean)
data$sd_bs_h2=apply(data[,c('WD_0712-WD_0718','WD_0718-WD_0720','WD_0720-WD_0727')],MARGIN=1,sd)

cor(data$mean_snp_h2,data$mean_bs_h2,use="complete.obs")
#[1] -0.0003982029

all_plots=list()

# WD_0712 SNP h2 and WD_0712-WD0718 BS h2
p1=ggplot(data,aes(x=WD_0712_snp_h2,y=`WD_0712-WD_0718`)) + geom_point() 
all_plots[[1]]=p1

# WD_0718 SNP h2 and WD_0712-WD0718 BS h2
p1=ggplot(data,aes(x=WD_0718_snp_h2,y=`WD_0712-WD_0718`)) + geom_point() 
all_plots[[2]]=p1

# WD_0712 SNP h2 and WD_0712-WD0720 BS h2
p1=ggplot(data,aes(x=WD_0712_snp_h2,y=`WD_0712-WD_0720`)) + geom_point() 
all_plots[[3]]=p1

# WD_0720 SNP h2 and WD_0712-WD0720 BS h2
p1=ggplot(data,aes(x=WD_0720_snp_h2,y=`WD_0712-WD_0720`)) + geom_point() 
all_plots[[4]]=p1

# WD_0712 SNP h2 and WD_0712-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0712_snp_h2,y=`WD_0712-WD_0727`)) + geom_point() 
all_plots[[5]]=p1

# WD_0727 SNP h2 and WD_0712-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0727_snp_h2,y=`WD_0712-WD_0727`)) + geom_point() 
all_plots[[6]]=p1

# WD_0718 SNP h2 and WD_0718-WD0720 BS h2
p1=ggplot(data,aes(x=WD_0718_snp_h2,y=`WD_0718-WD_0720`)) + geom_point() 
all_plots[[7]]=p1

# WD_0720 SNP h2 and WD_0718-WD0720 BS h2
p1=ggplot(data,aes(x=WD_0720_snp_h2,y=`WD_0718-WD_0720`)) + geom_point() 
all_plots[[8]]=p1

# WD_0718 SNP h2 and WD_0718-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0718_snp_h2,y=`WD_0718-WD_0727`)) + geom_point() 
all_plots[[9]]=p1

# WD_0727 SNP h2 and WD_0718-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0727_snp_h2,y=`WD_0718-WD_0727`)) + geom_point() 
all_plots[[10]]=p1

# WD_0720 SNP h2 and WD_0720-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0720_snp_h2,y=`WD_0720-WD_0727`)) + geom_point() 
all_plots[[11]]=p1

# WD_0727 SNP h2 and WD_0720-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0727_snp_h2,y=`WD_0720-WD_0727`)) + geom_point() 
all_plots[[12]]=p1

#mena snp h2 vs mean bs h2
p1=ggplot(data,aes(x=mean_snp_h2,y=mean_bs_h2)) + geom_point() 
all_plots[[13]]=p1

pdf('images/h2_snp_vs_bs_resids.pdf')
for(i in 1:13){
	print(all_plots[[i]])
}
dev.off()

fwrite(data,'eqtl/data/h2_correlations.txt',row.names=F,quote=F,sep='\t')

# BS h2 and resid SNP h2

all_plots=list()
# WD_0712 SNP resid h2 and WD_0712-WD0718 BS h2
p1=ggplot(data,aes(x=WD_0712_resid_snp_h2,y=`WD_0712-WD_0718`)) + geom_point() 
all_plots[[1]]=p1

# WD_0712 SNP resid h2 and WD_0712-WD0720 BS h2
p1=ggplot(data,aes(x=WD_0712_resid_snp_h2,y=`WD_0712-WD_0720`)) + geom_point() 
all_plots[[2]]=p1

# WD_0712 SNP h2 and WD_0712-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0712_resid_snp_h2,y=`WD_0712-WD_0727`)) + geom_point() 
all_plots[[3]]=p1

pdf('images/h2_snp_resids_vs_bs_resids.pdf')
for(i in 1:3){
	print(all_plots[[i]])
}
dev.off()


library(psych);
# r2: vector of correlation for each gene
# n_samp: number of samples

#49 "WD_0712-WD_0718"
#42 "WD_0712-WD_0720"
#24 "WD_0712-WD_0727"
#159  "WD_0718-WD_0720"
#108  "WD_0718-WD_0727"
#201 "WD_0720-WD_0727"




pdf('images/exp_correlations.pdf')
r2=data$`WD_0712-WD_0718`
n_samp=49
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0712-WD_0720`
n_samp=42
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0712-WD_0727`
n_samp=24
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0718-WD_0720`
n_samp=159
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0718-WD_0727`
n_samp=108
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0720-WD_0727`
n_samp=201
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)
dev.off()

# What is the correlation of heritabilities for with and without weights?

for(time in times){
	h20=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
	h2p=fread(sprintf('eqtl/data/%s_SNP_resid_h2.txt',time),data.table=F)
	h20$h2p=h2p[match(h20$gene,h2p$gene),]$h2
	print(cor(h20$h2,h20$h2p,use="complete.obs")**2)
}
#WD_0712 r2=0.1136042


