#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")

allgenes=c()
for(time in times){
	norm=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
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


data=as.data.frame(matrix(nrow=length(allgenes),ncol=11))
names(data)=c('gene','WD_0712-WD_0718','WD_0712-WD_0720','WD_0712-WD_0727','WD_0718-WD_0720','WD_0718-WD_0727','WD_0720-WD_0727','WD_0712_snp_h2','WD_0718_snp_h2',"WD_0720_snp_h2","WD_0727_snp_h2")
data$gene=allgenes

for(time in times){
	column=paste0(time,'_snp_h2')
	#snp=fread(sprintf('eqtl/data/%s_SNP_resid_h2.txt',time),data.table=F)

	snp=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
	data[,column]=snp[match(data$gene,snp$gene),]$h2
}

data$mean_snp_h2=apply(data[,c('WD_0712_snp_h2','WD_0718_snp_h2',"WD_0720_snp_h2","WD_0727_snp_h2")],MARGIN=1,mean(na.rm=T))
data$sd_snp_h2=apply(data[,c('WD_0712_snp_h2','WD_0718_snp_h2',"WD_0720_snp_h2","WD_0727_snp_h2")],MARGIN=1,sd(na.rm=T))

data=fread('eqtl/data/h2_correlations_2.txt',data.table=F)

for(time in times){
	column=paste0(time,'_resid_snp_h2')
	snp=fread(sprintf('eqtl/data/%s_SNP_resid_h2_2.txt',time),data.table=F)

	#snp=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
	data[,column]=snp[match(data$gene,snp$gene),]$h2
}

data$mean_resid_snp_h2=apply(data[,c('WD_0712_resid_snp_h2','WD_0718_resid_snp_h2',"WD_0720_resid_snp_h2","WD_0727_resid_snp_h2")],MARGIN=1,function(x) mean(x,na.rm=T))
data$sd_resid_snp_h2=apply(data[,c('WD_0712_resid_snp_h2','WD_0718_resid_snp_h2',"WD_0720_resid_snp_h2","WD_0727_resid_snp_h2")],MARGIN=1,function(x) sd(x,na.rm=T))

mincutoff=1e-3

kept=data[data$mean_bs_h2>mincutoff & data$mean_resid_snp_h2>mincutoff,]
#18573

for(i in 1:3){
	time1=times[i]
	print(time1)
	norm1=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
	pc1=fread(sprintf('eqtl/normalized/%s_PCA_covariates_2.txt',time1),data.table=F)
	rownames(norm1)=norm1$V1
	for(j in (i+1):4){
		time2=times[j]
		print(time2)
		norm2=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time2),data.table=F)
		rownames(norm2)=norm2$V1
		pc2=fread(sprintf('eqtl/normalized/%s_PCA_covariates_2.txt',time2),data.table=F)
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

data$mean_bs_h2=apply(data[,c('WD_0712-WD_0718','WD_0718-WD_0720','WD_0720-WD_0727')],MARGIN=1,function(x) mean(x,na.rm=T))
data$sd_bs_h2=apply(data[,c('WD_0712-WD_0718','WD_0718-WD_0720','WD_0720-WD_0727')],MARGIN=1,function(x) sd(x,na.rm=T))

cor(data$mean_resid_snp_h2,data$mean_bs_h2,use="complete.obs")
#[1] 0.5790619


summary(data$mean_bs_h2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0121  0.0263  0.0447  0.0539  0.4888    2978

summary(data$`WD_0712-WD_0718`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0054  0.0242  0.0558  0.0713  0.7259    2806

summary(data$`WD_0712-WD_0720`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0050  0.0231  0.0518  0.0686  0.6800    2850

summary(data$`WD_0712-WD_0727`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0071  0.0311  0.0704  0.0894  0.8720    2979

summary(data$`WD_0718-WD_0720`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0034  0.0150  0.0318  0.0415  0.3822    2087 

summary(data$`WD_0718-WD_0727`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0068  0.0281  0.0557  0.0725  0.6478    2343

summary(data$`WD_0720-WD_0727`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0052  0.0206  0.0487  0.0592  0.5676    1691


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

pdf('images/h2_snp_vs_bs_resids_2.pdf')
for(i in 1:13){
	print(all_plots[[i]])
}
dev.off()

fwrite(data,'eqtl/data/h2_correlations_2.txt',row.names=F,quote=F,sep='\t')

# BS h2 and resid SNP h2

all_plots=list()
# WD_0712 SNP resid h2 and WD_0712-WD0718 BS h2
p1=ggplot(data,aes(x=WD_0712_resid_snp_h2,y=`WD_0712-WD_0718`)) + geom_point() 
all_plots[[1]]=p1

# WD_0718 SNP resid h2 and WD_0712-WD0718 BS h2
p1=ggplot(data,aes(x=WD_0718_resid_snp_h2,y=`WD_0712-WD_0718`)) + geom_point() 
all_plots[[2]]=p1

# WD_0712 SNP resid h2 and WD_0712-WD0720 BS h2
p1=ggplot(data,aes(x=WD_0712_resid_snp_h2,y=`WD_0712-WD_0720`)) + geom_point() 
all_plots[[3]]=p1

# WD_0720 SNP resid h2 and WD_0712-WD0720 BS h2
p1=ggplot(data,aes(x=WD_0720_resid_snp_h2,y=`WD_0712-WD_0720`)) + geom_point() 
all_plots[[4]]=p1

# WD_0712 SNP h2 and WD_0712-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0712_resid_snp_h2,y=`WD_0712-WD_0727`)) + geom_point() 
all_plots[[5]]=p1

# WD_0727 SNP h2 and WD_0712-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0727_resid_snp_h2,y=`WD_0712-WD_0727`)) + geom_point() 
all_plots[[6]]=p1

# WD_0718 SNP resid h2 and WD_0718-WD0720 BS h2
p1=ggplot(data,aes(x=WD_0718_resid_snp_h2,y=`WD_0718-WD_0720`)) + geom_point() 
all_plots[[7]]=p1

# WD_0720 SNP resid h2 and WD_0718-WD0720 BS h2
p1=ggplot(data,aes(x=WD_0720_resid_snp_h2,y=`WD_0718-WD_0720`)) + geom_point() 
all_plots[[8]]=p1

# WD_0718 SNP resid h2 and WD_0718-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0718_resid_snp_h2,y=`WD_0718-WD_0727`)) + geom_point() 
all_plots[[9]]=p1

# WD_0727 SNP resid h2 and WD_0718-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0727_resid_snp_h2,y=`WD_0718-WD_0727`)) + geom_point() 
all_plots[[10]]=p1

# WD_0720 SNP resid h2 and WD_0720-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0720_resid_snp_h2,y=`WD_0720-WD_0727`)) + geom_point() 
all_plots[[11]]=p1

# WD_0727 SNP resid h2 and WD_0720-WD0727 BS h2
p1=ggplot(data,aes(x=WD_0727_resid_snp_h2,y=`WD_0720-WD_0727`)) + geom_point() 
all_plots[[12]]=p1

# mean SNP resid h2 and mean BS h2
p1=ggplot(data,aes(x=mean_resid_snp_h2,y=`mean_bs_h2`)) + geom_point() 
all_plots[[13]]=p1

pdf('images/h2_snp_resids_vs_bs_resids_2.pdf')
for(i in 1:13){
	print(all_plots[[i]])
}
dev.off()


library(psych)
# r2: vector of correlation for each gene
# n_samp: number of samples

#49 "WD_0712-WD_0718" 40
#42 "WD_0712-WD_0720" 37
#24 "WD_0712-WD_0727" 31
#159  "WD_0718-WD_0720" 151
#108  "WD_0718-WD_0727" 111
#201 "WD_0720-WD_0727" 203




pdf('images/exp_correlations_2.pdf')
r2=data$`WD_0712-WD_0718`
n_samp=40 #49
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0712-WD_0720`
n_samp=37   #42
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0712-WD_0727`
n_samp=31 #24
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0718-WD_0720`
n_samp=151 #159
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0718-WD_0727`
n_samp=111 #108
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0720-WD_0727`
n_samp=203 #201
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)
dev.off()

# What is the correlation of heritabilities for with and without weights?

for(time in times){
	h20=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
	h2p=fread(sprintf('eqtl/data/%s_SNP_resid_h2_2.txt',time),data.table=F)
	h20$h2p=h2p[match(h20$gene,h2p$gene),]$h2
	print(cor(h20$h2,h20$h2p,use="complete.obs")**2)
}
#WD_0712 r2=0.1136042


