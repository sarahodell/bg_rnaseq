#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library(psych)

times=c("WD_0712","WD_0718","WD_0720","WD_0727")

allgenes=c()
allinds=c()
for(time in times){
	norm=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	allgenes=union(allgenes,names(norm)[-1])
	allinds=union(allinds,norm$V1)
}
#19554 now 18950
# 278 individuals
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

calc_cor2=function(row){
	d1=unlist(norm1[row,])
	#names(d1)=c('exp1')
	#d1$PC1=pc1[match(names(d1),pc1$V1),]$PC1
    #d1$PC2=pc1[match(names(d1)[1],pc1$V1),]$PC2
    #d1$PC3=pc1[match(names(d1)[1],pc1$V1),]$PC3
	#m1=lm(exp1 ~ PC1 + PC2 + PC3,d1)
	#resid1=summary(m1)$residuals
	
	d2=unlist(norm2[row,])
	#names(d2)=c('exp2')
	#d2$PC1=pc1[match(rownames(d2),pc1$V1),]$PC1
    #d2$PC2=pc1[match(rownames(d2),pc1$V1),]$PC2
    #d2$PC3=pc1[match(rownames(d2),pc1$V1),]$PC3
	#m2=lm(exp2 ~ PC1 + PC2 + PC3,d2)
	#resid2=summary(m2)$residuals
	r2=cor(d1,d2,use="complete.obs")**2
	return(r2)
}


data=as.data.frame(matrix(nrow=length(allgenes),ncol=11))
names(data)=c('gene','WD_0712-WD_0718','WD_0712-WD_0720','WD_0712-WD_0727','WD_0718-WD_0720','WD_0718-WD_0727','WD_0720-WD_0727','WD_0712_resid_snp_h2','WD_0718_resid_snp_h2',"WD_0720_resid_snp_h2","WD_0727_resid_snp_h2")
data$gene=allgenes

#for(time in times){
#	column=paste0(time,'_snp_h2')
#	#snp=fread(sprintf('eqtl/data/%s_SNP_resid_h2.txt',time),data.table=F)
#
#	snp=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
#	data[,column]=snp[match(data$gene,snp$gene),]$h2
#}

data$mean_snp_h2=apply(data[,c('WD_0712_snp_h2','WD_0718_snp_h2',"WD_0720_snp_h2","WD_0727_snp_h2")],MARGIN=1,mean(na.rm=T))
data$sd_snp_h2=apply(data[,c('WD_0712_snp_h2','WD_0718_snp_h2',"WD_0720_snp_h2","WD_0727_snp_h2")],MARGIN=1,sd(na.rm=T))

data=fread('eqtl/data/gene_h2_correlations_2.txt',data.table=F)

data2=fread('eqtl/data/ind_h2_correlations_2.txt',data.table=F)

for(time in times){
	column=paste0(time,'_resid_snp_h2')
	snp=fread(sprintf('eqtl/data/%s_SNP_resid_h2_2.txt',time),data.table=F)

	#snp=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
	data[,column]=snp[match(data$gene,snp$gene),]$h2
}

data$mean_resid_snp_h2=apply(data[,c('WD_0712_resid_snp_h2','WD_0718_resid_snp_h2',"WD_0720_resid_snp_h2","WD_0727_resid_snp_h2")],MARGIN=1,function(x) mean(x,na.rm=T))
data$sd_resid_snp_h2=apply(data[,c('WD_0712_resid_snp_h2','WD_0718_resid_snp_h2',"WD_0720_resid_snp_h2","WD_0727_resid_snp_h2")],MARGIN=1,function(x) sd(x,na.rm=T))

data2=as.data.frame(matrix(nrow=length(allinds),ncol=7))
names(data2)=c('ind','WD_0712-WD_0718','WD_0712-WD_0720','WD_0712-WD_0727','WD_0718-WD_0720','WD_0718-WD_0727','WD_0720-WD_0727')
data2$ind=allinds

#mincutoff=1e-3

#kept=data[data$mean_bs_h2>mincutoff & data$mean_resid_snp_h2>mincutoff,]
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
			print(length(sharedg))
			#print(length(sharedg))
			norm1=norm1[,sharedg]
			norm2=norm2[,sharedg]
			#r2=sapply(seq(1,length(sharedg)),function(x) calc_cor(x))
			r2_ind=sapply(seq(1,length(sharedi)),function(x) calc_cor2(x))
			#subdata=data.frame(gene=sharedg,r2=r2,stringsAsFactors=F)
			subdata2=data.frame(ind=sharedi,r2=r2_ind,stringsAsFactors=F)
			column=sprintf('%s-%s',time1,time2)
			print(column)
			#data[,column]=subdata[match(data$gene,subdata$gene),]$r2
			data2[,column]=subdata2[match(data2$ind,subdata2$ind),]$r2
		}
	}
}

#data$mean_bs_h2=apply(data[,c('WD_0712-WD_0718','WD_0712-WD_0720','WD_0712-WD_0727','WD_0718-WD_0720','WD_0718-WD_0727','WD_0720-WD_0727')],MARGIN=1,function(x) mean(x,na.rm=T))
#data$sd_bs_h2=apply(data[,c('WD_0712-WD_0718','WD_0712-WD_0720','WD_0712-WD_0727','WD_0718-WD_0720','WD_0718-WD_0727','WD_0720-WD_0727')],MARGIN=1,function(x) sd(x,na.rm=T))

#fwrite(data,'eqtl/data/gene_h2_correlations_2.txt',row.names=F,quote=F,sep='\t')

data2$mean_bs_h2=apply(data2[,c('WD_0712-WD_0718','WD_0712-WD_0720','WD_0712-WD_0727','WD_0718-WD_0720','WD_0718-WD_0727','WD_0720-WD_0727')],MARGIN=1,function(x) mean(x,na.rm=T))
data2$sd_bs_h2=apply(data2[,c('WD_0712-WD_0718','WD_0712-WD_0720','WD_0712-WD_0727','WD_0718-WD_0720','WD_0718-WD_0727','WD_0720-WD_0727')],MARGIN=1,function(x) sd(x,na.rm=T))

fwrite(data2,'eqtl/data/ind_h2_correlations_2.txt',row.names=F,quote=F,sep='\t')

#summary(data2$mean_bs_h2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.5127  0.7558  0.8029  0.7858  0.8287  0.8898      82 

pdf('images/exp_ind_correlations_2.pdf')
r2=data2$`WD_0712-WD_0718`[!is.na(data2$`WD_0712-WD_0718`)]
n_samp=15926
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data2$`WD_0712-WD_0720`[!is.na(data2$`WD_0712-WD_0720`)]
n_samp=15891   
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data2$`WD_0712-WD_0727`[!is.na(data2$`WD_0712-WD_0727`)]
n_samp=15803 
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data2$`WD_0718-WD_0720`[!is.na(data2$`WD_0718-WD_0720`)]
n_samp=17188
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data2$`WD_0718-WD_0727`[!is.na(data2$`WD_0718-WD_0727`)]
n_samp=16932
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data2$`WD_0720-WD_0727`[!is.na(data2$`WD_0720-WD_0727`)]
n_samp=17487
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)
dev.off()



cor(data$mean_resid_snp_h2,data$mean_bs_h2,use="complete.obs")
#[1] 0.5790619


summary(data$mean_bs_h2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0302  0.0570  0.1052  0.1116  0.8562    1174 

summary(data$`WD_0712-WD_0718`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.0000  0.0071  0.0307  0.0768  0.0936  0.8764    3024

summary(data$`WD_0712-WD_0720`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0084  0.0356  0.0863  0.1051  0.9338    3059

summary(data$`WD_0712-WD_0727`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0112  0.0484  0.1022  0.1364  0.9450    3147 

summary(data$`WD_0718-WD_0720`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0132  0.0530  0.1154  0.1462  0.8875    1762

summary(data$`WD_0718-WD_0727`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0127  0.0514  0.1111  0.1396  0.8865    2018

summary(data$`WD_0720-WD_0727`)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0130  0.0491  0.1102  0.1320  0.8694    1463


all_plots=list()

# WD_0712 SNP h2 and WD_0712-WD0718 BS h2
r=cor(data$WD_0712_resid_snp_h2,data$`WD_0712-WD_0718`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0712_resid_snp_h2,y=`WD_0712-WD_0718`)) + geom_point() + ggtitle(sprintf("WD_0712 h^2 and H^2 WD_0712-WD_0718, r=%.2f",r))
all_plots[[1]]=p1

# WD_0718 SNP h2 and WD_0712-WD0718 BS h2
r=cor(data$WD_0718_resid_snp_h2,data$`WD_0712-WD_0718`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0718_resid_snp_h2,y=`WD_0712-WD_0718`)) + geom_point()  + ggtitle(sprintf("WD_0718 h^2 and H^2 WD_0712-WD_0718, r=%.2f",r))
all_plots[[2]]=p1

# WD_0712 SNP h2 and WD_0712-WD0720 BS h2
r=cor(data$WD_0712_resid_snp_h2,data$`WD_0712-WD_0720`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0712_resid_snp_h2,y=`WD_0712-WD_0720`)) + geom_point() + ggtitle(sprintf("WD_0712 h^2 and H^2 WD_0712-WD_0720, r=%.2f",r))
all_plots[[3]]=p1

# WD_0720 SNP h2 and WD_0712-WD0720 BS h2
r=cor(data$WD_0720_resid_snp_h2,data$`WD_0712-WD_0720`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0720_resid_snp_h2,y=`WD_0712-WD_0720`)) + geom_point() + ggtitle(sprintf("WD_0720 h^2 and H^2 WD_0712-WD_0720, r=%.2f",r))
all_plots[[4]]=p1

# WD_0712 SNP h2 and WD_0712-WD0727 BS h2
r=cor(data$WD_0712_resid_snp_h2,data$`WD_0712-WD_0727`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0712_resid_snp_h2,y=`WD_0712-WD_0727`)) + geom_point() + ggtitle(sprintf("WD_0712 h^2 and H^2 WD_0712-WD_0727, r=%.2f",r))
all_plots[[5]]=p1

# WD_0727 SNP h2 and WD_0712-WD0727 BS h2
r=cor(data$WD_0727_resid_snp_h2,data$`WD_0712-WD_0727`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0727_resid_snp_h2,y=`WD_0712-WD_0727`)) + geom_point() + ggtitle(sprintf("WD_0727 h^2 and H^2 WD_0712-WD_0727, r=%.2f",r))
all_plots[[6]]=p1

# WD_0718 SNP h2 and WD_0718-WD0720 BS h2
r=cor(data$WD_0718_resid_snp_h2,data$`WD_0718-WD_0720`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0718_resid_snp_h2,y=`WD_0718-WD_0720`)) + geom_point() + ggtitle(sprintf("WD_0718 h^2 and H^2 WD_0718-WD_0720, r=%.2f",r))
all_plots[[7]]=p1

# WD_0720 SNP h2 and WD_0718-WD0720 BS h2
r=cor(data$WD_0720_resid_snp_h2,data$`WD_0718-WD_0720`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0720_resid_snp_h2,y=`WD_0718-WD_0720`)) + geom_point() + ggtitle(sprintf("WD_0720 h^2 and H^2 WD_0718-WD_0720, r=%.2f",r))
all_plots[[8]]=p1

# WD_0718 SNP h2 and WD_0718-WD0727 BS h2
r=cor(data$WD_0718_resid_snp_h2,data$`WD_0718-WD_0727`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0718_resid_snp_h2,y=`WD_0718-WD_0727`)) + geom_point() + ggtitle(sprintf("WD_0718 h^2 and H^2 WD_0718-WD_0727, r=%.2f",r))
all_plots[[9]]=p1

# WD_0727 SNP h2 and WD_0718-WD0727 BS h2
r=cor(data$WD_0727_resid_snp_h2,data$`WD_0718-WD_0727`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0727_resid_snp_h2,y=`WD_0718-WD_0727`)) + geom_point() + ggtitle(sprintf("WD_0727 h^2 and H^2 WD_0718-WD_0727, r=%.2f",r))
all_plots[[10]]=p1

# WD_0720 SNP h2 and WD_0720-WD0727 BS h2
r=cor(data$WD_0720_resid_snp_h2,data$`WD_0720-WD_0727`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0720_resid_snp_h2,y=`WD_0720-WD_0727`)) + geom_point() + ggtitle(sprintf("WD_0720 h^2 and H^2 WD_0720-WD_0727, r=%.2f",r))
all_plots[[11]]=p1

# WD_0727 SNP h2 and WD_0720-WD0727 BS h2
r=cor(data$WD_0727_resid_snp_h2,data$`WD_0720-WD_0727`,use="complete.obs")
p1=ggplot(data,aes(x=WD_0727_resid_snp_h2,y=`WD_0720-WD_0727`)) + geom_point() + ggtitle(sprintf("WD_0727 h^2 and H^2 WD_0720-WD_0727, r=%.2f",r))
all_plots[[12]]=p1

#mena snp h2 vs mean bs h2
r=cor(data$mean_resid_snp_h2,data$mean_bs_h2,use="complete.obs")
p1=ggplot(data,aes(x=mean_resid_snp_h2,y=mean_bs_h2)) + geom_point() + ggtitle(sprintf("mean h^2 and mean H^2, r=%.2f",r))
all_plots[[13]]=p1

pdf('images/h2_snp_vs_bs_resids_2.pdf')
for(i in 1:13){
	print(all_plots[[i]])
}
dev.off()



# r2: vector of correlation for each gene
# n_samp: number of samples

# "WD_0712-WD_0718" 37
# "WD_0712-WD_0720" 32
# "WD_0712-WD_0727" 22
# "WD_0718-WD_0720" 110
# "WD_0718-WD_0727" 71
# "WD_0720-WD_0727" 152

library(psych)


pdf('images/exp_correlations_2.pdf')
r2=data$`WD_0712-WD_0718`
n_samp=37 #49
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0712-WD_0720`
n_samp=32   #42
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0712-WD_0727`
n_samp=22 #24
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0718-WD_0720`
n_samp=110 #159
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0718-WD_0727`
n_samp=71 #108
qqplot(t2r(rnorm(length(r2)),n_samp)^2,r2);abline(0,1)

r2=data$`WD_0720-WD_0727`
n_samp=203 #152
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


