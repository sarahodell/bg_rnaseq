#!/usr/bin/env Rscript

library('data.table')
library('stats')
library('ggplot2')

beta_merge=fread('eqtl/results/beta_merge_all_v3.txt',data.table=F)
beta_merge$pheno_f=as.factor(beta_merge$phenotype)
beta_merge$abs_beta_z=abs(beta_merge$beta_z)


# Extreme bins are >=3/<=-3 SD
bbreaks=c(0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,max(beta_merge$abs_beta_z)+0.1)
beta_merge$abs_bin=cut(beta_merge$abs_beta_z, breaks=bbreaks, label=F)

beta_merge=beta_merge[beta_merge$freq>=0.01,]
submerge=beta_merge[beta_merge$phenotype=="harvest_grain_moisture",]
avg_z_abs=submerge %>% group_by(abs_bin) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),n=length(freq),ngenes=length(unique(gene_founder)))

log_fit <- drm(freq ~ abs_beta_z, fct = L.4(), data = submerge)

log_fit <- drm(avg_freq ~ abs_bin, fct = L.4(),weights=avg_z_abs$ngenes,data = avg_z_abs)
log_logist_fit <- drm(avg_freq ~ abs_bin, fct = LL.4(),weights=avg_z_abs$ngenes,data = avg_z_abs)


png('HGM_test.png')
plot(log_fit)
plot(log_logist_fit, add=T, col = "red")
dev.off()


avg_z_abs=fread('paper_figures/pheno_abs_bin_model_data.txt',data.table=F)
avg_z_abs=avg_z_abs[avg_z_abs$ngenes>=100,]

phenos=unique(avg_z_abs$phenotype)
colors=palette()[2:7]
p=ggplot(avg_z_abs,aes(x=abs_bin,y=avg_freq,group=phenotype)) + geom_point(aes(color=phenotype,size=ngenes)) +
scale_color_manual(values=colors,breaks=phenos,labels=c("None","DTS","HGM","DTA","TKW","TPH"))

n=1
newpreds=c()
for(pheno in phenos){
	submerge=avg_z_abs[avg_z_abs$phenotype==pheno,]
	log_fit <- drm(avg_freq ~ abs_bin, fct = L.4(),weights=submerge$ngenes,data = submerge)
	newdata=submerge[,c('abs_bin'),drop=F]
	pm <- predict(log_fit, newdata=newdata, interval="confidence") 
    newdata$p <- pm[,1]
    newdata$pmin <- pm[,2]
    newdata$pmax <- pm[,3]
    newdata$phenotype=pheno
    newdata=rbind(newpreds,newdata)
    coln=colors[n]
    p=p+geom_ribbon(data=newdata, aes(x=abs_bin, y=p, ymin=pmin, ymax=pmax),fill=coln,alpha=0.2) +
    geom_line(data=newdata, aes(x=abs_bin, y=p),color=coln)
    n=n+1

}	

png('paper_figures/log_drm_test.png')
print(p)
dev.off()


avg_z_abs=fread('paper_figures/pheno_abs_bin_model_data.txt',data.table=F)
avg_z_abs=avg_z_abs[avg_z_abs$ngenes>=100,]
phenos=unique(avg_z_abs$phenotype)
colors=palette()[2:7]
p=ggplot(avg_z_abs,aes(x=abs_bin,y=avg_freq,group=phenotype)) + geom_point(aes(color=phenotype,size=ngenes)) +
scale_color_manual(values=colors,breaks=phenos,labels=c("None","DTS","HGM","DTA","TKW","TPH"))


n=1
newpreds=c()
for(pheno in phenos){
	submerge=avg_z_abs[avg_z_abs$phenotype==pheno,]
	log_fit <- drm(avg_freq ~ abs_bin, fct = LL.4(),weights=submerge$ngenes,data = submerge)
	newdata=submerge[,c('abs_bin'),drop=F]
	pm <- predict(log_fit, newdata=newdata, interval="confidence") 
    newdata$p <- pm[,1]
    newdata$pmin <- pm[,2]
    newdata$pmax <- pm[,3]
    newdata$phenotype=pheno
    newdata=rbind(newpreds,newdata)
    coln=colors[n]
    p=p+geom_ribbon(data=newdata, aes(x=abs_bin, y=p, ymin=pmin, ymax=pmax),fill=coln,alpha=0.2) +
    geom_line(data=newdata, aes(x=abs_bin, y=p),color=coln)
    n=n+1

}	

png('paper_figures/log_logist_drm_test.png')
print(p)
dev.off()