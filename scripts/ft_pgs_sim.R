#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')

allrep=c()
for(r in 1:100){
	rep=fread(sprintf('QTL/pgs_sim/Founder_ALL_DTA_rep%.0f_polygenic_scores_v3.txt',r),data.table=F)
	rep$rep=r
	allrep=rbind(allrep,rep)
}
allrep=as.data.frame(allrep,stringsAsFactors=F)

repvar=allrep %>% group_by(rep) %>% summarize(var=var(dta))

true=fread('QTL/Founder_ALL_DTA_polygenic_scores.txt',data.table=F)
truevar=var(true$dta)
truemean=mean(true$dta)
varquant=quantile(repvar$var,0.975)

p1=ggplot(repvar,aes(x=var)) + geom_density() + geom_vline(xintercept=truevar,color='red') + 
xlab("Variance in FT PGS") 

png('QTL/images/FT_pgs_sim_variance_v3.png')
print(p1)
dev.off()

repmean=allrep %>% group_by(rep) %>% summarize(mean=mean(dta))

meanquant=quantile(repmean$mean,0.975)
# mean is higher than 89.8% of reps

p2=ggplot(repmean,aes(x=mean)) + geom_density() + geom_vline(xintercept=truemean,color='red') + 
xlab("Avg FT PGS") 

png('QTL/images/FT_pgs_sim_mean_v3.png')
print(p2)
dev.off()

true$rep='real'
#allrep=rbind(allrep,true)

p3=ggplot(allrep,aes(x=dta,group=as.factor(rep))) + geom_density(color='black',alpha=0.4) +
geom_density(data=true,aes(x=dta),color='red')

png('QTL/images/FT_pgs_sim_data_v3.png')
print(p3)
dev.off()