#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('reshape2')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)
has_mite=which(mite==T,mite)
no_mite=which(mite==F,mite)

time="WD_0727"

fit=readRDS(sprintf('limma_results/%s_founder_model.rds',time))
#fit=readRDS('vgt1_founder_model.rds')
sd=fit$stdev.unscaled
foldchange=fit$coefficients
zcn8="Zm00001d010752"
rap27_1="Zm00001d010987"
rap27_2="Zm00001d010988"
mads69="Zm00001d042315"

allcoln=paste0('founder_f',founders)
coln=allcoln[allcoln %in% colnames(foldchange)]
has_mite=mite[allcoln %in% colnames(foldchange)]
present=founders[allcoln %in% colnames(foldchange)]

# ZmRap2.7

tmp3=foldchange[rap27_1,coln]
sd3=sd[rap27_1,coln]
melt1=melt(tmp3)
rownames(melt1)=present
sd_melt1=melt(sd3)
names(sd_melt1)=c('stdev')
rownames(sd_melt1)=present
melt1=cbind(melt1,sd_melt1$stdev)
melt1$mite=has_mite
melt1$founder=present
melt1=melt1[order(melt1$value),]
rownames(melt1)=seq(1,nrow(melt1))
melt1=as.data.frame(melt1,stringsAsFactors=F)
melt1$founder_f=factor(melt1$founder,levels=c(melt1$founder))
names(melt1)=c('value','stdev','mite','founder','founder_f')


png(sprintf('limma_results/images/%s_Rap27_1_WLS.png',time),width=800,heigh=650)
print(ggplot(melt1,aes(x=founder_f,y=value,color=mite))+geom_point()+geom_errorbar(aes(ymax=value+stdev,ymin=value-stdev),position="dodge") +
 xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle(sprintf("%s Expression",rap27_1)))
dev.off()

tmp3=foldchange[rap27_2,coln]
sd3=sd[rap27_2,coln]
melt1=melt(tmp3)
rownames(melt1)=present
sd_melt1=melt(sd3)
names(sd_melt1)=c('stdev')
rownames(sd_melt1)=present
melt1=cbind(melt1,sd_melt1$stdev)
melt1$mite=has_mite
melt1$founder=present
melt1=melt1[order(melt1$value),]
rownames(melt1)=seq(1,nrow(melt1))
melt1=as.data.frame(melt1,stringsAsFactors=F)
melt1$founder_f=factor(melt1$founder,levels=c(melt1$founder))
names(melt1)=c('value','stdev','mite','founder','founder_f')


png(sprintf('limma_results/images/%s_Rap27_2_WLS.png',time),width=800,heigh=650)
print(ggplot(melt1,aes(x=founder_f,y=value,color=mite))+geom_point()+geom_errorbar(aes(ymax=value+stdev,ymin=value-stdev),position="dodge") +
 xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle(sprintf("%s Expression",rap27_2)))
dev.off()


#ZCN8

tmp1=foldchange[zcn8,coln]
sd1=sd[zcn8,coln]
melt1=melt(tmp1)
rownames(melt1)=present
sd_melt1=melt(sd1)
names(sd_melt1)=c('stdev')
rownames(sd_melt1)=present
melt1=cbind(melt1,sd_melt1$stdev)
melt1$mite=has_mite
melt1$founder=present
melt1=melt1[order(melt1$value),]
rownames(melt1)=seq(1,nrow(melt1))
melt1=as.data.frame(melt1,stringsAsFactors=F)
melt1$founder_f=factor(melt1$founder,levels=c(melt1$founder))
names(melt1)=c('value','stdev','mite','founder','founder_f')


png(sprintf('limma_results/images/%s_ZCN8_WLS.png',time),width=800,heigh=650)
print(ggplot(melt1,aes(x=founder_f,y=value,color=mite))+geom_point()+geom_errorbar(aes(ymax=value+stdev,ymin=value-stdev),position="dodge") + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZCN8 Expression"))
dev.off()




# ZmMADS69

tmp2=foldchange[mads69,coln]
sd1=sd[zcn8,coln]
melt1=melt(tmp2)
rownames(melt1)=present
sd_melt1=melt(sd1)
names(sd_melt1)=c('stdev')
rownames(sd_melt1)=present
melt1=cbind(melt1,sd_melt1$stdev)
melt1$mite=has_mite
melt1$founder=present
melt1=melt1[order(melt1$value),]
rownames(melt1)=seq(1,nrow(melt1))
melt1=as.data.frame(melt1,stringsAsFactors=F)
melt1$founder_f=factor(melt1$founder,levels=c(melt1$founder))
names(melt1)=c('value','stdev','mite','founder','founder_f')

png(sprintf('limma_results/images/%s_ZmMADS69_WLS.png',time),width=800,heigh=650)
print(ggplot(melt1,aes(x=founder_f,y=value,color=mite))+geom_point()+geom_errorbar(aes(ymax=value+stdev,ymin=value-stdev),position="dodge") + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZmMADS69 Expression"))
dev.off()


### Test genes
testgene="Zm00001d010804"

testgene="Zm00001d010946"

testgene="Zm00001d011152"

testgene="Zm00001d011197"

tmp2=foldchange[testgene,coln]
sd1=sd[zcn8,coln]
melt1=melt(tmp2)
rownames(melt1)=present
sd_melt1=melt(sd1)
names(sd_melt1)=c('stdev')
rownames(sd_melt1)=present
melt1=cbind(melt1,sd_melt1$stdev)
melt1$mite=has_mite
melt1$founder=present
melt1=melt1[order(melt1$value),]
rownames(melt1)=seq(1,nrow(melt1))
melt1=as.data.frame(melt1,stringsAsFactors=F)
melt1$founder_f=factor(melt1$founder,levels=c(melt1$founder))
names(melt1)=c('value','stdev','mite','founder','founder_f')

png(sprintf('limma_results/images/%s_%s_WLS.png',time,testgene),width=800,heigh=650)
print(ggplot(melt1,aes(x=founder_f,y=value,color=mite))+geom_point()+geom_errorbar(aes(ymax=value+stdev,ymin=value-stdev),position="dodge") +
 xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle(sprintf("%s Expression",testgene)))
dev.off()

## MITE vs. NO_MITE

fit=readRDS(sprintf('limma_results/%s_MITE_full_model.rds',time))
zcn8="Zm00001d010752"
rap27="Zm00001d010987"
mads69="Zm00001d042315"

#fit=readRDS('vgt1_founder_model.rds')
sd=fit$stdev.unscaled
foldchange=fit$coefficients

#ZCN8

tmp2=foldchange[zcn8,]
melt2=melt(tmp2)
melt2$variable=rownames(melt2)
#melt2$mite=c(T,T,F,F,T,T,T,T,T,T,T,T,T,F,F,T)
#melt2=melt2[order(melt2$value),]
#melt2=as.data.frame(melt2,stringsAsFactors=F)
#rownames(melt2)=seq(1,dim(melt2)[1])
#melt2$founder_f=factor(melt2$variable,levels=melt2$variable)

png(sprintf('limma_results/images/%s_ZCN8_MITE_logFC.png',time),width=800,heigh=650)
print(ggplot(melt2,aes(x=variable,y=value))+geom_boxplot() + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZCN8 Expression in Lines With and Without MITE at vgt1"))
dev.off()


# ZmRap2.7

tmp3=foldchange[rap27,]
melt3=melt(tmp3)
melt3$variable=rownames(melt3)

#melt3$mite=c(T,T,F,F,T,T,T,T,T,T,T,T,T,F,F,T)
#melt3=melt3[order(melt3$value),]
#melt3=as.data.frame(melt3,stringsAsFactors=F)
#rownames(melt3)=seq(1,dim(melt3)[1])
#melt3$founder_f=factor(melt3$variable,levels=melt3$variable)

png(sprintf('limma_results/images/%s_ZmRap2.7_MITE_logFC.png',time),width=800,heigh=650)
print(ggplot(melt3,aes(x=variable,y=value))+geom_boxplot() + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZmRap2.7 Expression in Lines With and Without MITE at vgt1"))
dev.off()

# ZmMADS69

tmp1=foldchange[mads69,]
melt1=melt(tmp1)
melt1$variable=rownames(melt1)

#melt1$mite=c(T,T,F,F,T,T,T,T,T,T,T,T,T,F,F,T)
#melt1=melt1[order(melt1$value),]
#melt1=as.data.frame(melt1,stringsAsFactors=F)
#rownames(melt1)=seq(1,dim(melt1)[1])
#melt1$founder_f=factor(melt1$variable,levels=melt1$variable)

png(sprintf('limma_results/images/%s_ZmMADS69_MITE_logFC.png',time),width=800,heigh=650)
print(ggplot(melt1,aes(x=variable,y=value))+geom_boxplot() + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZmMADS69 Expression in Lines With and Without MITE at vgt1"))
dev.off()


#degs=fread('NO_MITE_v_MITE_DEGs.txt')
