#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

mite=c(T,F,T,T,F,T,T,T,T,T,T,F,T,T,T,F)
has_mite=which(mite==T,mite)
no_mite=which(mite==F,mite)

fit=readRDS('vgt1_founder_model.rds')
sd=fit$stdev.unscaled
foldchange=fit$coefficients
zcn8="Zm00001d010752"
rap27="Zm00001d010987"
mads69="Zm00001d042315"

#ZCN8

tmp1=foldchange[rownames(foldchange)==zcn8,]
sd1=sd[rownames(sd)==zcn8,]
melt1=melt(tmp1)
rownames(melt1)=c("A632_usa","A654_inra","B73_inra","B96",
"C103_inra","CO255_inra","D105_inra","DK63","EP1_inra",
"F492","FV2_inra","FV252_inra","ND245","OH43_inra","VA85","W117_inra")
sd_melt1=melt(sd1)
names(sd_melt1)=c('stdev')
rownames(sd_melt1)=c("A632_usa","A654_inra","B73_inra","B96",
"C103_inra","CO255_inra","D105_inra","DK63","EP1_inra",
"F492","FV2_inra","FV252_inra","ND245","OH43_inra","VA85","W117_inra")
melt1=cbind(melt1,sd_melt1$stdev)
melt1$mite=c(T,T,F,F,T,T,T,T,T,T,T,T,T,F,F,T)
melt1=melt1[order(melt1$value),]
melt1=as.data.frame(melt1,stringsAsFactors=F)
melt1$founder_f=factor(rownames(melt1),levels=rownames(melt1))

rownames(melt1)=seq(1,dim(melt1)[1])
names(melt1)=c('value','stdev','mite','founder_f')


png('ZCN8_WLS.png',width=800,heigh=650)
print(ggplot(melt1,aes(x=founder_f,y=value,color=mite))+geom_errorbar(aes(ymax=value+stdev,ymin=value-stdev),position="dodge") + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZCN8 Expression"))
dev.off()


# ZmRap2.7

tmp3=foldchange[rownames(foldchange)==rap27,]
sd3=sd[rownames(sd)==rap27,]
melt3=melt(tmp3)
rownames(melt3)=c("A632_usa","A654_inra","B73_inra","B96",
"C103_inra","CO255_inra","D105_inra","DK63","EP1_inra",
"F492","FV2_inra","FV252_inra","ND245","OH43_inra","VA85","W117_inra")
sd_melt3=melt(sd3)
names(sd_melt3)=c('stdev')
rownames(sd_melt3)=c("A632_usa","A654_inra","B73_inra","B96",
"C103_inra","CO255_inra","D105_inra","DK63","EP1_inra",
"F492","FV2_inra","FV252_inra","ND245","OH43_inra","VA85","W117_inra")
melt3=cbind(melt3,sd_melt3$stdev)
melt3$mite=c(T,T,F,F,T,T,T,T,T,T,T,T,T,F,F,T)
melt3=melt3[order(melt3$value),]
melt3=as.data.frame(melt3,stringsAsFactors=F)
melt3$founder_f=factor(rownames(melt3),levels=rownames(melt3))

rownames(melt3)=seq(1,dim(melt3)[1])
names(melt3)=c('value','stdev','mite','founder_f')


png('ZmRap2.7_WLS.png',width=800,heigh=650)
print(ggplot(melt3,aes(x=founder_f,y=value,color=mite))+geom_errorbar(aes(ymax=value+stdev,ymin=value-stdev),position="dodge") + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZmRap2.7 Expression"))
dev.off()

# ZmMADS69

tmp2=foldchange[rownames(foldchange)==mads69,]
sd2=sd[rownames(sd)==mads69,]
melt2=melt(tmp2)
rownames(melt2)=c("A632_usa","A654_inra","B73_inra","B96",
"C103_inra","CO255_inra","D105_inra","DK63","EP1_inra",
"F492","FV2_inra","FV252_inra","ND245","OH43_inra","VA85","W117_inra")
sd_melt2=melt(sd2)
names(sd_melt2)=c('stdev')
rownames(sd_melt2)=c("A632_usa","A654_inra","B73_inra","B96",
"C103_inra","CO255_inra","D105_inra","DK63","EP1_inra",
"F492","FV2_inra","FV252_inra","ND245","OH43_inra","VA85","W117_inra")
melt2=cbind(melt2,sd_melt2$stdev)
melt2$mite=c(T,T,F,F,T,T,T,T,T,T,T,T,T,F,F,T)
melt2=melt2[order(melt2$value),]
melt2=as.data.frame(melt2,stringsAsFactors=F)
melt2$founder_f=factor(rownames(melt2),levels=rownames(melt2))

rownames(melt2)=seq(1,dim(melt2)[1])
names(melt2)=c('value','stdev','mite','founder_f')


png('ZmMADS69_WLS.png',width=800,heigh=650)
print(ggplot(melt2,aes(x=founder_f,y=value,color=mite))+geom_errorbar(aes(ymax=value+stdev,ymin=value-stdev),position="dodge") + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZmMADS69 Expression"))
dev.off()





## MITE vs. NO_MITE

foldchange=fread('batch_1_normalized_coefficients_MITE.txt')
zcn8="Zm00001d010752"
rap27="Zm00001d010987"
mads69="Zm00001d042315"

#ZCN8

tmp2=foldchange[foldchange$V1==zcn8,]
melt2=melt(tmp2)
#melt2$mite=c(T,T,F,F,T,T,T,T,T,T,T,T,T,F,F,T)
#melt2=melt2[order(melt2$value),]
#melt2=as.data.frame(melt2,stringsAsFactors=F)
#rownames(melt2)=seq(1,dim(melt2)[1])
#melt2$founder_f=factor(melt2$variable,levels=melt2$variable)

png('ZCN8_MITE_logFC.png',width=800,heigh=650)
print(ggplot(melt2,aes(x=variable,y=value))+geom_boxplot() + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZCN8 Expression in Lines With and Without MITE at vgt1"))
dev.off()


# ZmRap2.7

tmp3=foldchange[foldchange$V1==rap27,]
melt3=melt(tmp3)
#melt3$mite=c(T,T,F,F,T,T,T,T,T,T,T,T,T,F,F,T)
#melt3=melt3[order(melt3$value),]
#melt3=as.data.frame(melt3,stringsAsFactors=F)
#rownames(melt3)=seq(1,dim(melt3)[1])
#melt3$founder_f=factor(melt3$variable,levels=melt3$variable)

png('ZmRap2.7_MITE_logFC.png',width=800,heigh=650)
print(ggplot(melt3,aes(x=variable,y=value))+geom_boxplot() + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZmRap2.7 Expression in Lines With and Without MITE at vgt1"))
dev.off()

# ZmMADS69

tmp1=foldchange[foldchange$V1==mads69,]
melt1=melt(tmp1)
#melt1$mite=c(T,T,F,F,T,T,T,T,T,T,T,T,T,F,F,T)
#melt1=melt1[order(melt1$value),]
#melt1=as.data.frame(melt1,stringsAsFactors=F)
#rownames(melt1)=seq(1,dim(melt1)[1])
#melt1$founder_f=factor(melt1$variable,levels=melt1$variable)

png('ZmMADS69_MITE_logFC.png',width=800,heigh=650)
print(ggplot(melt1,aes(x=variable,y=value))+geom_boxplot() + xlab("Founder") + ylab("Expression Weighted Least Squares") + ggtitle("ZmMADS69 Expression in Lines With and Without MITE at vgt1"))
dev.off()


degs=fread('NO_MITE_v_MITE_DEGs.txt')
