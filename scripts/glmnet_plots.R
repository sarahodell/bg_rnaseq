#!/usr/bin/env Rscript

library('ggplot2')
library('data.table')
library('reshape2')
library('cowplot')
library('grid')

summary=fread('QTT/glmnet_summary_stats.txt',data.table=F)

#summary[summary$]


theme_set(theme_classic())
#theme_update(text=element_text(family="Times"))
theme_update(plot.caption = element_text(hjust = 0))
theme_update(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
theme_update(plot.title = element_text(size=20),axis.title=element_text(size=20))
theme_update(panel.background=element_blank())
theme_update(plot.caption=element_text(size=20))

phenos=c("tkw_15","grain_yield_15")

envs=c("ALL","EXP_STPAUL_2017_WD")

# Which ones have test R2 > 0
# tkw_15 ALL T12 0.0258
# tkw_15 EXP_STPAUL_2017_WD T20 0.0360
# grain_yield_15 ALL T27 0.0467


times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(env in envs){
	for(pheno in phenos){
		pred=fread(sprintf('QTT/%s_%s_zscore_fitted_values2.txt',pheno,env),data.table=F)
		predmelt=reshape2::melt(pred,c("ID","pheno"))
		predmelt=predmelt[complete.cases(predmelt),]
		plotlist=list()
		count=1
		for(time1 in times){
			submelt=predmelt[predmelt$variable==time1,]
			rs=cor(submelt$pheno,submelt$value)
			 
			grob <- grobTree(textGrob(sprintf("r=%.2f",rs), x=0.80,  y=0.90, hjust=-0.1,
  			gp=gpar(col="black", fontsize=13)))
			# Plot
			p1=ggplot(aes(x=value,y=pheno),data=submelt) + geom_point() + 
		 	geom_smooth(method="lm",formula=y~x,color='darkred') +
		 	#annotate("text", x = 70, y=100, label = sprintf("r=%.2f",rs),size=6) + 
			xlab("Fitted Values") + ylab("True Values") +
			annotation_custom(grob)
			#facet_grid(cols=vars(variable)) +
			#ggtitle(sprintf("%s",time))
		
			plotlist[[count]]=p1
			count=count+1
		}
		
		prow=plot_grid(plotlist=plotlist,nrow=1,ncol=4,labels=c("T12","T18","T20","T27"))
		png(sprintf("QTT/images/%s_%s_QTT_Z_score_RR_pred.png",pheno,env),width=1600,height=600)
		print(prow)
		dev.off()
		
	}
}

# Total Exp
# Which ones have test R2 > 0
# tkw_15 EXP_STPAUL_2017_WD T20 0.00100
# grain_yield_15 ALL T27 0.00269
# grain_yield_15 EXP_STPAUL_2017_WD T20 0.0403
# grain_yield_15 EXP_STPAUL_2017_WD T27 0.0329


for(env in envs){
	for(pheno in phenos){
		pred=fread(sprintf('QTT/%s_%s_totexp_fitted_values2.txt',pheno,env),data.table=F)
		predmelt=reshape2::melt(pred,c("ID","pheno"))
		predmelt=predmelt[complete.cases(predmelt),]
		plotlist=list()
		count=1
		for(time1 in times){
			submelt=predmelt[predmelt$variable==time1,]
			rs=cor(submelt$pheno,submelt$value)
			 
			grob <- grobTree(textGrob(sprintf("r=%.2f",rs), x=0.80,  y=0.90, hjust=-0.1,
  			gp=gpar(col="black", fontsize=13)))
			# Plot
			p1=ggplot(aes(x=value,y=pheno),data=submelt) + geom_point() + 
		 	geom_smooth(method="lm",formula=y~x,color='darkred') +
		 	#annotate("text", x = 70, y=100, label = sprintf("r=%.2f",rs),size=6) + 
			xlab("Fitted Values") + ylab("True Values") +
			annotation_custom(grob)
			#facet_grid(cols=vars(variable)) +
			#ggtitle(sprintf("%s",time))
		
			plotlist[[count]]=p1
			count=count+1
		}
		
		prow=plot_grid(plotlist=plotlist,nrow=1,ncol=4,labels=c("T12","T18","T20","T27"))
		png(sprintf("QTT/images/%s_%s_QTT_totexp_RR_pred.png",pheno,env),width=1600,height=600)
		print(prow)
		dev.off()
		
	}
}