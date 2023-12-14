#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('glmnet')

library('lme4')
library('lmerTest')
library('plyr')
library('readr')
library('caret')
library('repr')
library('dplyr')
library('grid')
library('cowplot')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")

# chroms=1:10
# totalrares=c()
# for(time1 in times){
# 	for(chr in chroms){
# 		allrares=fread(sprintf('eqtl/results/%s_%s_5kb_rare_counts_all_v3.txt',time1,chr),data.table=F)
# 		#allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
# 		#allrares=allrares[allrares$max_f!="B73_inra",]
# 		totalrares=rbind(totalrares,allrares)
# 	}
# }
# totalrares=as.data.frame(totalrares,stringsAsFactors=F)
# totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
# totalrares$gene_time_id=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$ID)
# 
# # first I need to break up by gene_time_founder
# totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)
# 
# #totalrares=totalrares[,-1]
# #totalrares=rbind(totalrares,totallow)
# #fwrite(totalrares,'eqtl/results/all_rare_counts_max_f_all_exp.txt',row.names=F,quote=F,sep='\t')
# 
# totalrares=totalrares[!is.na(totalrares$max_f),]
# totalrares=totalrares[totalrares$max_f!="",]
# totalrares=totalrares[!is.infinite(totalrares$beta),]
# totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))
# 
# beta_z=totf%>% group_by(gene_time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
# totalrares$beta_z=beta_z[match(totalrares$gene_time_founder,beta_z$gene_time_founder),]$beta_z
# 
# totalrares=totalrares[totalrares$max_f!="B73_inra",]
# totalrares=totalrares[!is.na(totalrares$max_f),]
# 
# fwrite(totalrares,'eqtl/results/rare_counts_max_f_v3_all.txt',row.names=F,quote=F,sep='\t')

log_inverse=function(x){return(2^x)}

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
envs=c("ALL","EXP_STPAUL_2017_WD")
phenos=c("grain_yield_15","tkw_15")
result_table=c()
totalrares=fread('eqtl/results/rare_counts_max_f_v3_all.txt',data.table=F)
totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]
totalrares=totalrares[!is.na(totalrares$beta_z),]
totalrares=totalrares[!is.infinite(totalrares$beta_z),]
for(time in times){
	subrares=totalrares[totalrares$time==time,]
	for(pheno in phenos){
		for(env in envs){
			phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
			phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
			rownames(phenotypes)=phenotypes$Genotype_code
			inds=rownames(phenotypes)
			phenotypes=phenotypes[inds,pheno,drop=F]
			phenotypes=phenotypes[!is.na(phenotypes[,pheno]),,drop=F]
			inds=rownames(phenotypes)
			# Total Expression (D_i)
			exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
			rownames(exp)=exp$V1
			exp=exp[,-1]
			inter=intersect(inds,rownames(exp))
			y=phenotypes[inter,,drop=F]
			exp=exp[inter,]
			unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
			rownames(unlog)=rownames(exp)
			colnames(unlog)=colnames(exp)
			avg_exp = apply(unlog,MARGIN=2,mean)
			names(avg_exp)=names(unlog)
			avg_exp=sort(avg_exp,decreasing=TRUE)
			unlog=unlog[inter,names(avg_exp)]
			avg_logexp=log2(avg_exp)
			dev=sapply(seq(1,length(avg_exp)),function(x) (unlog[,x]-avg_exp[x])**2)
			rownames(dev)=rownames(unlog)
			colnames(dev)=colnames(unlog)
			di=as.data.frame(apply(dev,MARGIN=1,function(x) log10(mean(x,na.rm=T))))
			df=cbind(y,di)
			names(df)=c('pheno','di_score')
			# Fit model
			di_model=lm(pheno ~ I(di_score),df)
			di_anova=anova(di_model)
			pvalue=di_anova[1,5]
			di_summary=summary(di_model)
			n=nrow(df)
			coeff=di_summary$coefficients[2,1]
			fittedvals=di_summary$coefficients[1,1] + coeff*df$di_score
			df$fitted=fittedvals
			fwrite(df,sprintf('QTT/%s_%s_%s_cumul_totexp_fitted_vals_all.txt',pheno,env,time),row.names=F,quote=F,sep='\t')
			r=cor(df$pheno,fittedvals)
			line1=data.frame(phenotype=pheno,environment=env,method="D_i",timepoint=time,pvalue=pvalue,coefficient=coeff,r=r,n=n,stringsAsFactors=F)
			result_table=rbind(result_table,line1)
			# Z-score expression
			dzi=subrares %>% group_by(ID) %>% reframe(dzi_score=mean(beta_z^2,na.rm=T))
			dzi=as.data.frame(dzi,stringsAsFactors=F)
			inter=intersect(inds,dzi$ID)
			y=phenotypes[inter,,drop=F]
			rownames(dzi)=dzi$ID
			dzi=dzi[inter,]
			df=cbind(dzi,y)
			names(df)=c('ID','dzi_score','pheno')
			# Fit model
			dzi_model=lm(pheno ~ I(dzi_score),df)
			dzi_anova=anova(dzi_model)
			pvalue=dzi_anova[1,5]
			dzi_summary=summary(dzi_model)
			n=nrow(df)
			coeff=dzi_summary$coefficients[2,1]
			fittedvals=dzi_summary$coefficients[1,1] + coeff*df$dzi_score
			df$fitted=fittedvals
			fwrite(df,sprintf('QTT/%s_%s_%s_cumul_zscore_fitted_vals_all.txt',pheno,env,time),row.names=F,quote=F,sep='\t')
			r=cor(df$pheno,fittedvals)
			line2=data.frame(phenotype=pheno,environment=env,method="D_zi",timepoint=time,pvalue=pvalue,coefficient=coeff,r=r,n=n,stringsAsFactors=F)
			result_table=rbind(result_table,line2)
		}
	}
}



result_table=as.data.frame(result_table,stringsAsFactors=F)
result_table$p_adjust=p.adjust(result_table$pvalue,method='fdr')

#result_table$p_adjusted=unlist(sapply(seq(1,nrow(result_table)),function(x) min((result_table$pvalue[x]*nrow(result_table)),1)))

fwrite(result_table,'QTT/cumul_dys_regression_stats_all.txt',row.names=F,quote=F,sep='\t')


#### Make plots

# Total Expression
result_table=fread('QTT/cumul_dys_regression_stats_all.txt',data.table=F)

for(pheno in phenos){
	for(env in envs){
		plotlist=list()
		count=1
		for(time1 in times){
			df=fread(sprintf('QTT/%s_%s_%s_cumul_totexp_fitted_vals.txt',pheno,env,time1),data.table=F)
			r=cor(df$pheno,df$fitted)
			di_model=lm(pheno ~ I(di_score),df)
			di_anova=anova(di_model)
			pvalue=di_anova[1,5]
			grob <- grobTree(textGrob(sprintf("r=%.2f, p-value=%.2f",r,pvalue), x=0.65,  y=0.95,
			gp=gpar(col="black", fontsize=13)))
			p1=ggplot(df,aes(x=fitted,y=pheno)) + geom_point() + geom_smooth(method="lm",formula=y~x,color='black') + 
			xlab(paste0(expression(D[i])," Fitted Values")) + ylab("True Values") + 
			annotation_custom(grob)
			plotlist[[count]]=p1
			count=count+1
		}
		prow=plot_grid(plotlist=plotlist,nrow=1,ncol=4,labels=c("T12","T18","T20","T27"))
		png(sprintf("QTT/images/%s_%s_cumul_totexp_pred_all.png",pheno,env),width=1600,height=600)
		print(prow)
		dev.off()
	}
}

for(pheno in phenos){
	for(env in envs){
		plotlist=list()
		count=1
		for(time1 in times){
			df=fread(sprintf('QTT/%s_%s_%s_cumul_totexp_fitted_vals.txt',pheno,env,time1),data.table=F)
			r=cor(df$pheno,df$fitted)
			di_model=lm(pheno ~ I(di_score),df)
			di_anova=anova(di_model)
			pvalue=di_anova[1,5]
			di_summary=summary(di_model)
			n=nrow(df)
			coeff=di_summary$coefficients[2,1]
			if(coeff<0){
				r=r*(-1)
			}
			grob <- grobTree(textGrob(sprintf("r=%.2f, p-value=%.2f",r,pvalue), x=0.65,  y=0.95,
			gp=gpar(col="black", fontsize=13)))
			p1=ggplot(df,aes(x=di_score,y=pheno)) + geom_point() + geom_smooth(method="lm",formula=y~x,color='black') + 
			xlab(paste0(expression(D[i])," Score")) + ylab("True Values") + 
			annotation_custom(grob)
			plotlist[[count]]=p1
			count=count+1
		}
		prow=plot_grid(plotlist=plotlist,nrow=1,ncol=4,labels=c("T12","T18","T20","T27"))
		png(sprintf("QTT/images/%s_%s_cumul_totexp_lm_all.png",pheno,env),width=1600,height=600)
		print(prow)
		dev.off()
	}
}

# Z-scores
for(pheno in phenos){
	for(env in envs){
		plotlist=list()
		count=1
		for(time1 in times){
			df=fread(sprintf('QTT/%s_%s_%s_cumul_zscore_fitted_vals_all.txt',pheno,env,time1),data.table=F)
			r=cor(df$pheno,df$fitted)
			dzi_model=lm(pheno ~ I(dzi_score),df)
			dzi_anova=anova(dzi_model)
			pvalue=dzi_anova[1,5]
			grob <- grobTree(textGrob(sprintf("r=%.2f, p-value=%.2f",r,pvalue), x=0.65,  y=0.95,
			gp=gpar(col="black", fontsize=13)))
			p1=ggplot(df,aes(x=fitted,y=pheno)) + geom_point() + geom_smooth(method="lm",formula=y~x,color='black') + 
			xlab(paste0(expression(D[zi])," Fitted Values")) + ylab("True Values") + 
			annotation_custom(grob)
			plotlist[[count]]=p1
			count=count+1
		}
		prow=plot_grid(plotlist=plotlist,nrow=1,ncol=4,labels=c("T12","T18","T20","T27"))
		png(sprintf("QTT/images/%s_%s_cumul_zscore_pred_all.png",pheno,env),width=1600,height=600)
		print(prow)
		dev.off()
	}
}

for(pheno in phenos){
	for(env in envs){
		plotlist=list()
		count=1
		for(time1 in times){
			df=fread(sprintf('QTT/%s_%s_%s_cumul_zscore_fitted_vals_all.txt',pheno,env,time1),data.table=F)
			r=cor(df$pheno,df$fitted)
			dzi_model=lm(pheno ~ I(dzi_score),df)
			dzi_anova=anova(dzi_model)
			pvalue=dzi_anova[1,5]
			dzi_summary=summary(dzi_model)
			n=nrow(df)
			coeff=dzi_summary$coefficients[2,1]
			if(coeff<0){
				r=r*(-1)
			}
			grob <- grobTree(textGrob(sprintf("r=%.2f, p-value=%.2f",r,pvalue), x=0.65,  y=0.95,
			gp=gpar(col="black", fontsize=13)))
			p1=ggplot(df,aes(x=dzi_score,y=pheno)) + geom_point() + geom_smooth(method="lm",formula=y~x,color='black') + 
			xlab(paste0(expression(D[zi])," Score")) + ylab("True Values") + 
			annotation_custom(grob)
			plotlist[[count]]=p1
			count=count+1
		}
		prow=plot_grid(plotlist=plotlist,nrow=1,ncol=4,labels=c("T12","T18","T20","T27"))
		png(sprintf("QTT/images/%s_%s_cumul_zscore_lm_all.png",pheno,env),width=1600,height=600)
		print(prow)
		dev.off()
	}
}



########### Using sums instead of average
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
envs=c("ALL","EXP_STPAUL_2017_WD")
phenos=c("grain_yield_15","tkw_15")
result_table=c()
totalrares=fread('eqtl/results/beta_merge_all_v3.txt',data.table=F)
totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]
totalrares=totalrares[!is.na(totalrares$beta_z),]
totalrares=totalrares[!is.infinite(totalrares$beta_z),]
for(time in times){
	#totalrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f_all.txt',time),data.table=F)
	#totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
	#totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)
	#totalrares=totalrares[!is.na(totalrares$max_f),]
	#totalrares=totalrares[totalrares$max_f!="",]
	#totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),gene_time=unique(gene_time),max_f=unique(max_f))
	#beta_z=totf%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
	#beta_z=as.data.frame(beta_z,stringsAsFactors=F)
	#beta_z=beta_z[!is.na(beta_z$beta_z),]
	#totalrares$beta_z=beta_z[match(totalrares$gene_time_founder,beta_z$gene_time_founder),]$beta_z
	for(pheno in phenos){
		for(env in envs){
			phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
			phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
			rownames(phenotypes)=phenotypes$Genotype_code
			inds=rownames(phenotypes)
			phenotypes=phenotypes[inds,pheno,drop=F]
			phenotypes=phenotypes[!is.na(phenotypes[,pheno]),,drop=F]
			inds=rownames(phenotypes)
			# Total Expression (D_i)
			exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
			rownames(exp)=exp$V1
			exp=exp[,-1]
			inter=intersect(inds,rownames(exp))
			y=phenotypes[inter,,drop=F]
			exp=exp[inter,]
			unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
			rownames(unlog)=rownames(exp)
			colnames(unlog)=colnames(exp)
			avg_exp = apply(unlog,MARGIN=2,mean)
			names(avg_exp)=names(unlog)
			avg_exp=sort(avg_exp,decreasing=TRUE)
			avg_exp=avg_exp[1:5000] # grap top 5000 highest expressed genes
			unlog=unlog[inter,names(avg_exp)]
			avg_logexp=log2(avg_exp)
			dev=sapply(seq(1,length(avg_exp)),function(x) (unlog[,x]-avg_exp[x])**2)
			rownames(dev)=rownames(unlog)
			colnames(dev)=colnames(unlog)
			di=as.data.frame(apply(dev,MARGIN=1,function(x) log10(sum(x,na.rm=T))))
			df=cbind(y,di)
			names(df)=c('pheno','di_score')
			# Fit model
			di_model=lm(pheno ~ I(di_score),df)
			di_anova=anova(di_model)
			pvalue=di_anova[1,5]
			di_summary=summary(di_model)
			n=nrow(df)
			coeff=di_summary$coefficients[2,1]
			fittedvals=di_summary$coefficients[1,1] + coeff*df$di_score
			df$fitted=fittedvals
			fwrite(df,sprintf('QTT/%s_%s_%s_cumul_sum_totexp_fitted_vals_all.txt',pheno,env,time),row.names=F,quote=F,sep='\t')
			r=cor(df$pheno,fittedvals)
			line1=data.frame(phenotype=pheno,environment=env,method="D_i",timepoint=time,pvalue=pvalue,coefficient=coeff,r=r,n=n,stringsAsFactors=F)
			result_table=rbind(result_table,line1)
			# Z-score expression
			dzi=totalrares %>% group_by(ID) %>% reframe(dzi_score=log10(sum(beta_z^2,na.rm=T)))
			dzi=as.data.frame(dzi,stringsAsFactors=F)
			inter=intersect(inds,dzi$ID)
			y=phenotypes[inter,,drop=F]
			rownames(dzi)=dzi$ID
			dzi=dzi[inter,]
			df=cbind(dzi,y)
			names(df)=c('ID','dzi_score','pheno')
			# Fit model
			dzi_model=lm(pheno ~ I(dzi_score),df)
			dzi_anova=anova(dzi_model)
			pvalue=dzi_anova[1,5]
			dzi_summary=summary(dzi_model)
			n=nrow(df)
			coeff=dzi_summary$coefficients[2,1]
			fittedvals=dzi_summary$coefficients[1,1] + coeff*df$dzi_score
			df$fitted=fittedvals
			fwrite(df,sprintf('QTT/%s_%s_%s_cumul_sum_zscore_fitted_vals_all.txt',pheno,env,time),row.names=F,quote=F,sep='\t')
			r=cor(df$pheno,fittedvals)
			line2=data.frame(phenotype=pheno,environment=env,method="D_zi",timepoint=time,pvalue=pvalue,coefficient=coeff,r=r,n=n,stringsAsFactors=F)
			result_table=rbind(result_table,line2)
		}
	}
}



result_table=as.data.frame(result_table,stringsAsFactors=F)
result_table$p_adjust=p.adjust(result_table$pvalue,method='fdr')

#result_table$p_adjusted=unlist(sapply(seq(1,nrow(result_table)),function(x) min((result_table$pvalue[x]*nrow(result_table)),1)))

fwrite(result_table,'QTT/cumul_sum_dys_regression_stats_all.txt',row.names=F,quote=F,sep='\t')
