#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
rep=as.numeric(args[[1]])
cores=as.numeric(args[[2]])

library('data.table')
library('dplyr')
library('ggplot2')
library('lme4')
library('lmerTest')
library('parallel')
library('MASS')
# Z scores of expression
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10

beta_merge=fread('eqtl/results/local_eqtl_effect_size_frequency_beta_z_all_genes.txt',data.table=F)
beta_merge=beta_merge[!is.na(beta_merge$beta_z),]
# Extreme bins are >=3/<=-3 SD
bbreaks=c(min(beta_merge$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_merge$beta_z)+0.1)
beta_merge$bin=cut(beta_merge$beta_z, breaks=bbreaks, label=F)
#beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
#beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

avg_z=beta_merge %>% group_by(bin) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)

ft_df=fread('eqtl/data/FT_genelist.txt',data.table=F)
high=beta_merge[abs(beta_merge$cor)>0.2,]
ft_corgenes=unique(high$Gene_ID)

# Is this more than we'd expect - for a random set of genes,
# Is there a different ratio of quadratic to linear? (More quadratic )

#ft_corgenes=intersect(ft_df$Gene_ID,ft_corgenes)
# 101 genes
#vhigh=high[high$Gene_ID %in% ft_corgenes,]



#vhigh=vhigh[vhigh$time.x=="WD_0712",]
vhigh=beta_merge[beta_merge$Gene_ID %in% ft_corgenes,]
# 593 genes
gene_times=unique(vhigh$gene_time)

rand=beta_merge[!(beta_merge$gene_time %in% gene_times),]

n=length(gene_times)
set.seed(rep)
n_reps=1:100
draw=sample(rand$gene_time,n)
#rep

rep_comp_df=c()
for(gt in draw){
	sub1=rand[rand$gene_time==gt,]
	gene=unique(sub1$Gene_ID)
	time1=unique(sub1$time.x)
	m0=lm(freq~1,sub1)
	m1=lm(freq~beta_z,sub1)
	# Linear model
	m1sum=summary(m1)
	m1_pval= m1sum$coefficients['beta_z',4]
	m1_slope= m1sum$coefficients['beta_z',1]
	# Is slope in same direction as correlation?
	bv_r=unique(sub1$cor)
	if(bv_r>0 & m1_slope>0){
		consistent=TRUE
	}else{
		consistent=FALSE
	}
	# Quadratic model
	m2=lm(freq~poly(beta_z,2),sub1)
	m2sum=summary(m2)
	m2_pval= m2sum$coefficients[3,4]
	m2_curve= m2sum$coefficients[3,1]
	#Model comparison
	modelcomp=anova(m0,m1,m2)
	modelcomp_res1=modelcomp$`Pr(>F)`[2]
	modelcomp_res2=modelcomp$`Pr(>F)`[3]
	if(modelcomp_res2<=0.05){
		if(m2_curve<0){
			result="p_quadratic"
		}else{
			result="n_quadratic"
		}
	}else if(modelcomp_res1<=0.05){
		result="linear"
	}else{
		result="neither"
	}
	line=data.frame(rep=rep,gene_time=gt,Gene_ID=gene,time=time1,m1_pvalue=m1_pval,m1_slope=m1_slope,bv_cor=bv_r,consistent=consistent,m2_pvalue=m2_pval,m2_curve=m2_curve,modelcomp_m1_pvalue=modelcomp_res1,modelcomp_m2_pvalue=modelcomp_res2,better=result)
	rep_comp_df=rbind(rep_comp_df,line)
}
	
rep_comp_df=as.data.frame(rep_comp_df,stringsAsFactors=F)
# 2,131


fwrite(rep_comp_df,sprintf('eqtl/results/div_dir_model_comparison_perm%s.txt',rep),row.names=F,quote=F,sep='\t')
