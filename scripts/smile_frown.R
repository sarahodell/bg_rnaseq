#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')

# Disregulation smile plot

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10

totalrares=c()
for(time1 in times){
	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
	allrares=allrares[allrares$max_f!="B73_inra",]
	totalrares=rbind(totalrares,allrares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
totalrares$gene_time_id=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$ID)

# first I need to break up by gene_time_founder
totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

#totalrares=totalrares[,-1]
#totalrares=rbind(totalrares,totallow)
#fwrite(totalrares,'eqtl/results/all_rare_counts_max_f_all_exp.txt',row.names=F,quote=F,sep='\t')

totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]
totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))

beta_z=totf%>% group_by(gene_time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))

#beta_z=totalrares%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
# Extreme bins are >=3/<=-3 SD
bbreaks=c(min(beta_z$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_z$beta_z)+0.1)
beta_z$bin=cut(beta_z$beta_z, breaks=bbreaks, label=F)
beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)
#avg_z=avg_z[!is.na(bin),]

fwrite(avg_z,'paper_figures/data/avg_rares_5kb_top5k.txt',row.names=F,quote=F,sep='\t')

p3=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/all_beta_z_bin_disregulation_smooth_smile.png')
print(p3)
dev.off()

p2=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(method = "lm",mapping=aes(weight=ngenes),formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))


png('eqtl/images/all_beta_z_bin_disregulation_smile.png')
print(p2)
dev.off()

quadratic_high=lm(avg_rares~poly(bin,2),avg_z)
anova(quadratic_high)
summary(quadratic_high)

##### Founder beta rank instead of z-score

# scale beta_rank by max rank

scale_rank=totf%>% group_by(gene_time) %>% mutate(rank_scaled=(beta_rank/max(beta_rank,na.rm=T)))
scale_rank=as.data.frame(scale_rank,stringsAsFactors=F)
scale_rank=scale_rank[!is.na(scale_rank$rank_scaled),]

bbreaks=seq(0.0620,1,length.out=17)
scale_rank$bin=cut(scale_rank$rank_scaled, breaks=bbreaks, label=F)
scale_rank$gene_founder=paste0(scale_rank$Gene_ID,'-',scale_rank$max_f)
avg_z=scale_rank %>% group_by(bin) %>% reframe(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)

fwrite(avg_z,'paper_figures/data/avg_rares_beta_rank_5kb_top5k.txt',row.names=F,quote=F,sep='\t')

p2=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(method = "lm",mapping=aes(weight=ngenes),formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))


png('eqtl/images/all_beta_z_bin_disregulation_smile.png')
print(p2)
dev.off()

##### Frown plot

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

fwrite(avg_z,'paper_figures/data/avg_founder_freq_frown.txt',row.names=F,quote=F,sep='\t')

quadratic_high=lm(avg_rares~poly(bin,2),avg_z)
anova(quadratic_high)
summary(quadratic_high)

# m2 curve 1.41938 pvalue=9.91e-09


p3=ggplot(aes(x=bin,y=avg_freq),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/all_gene_allele_frequency_bin_smooth_smile.png')
print(p3)
dev.off()



### Try splitting FT genes and non FT genes
rankscomp=fread('paper_figures/local_eQTL_QTL_ranks.txt',data.table=F)
ftcomp=rankscomp[rankscomp$phenotype %in% c("male_flowering_d6"),]
ftcomp=ftcomp[ftcomp$abs_r>=0.3,]
length(unique(ftcomp$Trait))

ftblup=ftcomp[ftcomp$env=="ALL",]

ft_df=fread('eqtl/data/FT_genelist.txt',data.table=F)

#beta_merge$ft=(beta_merge$Gene_ID %in% ft_df$Gene_ID | abs(beta_merge$cor)>0.2)
beta_merge$ft=(beta_merge$Gene_ID %in% ftblup$Trait)

# 101 genes that are in the gene list and correlated with FT

bbreaks=c(min(beta_merge$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_merge$beta_z)+0.1)
#beta_merge$bin=cut(beta_merge$beta_z, breaks=bbreaks, label=F)
#beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
#beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

ft_avg_z=beta_merge %>% group_by(bin,ft) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
ft_avg_z=as.data.frame(ft_avg_z,stringsAsFactors=F)
nbins=max(ft_avg_z$bin)

p3=ggplot(aes(x=bin,y=avg_freq),data=ft_avg_z) + geom_point(aes(size=ngenes,group=ft,color=ft)) + stat_smooth(mapping=aes(weight=ngenes,color=ft),linewidth = 1,method="lm",formula = y ~ x + I(x^2)) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45)) +
geom_hline(yintercept=0.0625,linetype='dotted',color='darkgrey')

png('eqtl/images/BLUP_high_cor_FT_nonFT_gene_allele_frequency.png')
print(p3)
dev.off()


theme_set(theme_classic())
#theme_update(text=element_text(family="Times"))
theme_update(plot.caption = element_text(hjust = 0))
theme_update(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
theme_update(plot.title = element_text(size=20),axis.title=element_text(size=20))
theme_update(panel.background=element_blank())
theme_update(plot.caption=element_text(size=20))
theme_update(legend.text=element_text(size=20))
#beta_merge$ft=(beta_merge$Gene_ID %in% ft_df$Gene_ID | abs(beta_merge$cor)>0.2)
beta_merge$ft=(beta_merge$Gene_ID %in% ftcomp$Trait)

# 101 genes that are in the gene list and correlated with FT

bbreaks=c(min(beta_merge$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_merge$beta_z)+0.1)
#beta_merge$bin=cut(beta_merge$beta_z, breaks=bbreaks, label=F)
#beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
#beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

ft_avg_z=beta_merge %>% group_by(bin,ft) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
ft_avg_z=as.data.frame(ft_avg_z,stringsAsFactors=F)
nbins=max(ft_avg_z$bin)

p3=ggplot(aes(x=bin,y=avg_freq),data=ft_avg_z) + geom_point(aes(size=ngenes,group=ft,color=ft)) + stat_smooth(mapping=aes(weight=ngenes,color=ft),linewidth = 1,method="lm",formula = y ~ x + I(x^2)) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
scale_color_hue(name="",labels=c("Non-Flowering Time Genes","Flowering Time Genes")) +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45)) + geom_hline(yintercept=0.0625,linetype='dotted',color='darkgrey')

png('eqtl/images/high_cor_FT_nonFT_gene_allele_frequency.png',width=1000,height=800)
print(p3)
dev.off()

quad_ft=lm(avg_freq ~ factor(ft) + poly(avg_rares,2),ft_avg_z)
