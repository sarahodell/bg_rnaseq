#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')
library('interactions')
library('lme4qtl')
library('lmerTest')
library('reshape2')

# Disregulation smile plot

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10

totalrares=c()
for(time1 in times){
	for(chr in chroms){
		allrares=fread(sprintf('eqtl/results/%s_%s_5kb_rare_counts_all_v3.txt',time1,chr),data.table=F)
		#allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
		allrares=allrares[allrares$max_f!="B73_inra",]
		totalrares=rbind(totalrares,allrares)
	}
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
totalrares=totalrares[!is.infinite(totalrares$beta),]
totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))

beta_z=totf%>% group_by(gene_time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))


alldf=fread('eqtl/data/founder_frequency.txt',data.table=F)
alldf$snp_founder=paste0(alldf$snp,'-',alldf$variable)

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

beta_merge=merge(beta_z,eqtl,by='gene_time')
beta_merge$snp_founder=paste0(beta_merge$X_ID,'-',beta_merge$max_f)
# 769 gene time, 233 genes
beta_merge$freq=alldf[match(beta_merge$snp_founder,alldf$snp_founder),]$value

fwrite(beta_merge,'eqtl/results/beta_merge_all_v3.txt',row.names=F,quote=F,sep='\t')

#beta_z=totalrares%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
#beta_merge=as.data.frame(beta_merge,stringsAsFactors=F)
#beta_z=beta_z[!is.na(beta_z$beta_z),]
# Extreme bins are >=3/<=-3 SD
bbreaks=c(min(beta_merge$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_merge$beta_z)+0.1)
beta_merge$bin=cut(beta_merge$beta_z, breaks=bbreaks, label=F)
#beta_merge$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
beta_merge$gene_founder=paste0(beta_merge$Gene_ID,'-',beta_merge$max_f)




avg_z=beta_merge %>% group_by(bin) %>% reframe(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)

fwrite(avg_z,'paper_figures/avg_rares_5kb_all.txt',row.names=F,quote=F,sep='\t')

quadratic_high2=lm(rare_count~beta_z + I(beta_z^2),beta_merge)
anova(quadratic_high2)
#summary(quadratic_high2)
#Response: rare_count
#                Df   Sum Sq Mean Sq F value    Pr(>F)    
#beta_z           1     2631  2631.1  150.26 < 2.2e-16 ***
#I(beta_z^2)      1    29693 29693.1 1695.72 < 2.2e-16 ***
#Residuals   888272 15554190    17.5                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


quadratic_high1=lm(rare_count~bin + I(bin^2),beta_merge)
anova(quadratic_high1)

newdata <- data.frame(bin=avg_z$bin)
#newdata$avg_rares <- predict(pheno_model2, newdata = newdata)

p <- predict(quadratic_high1, newdata = newdata, se.fit=TRUE,     
         interval="confidence")
newdata=cbind(newdata,p$fit)

fwrite(newdata,'paper_figures/bin_smile_predicted_all.txt',row.names=F,quote=F,sep='\t')
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


quadratic_high2=lm(rare_count~poly(bin,2),beta_z)
anova(quadratic_high2)
summary(quadratic_high2)

quadratic_high1=lm(rare_count~I(bin),beta_z)

quadratic_high3=lm(rare_count~poly(bin,3),beta_z)
anova(quadratic_high2,quadratic_high3)
anova(quadratic_high1,quadratic_high2)



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

scaled_high2=lm(rare_count~poly(rank_scaled,2),scale_rank)
anova(scaled_high2)
summary(scaled_high2)

scaled_high=lm(avg_rares~poly(bin,2),avg_z)
anova(scaled_high)
summary(scaled_high)


##### Frown plot

#beta_merge=fread('eqtl/results/local_eqtl_effect_size_frequency_beta_z_all_genes.txt',data.table=F)

beta_merge=fread('eqtl/results/beta_merge_all_v3.txt',data.table=F)
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

#quadratic_high=lm(avg_rares~poly(bin,2),avg_z)
#anova(quadratic_high)
#summary(quadratic_high)

quadratic_frown=lm(freq~beta_z + I(beta_z^2),beta_merge)
anova(quadratic_frown)

png('eqtl/images/frown_diagnostic_plots.png')
print(pls205_diagnostics(quadratic_frown))
dev.off()

quadratic_frown_bin=lm(freq~bin + I(bin^2),beta_merge)
anova(quadratic_frown_bin)

png('eqtl/images/frown_bin_diagnostic_plots.png')
print(pls205_diagnostics(quadratic_frown_bin))
dev.off()
#Response: freq
#                Df Sum Sq Mean Sq F value    Pr(>F)    
#beta_z           1   0.01 0.00916   19.88 8.248e-06 ***
#I(beta_z^2)      1   1.85 1.84938 4013.15 < 2.2e-16 ***
#Residuals   888272 409.34 0.00046  

#summary(quadratic_high2)

#quadratic_frown2=lmer(freq~ (1|X_ID) + beta_z + I(beta_z^2),beta_merge)
#anova(quadratic_frown2,ddf='K')


#submerge=beta_merge[beta_merge$time.x=="WD_0712",]

quadratic_frown3=lm(freq~ X_ID + bin + I(bin^2),submerge)
#anova(quadratic_frown3)

newdata <- data.frame(bin=avg_z$bin)
#newdata$avg_rares <- predict(pheno_model2, newdata = newdata)

p <- predict(quadratic_frown, newdata = newdata, se.fit=TRUE,     
         interval="confidence")
newdata=cbind(newdata,p$fit)

fwrite(newdata,'paper_figures/frown_predicted.txt',row.names=F,quote=F,sep='\t')


# m2 curve 1.41938 pvalue=9.91e-09


p3=ggplot(aes(x=bin,y=avg_freq),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/all_gene_allele_frequency_bin_smooth_smile.png')
print(p3)
dev.off()

beta_merge=fread('eqtl/results/beta_merge_all_v3.txt',data.table=F)
#beta_merge=fread('eqtl/results/local_eqtl_effect_size_frequency_beta_z_all_genes.txt',data.table=F)
rankscomp=fread('paper_figures/local_eQTL_QTL_ranks.txt',data.table=F)
rankscomp=rankscomp[rankscomp$abs_r>=0.3,]

#ftcomp=rankscomp[rankscomp$phenotype %in% c("male_flowering_d6"),]
#ftcomp=ftcomp[ftcomp$abs_r>=0.3,]
#beta_merge$ft=(beta_merge$Gene_ID %in% ftcomp$Trait)
# If at any timepoint a gene overlaps a QTL and is correlated with a trait, get the phenotype it is correlated with
beta_merge$phenotype="None"
gts=unique(rankscomp$Trait)
gts=intersect(gts,beta_merge$Trait)
for(i in gts){
  sub1=rankscomp[rankscomp$Trait==i,]
  if(nrow(sub1)>1){
    pheno=sub1[which.max(sub1$abs_r),]$phenotype
  }else{
    pheno=sub1$phenotype
  }
  beta_merge[beta_merge$Trait==i,]$phenotype=pheno
}
fwrite(beta_merge,'eqtl/results/beta_merge_all_v3.txt',row.names=F,quote=F,sep='\t')

#png('eqtl/images/quad_pheno_interaction_plot.png')
#print(interactions::interact_plot(pheno_model, pred = beta_z, modx = phenotype, data = beta_merge))
#dev.off(),row.names=F,quote=F,sep='\t')
### Try splitting FT genes and non FT genes

pheno_model2=lm(freq ~ phenotype*beta_z + phenotype*I(beta_z^2),beta_merge)
anova(pheno_model2)
#Analysis of Variance Table

#Response: freq
#                          Df Sum Sq Mean Sq   F value    Pr(>F)    
#phenotype                  5   0.41 0.08165  177.3827 < 2.2e-16 ***
#beta_z                     1   0.01 0.00916   19.9013 8.155e-06 ***
#I(beta_z^2)                1   1.85 1.84879 4016.2240 < 2.2e-16 ***
#phenotype:beta_z           5   0.01 0.00143    3.1141  0.008185 ** 
#phenotype:I(beta_z^2)      5   0.03 0.00699   15.1786 6.065e-15 ***
#Residuals             888257 408.89 0.00046                        
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#summary(pheno_model)



png('eqtl/images/quad_pheno_interaction_plot.png')
print(interactions::interact_plot(pheno_model2, pred = beta_z, modx = phenotype, data = beta_merge))
dev.off()


beta_merge$pheno_f=as.factor(beta_merge$phenotype)
beta_merge$abs_beta_z=abs(beta_merge$beta_z)

beta_z=beta_merge %>% group_by(gene_time) %>% mutate(abs_beta_z=abs((beta-mean(beta,na.rm=T))/sd(beta,na.rm=T)))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
# Extreme bins are >=3/<=-3 SD
bbreaks=c(0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,max(beta_z$abs_beta_z)+0.1)
beta_z$abs_bin=cut(beta_z$abs_beta_z, breaks=bbreaks, label=F)
beta_merge$abs_bin=beta_z[match(beta_merge$gene_time_founder,beta_z$gene_time_founder),]$abs_bin

avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)



pheno_model1=lm(freq ~ pheno_f + abs_bin,beta_merge)
anova(pheno_model1)
summary(pheno_model1)

library(lsmeans)

pheno_model2=lm(freq ~ abs_bin*pheno_f,beta_merge)
anova(pheno_model2)
#Analysis of Variance Table
#
#Response: freq
#                    Df Sum Sq Mean Sq  F value    Pr(>F)    
#abs_bin              1   1.93 1.93379 4191.397 < 2.2e-16 ***
#pheno_f              5   0.42 0.08438  182.887 < 2.2e-16 ***
#abs_bin:pheno_f      5   0.04 0.00709   15.369 3.837e-15 ***
#Residuals       920257 424.58 0.00046                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(pheno_model2)

pheno_model2$coefficients
m.lst <- lstrends(pheno_model2, "pheno_f", var="abs_bin")
pairs(m.lst)
#contrast                                          estimate           SE     df
# female_flowering_d6 - harvest_grain_moisture -0.0004464901 1.156657e-04 920257 **
# female_flowering_d6 - male_flowering_d6      -0.0004095976 9.994880e-05 920257 **
# female_flowering_d6 - None                   -0.0004808994 9.098437e-05 920257 **
# female_flowering_d6 - tkw_15                 -0.0001363689 1.136680e-04 920257
# female_flowering_d6 - total_plant_height     -0.0001151348 1.148497e-04 920257
# harvest_grain_moisture - male_flowering_d6    0.0000368925 8.387308e-05 920257
# harvest_grain_moisture - None                -0.0000344094 7.295948e-05 920257
# harvest_grain_moisture - tkw_15               0.0003101212 9.982655e-05 920257 **
# harvest_grain_moisture - total_plant_height   0.0003313553 1.011701e-04 920257 **
# male_flowering_d6 - None                     -0.0000713018 4.398053e-05 920257
# male_flowering_d6 - tkw_15                    0.0002732287 8.109587e-05 920257 **
# male_flowering_d6 - total_plant_height        0.0002944628 8.274413e-05 920257 **
# None - tkw_15                                 0.0003445305 6.974907e-05 920257 **
# None - total_plant_height                     0.0003657646 7.165880e-05 920257 **
# tkw_15 - total_plant_height                   0.0000212341 9.887991e-05 920257
# t.ratio p.value
#  -3.860  0.0016
#  -4.098  0.0006
#  -5.286  <.0001
#  -1.200  0.8372
#  -1.002  0.9173
#   0.440  0.9979
#  -0.472  0.9971
#   3.107  0.0233
#   3.275  0.0135
#  -1.621  0.5844
#   3.369  0.0098
#   3.559  0.0050
#   4.940  <.0001
#   5.104  <.0001
#   0.215  0.9999

# Different from None
# DTS, TKW, TPH

# Different from DTA
#TPH, TKW, DTS

# Different from HGM
# DTS, TKW, TPH

### Intercept  - Interaction term is true, so watch out!
library("emmeans")
emm = emmeans(pheno_model2, "pheno_f", at = list(abs_bin = 14))
pairs(emm)

#contrast                                      estimate       SE     df t.ratio
# contrast                                      estimate       SE     df t.ratio
# female_flowering_d6 - harvest_grain_moisture -0.007296 0.001001 920257  -7.288 **
# female_flowering_d6 - male_flowering_d6      -0.003831 0.000866 920257  -4.425 **
# female_flowering_d6 - None                   -0.006260 0.000788 920257  -7.949 **
# female_flowering_d6 - tkw_15                 -0.001115 0.000985 920257  -1.133 
# female_flowering_d6 - total_plant_height     -0.000718 0.000995 920257  -0.722
# harvest_grain_moisture - male_flowering_d6    0.003465 0.000726 920257   4.770 **
# harvest_grain_moisture - None                 0.001036 0.000631 920257   1.641
# harvest_grain_moisture - tkw_15               0.006181 0.000865 920257   7.147 **
# harvest_grain_moisture - total_plant_height   0.006578 0.000877 920257   7.503 **
# male_flowering_d6 - None                     -0.002429 0.000382 920257  -6.362 **
# male_flowering_d6 - tkw_15                    0.002716 0.000704 920257   3.860 **
# male_flowering_d6 - total_plant_height        0.003113 0.000718 920257   4.334 **
# None - tkw_15                                 0.005145 0.000605 920257   8.505 **
# None - total_plant_height                     0.005542 0.000622 920257   8.911 **
# tkw_15 - total_plant_height                   0.000397 0.000858 920257   0.463
# p.value
#  <.0001 
#  0.0001
#  <.0001
#  0.8679
#  0.9793
#  <.0001
#  0.5715
#  <.0001
#  <.0001
#  <.0001
#  0.0016
#  0.0002
#  <.0001
#  <.0001
#  0.9974

# Different from None At bin 12
# DTS, DTA, TKW, TPH

# Different from DTA
# DTS, HGM, TKW, TPH, None

# Different from HGM


newdata <- data.frame(abs_bin=pheno_abs_z$abs_bin, pheno_f=pheno_abs_z$pheno_f)
newdata$avg_freq <- predict(pheno_model2, newdata = newdata)

p <- predict(pheno_model2, newdata = newdata, se.fit=TRUE,     
         interval="confidence")
newdata=cbind(newdata,p$fit)

fwrite(newdata,'paper_figures/pheno_model_predictions.txt',row.names=F,quote=F,sep='\t')
#crit <- qt(p= 0.05/2, df = summary(pheno_model2)$df[2], lower.tail = FALSE)
#crit

# check that computing the half-length of the interval and 
# dividing it by the critical value gives the same result 
# as that reported by se.fit
#(p$fit[,"upr"] - p$fit[,"fit"])/crit
#newdata$Quadratic <- predict(calquad.1, newdata = newdata)
#library(reshape2)
#newdata <- melt(newdata, id.vars = "abs_bin", variable.name = "pheno_f")
#ggplot(calvarbyruno.1, aes(x = PAR, y = Nominal, weight=Nominal^calweight)) + 
#    geom_line(data = newdata, aes(x = value, colour = Model)) + 
#    geom_point()


#pheno_means=e
bbreaks=c(0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,max(beta_z$abs_beta_z)+0.1)
pheno_abs_z=beta_merge %>% group_by(abs_bin,pheno_f) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
pheno_abs_z=as.data.frame(pheno_abs_z,stringsAsFactors=F)
nbins=max(pheno_abs_z$abs_bin)
fwrite(pheno_abs_z,'paper_figures/data/avg_freq_by_abs_bin_pheno.txt',row.names=F,quote=F,sep='\t')




p4=ggplot(aes(x=abs_bin,y=avg_freq),data=pheno_abs_z) + geom_point(aes(size=ngenes,group=pheno_f,color=pheno_f)) +
 #stat_smooth(mapping=aes(weight=ngenes,color=pheno_f),linewidth = 1,method="lm",formula = y ~ x) +
 geom_line(data = newdata, aes(x = abs_bin, y=avg_freq, colour = pheno_f)) +
 geom_ribbon(data=newdata,aes(ymin=newdata$lwr, ymax=newdata$upr,fill=pheno_f), linetype=2, alpha=0.4) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
#scale_color_hue(name="",labels=c("Non-Flowering Time Genes","Flowering Time Genes")) +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("(0.0..0.25]","(0.25..0.5]","(0.5..0.75]","(0.75..1.0]","(1.0..1.25]","(1.25..1.5]","(1.5..1.75]","(1.75..2.0]","(2.0..2.25]","(2.25..2.5]","(2.5..2.75]","(2.75..3.0]",">3.0)"))) +
scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45)) + geom_hline(yintercept=0.0625,linetype='dotted',color='darkgrey')

png('paper_figures/frown_beta_z_by_phenotypes.png')
print(p4)
dev.off()





bbreaks=c(min(beta_merge$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_merge$beta_z)+0.1)



fwrite(pheno_avg_z,'paper_figures/data/avg_freq_by_bin_pheno.txt',row.names=F,quote=F,sep='\t')


pheno_avg_z=beta_merge %>% group_by(bin,pheno_f) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
pheno_avg_z=as.data.frame(pheno_avg_z,stringsAsFactors=F)
nbins=max(pheno_avg_z$bin)

quadratic_pheno=lm(freq~pheno_f + bin + I(bin^2),beta_merge)
anova(quadratic_pheno)
#summary(quadratic_high2)

newdata <- data.frame(bin=pheno_avg_z$bin,pheno_f=pheno_avg_z$pheno_f)
#newdata$avg_rares <- predict(pheno_model2, newdata = newdata)

p <- predict(quadratic_pheno, newdata = newdata, se.fit=TRUE,     
         interval="confidence")
newdata=cbind(newdata,p$fit)

fwrite(newdata,'paper_figures/pheno_frown_predicted.txt',row.names=F,quote=F,sep='\t')



p4=ggplot(aes(x=bin,y=avg_freq),data=pheno_avg_z) + geom_point(aes(size=ngenes,group=pheno_f,color=pheno_f)) + stat_smooth(mapping=aes(weight=ngenes,color=pheno_f),linewidth = 1,method="lm",formula = y ~ x + I(x^2)) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
#scale_color_hue(name="",labels=c("Non-Flowering Time Genes","Flowering Time Genes")) +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45)) + geom_hline(yintercept=0.0625,linetype='dotted',color='darkgrey')

png('paper_figures/frown_quad_beta_z_by_phenotypes.png')
print(p4)
dev.off()

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



#######################################
# Homozygous rare alleles
#######################################

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10

totalrares=c()
for(time1 in times){
	for(chr in chroms){
		allrares=fread(sprintf('eqtl/results/%s_%s_5kb_homozygous_rare_counts_all_v3.txt',time1,chr),data.table=F)
		#allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
		allrares=allrares[allrares$max_f!="B73_inra",]
		totalrares=rbind(totalrares,allrares)
	}
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
totalrares$gene_time_id=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$ID)

totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]
totalrares=totalrares[!is.infinite(totalrares$beta),]
totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))

beta_z=totf%>% group_by(gene_time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))


alldf=fread('eqtl/data/founder_frequency.txt',data.table=F)
alldf$snp_founder=paste0(alldf$snp,'-',alldf$variable)

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

beta_merge=merge(beta_z,eqtl,by='gene_time')
beta_merge$snp_founder=paste0(beta_merge$X_ID,'-',beta_merge$max_f)
beta_merge$freq=alldf[match(beta_merge$snp_founder,alldf$snp_founder),]$value
beta_merge$gene_founder=paste0(beta_merge$Trait,'_',beta_merge$max_f)
fwrite(beta_merge,'eqtl/results/beta_merge_all_homozygous_v3.txt',row.names=F,quote=F,sep='\t')

beta_merge=beta_merge[!is.na(beta_merge$beta_z),]
bbreaks=c(min(beta_merge$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_merge$beta_z)+0.1)
beta_merge$bin=cut(beta_merge$beta_z, breaks=bbreaks, label=F)

avg_z=beta_merge %>% group_by(bin) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)

fwrite(avg_z,'paper_figures/data/avg_homozygous_founder_freq_frown.txt',row.names=F,quote=F,sep='\t')

#quadratic_high=lm(avg_rares~poly(bin,2),avg_z)
#anova(quadratic_high)
#summary(quadratic_high)

quadratic_frown=lm(freq~beta_z + I(beta_z^2),beta_merge)
anova(quadratic_frown)
#Response: freq
#                Df Sum Sq Mean Sq F value    Pr(>F)    
#beta_z           1   0.01 0.00916   19.88 8.248e-06 ***
#I(beta_z^2)      1   1.85 1.84938 4013.15 < 2.2e-16 ***
#Residuals   888272 409.34 0.00046                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

quadratic_frown=lm(freq~bin + I(bin^2),beta_merge)
anova(quadratic_frown)

newdata <- data.frame(bin=avg_z$bin)
#newdata$avg_rares <- predict(pheno_model2, newdata = newdata)

p <- predict(quadratic_frown, newdata = newdata, se.fit=TRUE,     
         interval="confidence")
newdata=cbind(newdata,p$fit)

fwrite(newdata,'paper_figures/homozygous_frown_predicted.txt',row.names=F,quote=F,sep='\t')


# m2 curve 1.41938 pvalue=9.91e-09


p3=ggplot(aes(x=bin,y=avg_freq),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/all_homozygous_rares_frequency_bin_smooth_smile.png')
print(p3)
dev.off()


tmp=beta_merge[,c('gene_time_founder','gene_founder','bin','beta_z','freq','rare_count')]
names(tmp)=c('gene_time_founder','gene_founder','bin','beta_z','freq','homo_rare_count')
tmp$all_rare_count=beta_merge2[match(tmp$gene_time_founder,beta_merge2$gene_time_founder),]$rare_count

tmp$het_rare_count=tmp$all_rare_count-tmp$homo_rare_count

tmpmelt=reshape2::melt(tmp,c('gene_time_founder','gene_founder','bin','beta_z','freq'))

avg_z=tmpmelt %>% group_by(bin,variable) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(value),sd_rares=sd(value),n=length(value),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)

avg_z=avg_z[avg_z$variable!="all_rare_count",]

p3=ggplot(aes(x=bin,y=avg_freq,group=variable),data=avg_z) + geom_point(aes(size=ngenes,color=variable)) + stat_smooth(mapping=aes(weight=ngenes,color=variable,fill=variable),linewidth = 1,alpha=0.4) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Founder Allele Frequency") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/all_homo_het_rares_frequency_bin_smooth_smile.png')
print(p3)
dev.off()

totalrares=c()
for(chr in chroms){
	rares=fread(sprintf('eqtl/results/chr%s_5kb_homozygous_founder_block_rare_counts_all_v3.txt',chr),data.table=F)
	totalrares=rbind(totalrares,rares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$snp_founder=paste0(totalrares$snp,'-',totalrares$founder)
alldf=fread('eqtl/data/founder_frequency.txt',data.table=F)

alldf$snp_founder=paste0(alldf$snp,'-',alldf$variable)
totalrares$freq=alldf[match(totalrares$snp_founder,alldf$snp_founder),]$value

totalrares$het_rare_count=totalrares$total_rare_count-totalrares$homo_rare_count
totalrares=totalrares[totalrares$founder!="B73_inra",]

fwrite(totalrares,'eqtl/results/total_rares_by_founder_block.txt',row.names=F,quote=F,sep='\t')

totalrares=fread('eqtl/results/total_rares_by_founder_block.txt',data.table=F)
K=fread('../GridLMM/K_matrices/founder_K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
K=K[,-1]
K=as.matrix(K)

totalrares=totalrares[-68139,]

m0 <- relmatLmer(freq~ (1|founder) + homo_rare_count,data=totalrares, relmat = list(founder = K))
m1 <- relmatLmer(freq~ (1|founder) + homo_rare_count + het_rare_count,data=totalrares, relmat = list(founder = K))


Anova(m0,m1)

model1=lm(freq~homo_rare_count,data=totalrares)
#Analysis of Variance Table
#
#Response: freq
#                   Df Sum Sq   Mean Sq F value    Pr(>F)    
#homo_rare_count     1  0.008 0.0080398  13.023 0.0003078 ***
#Residuals       75454 46.581 0.0006173                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
model2=lm(freq~homo_rare_count+het_rare_count,data=totalrares)
#Analysis of Variance Table

#Response: freq
                   Df Sum Sq   Mean Sq F value    Pr(>F)    
#homo_rare_count     1  0.008 0.0080398 13.0232 0.0003079 ***
#het_rare_count      1  0.000 0.0003030  0.4908 0.4835924    
#Residuals       75453 46.580 0.0006173                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(model1,model2)
#Model 1: freq ~ homo_rare_count
#Model 2: freq ~ homo_rare_count + het_rare_count
#  Res.Df    RSS Df  Sum of Sq      F Pr(>F)
#1  75454 46.581                            
#2  75453 46.580  1 0.00030296 0.4908 0.4836



#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     6.247e-02  9.093e-05 686.956  < 2e-16 ***
#homo_rare_count 1.482e-05  4.107e-06   3.609 0.000308 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

