#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')
library('emmeans')
#library('lsmeans')


##### Frown plot
beta_merge=fread('eqtl/results/beta_merge_all_v3.txt',data.table=F)
#beta_merge=fread('paper_figures/data/local_eqtl_effect_size_frequency_beta_z_all_genes.txt',data.table=F)
#beta_merge=beta_merge[!is.na(beta_merge$beta_z),]


#pheno_model=lm(freq ~ factor(phenotype) + beta_z + I(beta_z^2),beta_merge)
#anova(pheno_model)
#summary(pheno_model)
#m.lst <- lstrends(pheno_model, "phenotype", var="I(beta_z^2)")

#pheno_model1=lm(freq ~ factor(phenotype) + bin + I(bin^2),beta_merge)
#anova(pheno_model1)
#summary(pheno_model1)

####### By BIN######
pheno_model1=lm(freq ~ phenotype*bin + phenotype*I(bin^2) ,beta_merge)

bbreaks=c(min(beta_merge$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_merge$beta_z)+0.1)
#beta_merge$abs_bin=cut(beta_merge$abs_beta_z, breaks=bbreaks, label=F)
#beta_merge$abs_bin=beta_z[match(beta_merge$gene_time_founder,beta_z$gene_time_founder),]$abs_bin

avg_z_bin=beta_merge %>% group_by(bin,phenotype) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),n=length(freq),ngenes=length(unique(gene_founder)))
avg_z_bin=as.data.frame(avg_z_bin,stringsAsFactors=F)
nbins=max(avg_z_bin$bin)

phenos=unique(beta_merge$phenotype)
newdata <- avg_z_bin[,c('bin','phenotype')]
newdata$avg_freq <- predict(pheno_model1, newdata = newdata)

p <- predict(pheno_model1, newdata = newdata, se.fit=TRUE,     
         interval="confidence")
newdata=cbind(newdata,p$fit)

fwrite(avg_z_bin,'paper_figures/pheno_bin_model_data.txt',row.names=F,quote=F,sep='\t')

fwrite(newdata,'paper_figures/pheno_bin_model_predictions.txt',row.names=F,quote=F,sep='\t')




####### By ABS BIN######
beta_merge$pheno_f=as.factor(beta_merge$phenotype)
beta_merge$abs_beta_z=abs(beta_merge$beta_z)


# Extreme bins are >=3/<=-3 SD
bbreaks=c(0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,max(beta_merge$abs_beta_z)+0.1)
beta_merge$abs_bin=cut(beta_merge$abs_beta_z, breaks=bbreaks, label=F)
#beta_merge$abs_bin=beta_z[match(beta_merge$gene_time_founder,beta_z$gene_time_founder),]$abs_bin

avg_z_abs=beta_merge %>% group_by(abs_bin,phenotype) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),n=length(freq),ngenes=length(unique(gene_founder)))
avg_z_abs=as.data.frame(avg_z_abs,stringsAsFactors=F)
nbins=max(avg_z_abs$abs_bin)


abs_pheno_model=lm(freq ~ phenotype*abs_bin,beta_merge)
anova(abs_pheno_model)

phenos=unique(beta_merge$phenotype)
newdata <- avg_z_abs[,c('abs_bin','phenotype')]
newdata$avg_freq <- predict(abs_pheno_model, newdata = newdata)

p <- predict(abs_pheno_model, newdata = newdata, se.fit=TRUE,     
         interval="confidence")
newdata=cbind(newdata,p$fit)


fwrite(avg_z_abs,'paper_figures/pheno_abs_bin_model_data.txt',row.names=F,quote=F,sep='\t')

fwrite(newdata,'paper_figures/pheno_abs_bin_model_predictions.txt',row.names=F,quote=F,sep='\t')


####### By Beta_Z ######
pheno_model2=lm(freq ~ phenotype*beta_z + phenotype*I(beta_z^2) ,beta_merge)
#Analysis of Variance Table
#
#Response: freq
#                          Df  Sum Sq Mean Sq  F value    Pr(>F)    
#phenotype                  5   0.125 0.02503  53.7479 < 2.2e-16 ***
#beta_z                     1   0.003 0.00309   6.6367  0.009991 ** 
#I(beta_z^2)                1   0.401 0.40097 861.1890 < 2.2e-16 ***
#phenotype:beta_z           5   0.003 0.00055   1.1721  0.320054    
#phenotype:I(beta_z^2)      5   0.016 0.00313   6.7157 2.891e-06 ***
#Residuals             256873 119.601 0.00047                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

poly1_by_pheno = emtrends(pheno_model2,specs = 'phenotype',var = 'beta_z')

poly2_by_pheno = emtrends(pheno_model2,specs = 'phenotype',var = 'I(beta_z^2)')
pheno_effects = contrast(poly2_by_pheno,
                           method = 'pairwise',
                           name = 'Pheno_effect')
#print(summary(pheno_effects,infer = c(T,T)))
beta_merge$bin_f<-cut(beta_merge$, 4, ordered = TRUE)
table(beta_merge$pheno_f)
beta_merge$pheno_f=factor(beta_merge$phenotype,levels=c(unique(beta_merge$phenotype)))
contrasts(beta_merge$pheno_f) = contr.poly(3)

summary(lm(write ~ readcat, hsb2))


sim_slopes(pheno_model2, pred = 'phenotype', modx = 'I(beta_z^2)',modx.values = c(-2, 0, 2), johnson_neyman = FALSE)

summary_table = summary(pheno_effects,infer = c(T,T),as.df = TRUE)

regrouped_pheno_effects = update(pheno_effects,by = 'Pheno_effect',adjust = 'none')
interaction_effects = contrast(regrouped_pheno_effects,
                                          'pairwise',
                                          name='Pheno_effects_on_poly2_effects')

print(summary(interaction_effects_vs_control,infer = c(T,T),level = 1-0.05/1))

# Not significantly different from one another
m.lst <- lstrends(pheno_model2, "phenotype", var="I(beta_z^2)")

model_0 = lm(freq  ~ poly(beta_z,2),data=beta_merge)

model_3 <- lm(freq ~ poly(beta_z, 2) * phenotype, data = beta_merge)

emmeans(pheno_model2, ~ phenotype | bin, at = list(bin = 1:14))

midbin_em=emmeans(pheno_model2, ~ phenotype | I(bin^2), at = list(bin = c(3,7,12)))
pairs(midbin_em)

anova(model_0,model_3)

p1=interactions::interact_plot(model_3, pred = beta_z, modx = phenotype, data = beta_merge)
p1=p1 +xlab("Founder Effect Z-score") + ylab("Frequency") +
scale_color_hue(name="Correlated Phenotype",labels=c("Days-To-Silking","Harvest Grain Moisture","Days-To-Anthesis","None","Thousand-Kernel Weight","Total Plant Height")) +
guides(linetype="none")

png('eqtl/images/quad_pheno_interaction_plot.png')
print(p1)
dev.off()


#### Are the poly2 terms for phenotype significantly different?

phenos=c("None","male_flowering_d6","female_flowering_d6","harvest_grain_moisture",
"total_plant_height","tkw_15")

#resdf1=c()
resdf2=c()

for(pheno in phenos){
	subdf=beta_merge[beta_merge$phenotype==pheno,]
	n=nrow(subdf)
	sub_model=lm(freq ~ beta_z + I(beta_z^2),subdf)
	#sub_model=lm(freq ~ bin + I(bin^2),subdf)

	#poly2 = lstrends(sub_model,var= 'I(beta_z^2)')

	anova(sub_model)
	res=summary(sub_model)
	#summary(pheno_model)
	#est1=res$coefficients[2,1]
	#s2=sd(subdf$freq)
	#smean=mean(subdf$freq)
	#se1=res$coefficients[2,2]
	#lwr1=est1 - (2*se1)
	#upr1=est1 + (2*se1)
	est2=res$coefficients[3,1]
	se2=res$coefficients[3,2]
	s_i=sqrt(se2^2 * n)
	#lwr2=est2 - (2*se2)
	#upr2=est2 + (2*se2)
	#sdp=sd(subdf$freq)
	#sem=sqrt(sdp**2/n)
	#line1=data.frame(phenotype=pheno,est=est1,se=se1,lwr=lwr1,upr=upr1,sd=sdp,sem=sem,n=n,stringsAsFactors=F)
	line2=data.frame(phenotype=pheno,est=est2,se=se2,s_i=s_i,n=n,stringsAsFactors=F)
	#resdf1=rbind(resdf1,line1)
	resdf2=rbind(resdf2,line2)
	
}
fwrite(resdf2,'eqtl/results/quad_pheno_pol2_all_estimates.txt',row.names=F,quote=F,sep='\t')

s2_pooled=0
for(i in 1:nrow(resdf2)){
	row1=resdf2[i,]
	#s2=row1$s_i^2
	n=row1$n
	sem=row1$se
	s2_pooled=s2_pooled+sem
}

df_full=nrow(beta_merge)-6

results=c()
for(i in 1:5){
	pheno_i=phenos[i]
	for(j in (i+1):6){
		pheno_j=phenos[j]
		ni=resdf2[resdf2$phenotype==pheno_i,]$n
		nj=resdf2[resdf2$phenotype==pheno_j,]$n
		
		est_i=resdf2[resdf2$phenotype==pheno_i,]$est
		est_j=resdf2[resdf2$phenotype==pheno_j,]$est
		SED=sqrt(s2_pooled/ni + s2_pooled/nj)
		t_c <- 4.03
		diff=est_i-est_j
		CI_l = diff - t_c*SED
		CI_h= diff + t_c*SED
		
		line1=data.frame(pheno_1=pheno_i,pheno_2=pheno_j,delta=diff,sed=SED,lwr=CI_l,upr=CI_h,stringsAsFactors=F)
		results=rbind(results,line1)
	}
}
results$tukey_t <- (results$delta)/results$sed

fwrite(results,'eqtl/results/quad_pheno_tukey_all_results.txt',row.names=F,quote=F,sep='\t')
#p-values=pt(-abs(results$tukey_t),df=pmin(num1,num2)-1)
#pt(t,df=pmin(num1,num2)-1)


#m.lst <- lstrends(pheno_model, "phenotype", var="I(beta_z^2)")
#pairs(m.lst)

beta_merge$pheno_f=as.factor(beta_merge$phenotype)
beta_merge$abs_beta_z=abs(beta_merge$beta_z)


# Extreme bins are >=3/<=-3 SD
bbreaks=c(0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,max(beta_merge$abs_beta_z)+0.1)
beta_merge$abs_bin=cut(beta_merge$abs_beta_z, breaks=bbreaks, label=F)
#beta_merge$abs_bin=beta_z[match(beta_merge$gene_time_founder,beta_z$gene_time_founder),]$abs_bin

avg_z=beta_merge %>% group_by(abs_bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$abs_bin)


pheno_model2=lm(freq ~ phenotype*abs_beta_z,beta_merge)
anova(pheno_model2)
#Analysis of Variance Table
#
#Response: freq
#                         Df Sum Sq Mean Sq  F value    Pr(>F)    
#phenotype                 5   0.41 0.08165  177.344 < 2.2e-16 ***
#abs_beta_z                1   1.78 1.77646 3858.235 < 2.2e-16 ***
#phenotype:abs_beta_z      5   0.03 0.00617   13.406 4.256e-13 ***
#Residuals            888263 408.99 0.00046                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(pheno_model2)



pheno_model2$coefficients
m.lst <- emtrends(pheno_model2, "phenotype", var="abs_beta_z")
pairs(m.lst)
#contrast                                      estimate       SE     df t.ratio
# female_flowering_d6 - harvest_grain_moisture -0.001665 0.000469 888263  -3.549 *
# female_flowering_d6 - male_flowering_d6      -0.001591 0.000406 888263  -3.913 *
# female_flowering_d6 - None                   -0.001848 0.000370 888263  -4.992 *
# female_flowering_d6 - tkw_15                 -0.000580 0.000463 888263  -1.252
# female_flowering_d6 - total_plant_height     -0.000428 0.000467 888263  -0.917
# harvest_grain_moisture - male_flowering_d6    0.000074 0.000339 888263   0.218 
# harvest_grain_moisture - None                -0.000183 0.000295 888263  -0.621 
# harvest_grain_moisture - tkw_15               0.001084 0.000406 888263   2.674
# harvest_grain_moisture - total_plant_height   0.001237 0.000410 888263   3.019 *
# male_flowering_d6 - None                     -0.000257 0.000179 888263  -1.438 
# male_flowering_d6 - tkw_15                    0.001010 0.000331 888263   3.051 *
# male_flowering_d6 - total_plant_height        0.001163 0.000336 888263   3.458 *
# None - tkw_15                                 0.001267 0.000285 888263   4.441 *
# None - total_plant_height                     0.001420 0.000291 888263   4.876 *
# tkw_15 - total_plant_height                   0.000152 0.000403 888263   0.377
# p.value
#  0.0052
#  0.0013
#  <.0001
#  0.8109
#  0.9424
#  0.9999
#  0.9895
#  0.0805
#  0.0305
#  0.7036
#  0.0277
#  0.0072
#  0.0001
#  <.0001
#  0.9990

# Different from None
# DTS, TKW, TPH

# Different from DTA
#TPH, TKW, DTS

# Different from HGM
# DTS, TKW, TPH
beta_merge=fread('eqtl/results/beta_merge_v3.txt',data.table=F)

#pheno_model3=lm(freq ~ phenotype*logit(abs_beta_z),beta_merge)
#anova(pheno_model3)
beta_merge=beta_merge[beta_merge$freq>=0.01,]
beta_merge$abs_beta_z=abs(beta_merge$beta_z)
for(time1 in times){
	submerge=beta_merge[beta_merge$time.x==time1,]
	p1=ggplot(data=submerge,aes(x=abs_beta_z,y=freq)) + geom_point(size=1) + 
	geom_smooth(mapping=aes(color=phenotype,fill=phenotype),method="lm",formula=y~x,linewidth = 1,alpha=0.3) +
	facet_wrap(.~phenotype) +
	xlab("Absolute Z-Score") + ylab("Founder Allele Frequency")

	png(sprintf('paper_figures/%s_top5k_abs_beta_z_by_pheno.png',time1))
	print(p1)
	dev.off()

}



p1=ggplot(data=beta_merge,aes(x=abs_beta_z,y=freq)) + geom_point(size=1) + 
geom_smooth(mapping=aes(color=phenotype,fill=phenotype),method="lm",formula=y~x,linewidth = 1,alpha=0.3) +
facet_wrap(.~phenotype) +
xlab("Absolute Z-Score") + ylab("Founder Allele Frequency")

png('paper_figures/abs_beta_z_by_pheno_top5k.png')
print(p1)
dev.off()


pheno_merge=beta_merge[beta_merge$phenotype!="None",]
snps=unique(pheno_merge$X_ID)
non_pheno=beta_merge[beta_merge$phenotype=="None",]
non_pheno=non_pheno[!(non_pheno$X_ID %in% snps),]

beta_merge2=rbind(pheno_merge,non_pheno)

p1=ggplot(data=beta_merge2,aes(x=abs_beta_z,y=freq)) + geom_point(size=1) + 
geom_smooth(mapping=aes(color=phenotype,fill=phenotype),method="lm",formula=y~x,linewidth = 1,alpha=0.3) +
facet_wrap(.~phenotype) +
xlab("Absolute Z-Score") + ylab("Founder Allele Frequency")

png('paper_figures/abs_beta_z_by_pheno_top5k_ld_drop.png')
print(p1)
dev.off()

p1=ggplot(data=beta_merge2,aes(x=abs_beta_z,y=freq)) + geom_point(size=1) + 
geom_smooth(mapping=aes(color=phenotype,fill=phenotype),linewidth = 1,alpha=0.3) +
facet_wrap(.~phenotype) +
xlab("Absolute Z-Score") + ylab("Founder Allele Frequency")

png('paper_figures/abs_beta_z_by_pheno_top5k_ld_drop_smooth.png')
print(p1)
dev.off()

pheno_model2=lm(freq ~ phenotype*abs_beta_z,beta_merge)
anova(pheno_model2)
p1=interactions::interact_plot(pheno_model2, pred = abs_beta_z, modx = phenotype, data = beta_merge)
p1=p1 +xlab("Absolute Founder Effect Z-score") + ylab("Frequency") +
scale_color_hue(name="Correlated Phenotype",labels=c("Days-To-Silking","Harvest Grain Moisture","Days-To-Anthesis","None","Thousand-Kernel Weight","Total Plant Height")) +
guides(linetype="none")

png('eqtl/images/abs_beta_z_pheno_interaction_plot.png')
print(p1)
dev.off()

newdata <- as.data.frame(expand.grid(abs_bin=bbreaks, phenotype=phenos))
newdata$avg_freq <- predict(pheno_model2, newdata = newdata)

p <- predict(pheno_model2, newdata = newdata, se.fit=TRUE,     
         interval="confidence")
newdata=cbind(newdata,p$fit)

fwrite(newdata,'paper_figures/pheno_model_predictions.txt',row.names=F,quote=F,sep='\t')

pheno_avg_z=beta_merge %>% group_by(phenotype,abs_bin) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),n=length(freq),ngenes=length(unique(gene_founder)))
pheno_avg_z=as.data.frame(pheno_avg_z,stringsAsFactors=F)
nbins=max(pheno_avg_z$abs_bin)
fwrite(pheno_avg_z,'paper_figures/pheno_abs_beta_z_bins.txt',row.names=F,quote=F,sep='\t')

newdata <- data.frame(abs_bin=pheno_avg_z$abs_bin, phenotype=pheno_avg_z$phenotype)
newdata$avg_freq <- predict(pheno_model2, newdata = newdata)

p <- predict(pheno_model2, newdata = newdata, se.fit=TRUE,     
         interval="confidence")
newdata=cbind(newdata,p$fit)

fwrite(newdata,'paper_figures/pheno_abs_model_predictions.txt',row.names=F,quote=F,sep='\t')
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





pheno_model2=lm(freq ~ phenotype*bin + phenotype*I(bin^2) ,beta_merge)

pheno_avg_z=beta_merge %>% group_by(bin,phenotype) %>% reframe(avg_freq=mean(freq),sd_freq=sd(freq),avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
pheno_avg_z=as.data.frame(pheno_avg_z,stringsAsFactors=F)
nbins=max(pheno_avg_z$bin)
fwrite(pheno_avg_z,'paper_figures/data/avg_freq_by_bin_pheno.txt',row.names=F,quote=F,sep='\t')



newdata <- data.frame(bin=pheno_avg_z$bin,phenotype=pheno_avg_z$phenotype)
#newdata$avg_rares <- predict(pheno_model2, newdata = newdata)

p <- predict(pheno_model2, newdata = newdata, se.fit=TRUE,     
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
