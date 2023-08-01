#!/usr/bin/env Rscript

library('data.table')
library('qqman')
library('qvalue')
library('dplyr')
library('ggplot2')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
df=c()
for(time in times){
	factors=fread(sprintf('eqtl/trans/%s_factors_FIXED.txt',time),header=F,data.table=F)
	for(f in factors$V1){
		for(c in 1:10){
  			d=fread(sprintf('eqtl/trans/results/%s_c%.0f_%s_trans_results_FIXED.txt',time,c,f),data.table=F)
 			pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  			d$time=time
  			d$factor=f
  			d$CHR=c
 			d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
 			d=d[,c('Trait','X_ID','p_value_ML','time','factor','CHR','BP')]
  			df=rbind(df,d)
		}
	}
}

df=as.data.frame(df)
df=df[!is.na(df$p_value_ML),]
df=df[order(df$p_value_ML),]
rownames(df)=seq(1,nrow(df))
qvalues=qvalue(df$p_value_ML,fdr.level=0.05)
summary(qvalues)

write.qvalue(qvalues,'eqtl/results/factor_eqtl_qvalues.txt',sep='\t')

#df$p_adjusted=p.adjust(df$p_value_ML,method='fdr')
df$p_adjusted=qvalues$qvalues
df$value=-log10(df$p_adjusted)


threshold=-log10(0.05)
print(threshold)

sig=df[df$value>=threshold,]
sig=sig[order(sig$CHR,sig$BP),]
fwrite(sig,'eqtl/results/all_factor_trans_eqtl_fdr_hits_FIXED.txt',row.names=F,quote=F,sep='\t')

png('eqtl/images/all_factor_fdr_qqplot_FIXED.png')
print(qqman::qq(df$p_value_ML))
dev.off()

pdf('eqtl/trans/images/factor_transeqtl_qvalue_diagnostics.pdf')
hist(qvalues)
plot(qvalues)
dev.off()


founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
greypalette=gray.colors(5)


gg.manhattan2 <- function(fdf, threshold, col, ylims,bounds){
  df.tmp <- fdf %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  #select(-chr_len) %>%
  left_join(fdf, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=as.numeric(BP+tot)) #%>%
  #gather(key="", value="value", -BP,-X_ID,-CHR,-BPcum,-tot,-Trait)
  df.tmp=df.tmp[!is.na(df.tmp$value),]
  #df.tmp$value=as.numeric(df.tmp$value)
  df.tmp$sig=df.tmp$value > threshold
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  rownames(df.tmp)=seq(1,nrow(df.tmp))
  ggplot(df.tmp, aes(x=BPcum, y=value)) +
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22)) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + 
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    geom_hline(yintercept = threshold,linetype="dashed") +
    geom_point(data=subset(df.tmp, sig==T), color="coral2", size=5) +
    theme_classic() +
    theme(
      text = element_text(size=30),
      plot.title = element_text(hjust = 0.5),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) + guides(color="none")
}

times=c("WD_0718")
for(time in times){
	subdf=df[df$time==time,]
	sig_f=sig[sig$time==time,]
	sig_factors=unique(sig_f$factor)
	threshold=-log10(0.05)
	for(factor in sig_factors){
		subdf2=subdf[subdf$factor==factor,]
		ymax=round(max(subdf$value))+1
		title=sprintf("%s trans-eQTL at timepoint %s",factor,time)
		a2<-gg.manhattan2(fdf=subdf2,threshold,
             col=greypalette,
             ylims=c(0,ymax)) + labs(caption = title)

		png(sprintf('eqtl/trans/images/%s_trans_eQTL_weights_fdr_manhattan_%s_FIXED.png',factor,time),width=2000,height=1500)
		print(a2)
		dev.off()
	}
	
}




############ Residuals ############

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
df=c()
for(time in times){
	factors=fread(sprintf('eqtl/trans/%s_residuals_factors_FIXED.txt',time),header=F,data.table=F)
	for(f in factors$V1){
		for(c in 1:10){
  			d=fread(sprintf('eqtl/trans/results/%s_residuals_c%.0f_%s_trans_results_FIXED.txt',time,c,f),data.table=F)
 			pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  			d$time=time
  			d$factor=f
  			d$CHR=c
 			d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
 			d=d[,c('Trait','X_ID','p_value_ML','time','factor','CHR','BP')]
  			df=rbind(df,d)
		}
	}
}

df=as.data.frame(df)
df=df[!is.na(df$p_value_ML),]
df=df[order(df$p_value_ML),]
rownames(df)=seq(1,nrow(df))
qvalues=qvalue(df$p_value_ML,fdr.level=0.05)
summary(qvalues)

write.qvalue(qvalues,'eqtl/results/residuals_factor_eqtl_qvalues.txt',sep='\t')

#df$p_adjusted=p.adjust(df$p_value_ML,method='fdr')
df$p_adjusted=qvalues$qvalues
df$value=-log10(df$p_adjusted)


threshold=-log10(0.05)
print(threshold)

sig=df[df$value>=threshold,]
sig=sig[order(sig$CHR,sig$BP),]
fwrite(sig,'eqtl/results/all_residuals_factor_trans_eqtl_fdr_hits_FIXED.txt',row.names=F,quote=F,sep='\t')

png('eqtl/images/all_residuals_factor_fdr_qqplot_FIXED.png')
print(qqman::qq(df$p_value_ML))
dev.off()

pdf('eqtl/trans/images/residuals_factor_transeqtl_qvalue_diagnostics.pdf')
hist(qvalues)
plot(qvalues)
dev.off()


founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
greypalette=gray.colors(5)


gg.manhattan2 <- function(fdf, threshold, col, ylims,bounds){
  df.tmp <- fdf %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  #select(-chr_len) %>%
  left_join(fdf, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=as.numeric(BP+tot)) #%>%
  #gather(key="", value="value", -BP,-X_ID,-CHR,-BPcum,-tot,-Trait)
  df.tmp=df.tmp[!is.na(df.tmp$value),]
  #df.tmp$value=as.numeric(df.tmp$value)
  df.tmp$sig=df.tmp$value > threshold
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  rownames(df.tmp)=seq(1,nrow(df.tmp))
  ggplot(df.tmp, aes(x=BPcum, y=value)) +
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22)) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + 
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    geom_hline(yintercept = threshold,linetype="dashed") +
    geom_point(data=subset(df.tmp, sig==T), color="coral2", size=5) +
    theme_classic() +
    theme(
      text = element_text(size=30),
      plot.title = element_text(hjust = 0.5),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) + guides(color="none")
}

times=c("WD_0712","WD_0718","WD_0720")
for(time in times){
	subdf=df[df$time==time,]
	sig_f=sig[sig$time==time,]
	sig_factors=unique(sig_f$factor)
	threshold=-log10(0.05)
	for(factor in sig_factors){
		subdf2=subdf[subdf$factor==factor,]
		ymax=round(max(subdf$value))+1
		title=sprintf("%s Residuals trans-eQTL at timepoint %s",factor,time)
		a2<-gg.manhattan2(fdf=subdf2,threshold,
             col=greypalette,
             ylims=c(0,ymax)) + labs(caption = title)

		png(sprintf('eqtl/trans/images/%s_residuals_trans_eQTL_weights_fdr_manhattan_%s_FIXED.png',factor,time),width=2000,height=1500)
		print(a2)
		dev.off()
	}
	
}


