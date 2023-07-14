#!/usr/bin/env Rscript

library('ggplot2')
library('data.table')
library('reshape2')
library('tibble')
library('dplyr')
library('tidyr')
library('qqman')
#library('cowplot')
# Libraries ====
library('readr')
library('ggrepel')
library('RColorBrewer')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
df=c()
for(time in times){
	for(c in 1:10){
  		d=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_weights_results_FIXED.txt',time,c))
  		d$time=time
  		pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  		d$CHR=c
  		d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
  		d=d[,c('Trait','X_ID','p_value_ML','CHR','BP','time')]
  		df=rbind(df,d)
	}
}
df=as.data.frame(df,stringsAsFactors=F)

df=df[!is.na(df$p_value_ML),]
df=df[order(df$p_value_ML),]
rownames(df)=seq(1,nrow(df))
df$p_adjusted=p.adjust(df$p_value_ML,method='fdr')
#df$value=-log10(df$p_adjusted)
df$value=-log10(df$p_value_ML)

#threshold=-log10(0.05)
adjust=nrow(df)
threshold=-log10(0.01/adjust)
sig=df[df$value>=threshold,]
sig=sig[order(sig$CHR,sig$BP),]

png('eqtl/images/all_cis_qqplot.png')
print(qqman::qq(df$p_value_ML))
dev.off()


#fwrite(sig,'eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',row.names=F,quote=F,sep='\t')
fwrite(sig,'eqtl/results/all_cis_eQTL_weights_1bonf_hits_FIXED.txt',row.names=F,quote=F,sep='\t')


founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
greypalette=gray.colors(5)


gg.manhattan2 <- function(fdf, threshold, col, ylims,bounds){
  df.tmp <- fdf %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
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

for(time in times){
	subdf=df[df$time==time,]
	ymax=round(max(subdf$value))+1
	#threshold=-log10(0.05)

	title=sprintf("cis-eQTL at timepoint %s",time)
	a2<-gg.manhattan2(fdf=subdf,threshold,
             col=greypalette,
             ylims=c(0,ymax)) + labs(caption = title)

	#png(sprintf('eqtl/cis/images/ciseQTL_weights_fdr_manhattan_%s_FIXED.png',time),width=2000,height=1500)

	png(sprintf('eqtl/cis/images/ciseQTL_weights_1bonf_manhattan_%s_FIXED.png',time),width=2000,height=1500)
	print(a2)
	dev.off()
}


