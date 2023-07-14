args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])


library('ggplot2')
library('data.table')
library('reshape2')
library('tibble')
library('dplyr')
library('tidyr')
#library('cowplot')
# Libraries ====
library('readr')
library('ggrepel')
library('RColorBrewer')
#library('qqman')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
greypalette=gray.colors(5)

#pcol=list("asi_P"=mypalette[7],"female_flowering_d6_P" = mypalette[1],"grain_yield_15_P"=mypalette[2],"harvest_grain_moisture_P"=mypalette[3],"male_flowering_d6_P"=mypalette[4],"tkw_15_P"=mypalette[5],"total_plant_height_P"=mypalette[6])

df=c()
for(c in 1:10){
  d=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_weights_results_FIXED.txt',time,c))
  #d=d[complete.cases(d),]
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  d$CHR=c
  d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
  df=rbind(df,d)
}

#times=c("WD_0712","WD_0718","WD_0720","WD_0727")
#alldf=c()
#for(time in times){
#	df=c()
#	for(c in 1:10){
 # 		d=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_fkeep_results.txt',time,c))
  #		#d=d[complete.cases(d),]
#  		pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
#  		d$CHR=c
#  		d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
#  		df=rbind(df,d)#
#	}
#	alldf=rbind(alldf,df)
#}

#alldf=alldf[!is.na(alldf$p_value_ML),]

#png(sprintf('eqtl/cis/images/%s_fkeep_qqplot.png',time))
#qqman::qq(df$p_value_ML)
#dev.off()

df=df[!is.na(df$p_value_ML),]

#order=match(df[order(df$p_value_ML),]$Gene,df$Gene)
#df=df[order(df$p_value_ML),]
#rownames(df)=seq(1,nrow(df))
#p_adjusted=p.adjust(df$p_value_ML,method='fdr')
#df$p_adjusted=p_adjusted
df$value=-log10(df$p_value_ML)
#df$value=-log10(df$p_adjusted)

#png(sprintf('eqtl/cis/images/%s_qqplot.png',time))
#qqman::qq(df$p_value_ML)
#dev.off()

df=df[,c('Trait','X_ID','p_value_ML','CHR','BP','value')]
names(df)=c('Gene','SNP','p_value_ML','CHR','BP','value')

df2=df[,c('Gene','CHR','BP','SNP','value')]
#qq=qqPlot(df$p_value_ML)
df2=df2[order(df$CHR,df$BP),]

gg.manhattan2 <- function(df, threshold, col, ylims,bounds){
  # format df
  df.tmp <- df %>%

    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=as.numeric(BP+tot)) %>%
    gather(key, value, -BP,-SNP,-CHR,-BPcum,-tot,-Gene)

    # Add highlight and annotation information
    #mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) #%>%
    #mutate( is_annotate=ifelse(P < threshold, "yes", "no"))

  df.tmp$sig=df.tmp$value > threshold
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  df.tmp=df.tmp[!is.na(df.tmp$value),]
  #df.tmp=df.tmp[df.tmp$value >= quantile(df.tmp$value,0.3),]
  rownames(df.tmp)=seq(1,nrow(df.tmp))
  ggplot(df.tmp, aes(x=BPcum, y=value)) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +

    scale_color_manual(values = rep(col, 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis

    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +

    # add genome-wide sig and sugg lines
    geom_hline(yintercept = threshold,linetype="dashed") +
    #geom_hline(yintercept = -log10(sugg), linetype="dashed") +

    # Add highlighted points
    geom_point(data=subset(df.tmp, sig==T), color="coral2", size=5) +
#    geom_label(data=subset(df.tmp,sig==T,aes(label=Gene),position = position_dodge(0.9),vjust = 0) +
    #geom_text(data=subset(df.tmp, sig==T),
    #        aes(label=Gene),,position = position_dodge(1.0),vjust = 0,size=16) +
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +

    # Custom the theme:
    theme_classic() +
    theme(
      text = element_text(size=30),
      plot.title = element_text(hjust = 0.5),
      #legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + guides(color="none")
}

ymax=round(max(df2$value))+1
#threshold=-log10(0.05)
adjust=31879
threshold=-log10(0.05/adjust)
title=sprintf("cis-eQTL at timepoint %s",time)
a2<-gg.manhattan2(df2,threshold,
             col=greypalette,
             ylims=c(0,ymax)) + labs(caption = title)



png(sprintf('eqtl/cis/images/ciseQTL_weights_manhattan_%s_FIXED.png',time),width=2000,height=1500)
print(a2)
dev.off()

#sigs=df[df$value>=threshold,]
#fwrite(sigs,sprintf('eqtl/results/%s_cis_eQTL_weights_hits.txt',time),row.names=F,quote=F,sep='\t')

df=df[order(df$p_value_ML),]
rownames(df)=seq(1,nrow(df))
#p_adjusted=p.adjust(df$p_value_ML,method='fdr')
df$p_adjusted=p.adjust(df$p_value_ML,method='fdr')
#df$value=-log10(df$p_value_ML)
df$value=-log10(df$p_adjusted)
df2=df[,c('Gene','CHR','BP','SNP','value')]
#qq=qqPlot(df$p_value_ML)
df2=df2[order(df$CHR,df$BP),]


ymax=round(max(df2$value))+1
threshold=-log10(0.05/4)
#adjust=31879
#threshold=-log10(0.05/adjust)
title=sprintf("cis-eQTL at timepoint %s",time)
a2<-gg.manhattan2(df2,threshold,
             col=greypalette,
             ylims=c(0,ymax)) + labs(caption = title)



png(sprintf('eqtl/cis/images/ciseQTL_weights_fdr_manhattan_%s_FIXED.png',time),width=2000,height=1500)
print(a2)
dev.off()

#sigs=df[df$value>=threshold,]
#fwrite(sigs,sprintf('eqtl/results/%s_cis_eQTL_weights_fdr_hits.txt',time),row.names=F,quote=F,sep='\t')
