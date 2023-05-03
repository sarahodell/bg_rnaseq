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

gg.manhattan2=function(df, threshold, col, ylims,bounds){
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
    #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    # Add label using ggrepel to avoid overlapping

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

#pcol=list("asi_P"=mypalette[7],"female_flowering_d6_P" = mypalette[1],"grain_yield_15_P"=mypalette[2],"harvest_grain_moisture_P"=mypalette[3],"male_flowering_d6_P"=mypalette[4],"tkw_15_P"=mypalette[5],"total_plant_height_P"=mypalette[6])

#threshtable=fread('eqtl/trans/threshold_0.05_table.txt',data.table=F)
#threshold=threshtable[threshtable$time==time & threshtable$factor==factor,]$threshold
#print(threshold)

#threshold=-log10((10^-threshold)/adjust)
df=c()
for(c in 1:10){
	print(c)
	dlist=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%.0f_c%s_fkeep_results.rds',c,time))
	d=as.data.frame(do.call(rbind, dlist))
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
	d$CHR=c
	d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
	df=rbind(df,d)
}
genes=unique(df$Trait)
df=df[!is.na(df$p_value_ML),]
df=df[!is.infinite(df$ML_logLik),]
n_m=4716
adjust=length(genes)*n_m*4
threshold=-log10(0.05/adjust)
print(threshold)


#factors=unique(full_df$Trait)
#factor=factors[f]
#factor=paste0('Factor',factors[f])
#df=full_df[full_df$Trait==factor,]


allsigs=c()
for(g in genes){
	sub=df[df$Trait==g,]
	sub=sub[order(sub$p_value_ML),]
	rownames(sub)=seq(1,nrow(sub))
	p_adjusted=p.adjust(sub$p_value_ML,method='fdr')
	#df$value=-log10(df$p_value_ML)
	sub$value=-log10(p_adjusted)

	sub=sub[,c('Trait','X_ID','p_value_ML','CHR','BP','value')]
	names(sub)=c('Gene','SNP','p_value_ML','CHR','BP','value')

	subdf=sub[,c('Gene','CHR','BP','SNP','value')]
	#qq=qqPlot(df$p_value_ML)
	subdf=subdf[order(subdf$CHR,subdf$BP),]
	subdf=subdf[,c('Gene','CHR','BP','SNP','value')]
	#qq=qqPlot(df$p_value_ML)
	subdf=subdf[order(subdf$CHR,subdf$BP),]

	ymax=round(max(subdf$value))+1
	#threshold=-log10(0.05)
	#threshold=-log10(0.05/(nrow(df)))

	title=sprintf("trans-eQTL for %s at timepoint %s",g,time)
	#title=sprintf("Residual trans-eQTL for %s at timepoint %s",factor,time)
	a2<-gg.manhattan2(subdf,threshold,
             col=greypalette,
             ylims=c(0,ymax)) + labs(caption = title)


	sigs2=subdf[subdf$value>=threshold,]

	if(dim(sigs2)[1]!=0){#png(sprintf('eqtl/trans/images/%s_pheno_residuals_trans_eQTL_manhattan_%s.png',factor,time),width=2000,height=1500)
  		png(sprintf('eqtl/trans/images/%s_%s_trans_eQTL_scan_manhattan.png',time,gene),width=2000,height=1500)
  		print(a2)
  		dev.off()
  		allsigs=rbind(allsigs,sigs2)
  		#fwrite(sigs2,sprintf('eqtl/results/%s_pheno_residuals_trans_%s_eQTL_hits.txt',factor,time),row.names=F,quote=F,sep='\t')
	}

}
allsigs=as.data.frame(allsigs,stringsAsFactors=F)
fwrite(allsigs,sprintf('eqtl/results/%s_trans_eQTL_scan_hits.txt',time),row.names=F,quote=F,sep='\t')

#order=match(df[order(df$p_value_ML),]$Gene,df$Gene)



#ymax=round(max(subdf$value))+1
#threshold=-log10(0.05)
#threshold=-log10(0.05/nrow(df))
#title=sprintf("trans-eQTL for %s at timepoint %s",factor,time)
#a2<-gg.manhattan2(subdf,threshold,
#             col=greypalette,
#             ylims=c(0,ymax)) + labs(caption = title)


#sigs1=subdf[subdf$value>=threshold,]

#if(dim(sigs1)[1]!=0){
#  png(sprintf('eqtl/trans/images/%s_trans_eQTL_manhattan_fdr_%s.png',factor,time),width=2000,height=1500)
#  print(a2)
#  dev.off()
#  fwrite(sigs1,sprintf('eqtl/results/%s_trans_%s_eQTL_hits.txt',factor,time),row.names=F,quote=F,sep='\t')
#}

# Using Bonferroni instead of 5% FDR
#df$value=-log10(df$p_value_ML)

#png(sprintf('eqtl/trans/images/%s_%s_trans_qqplot.png',time,factor))
#qqman::qq(df$p_value_ML)
#dev.off()


