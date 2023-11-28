#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')
library('cowplot')
#library('png')

##### Factor-eQTL manhattan plot
library('reshape2')
library('tibble')
library('tidyr')
#library('cowplot')
# Libraries ====
library('readr')
library('ggrepel')
library('RColorBrewer')
library('qqman')


time="WD_0720"
factor="Factor17"
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

full_nrow=433872
resid_nrow=358416

mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
greypalette=gray.colors(5)

#pcol=list("asi_P"=mypalette[7],"female_flowering_d6_P" = mypalette[1],"grain_yield_15_P"=mypalette[2],"harvest_grain_moisture_P"=mypalette[3],"male_flowering_d6_P"=mypalette[4],"tkw_15_P"=mypalette[5],"total_plant_height_P"=mypalette[6])
qtl=fread(sprintf('eqtl/results/%s_pheno_trans_%s_eQTL_hits.txt',factor,time),data.table=F)
#qtl=fread(sprintf('eqtl/results/%s_pheno_residuals_trans_%s_eQTL_hits.txt',factor,time),data.table=F)

qtl=qtl[,c('Factor','CHR','BP','SNP','new_value')]
names(qtl)=c('Factor','CHR','BP','SNP','value')
#qq=qqPlot(df$p_value_ML)

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
    gather(key, value, -BP,-SNP,-CHR,-BPcum,-tot,-Factor)

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

ymax=round(max(qtl$value))+1
#threshold=-log10(0.05)
resid_threshold=-log10(0.05/resid_nrow)
full_threshold=-log10(0.05/full_nrow)

#threshold=resid_threshold
threshold=full_threshold
title=sprintf("trans-eQTL for %s at timepoint %s",factor,time)
#title=sprintf("Residual trans-eQTL for %s at timepoint %s",factor,time)
a2<-gg.manhattan2(qtl,threshold,
             col=greypalette,
             ylims=c(0,ymax)) + labs(caption = title)


sigs2=qtl[qtl$new_value>=threshold,]

if(dim(sigs2)[1]!=0){
  png(sprintf('eqtl/images/%s_pheno_trans_eQTL_manhattan_%s.png',factor,time),width=2000,height=1500)
  #png(sprintf('eqtl/images/%s_pheno_residuals_trans_eQTL_manhattan_%s.png',factor,time),width=2000,height=1500)
  print(a2)
  dev.off()
  #fwrite(sigs2,sprintf('eqtl/results/%s_pheno_trans_%s_eQTL_hits.txt',factor,time),row.names=F,quote=F,sep='\t')
  #fwrite(sigs2,sprintf('eqtl/results/%s_pheno_residuals_trans_%s_eQTL_hits.txt',factor,time),row.names=F,quote=F,sep='\t')
}




#### Factor eQTL GO terms
factoreqtl=fread('eqtl/results/all_factor_fdr_SIs_FIXED.txt',data.table=F)
go=fread('GO/MegaLMM_0.10_GOSeq_enriched_FIXED.txt')

# WE
for(i in 1:nrow(factoreqtl)){
	row1=factoreqtl[i,]
	time1=row1$time
	factor1=row1$factor
	subgo=go[go$time==time1 & go$factor==factor1,]
	subgo=subgo[!is.na(subgo$term) & subgo$term!="",]
	subgo=subgo[!(subgo$term %in% c("biological process","molecular function"))]
	print(sprintf("%s %s has %.0f GO terms",time1,factor1,nrow(subgo)))
	
	subgo=go[go$time==time1 & go$factor==factor1,c('category','over_represented_pvalue')]
	fwrite(subgo,sprintf('GO/MegaLMM_%s_%s_enriched_GO_terms.txt',time1,factor1),row.names=F,quote=F,sep='\t')
}

#[1] "WD_0718 Factor4 has 99 GO terms"
#[1] "WD_0712 Factor2 has 6 GO terms"
#[1] "WD_0720 Factor17 has 53 GO terms"
#[1] "WD_0712 Factor12 has 0 GO terms"

# DRE
factoreqtl=fread('eqtl/results/all_residual_factor_fdr_SIs_FIXED.txt',data.table=F)
go=fread('GO/MegaLMM_residuals_0.10_GOSeq_enriched_FIXED.txt')

for(i in 1:nrow(factoreqtl)){
	row1=factoreqtl[i,]
	time1=row1$time
	factor1=row1$factor
	subgo=go[go$time==time1 & go$factor==factor1,]
	subgo=subgo[!is.na(subgo$term) & subgo$term!="",]
	subgo=subgo[!(subgo$term %in% c("biological_process","molecular_function"))]
	print(head(subgo[subgo$ontology=="BP",]))
	print(sprintf("%s %s has %.0f GO terms",time1,factor1,nrow(subgo)))
	#subgo=go[go$time==time1 & go$factor==factor1,c('category','over_represented_pvalue')]

	#fwrite(subgo,sprintf('GO/MegaLMM_residuals_%s_%s_enriched_GO_terms.txt',time1,factor1),row.names=F,quote=F,sep='\t')
}

#[1] "WD_0718 Factor4 has 134 GO terms"
#[1] "WD_0718 Factor18 has 25 GO terms"
#[1] "WD_0718 Factor7 has 196 GO terms"
#[1] "WD_0712 Factor2 has 16 GO terms"
#[1] "WD_0712 Factor5 has 2 GO terms"
#[1] "WD_0727 Factor12 has 3 GO terms"
#[1] "WD_0727 Factor17 has 8 GO terms"
#[1] "WD_0712 Factor23 has 10 GO terms"
#[1] "WD_0720 Factor19 has 75 GO terms"
#[1] "WD_0720 Factor2 has 423 GO terms"
#[1] "WD_0718 Factor12 has 60 GO terms"
#[1] "WD_0712 Factor9 has 24 GO terms"
#[1] "WD_0712 Factor23 has 10 GO terms"
#[1] "WD_0712 Factor23 has 10 GO terms"
#[1] "WD_0720 Factor17 has 144 GO terms"
#[1] "WD_0727 Factor19 has 7 GO terms"

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

# one of them is WD_0720_Factor17_revigo_FIXED.png
#pb=readPNG('paper_figures/WD_0720_Factor17_revigo_FIXED.png')
#pb<-ggdraw() + draw_image('paper_figures/WD_0720_Factor17_revigo_FIXED.png')


# DRE T20 Factor 17
prop_var=fread('MegaLMM/MegaLMM_residuals_WD_0720_prop_variance_FIXED.txt',data.table=F)
rownames(prop_var)=prop_var$V1

dre_f19=prop_var[prop_var$Factor19>=0.1,c('V1','Factor19')]
dre_f19=merge(dre_f19,genetable,by.x="V1",by.y="Gene_ID")

env1=as.data.table(factoreqtl[9,])
env2=as.data.table(dre_f19)
setkey(env2,CHROM,START,END)
comp2=foverlaps(env1,env2,by.x=c("CHR","left_bound_bp","alt_right_bound_bp"),by.y=c("CHROM","START","END"))

#####
# 579 shared genes between the 2
factoreqtl=fread('eqtl/results/all_factor_fdr_SIs_FIXED.txt',data.table=F)
factoreqtl$factor_time=paste0(factoreqtl$time,'-',factoreqtl$factor)
factoreqtls=unique(factoreqtl$factor_time)

time1="WD_0720"

dre_prop_var=fread('MegaLMM/MegaLMM_WD_0720_residuals_all_F_means_FIXED.txt',data.table=F)
rownames(dre_prop_var)=dre_prop_var$V1
we_prop_var=fread('MegaLMM/MegaLMM_WD_0720_all_F_means_FIXED.txt',data.table=F)
rownames(we_prop_var)=we_prop_var$V1
ginter=intersect(rownames(we_prop_var),rownames(dre_prop_var))
we_prop_var=we_prop_var[ginter,]
dre_prop_var=dre_prop_var[ginter,]

prop_var=fread('MegaLMM/MegaLMM_WD_0720_prop_variance_FIXED.txt',data.table=F)
f17=prop_var[prop_var$Factor17>=0.1,c('V1','Factor17')]
f17=merge(f17,genetable,by.x="V1",by.y="Gene_ID")

env1=as.data.table(factoreqtl[9,])
env2=as.data.table(f17)
setkey(env2,CHROM,START,END)
comp2=foverlaps(env1,env2,by.x=c("CHR","left_bound_bp","alt_right_bound_bp"),by.y=c("CHROM","START","END"))


env1=as.data.table(factoreqtl[3,])
env2=as.data.table(f17)
setkey(env2,CHROM,START,END)
comp=foverlaps(env1,env2,by.x=c("CHR","left_bound_bp","alt_right_bound_bp"),by.y=c("CHROM","START","END"))
#Zm00001d041620 0.9106396 CTP synthase

axisdf=fread('eqtl/data/chromosome_axis.txt',data.table=F)
cumtot=2105119857
############### 
#Plot out position by prop_var
f="WD_0720-Factor17"
sub=factoreqtl[factoreqtl$factor_time==f,]
time=unique(sub$time)
factor=unique(sub$factor)
	
prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance_FIXED.txt',time),data.table=F)
df=prop_var[prop_var[,factor]>0.1,c('V1',factor)]
names(df)=c('Gene_ID','prop_var')
df=merge(df,genetable,by="Gene_ID")
df$midgene=round(df$START + (df$END-df$START)/2)
df$tot=axisdf[match(df$CHROM,axisdf$gene_chr),]$tot
df$midgene_cum=df$midgene + df$tot 
	
sub$tot=axisdf[match(sub$CHR,axisdf$gene_chr),]$tot
sub$bp_cum=sub$BP + sub$tot 
 
ngenes=nrow(df)
print(sprintf('%.0f Genes loaded on %s %s',ngenes,time, factor))


# Add shaded regions for qHGM3_2 SI
qleft=133491494 + 601388623
qright=163726450 + 601388623

fleft=121952745  +601388623
fright=139134048 +601388623
# Add shaded region around factor-eQTL support interval

df.tmp3 <- axisdf %>%
    left_join(df, ., by=c("CHROM"="gene_chr")) %>%
    arrange(CHROM, midgene) %>%
    mutate(midgene_cum=as.numeric(midgene+chr_start))



axisdf$chr_start=as.numeric(axisdf$chr_start)
axisdf$chr_end=as.numeric(axisdf$chr_end)
axisdf$center_shift=as.numeric(axisdf$center_shift)
minor_breaks=sort(c(0,axisdf$chr_start[-1],axisdf$chr_end))
cumtot=axisdf[axisdf$gene_chr==10,]$chr_end
cumtot=2330119857

fwrite(df.tmp3,'paper_figures/data/T20_Factor17_gene_prop_var.txt',row.names=F,quote=F,sep='\t')

pa=ggplot(data=df.tmp3,aes(x=midgene_cum,y=prop_var)) +
    scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center_shift,minor_breaks=minor_breaks,limits=c(0, cumtot)) +
    ylim(0.1,1) +
    #scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center_shift,minor_breaks=minor_breaks,limits=c(0.1, 1)) +
    #geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
	geom_point(aes(color=prop_var),size=1) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
	xlab("Position") + ylab("Proportion Variance") +
	ggtitle(sprintf("Gene Loadings on T20 %s, n=%.0f",factor,ngenes)) +
	theme_classic() +
    theme(panel.grid.minor.x=element_line(colour="darkgrey"),panel.grid.major=element_blank(),panel.grid.minor.y=element_blank())
 	
 #for(i in 1:nrow(sub)){
 #	row=sub[i,]
BP_cum=127983917+601388623 
 #}
 pa=pa + annotate("rect",xmin=qleft,xmax=qright,ymin=0.1,ymax=1,alpha=0.5,fill='yellow') + 
 #geom_vline(xintercept=BP_cum,color="coral") +
 annotate("rect",xmin=fleft,xmax=fright,ymin=0.1,ymax=1,alpha=0.5,fill='coral')
 
pdf(sprintf('eqtl/images/%s_%s_plot.pdf',time,factor))
print(pa)
dev.off()





id="qHGM3_2"
pheno="harvest_grain_moisture"
env="ALL"
qsnp="PZE-103084819"

time1="WD_0720"
factor1="Factor17"
esnp="AX-90826809"
#-0.6890176 T20 Factor 17

effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr2,pheno,env),data.table=F)
effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
effect_size=unlist(effect_size[,c(6:21)])
wn=which(!is.na(effect_size))[1]
effect_size[-wn]=effect_size[-wn]+effect_size[wn]

results=fread(sprintf('eqtl/trans/results/%s_c%s_%s_trans_results_FIXED.txt',time1,chr2,factor1),data.table=F)
result=results[results$X_ID==esnp,]
betas=unlist(result[,founders])
wn=which(!is.na(betas))[1]
betas[-wn]=betas[-wn]+betas[wn]

r=cor(betas,effect_size,use="complete.obs")

tmpdf=data.frame(ebeta=betas,qbeta=effect_size,founder=founders,stringsAsFactors=F)

fwrite(tmpdf,'paper_figures/data/T20_Factor17_qHGM3_2_effect_sizes.txt',row.names=F,quote=F,sep='\t')


colorcodes=fread('metadata/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
founderlabels=c("B73","A632","CO255","F252","OH43","A654","F2","C103","EP1","D105","W117","B96","DK63","F492","ND245","VA85")
colorcodes=colorcodes[founders,]


pc=ggplot(aes(x=ebeta,y=qbeta),data=tmpdf) + geom_point(aes(color=founder)) + xlab("eQTL Effect Size") +
ggtitle(label=sprintf("r=%.2f",r)) + theme_classic() +
scale_color_manual(values=colorcodes$hex_color,labels=founderlabels) +
ylab("QTL Effect Size") 


prow <- plot_grid(
  pa + theme(legend.position="none"),
  pb,
  pc,
  #c2 + theme(legend.position="none",axis.ticks=element_blank()),
  #  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1,
  ncol=3
)

png('paper_figures/figure2.png')
print(prow)
dev.off()




### WE

factor_df=c()
genes_per=c()
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time1 in times){
	#fval=fread(sprintf('MegaLMM/MegaLMM_%s_all_F_means_FIXED.txt',time1),data.table=F)
	prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance_FIXED.txt',time1),data.table=F)
	factors=names(prop_var)[-1]
	for(f in factors){
		ngenes=nrow(prop_var[prop_var[,f]>=0.1,])
		line=data.frame(time=time1,factor=f,ngenes=ngenes,stringsAsFactors=F)
		factor_df=rbind(factor_df,line)
	}
	rownames(prop_var)=prop_var$V1
	prop_var=prop_var[,-1]
	nper=apply(prop_var,MARGIN=1,function(x) sum(x>=0.1))
	line2=data.frame(time=time1,mean=mean(nper),median=median(nper),stringsAsFactors=F)
	genes_per=rbind(genes_per,line2)
}
summary(factor_df$ngenes)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   61.0   515.2  1352.5  2540.1  3933.0  9797.0 
fwrite(factor_df,'paper_figures/data/WE_factor_gene_count.txt',row.names=F,quote=F,sep='\t')

p1=ggplot(data=factor_df,aes(x=ngenes)) + geom_histogram(bins=20) + theme_classic()+
xlab("Number of Genes Per Factor") + ylab("Frequency") + geom_vline(xintercept=mean(factor_df$ngenes))


png('paper_figures/WE_genes_per_factor.png')
print(p1)
dev.off()

## DRE

factor_df2=c()
genes_per2=c()
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time1 in times){
	#fval=fread(sprintf('MegaLMM/MegaLMM_%s_all_F_means_FIXED.txt',time1),data.table=F)
	prop_var=fread(sprintf('MegaLMM/MegaLMM_residuals_%s_prop_variance_FIXED.txt',time1),data.table=F)
	factors=names(prop_var)[-1]
	for(f in factors){
		ngenes=nrow(prop_var[prop_var[,f]>=0.1,])
		line=data.frame(time=time1,factor=f,ngenes=ngenes,stringsAsFactors=F)
		factor_df2=rbind(factor_df2,line)
	}
	rownames(prop_var)=prop_var$V1
	prop_var=prop_var[,-1]
	nper=apply(prop_var,MARGIN=1,function(x) sum(x>=0.1))
	line2=data.frame(time=time1,mean=mean(nper),median=median(nper),stringsAsFactors=F)
	genes_per2=rbind(genes_per2,line2)
}

summary(factor_df2$ngenes)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   12.0   508.5  2186.0  2462.6  3869.5  7958.0 

fwrite(factor_df2,'paper_figures/data/DRE_factor_gene_count.txt',row.names=F,quote=F,sep='\t')


p1=ggplot(data=factor_df2,aes(x=ngenes)) + geom_histogram(bins=20) + theme_classic()+
xlab("Number of Genes Per Factor") + ylab("Frequency") + geom_vline(xintercept=mean(factor_df2$ngenes))


png('paper_figures/DRE_genes_per_factor.png')
print(p1)
dev.off()



prop_var=fread('MegaLMM')

# Are any of the qDTA3_2 candidate genes in T20 Factor 17
prop_var=fread('MegaLMM/MegaLMM_WD_0720_prop_variance_FIXED.txt',data.table=F)

f17genes=prop_var[prop_var$Factor17>=0.1,c('V1','Factor17')]
cand=fread('QTT/sig_candidate_genes.txt',data.table=F)
cgenes=cand[cand$ID=="qDTA3_2",]$Trait

sum(cgenes %in% f17genes$V1)


# DRE
dre=fread('eqtl/results/all_residual_factor_fdr_SIs_FIXED.txt',data.table=F)
