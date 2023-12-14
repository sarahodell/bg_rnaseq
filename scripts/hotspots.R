#!/usr/bin/env Rscript

library('ggplot2')
library('data.table')
library('dplyr')

times=c('WD_0712',"WD_0718","WD_0720","WD_0727")
allgenes=c()
for(time in times){
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	genes=names(exp)[-1]
	allgenes=c(allgenes,genes)
}


eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

otrans=fread('eqtl/results/all_trans_fdr_hits_FIXED.txt',data.table=F)
otrans$gene_time_snp=paste0(otrans$Trait,'-',otrans$time,'-',otrans$X_ID)


trans=fread('eqtl/results/all_trans_fdr_SIs_ld_FIXED.txt',data.table=F)
trans$p_value_ML=otrans[match(trans$gene_time_snp,otrans$gene_time_snp),]$p_value_ML
fwrite(trans,'eqtl/results/all_trans_fdr_SIs_ld_FIXED.txt',row.names=F,quote=F,sep='\t')

lgenes=unique(eqtl[eqtl$X_ID=="AX-91102858",]$Trait)
#  "Zm00001d010974" "Zm00001d010975" in FT gene list
ft_genelist=fread('eqtl/data/FT_genelist.txt',data.table=F)
time1="WD_0727"
tmp=trans[trans$time==time1,]

hotspot=trans[trans$time=="WD_0727" & trans$SNP=="AX-91102858",]
# Zm00001d011457 is a distal eQTL in T12 and T20
# All others (of 101) are unique to the timepoint
#  time        n
#  <chr>   <int>
#1 WD_0712    21
#2 WD_0718     7
#3 WD_0720    12
#4 WD_0727    62

prop_var=fread('MegaLMM/MegaLMM_WD_0727_prop_variance_FIXED.txt',data.table=F)
rownames(prop_var)=prop_var$V1
tgenes=unique(hotspot$gene)
allgenes=c(tgenes,lgenes)
prop_var=prop_var[allgenes,]
apply(prop_var[,-1], MARGIN=2,function(x) sum(x>0.1))

f4=prop_var[prop_var$Factor4>=0.1,c('V1','Factor4')]

#"Zm00001d010975" in F4

time1="WD_0712"
hotspot=trans[trans$time==time1 & trans$SNP=="AX-91102858",]

lspot=eqtl[eqtl$time==time1 & eqtl$X_ID=="AX-91102858", ]
lgene="Zm00001d010974"
lspot=lspot[,c("Trait","time","CHR","X_ID","BP",'p_value_ML','value')]
names(lspot)=c("gene","time","CHR","SNP","BP",'p_value_ML',"log10qvalue")
lspot$class="local"

hotspot=hotspot[,c("gene","time","CHR","SNP","BP","p_value_ML","value")]
names(hotspot)=c("gene","time","CHR","SNP","BP",'p_value_ML',"log10qvalue")

hotspot$class="distal"


hotspot=rbind(lspot,hotspot)
hotspot$value=-log10(hotspot$p_value_ML)

prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance_FIXED.txt',time1),data.table=F)
rownames(prop_var)=prop_var$V1
tgenes=unique(hotspot$gene)
#allgenes=c(tgenes,lgenes)
prop_var=prop_var[tgenes,]
apply(prop_var[,-1], MARGIN=2,function(x) sum(x>0.1))

# Factor 3 and Factor 8 

hotspot$f1_prop_var=prop_var[match(hotspot$gene,prop_var$V1),]$Factor1
hotspot$f6_prop_var=prop_var[match(hotspot$gene,prop_var$V1),]$Factor6
hotspot$f10_prop_var=prop_var[match(hotspot$gene,prop_var$V1),]$Factor10

hotspot$f12_prop_var=prop_var[match(hotspot$gene,prop_var$V1),]$Factor12

#hotspot$f3_prop_var=prop_var[match(hotspot$gene,prop_var$V1),]$Factor3
#hotspot$f8_prop_var=prop_var[match(hotspot$gene,prop_var$V1),]$Factor8

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

hotspot=merge(hotspot,genetable,by.x="gene",by.y="Gene_ID")

hotspot$mid_gene=apply(hotspot[,c('START','END')],1,function(x) round(mean(x)))

#Chr1: 92215734..92221005)#


localr=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
pei="male_flowering_d6-ALL-qDTA8"
#="male_flowering_d6-STPAUL_2017_WD-qDTA8"
#"female_flowering_d6-EXP_STPAUL_2017_WD-qDTS8"

distalr=fread('QTT/QTL_trans_eQTL_interval_overlap_cor.txt',data.table=F)

rs=localr[localr$Trait==lspot$gene & localr$time==lspot$time & localr$X_ID==lspot$SNP & localr$pheno_env_ID==pei,]$r
hotspot$gene_time_snp=paste0(hotspot$gene,'-',hotspot$time,'-',hotspot$SNP)
subdistal=distalr[distalr$pheno_env_ID==pei,]
hotspot$r=subdistal[match(hotspot$gene_time_snp,subdistal$gts),]$r
hotspot[hotspot$class=="local",]$r=rs

fwrite(hotspot,sprintf('QTT/%s_qDTA8_hotspot.txt',time1),row.names=F,quote=F,sep='\t')

df.tmp4 <- axisdf %>%
    left_join(hotspot, ., by=c("CHROM"="gene_chr")) %>%
    arrange(CHROM, mid_gene) %>%
    mutate(midgene_cum=as.numeric(mid_gene+tot))


theme_set(theme_classic())
#theme_update(text=element_text(family="Times"))
theme_update(plot.caption = element_text(hjust = 0))
theme_update(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
theme_update(plot.title = element_text(size=20),axis.title=element_text(size=20))
theme_update(panel.background=element_blank())
theme_update(plot.caption=element_text(size=20))

p1=ggplot(aes(x=midgene_cum,y=f12_prop_var),data=df.tmp4) +
 geom_point(aes(color=r,shape=class),size=3) +
 scale_color_gradientn(colors=c("red","grey","blue"),limits=c(-1,1), breaks=seq(-1,1,by=0.5)) + 
 xlab("Chromosome") +
ylab("T12-F12 Proportion Variance") + 
scale_x_continuous( label = axisdf$gene_chr, breaks= axisdf$center,minor_breaks=axisdf$tot) +
theme(panel.grid.minor.x=element_line(colour="darkgrey"),panel.grid.minor.y=element_blank())


png(sprintf('eqtl/images/%s_chr8_hotspot_BLUP_DTA8_cor.png',time1),width=1500,height=1000)
print(p1)
dev.off()
#######
rap27="Zm00001d010987"
# distal in T20 on chr1

#rap27_2="Zm00001d010988"
zcn8="Zm00001d010752"
# distal in T18, chr 2
# distal in T20 chr 1, 2, 5
# distal in T27 chr 10

mads69="Zm00001d042315"
# T12 chr 1, 7
# T20 chr 3
# T27 chr 2, 3, and 9


########


f4genes=tgenes[tgenes %in% f4$V1]
local=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
sub1=local[local$Trait %in% tgenes,]
#Factor2  Factor3  Factor4  Factor5  Factor6  Factor7  Factor8  Factor9 
#       8       21       42        4        6        1        1       16 
#Factor10 Factor11 Factor12 Factor13 
#      13       38       13        1 

# WD_0727 factor 4 go terms  - 223 floral development!


prop_var=fread('MegaLMM/MegaLMM_residuals_WD_0727_prop_variance_FIXED.txt',data.table=F)
rownames(prop_var)=prop_var$V1
tgenes=unique(hotspot$gene)
prop_var=prop_var[tgenes,]
apply(prop_var[,-1], MARGIN=2,function(x) sum(x>0.1))
#Factor1  Factor2  Factor3  Factor4  Factor5  Factor6  Factor7  Factor8 
#       4       15       24       15       15        5        7        7 
# Factor9 Factor10 Factor11 Factor12 Factor13 Factor14 Factor15 Factor16 
#       8       28       17        0        0        3        4        0 
#Factor17 Factor18 Factor19 
#       2        0        0 

bysnp=trans %>% group_by(time,SNP) %>% reframe(ngenes=length(unique(gene)))

# 18252 time-SNPs
summary(bysnp$ngenes)
# Just WD_0712
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   3.000   5.000   6.291   8.000  45.000

# ALL times
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   3.000   6.000   6.669   9.000  62.000

hotspots=bysnp[bysnp$ngenes>1,]
# 4699 SNPs, 17036 time-SNPs

fwrite(hotspots,'paper_figures/data/trans_hotspot_SNP_counts.txt',row.names=F,quote=F,sep='\t')

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

#8 AX-91102858 135530058 135656887  - 126829bp window

sub=trans[trans$SNP=="AX-91102858" & trans$time=="WD_0727",]
cand=fread('QTT/sig_candidate_genes.txt',data.table=F)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

sub=merge(sub,genetable,by.x="gene",by.y="Gene_ID")


# GO IDs of these genes


library('goseq')
library('rols')
time1="WD_0727"
genetable$LENGTH=genetable$END-genetable$START
genetable=genetable[order(genetable$CHROM,genetable$START),]
rownames(genetable)=seq(1,nrow(genetable))

annotation=fread('GO/maize.B73.AGPv4.aggregate.withancestors.csv',data.table=F)
annotation=merge(annotation,genetable,by.x="Gene",by.y="Gene_ID")
annotation=annotation[order(annotation$CHROM,annotation$START),]

log_inverse=function(x){
  	return(2^x)
}

genes=unique(annotation$Gene)
	
exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
rownames(exp)=exp$V1
exp=exp[,-1]


unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
avg_exp = apply(unlog,2,mean)
names(avg_exp)=names(exp)
avg_exp[avg_exp<1]=0
# re-log the input data
avg_logexp=log2(avg_exp)
avg_logexp[is.infinite(avg_logexp)]=0
	
#all_genes=intersect(unique(genetable$Gene_ID),unique(annotation$Gene))


pdf(sprintf('images/%s_nullp_plots.pdf',time1))

genes=names(avg_logexp)

test=unique(sub$gene)
deg=ifelse(genes %in% test,1,0)
names(deg)=genes
inter=intersect(test,genes[which(genes %in% test)])
test=inter
go=nullp(DEgenes=deg,bias.data=avg_logexp,plot.fit=T)
GO.samp=goseq(go,gene2cat=annotation[,c('Gene','GO')],use_genes_without_cat=TRUE,test.cats=c("GO:BP"))
GO.samp=as.data.frame(GO.samp,stringsAsFactors=F)
GO.samp$time=time1

dev.off()



#fwrite(all_go,sprintf('eqtl/results/chr8_trans_hotspot_GO_terms.txt',cutoff),row.names=F,quote=F,sep='\t')

GO.samp$adjusted_p=p.adjust(GO.samp$over_represented_pvalue,method="fdr")
enriched_GO=GO.samp[GO.samp$adjusted_p<=0.05,]
fwrite(enriched_GO,'eqtl/results/WD_0727_chr8_trans_hotspot_GO_terms.txt',row.names=F,quote=F,sep='\t')

# All chloroplast & thylakoid/plastid localized, no BP or MF terms

# Where are they located?
# Are these genes just in interchromosomal LD with this region?

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

eqtl[eqtl$X_ID=="AX-91102858",]
# Zm00001d010974 Zm00001d010975 Zm00001d010976
# Which of these has more genes with correlated effect sizes with trans-eQTL?
