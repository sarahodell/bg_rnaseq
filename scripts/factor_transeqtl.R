#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')



factoreqtl=fread('eqtl/results/all_factor_fdr_peaks.txt',data.table=F)
factordf=fread('eqtl/results/all_factor_trans_eqtl_fdr_genes.txt',data.table=F)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
axisdf=fread('eqtl/data/chromosome_axis.txt',data.table=F)
cumtot=2105119857

# For each factor
#		- How many genes are loaded on them
#		- Where are the genes located in the genome
#		- What prop variance is explained by the factor for each gene
#		- How much variation in the factor is explained by the factor eQTL
#		- Are there trans-eQTL or cis-eQTL that overlap with the factor-eQTL variant? 
#		- Are there trans-eQTL or cis-eQTL genes that load on the factor?

factoreqtl$time_factor=paste0(factoreqtl$time,'-',factoreqtl$Trait)
tf=unique(factoreqtl$time_factor)
for(f in tf){
	sub=factoreqtl[factoreqtl$time_factor==f,]
	time=unique(sub$time)
	factor=unique(sub$factor)
	df=factordf[factordf$time==time & factordf$Trait==factor,]
	df$midgene=round(df$gene_start + (df$gene_end-df$gene_start)/2)
	df$tot=axisdf[match(df$gene_chr,axisdf$gene_chr),]$tot
	df$tot.x=axisdf[match(df$CHR,axisdf$gene_chr),]$tot
	BP_cum=unique(df$BP) + unique(df$tot.x)
	df$midgene_cum=df$midgene + df$tot 
	prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance.txt',time),data.table=F)
	print(time)
	print(factor)
	fgenes=prop_var[prop_var[,factor]>0.1,]$V1
	ngenes=length(fgenes)
	print(sprintf('%.0f Genes loaded on %s',ngenes,factor))
	
	p=ggplot(df,aes(x=midgene_cum,y=prop_var)) +
	geom_vline(xintercept=BP_cum,color="coral") +
	scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
	ylim(0.1,1) +
    #scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
    #geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
	geom_point(aes(color=prop_var)) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
	xlab("Position") + ylab("Proportion Variance") +
	ggtitle(sprintf("Gene Loadings on %s %s, n=%.0f",time,factor,ngenes)) +
	theme_classic() +
    theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_line(colour="black"))
 
	
	pdf(sprintf('eqtl/images/%s_%s_plot.pdf',time,factor))
	print(p)
	dev.off()
	# Plot x-axis location
	# y-axis prop_var
	# vhline of location of factor eQTL

}

#[1] "WD_0727"
#[1] "Factor22"
#[1] "6 Genes loaded on Factor22"

#[1] "WD_0712"
#[1] "Factor1"
#[1] "5030 Genes loaded on Factor1"

#[1] "WD_0718"
#[1] "Factor16"
#[1] "76 Genes loaded on Factor16"

#[1] "WD_0720"
#[1] "Factor9"
#[1] "2725 Genes loaded on Factor9"
#[1] "WD_0727"

#[1] "Factor14"
#[1] "955 Genes loaded on Factor14"

#[1] "WD_0727"
#[1] "Factor2"
#[1] "2103 Genes loaded on Factor2"

df=fread('eqtl/results/all_factor_fdr_peaks.txt',data.table=F)

qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
fqtl=qtl[qtl$Method=="Founder_probs",]
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

df$block_start=all_founder_blocks[match(df$X_ID,all_founder_blocks$focal_snp),]$start
df$block_end=all_founder_blocks[match(df$X_ID,all_founder_blocks$focal_snp),]$end

fqtl$block_start=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$start
fqtl$block_end=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$end


# Overlap of factor eQTL and QTL_F recombination blocks
env1=df
env1=as.data.table(env1)
env2=as.data.table(fqtl)
setkey(env2,Chromosome,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','block_start','block_end'),nomatch=NULL)
#### Factor 2 WD_0727 overlaps with qDTA8 in NERAC_2016_WD

# Overlap of factor eQTL and any QTL support interval
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
# WD_0727 Factor 2 overlaps with qDTA8/qDTS8 [1] "NERAC_2016_WD"     "GRANEROS_2015_OPT" "ALL"              
#[4] "STPAUL_2017_WD"    "BLOIS_2017_OPT"    "SZEGED_2017_OPT"

# Overlap of factor eQTL and any 10% QTL support interval
qtl10=fread('../GridLMM/Biogemma_10p_QTL.csv',data.table=F)
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl10)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
###Factor 2 WD_0727 eQTL hit overlaps with qTPH8 in GRANEROS_2015_OPT

# Overlap of factor eQTL and St.Paul recomb blocks
#DTA
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 4 SNPs for male flowering time in EXP_STPAUL overlap with Factor 2 WD_0727

#Just DTA peaks
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison5=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# None

#DTS
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 4 SNPs for female flowering time in EXP_STPAUL overlap with Factor 2 WD_0727

#Just DTS peaks
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison5=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 1 SNPs for female flowering time in EXP_STPAUL overlap with Factor 2 WD_0727

#tkw_15
qtl2=fread('QTL/tkw_15_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

### Overlap of factor eQTL with cis-eQTL?
ciseqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits.txt',data.table=F)
ciseqtl$block_start=all_founder_blocks[match(ciseqtl$X_ID,all_founder_blocks$focal_snp),]$start
ciseqtl$block_end=all_founder_blocks[match(ciseqtl$X_ID,all_founder_blocks$focal_snp),]$end
env1=df
env1=as.data.table(env1)
env2=as.data.table(ciseqtl)
setkey(env2,CHR,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# No overlap with cis-eQTL

# cis and trans
trans=fread('eqtl/results/all_trans_fdr_hits.txt',data.table=F)
trans$block_start=all_founder_blocks[match(trans$X_ID,all_founder_blocks$focal_snp),]$start
trans$block_end=all_founder_blocks[match(trans$X_ID,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(ciseqtl)
setkey(env2,CHR,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 8 overlaps 
# cis WD_0712 Zm00001d031961 overlaps with trans Zm00001d047592 WD_0727 and Zm00001d035217 WD_0727
# cis WD_0720 Zm00001d025017 overlaps with WD_0720 Zm00001d043918
# cis WD_0718 Zm00001d044657 overlaps with WD_0720 Zm00001d053411, WD_0718 Zm00001d017387, and WD_0718 Zm00001d019724
# cis WD_0718 Zm00001d032099 overlaps with WD_0720 Zm00001d046714 and WD_0712 Zm00001d031961

# With peaks - it is just  one cis-eQTL Zm00001d032099 (WD0718) overlaps with trans-eQTL for Zm00001d046714 (WD0720)

# trans and factor
env1=df
env1=as.data.table(env1)
env2=as.data.table(trans)
setkey(env2,CHR,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 97 overlapping SNPs for trans and factor eQTL 
# Factor 1 WD_0712 overlaps with Zm00001d013230 (cytochrome P450), loaded on the factor

# Factor 2 WD_0727 overlaps with one trans-eQTL from WD_0712, Zm00001d035087, not loaded on the factor

# Factor 14 WD_0727 overlaps 21 trans-eQTL

# Factor 16 WD_0718 overlaps with 3 trans-eqTL
# WD_0727 trans-eQTL Zm00001d023472, Zm00001d042948; WD_0720 Zm00001d016783 - none of which are loaded on the Factor

# Factor 22 in WD0727 overlaps with two trans-eQTL: Zm00001d041712 and Zm00001d009969 - both are loaded on the factor


# cis eQTL and QTL
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(fqtl)
setkey(env2,Chromosome,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','block_start','block_end'),nomatch=NULL)
# None

# Overlap of cis-eQTL and any QTL support interval
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
#None

# Overlap of cis eQTL and any 10% QTL support interval
qtl10=fread('../GridLMM/Biogemma_10p_QTL.csv',data.table=F)
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(qtl10)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
#None

# Overlap of factor eQTL and St.Paul recomb blocks
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

#DTS
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

#tkw_15
qtl2=fread('QTL/tkw_15_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

# trans eQTL and QTL highest SNP recomb block
env1=trans
env1=as.data.table(env1)
env2=as.data.table(fqtl)
setkey(env2,Chromosome,block_start,block_end)
comparison1=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','block_start','block_end'),nomatch=NULL)
# 17 genes and 10 QTL
#  [1] "qDTA3_2" "qDTS3_2" "qTPH7"   "qDTA9"   "qHGM3_1" "qDTS8"   "qDTA8"  
# [8] "qHGM3_2" "qTKW7_1" "qDTS9"

# Overlap of trans-eQTL and any QTL support interval
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
#2504 overlaps with QTL support intervals
# 114 trans eQTL overlap with 20 QTL IDs

# Overlap of trans eQTL and any 10% QTL support interval
qtl10=fread('../GridLMM/Biogemma_10p_QTL.csv',data.table=F)
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl10)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison3=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
#1-29 overlaps with 10% QTL support intervals
# 170 trans-eQTL overlap with 17 QTL IDs

# Overlap of trans eQTL and St.Paul recomb blocks
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison4=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 33 overlaps, 14 genes, 6 SNPs
#qDTA3_2 and WD_0727 genes Zm00001d038732 and WD_0712 Zm00001d053060
#qDTA8_1 WD_0712 Zm00001d027535 (8:127Mb)
#qDTA8_2 (135Mb-137Mb) WD_0712 Zm00001d023869, WD_0720 Zm00001d017387
#qDTA8_2 (142-148Mb) WD_0718 Zm00001d027184, WD_0720 Zm00001d033138
# qDTA8_3 (150-151Mb) WD_0718 Zm00001d035087 and WD_0720 Zm00001d007286
#qDTA9  and WD_0712 Zm00001d046357, WD_0712 Zm00001d034175, WD_0718 Zm00001d028161, WD_0727 Zm00001d026703, and WD_0712 Zm00001d012087

#Just DTA peaks
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison5=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#qDTA3_2 and WD_0727 genes Zm00001d038732 and WD_0712 Zm00001d053060
#qDTA9  and WD_0712 Zm00001d046357, WD_0712 Zm00001d034175, WD_0727 Zm00001d026703, and WD_0712 Zm00001d012087
#qDTA8_2 (142-148Mb) WD_0718 Zm00001d027184, WD_0720 Zm00001d033138


#DTS
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# qDTS8_2 (146Mb) WD_0720 Zm00001d033138 and WD_0718 Zm00001d027184, and WD_0720 Zm00001d007286


#Just DTS peaks
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison5=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

#tkw_15
qtl2=fread('QTL/tkw_15_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# Zm00001d041608 overlaps with qTKW2


# What about just trans peaks, not all SNPs? Do QTL and trans-eQTL have the same peak SNPs?
trans=fread('eqtl/results/all_trans_fdr_peaks.txt',data.table=F)
trans$block_start=all_founder_blocks[match(trans$X_ID,all_founder_blocks$focal_snp),]$start
trans$block_end=all_founder_blocks[match(trans$X_ID,all_founder_blocks$focal_snp),]$end

#DTA
qtl2=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# WD_0712 Zm00001d053060 and qDTA3_2
# WD_0712 Zm00001d046357, WD_0712 Zm00001d034175, WD_0712 Zm00001d012087, WD_0727 Zm00001d026703 and qDTA9

#DTS
qtl2=fread('QTL/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
#None

#TKW
qtl2=fread('QTL/tkw_15_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=trans
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# None



# Calculate correlation between effect sizes for eQTL and phenotypes in overlapping regions
# First do this with ones that overlap for the peak

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


phenotypes=fread('phenotypes/EXP_STPAUL_2017_WD_phenotypes.csv',data.table=F)

# First, DTA and trans-eQTL peaks
pheno="male_flowering_d6"
cors=c()
ps=c()
for(i in 1:nrow(comparison)){
	row=comparison[i,]
	time=row$time
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
	gene=row$Trait
	snp=row$SNP
	chr=row$CHR
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])
	inter=intersect(exp$V1,inds)

	effect_sizes=fread(sprintf('QTL/Biogemma_chr%s_male_flowering_d6_x_EXP_STPAUL_2017_WD_vst_founderprobs.txt',chr),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
	
	results=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results.rds',time,chr))
	w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	results=results[[w]]
	results=results[results$X_ID==snp,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	ps=c(ps,p)
	
	df=data.frame(founder=founders,pheno=effect_size,eqts=betas,stringsAsFactors=F)
	p1=ggplot(df,aes(x=eqts,y=pheno)) + geom_point(aes(color=founder))
	png(sprintf('QTT/%s_%s_%s_%s_effect_sizes.png',time,gene,snp,pheno))
	print(p1)
	dev.off()
	#subpheno=phenotypes[phenotypes$ID %in% inter,c('ID','female_flowering_d6','male_flowering_d6','tkw_15')]
	#fdata$dta=subpheno[match(fdata$ID,subpheno$ID),]$male_flowering_d6
	#fdata$dts=subpheno[match(fdata$ID,subpheno$ID),]$female_flowering_d6

}

comparison$cor=cors
comparison$pvalue=ps

# trans Peaks and DTA
# Highest correlations are  with AX-91145110 for Zm00001d034175 (r=-0.5179,p=0.102) and Zm00001d046357 (r=-0.507,p=0.110)

# all trans hits and DTA




####### Fvalue for Factor 2 and DTA

# Factor 2 and DTS r=-0.2436067, p-value = 0.3816

cors=c()
ps=c()

for(i in 1:nrow(comparison)){
	row=comparison[i,]
	time=row$time
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
	gene=row$factor
	snp=row$SNP
	chr=row$CHR
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])
	inter=intersect(exp$V1,inds)
	
	effect_sizes=fread(sprintf('QTL/Biogemma_chr%s_female_flowering_d6_x_EXP_STPAUL_2017_WD_vst_founderprobs.txt',chr),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,founders])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
		
	results=fread(sprintf('eqtl/trans/results/%s_c%s_%s_trans_results.txt',time,chr,gene),data.table=F)
	#w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	#results=results[[w]]
	results=results[results$X_ID==snp,]
	betas=unlist(results[,founders])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	ps=c(ps,p)
}
comparison$cor=cors
comparison$pvalue=ps

# highest correlation with DTS and Factor 2 trans-eQTL is AX-91202104, r=-0.327, p.value=0.253


# DTA
cors=c()
ps=c()

for(i in 1:nrow(comparison)){
	row=comparison[i,]
	time=row$time
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
	gene=row$factor
	snp=row$SNP
	chr=row$CHR
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])
	inter=intersect(exp$V1,inds)
	
	effect_sizes=fread(sprintf('QTL/Biogemma_chr%s_male_flowering_d6_x_EXP_STPAUL_2017_WD_vst_founderprobs.txt',chr),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,founders])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
		
	results=fread(sprintf('eqtl/trans/results/%s_c%s_%s_trans_results.txt',time,chr,gene),data.table=F)
	#w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	#results=results[[w]]
	results=results[results$X_ID==snp,]
	betas=unlist(results[,founders])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	ps=c(ps,p)
}
comparison$cor=cors
comparison$pvalue=ps

# AX-9110697 r=-0.383, p-value=0.176 is strongest correlation

# How do I get a null expectation?

### Look more closely at WD_0727 Factor 2 3####
library('data.table')
library('ggplot2')

ld=fread('../stats/ld_decay/circos/ld_bundled_links_filtered.txt',data.table=F)
names(ld)=c('CHR_A','START_A','END_A','CHR_B','START_B','END_B','size')
chra=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(ld[x,]$CHR_A,'chr')[[1]][2]))
ld$CHR_A=chra
chrb=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(ld[x,]$CHR_B,'chr')[[1]][2]))
ld$CHR_B=chrb

pmap=c()
for(c in 1:10){
	p=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
	pmap=rbind(pmap,p)
}

snp="AX-91772415"
snp_pos=pmap[pmap$marker==snp,]$pos

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)


df=fread('eqtl/results/factor_transQTL_all.txt',data.table=F)
df=df[df$Factor=='Factor2',]
factor_groups=readRDS('MegaLMM/MegaLMM_WD_0727_factor_groups.rds')
fgenes=factor_groups[['Factor2']]$genes

genetable=genetable[genetable$Gene_ID %in% fgenes,]
 genetable %>% group_by(CHROM) %>% count
# A tibble: 10 × 2
# Groups:   CHROM [10]
#CHROM     n
#   <int> <int>
# 1     1   338
# 2     2   237
# 3     3   228
# 4     4   197
# 5     5   270
# 6     6   167
# 7     7   153
# 8     8   203
# 9     9   161
# 10    10   149


# grab rows where SNPs in factor eQTL are being compared to regions in genetanble

## plot out genes in factor relative to eQTL snp locations

## What regions are in high interchrom ld

## Group inds by founder and plot F-values

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

exp=fread('eqtl/normalized/WD_0727_voom_normalized_gene_counts_formatted.txt',data.table=F)
K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])

inter=intersect(exp$V1,inds)

X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
colnames(X) = founders
rownames(X) = inter




founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))
fdata=data.frame(ID=inter,founder=founder,stringsAsFactors=F)

f_all_means=fread('MegaLMM/MegaLMM_WD_0727_all_F_means.txt',data.table=F)

fdata$fvalue=f_all_means[match(fdata$ID,f_all_means$V1),]$Factor2

p1=ggplot(aes(x=founder,y=fvalue),data=fdata) + geom_boxplot() + geom_jitter() +
theme(axis.text.x = element_text(angle = 45, size=10,hjust=1))
png('images/WD_0727_Factor2_F_by_founder.png')
print(p1)
dev.off()

# Ind EB.10H.H.00025 may be a little bit of an outlier?

# Correlation with flowering time?


phenotypes=fread('phenotypes/EXP_STPAUL_2017_WD_phenotypes.csv',data.table=F)
#phenotypes$Genotype_code=gsub("-",".",phenotypes$Genotype_code)
subpheno=phenotypes[phenotypes$ID %in% inter,c('ID','female_flowering_d6','male_flowering_d6')]
fdata$dta=subpheno[match(fdata$ID,subpheno$ID),]$male_flowering_d6
fdata$dts=subpheno[match(fdata$ID,subpheno$ID),]$female_flowering_d6

cor(fdata$fvalue,fdata$dts)
#[1] -0.01771362
> cor(fdata$fvalue,fdata$dta)
#[1] -0.07511269
highest_SNPs=c("AX-91772402","AX-91772415","AX-91202104","AX-91107495")


results=fread('eqtl/trans/results/WD_0727_c8_Factor2_trans_results.txt',data.table=F)
results$log10p=-log10(results$p_value_ML)
results=results[results$X_ID %in% df$SNP,]
highest=results[which.max(results$log10p),]
betas=unlist(highest[,c(6,8:22)])
betas[-1]=betas[-1]+betas[1]

bvs=X %*% betas
fdata$bv=bvs
cor(fdata$bv,fdata$dta)
 #         [,1]
#[1,] 0.1572119

cor(fdata$dts,fdata$bv)
#          [,1]
#[1,] 0.1186966

pheno="male_flowering_d6"
env="STPAUL_2017_WD"
effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%s_%s_x_%s_founderprobs.rds',chr,pheno,env))
effect_sizes=unlist(unname(effect_sizes[effect_sizes$X_ID==snp,6:21]))
effect_sizes[-1]=effect_sizes[1] + effect_sizes[-1]

ft_bvs=X %*% effect_sizes
fdata$ft_bvs_stpaul=ft_bvs

cor(fdata$bv,fdata$ft_bvs)
# DTA BLUPS
#          [,1]
#[1,] 0.2387707

cor(fdata$bv,fdata$ft_bvs_stpaul)
#         [,1]
#[1,] 0.1786436

fdata=fdata[order(fdata$fvalue),]
rownames(fdata)=seq(1,nrow(fdata))
fdata$variable_f=factor(fdata$founder,levels=c(founders))

phenotypes=fread('phenotypes/phenotypes_BLUPS.csv',data.table=F)
fdata$dta_blup=phenotypes[match(fdata$ID,phenotypes$ID),'BLUP-male_flowering_d6']
cor(fdata$fvalue,fdata$dta_blup)
#[1] -0.0713566

phenotypes=fread('phenotypes/phenotypes_asi.csv',data.table=F)
phenotypes$Genotype_code=gsub("-",".",phenotypes$Genotype_code)
envs=c("BLOIS_2014_OPT","BLOIS_2017_OPT","GRANEROS_2015_OPT","NERAC_2016_WD","STPAUL_2017_WD","SZEGED_2017_OPT")
for(env in envs){
	subpheno=phenotypes[phenotypes$Loc.Year.Treat==env,]
	print(env)
	phenos=unlist(subpheno[match(fdata$ID,phenotypes$Genotype_code),]$male_flowering_d6)
	print(cor(fdata$fvalue,phenos,use='complete.obs'))
}
#[1] "BLOIS_2014_OPT"
#[1] -0.02306991
#[1] "BLOIS_2017_OPT"
#[1] -0.1415513
#[1] "GRANEROS_2015_OPT"
#[1] -0.02039115
#[1] "NERAC_2016_WD"
#[1] -0.09728809
#[1] "STPAUL_2017_WD"
#[1] 0.0223498
#[1] "SZEGED_2017_OPT"
#[1] -0.06438033
cor(fdata$fvalue,fdata$dta_blup)

colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder

p1=ggplot(aes(x=fvalue,y=ft_bvs),data=fdata)+geom_point(aes(color=founder)) +
scale_color_manual(values=colorcodes[levels(fdata$variable_f),]$hex_color,labels=levels(fdata$variable_f))+


png('images/WD_0727_Factor2_eQTL_by_DTA_BLUP.png')
print(p1)
dev.off()


p2=ggplot(aes(x=fvalue,y=dta),data=fdata)+geom_point(aes(color=founder)) +
scale_color_manual(values=colorcodes[levels(fdata$variable_f),]$hex_color,labels=levels(fdata$variable_f))+

png('images/WD_0727_Factor2_eQTL_by_DTA.png')
print(p2)
dev.off()

p3=ggplot(aes(x=bv,y=ft_bvs),data=fdata)+geom_point(aes(color=founder)) +
scale_color_manual(values=colorcodes[levels(fdata$variable_f),]$hex_color,labels=levels(fdata$variable_f))+
xlab("F Genetic Values") + ylab("DTA BLUP Breeding Values") + ggtitle(sprintf('Flowering Time by Factor 2 F for %s',snp))

png('images/WD_0727_Factor2_eQTL_Fbv_by_DTA_BLUP.png')
print(p3)
dev.off()


### What SNPs are nearest to the genes loading onto this factor
pmap$end=pmap$pos+1

env1=genetable
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(pmap)
#env2$end=env2$end-1
setkey(env2,chr,pos,end)
comparison=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('chr','pos','end'),nomatch=NULL)

# 1731 have SNPs within them

find_nearest_snp=function(row){
    index=which.min(abs(row$START-pmap[pmap$chr==row$CHROM,]$pos))
    return(pmap[index,]$marker)
}

nearest_snps=sapply(seq(1,nrow(genetable)),function(x) find_nearest_snp(genetable[x,]))
genetable$SNP=nearest_snps

#Instead of doing this, I need to find the founder recombination blocks that have those SNPs
# And use those instead of the nearest SNP
all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

env1=genetable
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(all_founder_blocks)
#env2$end=env2$end-1
setkey(env2,chr,start,end)
comparison=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('chr','start','end'),nomatch=NULL)

## Calculate breeding value of DTA using SNPs near the genes in Factor 2

bvs=c()

env="ALL"
pheno="male_flowering_d6"
for(chr in 1:10){
	subtable=genetable[genetable$CHROM==chr,]
	subcomp=comparison[comparison$CHROM==chr,]
	snps=unique(subcomp$focal_snp)
	
	effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',chr,pheno,env))
	
	es=effect_sizes[effect_sizes$X_ID %in% snps,]
	

	#effect_sizes[-1]=effect_sizes[1] + effect_sizes[-1]
	
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])

	inter=intersect(exp$V1,inds)

	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
}


#########


prop_var=fread('MegaLMM/MegaLMM_WD_0727_prop_variance.txt',data.table=F)
# How much variation is enough? 
subvar=prop_var[,c('V1','Factor2')]
cutoff=0.9
test=subvar[subvar$Factor2>=cutoff,]$V1

testlist=data.frame(genes=test,stringsAsFactors=F)
fwrite(testlist,sprintf('WD_0727_Factor2_%.1f_gene_list.txt',cutoff),row.names=F,col.names=F,quote=F,sep='\t')


######### GOSeq ########
pheno_factors=c('Factor2','Factor14')

#pheno_factors=unique(pheno_df$factor)
genes=names(avg_exp)
inter2=length(intersect(names(avg_exp),annotation$Gene))
#genelength=genetable[match(genes,genetable$Gene_ID),]$LENGTH
#names(genelength)=genes

fulllist=data.frame(genes=genes,stringsAsFactors=F)
fwrite(fulllist,'WD_0727_full_gene_list.txt',row.names=F,col.names=F,quote=F,sep='\t')


ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
# Enrichment of flowering time genes?
cutoff=0.1
test=subvar[subvar$Factor2>=cutoff,]$V1
gtable=genetable[genetable$Gene_ID %in% test,]


find_nearest_snp=function(row){
    index=which.min(abs(row$START-pmap[pmap$chr==row$CHROM,]$pos))
    return(pmap[index,]$marker)
}

nearest_snps=sapply(seq(1,nrow(gtable)),function(x) find_nearest_snp(gtable[x,]))
gtable$SNP=nearest_snps

snplist=gtable[,'SNP',drop=F]
fwrite(snplist,'WD_0727_Factor2_snplist.txt',row.names=F,col.names=F,sep='\t',quote=F)

####


snplist=fread('WD_0727_Factor2_snplist.txt',data.table=F,header=F)
snplist=snplist[-1,]
snp="AX-91772415"

lambda=fread('MegaLMM/MegaLMM_WD_0727_all_Lambda_means.txt',data.table=F)
ld=fread('Factor2_rsquared.ld',data.table=F)
ld=ld[ld$SNP_A==snp,]

ld2=ld[ld$SNP_B %in% snplist,]


##### WD_0727 trans-eQTL

eqtl=fread('eqtl/results/WD_0727_trans_eQTL_scan_hits.txt',data.table=F)

pmap=c()
for(c in 1:10){
	p=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
	pmap=rbind(pmap,p)
}

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

eqtl$Gene_CHR=genetable[match(eqtl$Gene,genetable$Gene_ID),]$CHROM
eqtl$Gene_START=genetable[match(eqtl$Gene,genetable$Gene_ID),]$START
eqtl$Gene_END=genetable[match(eqtl$Gene,genetable$Gene_ID),]$END

#[1] "Zm00001d041650" "Zm00001d020311" "Zm00001d009688" "Zm00001d045677"
#[5] "Zm00001d047592"

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

eqtl$block_start=all_founder_blocks[match(eqtl$SNP,all_founder_blocks$focal_snp),]$start
eqtl$block_end=all_founder_blocks[match(eqtl$SNP,all_founder_blocks$focal_snp),]$end

fqtl$block_start=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$start
fqtl$block_end=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$end


env1=eqtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(qtl)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

f14_27=fread('eqtl/results/Factor14_trans_WD_0727_eQTL_fkeep_hits.txt',data.table=F)

goi=c("Zm00001d041650","Zm00001d009688","Zm00001d045677")
factor_groups=readRDS('MegaLMM/MegaLMM_WD_0727_factor_groups2.rds')
f2genes=factor_groups[['Factor2']]$genes
f14genes=factor_groups[['Factor14']]$genes

# nearest SNPs for genes in Factor 14

find_nearest_snp=function(row){
    index=which.min(abs(row$START-pmap[pmap$chr==row$CHROM,]$pos))
    return(pmap[index,]$marker)
}

gtable=genetable[genetable$Gene_ID %in% f14genes,]

nearest_snps=sapply(seq(1,nrow(gtable)),function(x) find_nearest_snp(gtable[x,]))
gtable$SNP=nearest_snps

snplist=gtable[,'SNP',drop=F]
snplist=rbind(snplist,data.frame(SNP=unique(f14_27$SNP)))
snp2=unique(snplist$SNP)
snp2=data.frame(SNP=snp2,stringsAsFactors=F)

fwrite(snp2,'WD_0727_Factor14_snplist.txt',row.names=F,col.names=F,sep='\t',quote=F)

ld=fread('Factor14_rsquared.ld',data.table=F)
ld=ld[ld$SNP_A==snp,]

ld2=ld[ld$SNP_B %in% snplist$V1,]


#####
# Correlation between Fvalues and phenotypes
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
rownames(exp)=exp$V1
exp=exp[,-1]
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)

metadata=metadata[metadata$experiment==time,]


geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
exp=exp[,kept_genes]

genes=names(exp)

#genes=genes[1:5]
df=fread('eqtl/results/factor_transQTL_all.txt',data.table=F)
df %>% group_by(time,Factor) %>% filter(value==max(value)) %>% arrange(time,Factor,SNP)

# A tibble: 10 × 6
# Groups:   time, Factor [3]
#   Factor     CHR        BP SNP           value time   
#   <chr>    <int>     <int> <chr>         <dbl> <chr>  
# 1 Factor16     4   5535185 AX-90856708    2.90 WD_0718
# 2 Factor16     4   5414982 AX-91597662    2.90 WD_0718
# 3 Factor16     4   5663654 AX-91849936    2.90 WD_0718
# 4 Factor14     7 171058661 AX-91065358    3.95 WD_0727
# 5 Factor14     7 170989903 AX-91743057    3.95 WD_0727
# 6 Factor14     7 170886882 PZE-107118743  3.95 WD_0727
# 7 Factor2      8 152630937 AX-91107495    3.43 WD_0727
# 8 Factor2      8 150428041 AX-91202104    3.43 WD_0727
# 9 Factor2      8 150234596 AX-91772402    3.43 WD_0727
#10 Factor2      8 150347882 AX-91772415    3.43 WD_0727

time="WD_0727"
factor="Factor2"
snp="AX-91772402"
chr="8"

time="WD_0727"
factor="Factor14"
snp="AX-91743057"
chr="7"

time="WD_0718"
factor="Factor16"
snp="AX-91597662"
chr="4"

f_all_means=fread(sprintf('MegaLMM/MegaLMM_%s_all_F_means.txt',time),data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
inds=rownames(X_list[[1]])

inter=intersect(exp$V1,inds)

X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
colnames(X) = founders
rownames(X) = inter


founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))
fdata=data.frame(ID=inter,founder=founder,stringsAsFactors=F)


fdata$fvalue=f_all_means[match(fdata$ID,f_all_means$V1),factor]


find_qtts=function(e,p){
	pheno=phenotypes[phenotypes$Loc.Year.Treat==e,c('Genotype_code',p)]
	pheno$f=fdata[match(pheno$Genotype_code,fdata$ID),'fvalue']
	test=cor.test(pheno[,p],pheno$f,use="complete.obs")
	return(test)
}

envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:8]

all_qtts=data.frame(matrix(ncol=5,nrow=0))
names(all_qtts)=c('pheno','env','factor','r','pvalue')

for(p in phenos){
	for(e in envs){
		correlations=find_qtts(e,p)
		line=data.frame(pheno=p,env=e,factor=factor,r=unlist(correlations['estimate']),pvalue=unlist(correlations['p.value']),stringsAsFactors=F)
		all_qtts=rbind(all_qtts,line)
	}
}

ntests=48
threshold=0.05/3
all_qtts$padjust=p.adjust(all_qtts$pvalue,method='fdr')

# WD_0727 Factor 2
### total plant height was actually th emost correlated with f-values
# r of -0.17
# none of them are significant though
# Is there a way I can compare this to a null distribution rather than fdr correction?

# WD_0727 Factor 14
#### highest correlation is with male_flowering_d6 (0.15) and tkw_15 (-0.15)
## neither of them are significant after multiple test correction

# WD_0718 Factor16
#### highest correlation is with grain_yield_15 (r=0.18) and total_plant_height (r=0.16) in SZEGED_2017_OPT
# Neither are significant after fdr correction

# Raw F values and raw phenotypes are not significantly correlated
# What about phenotype genetic values?


fwrite(all_qtts,sprintf('eqtl/results/%s_%s_F_pheno_corrs.txt',time,factor),row.names=F,quote=F,sep='\t')


#### Do any factors have enrichment of FT genes?
time="WD_0727"
prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance.txt',time),data.table=F)
ft_genes=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
names(ft_genes)=c('CHR','START','END','GENE_ID')

allft=c()

factors=names(prop_var)[-1]
for(factor in factors){
	fgenes=prop_var[prop_var[,factor]>=0.1,]$V1
	tn=sum(fgenes %in% ft_genes$GENE_ID)
	ft_geneoverlap=data.frame(factor=factor,ngenes=length(fgenes),inft=tn,stringsAsFactors=F)
	allft=rbind(allft,ft_geneoverlap)
}

allgenes=prop_var$V1

null=c()
for(i in 1:nrow(allft)){
	row=allft[i,]
	ngenes=row$ngenes
	tgenes=row$inft
	alltn=c()
	for(j in 1:1000){
		draw=sample(allgenes,ngenes)
		tn=sum(draw %in% ft_genes$GENE_ID)
		alltn=c(alltn,tn)
	}
	adjust=68
	q5=quantile(alltn,1-(0.05/adjust))
	null=c(null,q5)

}
allft$null=null

ft_geneoverlap=c()
for(i in 1:length(factor_groups)){
  #inds=which(pheno_factors==f)
  #i=pheno_loc[inds]
  fgenes=factor_groups[[i]]$genes
  f=factor_groups[[i]]$factor
  ngenes=length(fgenes)
  match=intersect(fgenes,ft_genes$GENE_ID)
  if(f %in% ft_factors){
    inft=T
  }else{
    inft=F
  }
  ft_geneoverlap=rbind(ft_geneoverlap,c(f,length(match),ngenes,inft))
}
ft_geneoverlap=as.data.frame(ft_geneoverlap,stringsAsFactors=F)
names(ft_geneoverlap)=c('factor','overlap','ngenes','inft')
ft_geneoverlap$overlap=as.numeric(ft_geneoverlap$overlap)
ft_geneoverlap$ngenes=as.numeric(ft_geneoverlap$ngenes)
ft_geneoverlap$inft=factor(ft_geneoverlap$inft,levels=c('FALSE','TRUE'))
m1=lm(overlap~ngenes+inft,ft_geneoverlap)
anova(m1)



### What do Fvalues look like for Factor 14
time="WD_0727"
fvalues=fread('MegaLMM/MegaLMM_WD_0727_all_F_means.txt',data.table=F)
snp="PZE-107118743"
samples=fread(sprintf('eqtl/data/%s_samples.txt',time),data.table=F,header=F)
chr="7"

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder

#pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
X_list=lapply(X_list,function(x) x[samples$V1,])
X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
#frep2=apply(X,MARGIN=2,function(x) round(sum(x[x>0.75])))
#fkeep=founders[frep2>3]

certain=apply(X,MARGIN=1,function(x) sum(x>0.75)>0)
X=X[certain,]
founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))
fvalues=fvalues[,c('V1','Factor14')]
rownames(fvalues)=fvalues$V1
fvalues=fvalues[certain,]
fvalues$founder=founder
fvalues=fvalues[order(fvalues$Factor14),]
fvalues$ind_f=factor(fvalues$V1,levels=c(fvalues$V1))
fvalues$founder_f=factor(fvalues$founder,levels=c(founders))



p1=ggplot(aes(x=ind_f,y=Factor14),data=fvalues)+geom_point(aes(color=founder_f)) +
scale_color_manual(values=colorcodes[levels(fvalues$founder_f),]$hex_color,labels=levels(fvalues$founder_f))


png('images/WD_0727_Factor14_Fvalue_by_ind.png')
print(p1)
dev.off()