#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')
library('ggnewscale')


#cis=fread('eqtl/results/WD_0712_cis_eQTL_weights_fdr_hits.txt',data.table=F)
#cis$class="cis"
#cis$time="WD_0712"
#cis2=fread('eqtl/results/WD_0720_cis_eQTL_weights_fdr_hits.txt',data.table=F)
#cis2$class="cis"
#cis2$time="WD_0720"
#cis=rbind(cis,cis2)
#names(cis)=c("Trait","SNP","p_value_ML","CHR","BP","value","p_adjusted","class","time")
#cis=cis[,c("Trait","CHR","BP","SNP","value","class","time")]

#fwrite(cis,'eqtl/results/all_cis_eQTL_weights_fdr_hits.txt',row.names=F,quote=F,sep='\t')

#trans=fread('eqtl/results/WD_0727_trans_eQTL_weights_fdr_hits.txt',data.table=F)
#trans$class="trans"
#trans$time="WD_0727"

#trans2=fread('eqtl/results/WD_0720_trans_eQTL_weights_fdr_hits.txt',data.table=F)
#trans2$class="trans"
#trans2$time="WD_0720"
#trans=rbind(trans,trans2)

#trans3=fread('eqtl/results/WD_0718_trans_eQTL_weights_fdr_hits.txt',data.table=F)
#trans3$class="trans"
#trans3$time="WD_0718"
#trans=rbind(trans,trans3)

#trans4=fread('eqtl/results/WD_0712_trans_eQTL_weights_fdr_hits.txt',data.table=F)
#trans4$class="trans"
#trans4$time="WD_0712"
#trans=rbind(trans,trans4)

#names(trans)=c("Trait","CHR","BP","SNP","p_adjusted","value","class","time")
#trans=trans[,c("Trait","CHR","BP","SNP","value","class","time")]
#fwrite(trans,'eqtl/results/all_trans_eQTL_weights_fdr_hits.txt',row.names=F,quote=F,sep='\t')
cis=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits.txt',data.table=F)
trans=fread('eqtl/results/all_trans_fdr_peaks.txt',data.table=F)

factoreqtl=fread('eqtl/results/Factor2_trans_WD_0727_eQTL_fkeep_vst_fdr_hits.txt',data.table=F)
factoreqtl$class="factor"
factoreqtl$time="WD_0727"

factoreqtl2=fread('eqtl/results/Factor14_trans_WD_0727_eQTL_fkeep_vst_fdr_hits.txt',data.table=F)
factoreqtl2$class="factor"
factoreqtl2$time="WD_0727"
factoreqtl=rbind(factoreqtl,factoreqtl2)

factoreqtl2=fread('eqtl/results/Factor16_trans_WD_0718_eQTL_fkeep_vst_fdr_hits.txt',data.table=F)
factoreqtl2$class="factor"
factoreqtl2$time="WD_0718"
factoreqtl=rbind(factoreqtl,factoreqtl2)

names(factoreqtl)=c("Trait","CHR","BP","SNP","value","class","time")
factoreqtl=factoreqtl[,c("Trait","CHR","BP","SNP","value","class","time")]
fwrite(factoreqtl,'eqtl/results/all_factor_eQTL_fdr_hits.txt',row.names=F,quote=F,sep='\t')

allhits=rbind(cis,trans)
allhits=rbind(allhits,factoreqtl)

fwrite(allhits,'eqtl/results/all_eQTL_fdr_hits.txt',row.names=F,quote=F,sep='\t')

allhits=fread('eqtl/results/all_eQTL_fdr_hits.txt',data.table=F)

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)


allhits$block_start=all_founder_blocks[match(allhits$SNP,all_founder_blocks$focal_snp),]$start
allhits$block_end=all_founder_blocks[match(allhits$SNP,all_founder_blocks$focal_snp),]$end


f27_groups=readRDS('MegaLMM/MegaLMM_WD_0727_factor_groups.rds')
f2genes=f27_groups[['Factor2']]$genes
f14genes=f27_groups[['Factor14']]$genes

prop27=fread('MegaLMM/MegaLMM_WD_0727_prop_variance.txt',data.table=F)
f2genes=prop27[prop27$Factor2>=0.10,]$V1
f14genes=prop27[prop27$Factor14>=0.10,]$V1


f18_groups=readRDS('MegaLMM/MegaLMM_WD_0718_factor_groups.rds')
f16genes=f18_groups[['Factor16']]$genes

prop18=fread('MegaLMM/MegaLMM_WD_0718_prop_variance.txt',data.table=F)
f16genes=prop18[prop18$Factor16>=0.10,]$V1

igenes=intersect(f2genes,f14genes)
length(igenes)
# 368
igenes=intersect(igenes,f16genes)
length(igenes)
# 1 "Zm00001d028809"

ugenes=union(f2genes,f14genes)
ugenes=union(ugenes,f16genes)
gtable=genetable[genetable$Gene_ID %in% c(f2genes,f14genes,f16genes),]

### make table of genes grouped by factor and location of SNP and gene
sub2=factoreqtl[factoreqtl$Trait=="Factor2",]
gtable=genetable[genetable$Gene_ID %in% f2genes,]

factordf=c()
for(i in 1:nrow(gtable)){
	subdf=sub2
	row=gtable[i,]
	subdf$Gene=row$Gene_ID
	subdf$gene_chr=row$CHROM
	subdf$gene_start=row$START
	subdf$gene_end=row$END
	factordf=rbind(factordf,subdf)
}

factordf$prop_var=prop27[match(factordf$Gene,prop27$V1),]$Factor2

sub2=factoreqtl[factoreqtl$Trait=="Factor14",]
gtable=genetable[genetable$Gene_ID %in% f14genes,]
for(i in 1:nrow(gtable)){
	subdf=sub2
	row=gtable[i,]
	subdf$Gene=row$Gene_ID
	subdf$gene_chr=row$CHROM
	subdf$gene_start=row$START
	subdf$gene_end=row$END
	subdf$prop_var=prop27[match(subdf$Gene,prop27$V1),]$Factor14
	factordf=rbind(factordf,subdf)
}

sub2=factoreqtl[factoreqtl$Trait=="Factor16",]
gtable=genetable[genetable$Gene_ID %in% f16genes,]

for(i in 1:nrow(gtable)){
	subdf=sub2
	row=gtable[i,]
	subdf$Gene=row$Gene_ID
	subdf$gene_chr=row$CHROM
	subdf$gene_start=row$START
	subdf$gene_end=row$END
	subdf$prop_var=prop18[match(subdf$Gene,prop18$V1),]$Factor16
	factordf=rbind(factordf,subdf)
}

fwrite(factordf,'eqtl/results/all_factor_trans_eqtl_fdr_genes.txt',row.names=F,quote=F,sep='\t')

factordf=fread('eqtl/results/all_factor_trans_eqtl_fdr_genes.txt',data.table=F)

pmap=c()
for(c in 1:10){
	p=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
	pmap=rbind(pmap,p)
}
chr_ends=pmap %>% group_by(chr) %>% summarize(chr_len=max(pos))


factordf$block_start=all_founder_blocks[match(factordf$SNP,all_founder_blocks$focal_snp),]$start
factordf$block_end=all_founder_blocks[match(factordf$SNP,all_founder_blocks$focal_snp),]$end

factordf$mid_block=apply(factordf[,c('block_start','block_end')],1,function(x) round(mean(x)))
factordf$mid_gene=apply(factordf[,c('gene_start','gene_end')],1,function(x) round(mean(x)))

factordf$chr_lenA=chr_ends[match(factordf$CHR,chr_ends$chr),]$chr_len
factordf$chr_lenB=chr_ends[match(factordf$gene_chr,chr_ends$chr),]$chr_len


# Factor eQTL - get chromosome locations of midgene  
df.tmp2 <- factordf %>%

    # Compute chromosome size
    group_by(gene_chr) %>%
    summarise(chr_len=max(chr_lenB)) %>%

    # Calculate cumulative position of each chromosome
    mutate(totB=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(factordf, ., by=c("gene_chr"="gene_chr")) %>%

    # Add a cumulative position of each SNP
    arrange(gene_chr, mid_gene) %>%
    mutate(midgene_cum=as.numeric(mid_gene+totB))
    
# Factor eQTL - get chromosome locations of mid-block
df.tmp = df.tmp2 %>%
	# Compute chromosome size
    group_by(gene_chr) %>%
    summarise(chr_len=max(chr_lenB)) %>%

    # Calculate cumulative position of each chromosome
    mutate(totA=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(df.tmp2, ., by=c("CHR"="gene_chr")) %>%
	arrange(CHR, mid_block) %>%
    mutate(midblock_cum=as.numeric(mid_block+totA))


axisdf <- factordf %>%
	# Compute chromosome size
    group_by(gene_chr) %>%
    summarise(chr_len=max(chr_lenB)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len,chr_len=chr_len) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>% mutate(center=(chr_len/2)+tot)


cumtot=2105119857

p2=ggplot(data=df.tmp, aes(x=midblock_cum, y=midgene_cum)) +
	geom_point(aes(color=prop_var,shape=Trait),,size=1,alpha=0.4) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
    scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    labs(x = "eQTL variants by Chromosome",y="eQTL transcripts by Chromosome") #+ 
    #guides(shape="none")

pdf('eqtl/images/factor_0.1_eQTL_plot.pdf')
print(p2)
dev.off()

#df.tmp$Trait=paste0(df.tmp$time,'-',df.tmp$Trait)

#for(i in unique(df.tmp$Trait)){
#	subtmp=df.tmp[df.tmp$Trait==i,]
#	p2=ggplot(data=subtmp, aes(x=midblock_cum, y=midgene_cum)) +
#	geom_point(aes(color=prop_var),size=1,alpha=0.4) +
#	scale_color_gradient(low = "lightblue", high = "darkblue") +
#   scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
#  scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
#    labs(x = "eQTL variants by Chromosome",y="eQTL transcripts by Chromosome") + 
#    ggtitle(sprintf("%s eQTL",i))
#
#	pdf(sprintf('eqtl/images/%s_eQTL_plot.pdf',i))
#	print(p2)
#	dev.off()
#}




#p2=ggplot(df.tmp,aes(x=midblock_cum,y=midgene_cum)) + geom_point(color=Trait)

### just cis and trans plot first

sub=allhits[allhits$class %in% c('cis','trans'),]
sub$gene_chr=genetable[match(sub$Trait,genetable$Gene_ID),]$CHROM
sub$gene_start=genetable[match(sub$Trait,genetable$Gene_ID),]$START
sub$gene_end=genetable[match(sub$Trait,genetable$Gene_ID),]$END

sub$mid_block=apply(sub[,c('block_start','block_end')],1,function(x) round(mean(x)))

sub$mid_gene=apply(sub[,c('gene_start','gene_end')],1,function(x) round(mean(x)))



sub$chr_lenA=chr_ends[match(sub$CHR,chr_ends$chr),]$chr_len
sub$chr_lenB=chr_ends[match(sub$gene_chr,chr_ends$chr),]$chr_len

df.tmp3 <- axisdf %>%
    left_join(sub, ., by=c("CHR"="gene_chr")) %>%
    arrange(CHR, mid_block) %>%
    mutate(midblock_cum=as.numeric(mid_block+tot))
    
    
df.tmp4 <- axisdf %>%
    left_join(df.tmp3, ., by=c("gene_chr"="gene_chr")) %>%
    arrange(gene_chr, mid_gene) %>%
    mutate(midgene_cum=as.numeric(mid_gene+tot.y))


cumtot=2105119857

##### Paper figure

cistrans=ggplot(df.tmp4,aes(x=midblock_cum,y=midgene_cum)) +
	geom_abline(slope=1, intercept=0,alpha=0.5) +
	scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
    scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
    geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
    geom_point(aes(color=time),size=1) + 
    labs(x = "Position of eQTL variants",y="Position of eQTL Transcripts") +
    theme_classic() +
    theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_blank())
 
pdf('eqtl/images/cis_trans.pdf')
print(cistrans)
dev.off()

output=df.tmp4[,c('Trait','SNP','CHR','block_start','block_end','gene_chr','gene_start','gene_end','class','time')]
fwrite(output,'eqtl/results/all_cis_trans_fdr_hits.txt',row.names=F,quote=F,sep='\t')
times=c("WD_0712","WD_0718","WD_0720","WD_0727")

for(time in times){
	subtmp=df.tmp4[df.tmp4$time==time,]
	cistrans=ggplot(subtmp,aes(x=midblock_cum,y=midgene_cum)) +
	geom_abline(slope=1, intercept=0,alpha=0.5) +
	scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
    scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
    geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
 	geom_point(aes(color=class),size=1) + 
 	scale_color_manual(values = c("cis" = "red","trans" = "blue")) + 
    labs(x = "Position of eQTL variants",y="Position of eQTL Transcripts") +
    theme_classic() +
    theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_line(colour="black"))
 
	pdf(sprintf('eqtl/images/%s_cis_trans_fdr.pdf',time))
	print(cistrans)
	dev.off()

}

#######

mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
greypalette=gray.colors(6)

p3=p2 + new_scale_color() +
 geom_point(data=df.tmp4,aes(x=midblock_cum,y=midgene_cum,color=class,shape=time)) + 
pdf('eqtl/images/all_0.1_eqtl_plot.pdf')
print(p3)
dev.off()


cumtot=1377773285

p1=ggplot(data=df.tmp, aes(x=midblock_cum, y=midgene_cum)) +
	geom_point(aes(color=class,shape=time)) +
    scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    labs(x = "eQTL variants by Chromosome",y="eQTL transcripts by Chromosome")

png('eqtl/images/eQTL_plot.png',width=1000,height=1000)
print(p1)
dev.off()



p1=ggplot(aes(x=mid_block,y=mid_gene),data=sub) + geom_point(aes(color=class,shape=time)) + xlab("eQTL variant") + ylab("eQTL transcript")

png('eqtl/images/eQTL_plot.png')
print(p1)
dev.off()


#### Are there cis-eQTL in the factors?  NO
cisgenes=cis$Trait
cisgenes %in% factordf$Gene
#[1] FALSE FALSE

### Are there trans-eQTL in the factors? YES
transgenes=unique(trans$Trait)
transgenes[transgenes %in% factordf$Gene]
#[1] "Zm00001d035217" "Zm00001d047592" "Zm00001d022126"
tgoi=transgenes[transgenes %in% factordf$Gene]
factordf[factordf$Gene %in% tgoi,]

## Zm00001d022126 trans-eQTL in WD_0718
# Variant on chr 1, gene on chr7
# Loaded on Factor2, which has an eQTL on chr8, prop_var is 0.267.
#A Total of 4% percent of variation explained by all factors - how does this compare to h2 (0.010625 var explained)
# h2 Zm00001d022126 0.1078221 WD0727 (max h2 of 0.324 in WD_0712)
# Lambda value for Factor2 is -0.103
# also loaded on Factor 1(28%), and Factor 6 (19%)


## Zm00001d035217 trans-eQTL in WD_0727
# Variant on chr 5, gene on chr 6
# Loaded on Factor 14, which has an eQTL on chr7, prop_var is 0.838
#A Total of 36% percent of variation explained by all factors - how does this compare to h2 (0.3040358 var explained)
# Lambda value for Factor14 is -0.551
# h2 Zm00001d035217 0.1018871 WD_0727 (max h2 of 0.45 in WD_0718)
# Only loaded on Factor 14

## Zm00001d047592 trans-eQTL in WD_0727
# Variants on chr 5, gene on chr 9
# Loaded on Factor 14, which has an eQTL on chr 7, prop_var is 0.471
#A Total of 14% percent of variation explained by all factors - how does this compare to h2 (0.06674546 var explained)
# h2 Zm00001d047592 0.1641610
# Lambda value for Factor 14 is -0.258 (max h2 of 0.345 in WD_0718)
# also loaded on Factor 7 (36.5%)
# also loaded on Factor 7 lambda is -.227

# Zm00001d039372 trans-eQTL in WD_0712
# Variants on chr 5, gene on chr 3
# Loaded on Factor 2, which has an eQTL on chr8, prop_var is 0.112.
# A total of 8% percent of variation explained by all factors - how does this compare to h2
# h2 WD_0712 Zm00001d039372 0.4668825
# h2 WD_072 7Zm00001d039372 0.0769847
# Lambda value for Factor 2 is -0.09681375 (max h2 of 0.466 in WD_0712)
# also loaded on Factor 6 (14%) and 7 (25%)

#### Are trans-eQTL genes in inter-chromosomal LD with their variants?


### Are genes in factors in inter-chromosomal LD with themselves?



### Do any of these genes have correlated expression? Cis and trans-genes?


### Are any of these genes in recombination blocks for FT qtl?


### For Factor 2, it has a factor eQTL on chr8.
# 1 - what other regions on chr8 are in high LD? Are any of them near vgt1?

# 2 - what other regions are in inter-chromosomal LD? How does R2 correlate with

#interld=ld[ld$CHR_A!=ld$CHR_B,]
#f2=factordf[factordf$Trait=="Factor2",]

env1=interld
env1=as.data.table(env1)
env2=as.data.table(f2)
setkey(env2,gene_chr,gene_start,gene_end)
comparisonA=foverlaps(env1,env2,by.x=c('CHR_A','BP_A','BP_A_END'),by.y=c('gene_chr','gene_start','gene_end'),nomatch=NULL)
# It seems like 7 genes on chromosome 1 are in high inter-chromosomal LD with 

env1=interld
env1=as.data.table(env1)
env2=as.data.table(f2)
setkey(env2,gene_chr,gene_start,gene_end)
comparisonB=foverlaps(env1,env2,by.x=c('CHR_B','BP_B','BP_B_END'),by.y=c('gene_chr','gene_start','gene_end'),nomatch=NULL)

pmap=c()
for(c in 1:10){
	p=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
	pmap=rbind(pmap,p)
}

chr_ends=pmap %>% group_by(chr) %>% summarize(chr_len=max(pos))
axisdf <- chr_ends %>%
	# Compute chromosome size
    #group_by(chr) %>%
    mutate(tot=cumsum(chr_len)-chr_len,chr_len=chr_len) %>%
    mutate(tot=cumsum(chr_len)-chr_len) %>% mutate(center=(chr_len/2)+tot)
cumtot=2105119857

find_nearest_snp=function(row){
	pc=pmap[pmap$chr==row$gene_chr,]
    index=which.min(abs(row$gene_start-pc$pos))
    return(pc[index,]$marker)
}

factordf=fread('eqtl/results/all_factor_trans_eqtl_genes.txt',data.table=F)
factordf$nearest_SNP=sapply(seq(1,nrow(factordf)),function(x) find_nearest_snp(factordf[x,]))

trans=fread('eqtl/results/all_trans_eQTL_weights_hits.txt',data.table=F)
#Factor 2
ld=fread('MegaLMM/Factor2_rsquared.ld',data.table=F)
ld=ld[ld$R2>=0.5,]

f2=factordf[factordf$Trait=="Factor2",]

df=data.frame(snp=unique(f2$nearest_SNP))
if(!unique(f2$SNP) %in% df$snp){
	df=rbind(df,unique(f2$SNP))

}
fwrite(df,'MegaLMM/WD_0727_Factor2_snplist.txt',row.names=F,quote=F,sep=',',col.names=F)

ld$Gene_A=f2[match(ld$SNP_A,f2$nearest_SNP),]$Gene
ld=ld[!is.na(ld$Gene_A),]
ld$Gene_B=f2[match(ld$SNP_B,f2$nearest_SNP),]$Gene
ld=ld[!is.na(ld$Gene_A) & !is.na(ld$Gene_B),]

ld$chr_lenA=axisdf[match(ld$CHR_A,axisdf$chr),]$chr_len
ld$chr_lenB=axisdf[match(ld$CHR_B,axisdf$chr),]$chr_len
ld$totA=axisdf[match(ld$CHR_A,axisdf$chr),]$tot
ld$totB=axisdf[match(ld$CHR_B,axisdf$chr),]$tot
ld$BP_A_cum = ld$BP_A+ld$totA
ld$BP_B_cum = ld$BP_B+ld$totB

for(c in 1:10){
	subld=ld[ld$CHR_A==c & ld$CHR_B==c,]
	p1=ggplot(subld,aes(x=BP_A/1e6,y=BP_B/1e6)) + geom_point(aes(color=R2)) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
	#scale_x_continuous(label = axisdf$chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    #scale_y_continuous(label = axisdf$chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
	xlab('Chromosome (Mb)') + ylab("Chromosome (Mb)") +
	ggtitle(sprintf("Factor 2 Genes in High LD on Chromosome %.0f",c))

	pdf(sprintf('images/Factor2_c%.0f_ld.pdf',c))
	print(p1)
	dev.off()
}


ld=ld[ld$SNP_A != ld$SNP_B,]

p1=ggplot(ld,aes(x=BP_A_cum,y=BP_B_cum)) + geom_point(aes(color=R2)) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
	scale_x_continuous(label = axisdf$chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    scale_y_continuous(label = axisdf$chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
	xlab('Chromosome (Mb)') + ylab("Chromosome (Mb)") +
	ggtitle("Factor 2 Genes in High LD ")

pdf('images/Factor2_ld.pdf')
print(p1)
dev.off()



### Factor 14

ld=fread('MegaLMM/Factor14_rsquared.ld',data.table=F)
ld=ld[ld$R2>=0.5,]

f14=factordf[factordf$Trait=="Factor14",]

df=data.frame(snp=unique(f14$nearest_SNP))
#if(!unique(f14$SNP) %in% df$snp){
df=rbind(df,data.frame(snp=unique(f14$SNP)))

#}
fwrite(df,'MegaLMM/WD_0727_Factor14_snplist.txt',row.names=F,quote=F,sep=',',col.names=F)


ld$Gene_A=f14[match(ld$SNP_A,f14$nearest_SNP),]$Gene
ld=ld[!is.na(ld$Gene_A),]
ld$Gene_B=f14[match(ld$SNP_B,f14$nearest_SNP),]$Gene
ld=ld[!is.na(ld$Gene_A) & !is.na(ld$Gene_B),]

ld$chr_lenA=axisdf[match(ld$CHR_A,axisdf$chr),]$chr_len
ld$chr_lenB=axisdf[match(ld$CHR_B,axisdf$chr),]$chr_len
ld$totA=axisdf[match(ld$CHR_A,axisdf$chr),]$tot
ld$totB=axisdf[match(ld$CHR_B,axisdf$chr),]$tot
ld$BP_A_cum = ld$BP_A+ld$totA
ld$BP_B_cum = ld$BP_B+ld$totB


gene="Zm00001d035217" 

gene="Zm00001d047592"

p1=ggplot(ld,aes(x=BP_A_cum,y=BP_B_cum)) + geom_point(aes(color=R2)) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
	scale_x_continuous(label = axisdf$chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    scale_y_continuous(label = axisdf$chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
	xlab("Chromosome") + ylab("Chromosome") +
	ggtitle("Factor 14 Genes in High LD ")

pdf('images/Factor14_ld.pdf')
print(p1)
dev.off()

for(c in 1:10){
	subld=ld[ld$CHR_A==c & ld$CHR_B==c,]
	p1=ggplot(subld,aes(x=BP_A/1e6,y=BP_B/1e6)) + geom_point(aes(color=R2)) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
	#scale_x_continuous(label = axisdf$chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    #scale_y_continuous(label = axisdf$chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
	xlab('Chromosome (Mb)') + ylab("Chromosome (Mb)") +
	ggtitle(sprintf("Factor 14 Genes in High LD on Chromosome %.0f",c))

	pdf(sprintf('images/Factor14_c%.0f_ld.pdf',c))
	print(p1)
	dev.off()
}


####

ld=fread('MegaLMM/Factor14_rsquared.ld',data.table=F)
ld=ld[ld$R2>=0.5,]
f14_snps=unique(f14$SNP)
ld=ld[ld$SNP_A %in% f14_snps | ld$SNP_B %in% f14_snps,]


ld=fread('MegaLMM/Factor2_rsquared.ld',data.table=F)
ld=ld[ld$R2>=0.5,]
f2_snps=unique(f2$SNP)
ld=ld[ld$SNP_A %in% f2_snps | ld$SNP_B %in% f2_snps,]