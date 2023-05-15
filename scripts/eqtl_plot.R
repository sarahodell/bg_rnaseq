#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')


cis=fread('eqtl/results/WD_0712_cis_eQTL_fkeep_hits.txt',data.table=F)
cis$class="cis"
cis$time="WD_0712"
names(cis)=c("Trait","CHR","BP","SNP","value","class","time")
trans=fread('eqtl/results/WD_0727_trans_eQTL_scan_hits.txt',data.table=F)
trans$class="trans"
trans$time="WD_0727"
names(trans)=c("Trait","CHR","BP","SNP","value","class","time")

factoreqtl=fread('eqtl/results/factor_transQTL_all.txt',data.table=F)
factoreqtl$class="factor"
names(factoreqtl)=c("Trait","CHR","BP","SNP","value","time","class")
factoreqtl=factoreqtl[,c("Trait","CHR","BP","SNP","value","class","time")]

allhits=rbind(cis,trans)
allhits=rbind(allhits,factoreqtl)

fwrite(allhits,'eqtl/results/all_eQTL_hits.txt',row.names=F,quote=F,sep='\t')

allhits=fread('eqtl/results/all_eQTL_hits.txt',data.table=F)

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
f2genes=prop27[prop27$Factor2>=.25,]$V1
f14genes=prop27[prop27$Factor14>=.25,]$V1


f18_groups=readRDS('MegaLMM/MegaLMM_WD_0718_factor_groups.rds')
f16genes=f18_groups[['Factor16']]$genes

prop18=fread('MegaLMM/MegaLMM_WD_0718_prop_variance.txt',data.table=F)
f16genes=prop18[prop18$Factor16>=.25,]$V1


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

fwrite(factordf,'eqtl/results/factor_gene_table.txt',row.names=F,quote=F,sep='\t')

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


axisdf <- df.tmp %>%

    # Compute chromosome size
    group_by(gene_chr) %>%
    summarise(chr_len=max(chr_lenB)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>% mutate(center=(chr_len/2)+tot)


cumtot=2105119857

p2=ggplot(data=df.tmp, aes(x=midblock_cum, y=midgene_cum)) +
	geom_point(aes(color=prop_var,shape=Trait),,size=1,alpha=0.4) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
    scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    labs(x = "eQTL variants by Chromosome",y="eQTL transcripts by Chromosome") + 
    guides(shape="none")

pdf('eqtl/images/factor_0.2_eQTL_plot.pdf')
print(p2)
dev.off()

df.tmp$Trait=paste0(df.tmp$time,'-',df.tmp$Trait)

for(i in unique(df.tmp$Trait)){
	subtmp=df.tmp[df.tmp$Trait==i,]
	p2=ggplot(data=subtmp, aes(x=midblock_cum, y=midgene_cum)) +
	geom_point(aes(color=prop_var),size=1,alpha=0.4) +
	scale_color_gradient(low = "lightblue", high = "darkblue") +
    scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    labs(x = "eQTL variants by Chromosome",y="eQTL transcripts by Chromosome") + 
    ggtitle(sprintf("%s eQTL",i))

	pdf(sprintf('eqtl/images/%s_eQTL_plot.pdf',i))
	print(p2)
	dev.off()
}




#p2=ggplot(df.tmp,aes(x=midblock_cum,y=midgene_cum)) + geom_point(color=Trait)

### just cis and trans plot first

sub=allhits[allhits$class %in% c('cis','trans'),]
sub$gene_start=genetable[match(sub$Trait,genetable$Gene_ID),]$START
sub$gene_end=genetable[match(sub$Trait,genetable$Gene_ID),]$END

sub$mid_block=apply(sub[,c('block_start','block_end')],1,function(x) round(mean(x)))

sub$mid_gene=apply(sub[,c('gene_start','gene_end')],1,function(x) round(mean(x)))



sub$chr_len=chr_ends[match(sub$CHR,chr_ends$chr),]$chr_len


df.tmp3 <- df.tmp %>%

    # Compute chromosome size
    group_by(gene_chr) %>%
    summarise(chr_len=max(chr_lenB)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(sub, ., by=c("CHR"="gene_chr")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, mid_block) %>%
    mutate(midblock_cum=as.numeric(mid_block+tot),midgene_cum=as.numeric(mid_gene+tot))
    #gather(key, value, -mid_block,-mid_gene,-SNP,-BP,-CHR,-midblock_cum,-midgene_cum,-tot,-Trait,-class,-time)

    # Add highlight and annotation information
    #mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) #%>%
    #mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  # get chromosome center positions for x-axis


#axisdf <- sub %>%

    # Compute chromosome size
#    group_by(CHR) %>%
#    summarise(chr_len=max(chr_len)) %>%

    # Calculate cumulative position of each chromosome
#    mutate(tot=cumsum(chr_len)-chr_len) %>% mutate(center=(chr_len/2)+tot)
#yaxisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(midgene_cum) + min(midgene_cum) ) / 2 )

#df.tmp=df.tmp[!is.na(df.tmp$value),]


mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
greypalette=gray.colors(6)

p3=p2 + new_scale_color() + geom_point(data=df.tmp3,aes(x=midblock_cum,y=midgene_cum,color=class,shape=time))
pdf('eqtl/images/all_0.2_eqtl_plot.pdf')
print(p3)
dev.off()


cumtot=1377773285

p1=ggplot(data=df.tmp, aes(x=midblock_cum, y=midgene_cum)) +
	geom_point(aes(color=class,shape=time)) +
    scale_x_continuous(label = axisdf$CHR,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    scale_y_continuous(label = axisdf$CHR,breaks=axisdf$center,minor_breaks=axisdf$tot,limits = c(0,cumtot)) +
    labs(x = "eQTL variants by Chromosome",y="eQTL transcripts by Chromosome")

png('eqtl/images/eQTL_plot.png',width=1000,height=1000)
print(p1)
dev.off()



    # Show all points

    #scale_color_manual(values = rep(greypalette, 22 )) +
	#geom_rect(data=axisdf,aes(xmax = tot, 
    #            xmin = tot-chr_len, 
    #            ymin = tot-chr_len, 
    #            ymax = tot, 
    #            fill = factor(CHR)), alpha = 0.4) + 
  	#scale_fill_manual(values = greypalette) +
    # custom X axis:
    # add genome-wide sig and sugg lines
    #geom_hline(yintercept = threshold,linetype="dashed") +
    #geom_hline(yintercept = -log10(sugg), linetype="dashed") +

    # Add highlighted points
    #geom_point(data=subset(df.tmp, sig==T), color="coral2", size=5) +
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


p1=ggplot(aes(x=mid_block,y=mid_gene),data=sub) + geom_point(aes(color=class,shape=time)) + xlab("eQTL variant") + ylab("eQTL transcript")

png('eqtl/images/eQTL_plot.png')
print(p1)
dev.off()
