#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')
library('ggnewscale')


cis=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits.txt',data.table=F)
cis$class="cis"
cis=cis[,c('Trait','X_ID','CHR','BP','value','class','time')]


trans=fread('eqtl/results/all_trans_fdr_peaks.txt',data.table=F)
#trans=fread('eqtl/results/all_trans_eQTL_weights_fdr_hits.txt',data.table=F)
trans$class="trans"
#names(trans)=c('Trait','CHR','BP','X_ID','value','class','time')
trans=trans[,c('Trait','CHR','X_ID','BP','value','class','time')]

#factoreqtl=fread('eqtl/results/all_factor_trans_eqtl_fdr_hits',data.table=F)
factoreqtl=fread('eqtl/results/all_factor_fdr_peaks.txt',data.table=F)
factoreqtl$class="factor"
factoreqtl=factoreqtl[,c('Trait','X_ID','CHR','BP','value','class','time')]


allhits=rbind(cis,trans)
allhits=rbind(allhits,factoreqtl)
names(allhits)[2]="SNP"
fwrite(allhits,'eqtl/results/all_eQTL_fdr_peaks.txt',row.names=F,quote=F,sep='\t')

allhits=fread('eqtl/results/all_eQTL_fdr_peaks.txt',data.table=F)

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)


allhits$block_start=all_founder_blocks[match(allhits$SNP,all_founder_blocks$focal_snp),]$start
allhits$block_end=all_founder_blocks[match(allhits$SNP,all_founder_blocks$focal_snp),]$end


#f27_groups=readRDS('MegaLMM/MegaLMM_WD_0727_factor_groups.rds')
#f2genes=f27_groups[['Factor2']]$genes
#f14genes=f27_groups[['Factor14']]$genes

prop27=fread('MegaLMM/MegaLMM_WD_0727_prop_variance.txt',data.table=F)
f2genes=prop27[prop27$Factor2>=0.10,]$V1
f14genes=prop27[prop27$Factor14>=0.10,]$V1
f22genes=prop27[prop27$Factor22>=0.10,]$V1


prop20=fread('MegaLMM/MegaLMM_WD_0720_prop_variance.txt',data.table=F)
f9genes=prop20[prop20$Factor9>=0.10,]$V1

#f18_groups=readRDS('MegaLMM/MegaLMM_WD_0718_factor_groups.rds')
#f16genes=f18_groups[['Factor16']]$genes

prop18=fread('MegaLMM/MegaLMM_WD_0718_prop_variance.txt',data.table=F)
f16genes=prop18[prop18$Factor16>=0.10,]$V1

prop12=fread('MegaLMM/MegaLMM_WD_0712_prop_variance.txt',data.table=F)
f1genes=prop12[prop12$Factor1>=0.10,]$V1



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

sub2=factoreqtl[factoreqtl$Trait=="Factor22",]
#names(sub2)[2]="SNP"
gtable=genetable[genetable$Gene_ID %in% f22genes,]
for(i in 1:nrow(gtable)){
	subdf=sub2
	row=gtable[i,]
	subdf$Gene=row$Gene_ID
	subdf$gene_chr=row$CHROM
	subdf$gene_start=row$START
	subdf$gene_end=row$END
	subdf$prop_var=prop27[match(subdf$Gene,prop27$V1),]$Factor22
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

sub2=factoreqtl[factoreqtl$Trait=="Factor9",]
gtable=genetable[genetable$Gene_ID %in% f9genes,]

for(i in 1:nrow(gtable)){
	subdf=sub2
	row=gtable[i,]
	subdf$Gene=row$Gene_ID
	subdf$gene_chr=row$CHROM
	subdf$gene_start=row$START
	subdf$gene_end=row$END
	subdf$prop_var=prop20[match(subdf$Gene,prop20$V1),]$Factor9
	factordf=rbind(factordf,subdf)
}

sub2=factoreqtl[factoreqtl$Trait=="Factor1",]
gtable=genetable[genetable$Gene_ID %in% f1genes,]

for(i in 1:nrow(gtable)){
	subdf=sub2
	row=gtable[i,]
	subdf$Gene=row$Gene_ID
	subdf$gene_chr=row$CHROM
	subdf$gene_start=row$START
	subdf$gene_end=row$END
	subdf$prop_var=prop12[match(subdf$Gene,prop12$V1),]$Factor1
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
names(factordf)[2]="SNP"

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


#fwrite(axisdf,'eqtl/data/chromosome_axis.txt',row.names=F,quote=F,sep='\t')
axisdf=fread('eqtl/data/chromosome_axis.txt',data.table=F)

shift=2.5e7
#axisdf$tot_shift=axisdf$tot+shift
#axisdf$center_shift=axisdf$center+shift
#axisdf[axisdf$gene_chr==1,]$tot_shift=0
#axisdf[axisdf$gene_chr==1,]$center_shift=axisdf[axisdf$gene_chr==1,]$center

tot_shift=c(0)
center_shift=c(axisdf[axisdf$gene_chr==1,]$center)
cums=c()
for(i in 1:10){
	seqs=length(1:i)-1
	addon=sum(axisdf$chr_len[1:i])
	addon=addon + (shift*seqs)
	cums=c(cums,addon)
}
axisdf$len_shift=cums

cstart=c(0)
cend=c(306972432)
cshift=c(153486216)
for(i in 1:nrow(axisdf)){
	if(i !=1){
		row=axisdf[i,]
		prev=axisdf[(i-1),]
		prev_end=prev$len_shift
		current_start=prev_end+shift
		cstart=c(cstart,current_start)
		current_end=row$len_shift
		cend=c(cend,current_end)
		center_shift=current_start+((current_end-current_start)/2)
		cshift=c(cshift,center_shift)
	}
}
axisdf$center_shift=cshift
axisdf$chr_start=cstart
axisdf$chr_end=cend

### just cis and trans plot first

sub=allhits[allhits$class %in% c('cis','trans'),]
sub$gene_chr=genetable[match(sub$Trait,genetable$Gene_ID),]$CHROM
sub$gene_start=genetable[match(sub$Trait,genetable$Gene_ID),]$START
sub$gene_end=genetable[match(sub$Trait,genetable$Gene_ID),]$END

sub$mid_block=apply(sub[,c('block_start','block_end')],1,function(x) round(mean(x)))

sub$mid_gene=apply(sub[,c('gene_start','gene_end')],1,function(x) round(mean(x)))



sub$chr_lenA=axisdf[match(sub$CHR,axisdf$gene_chr),]$chr_len
sub$chr_lenB=axisdf[match(sub$gene_chr,axisdf$gene_chr),]$chr_len

df.tmp3 <- axisdf %>%
    left_join(sub, ., by=c("CHR"="gene_chr")) %>%
    arrange(CHR, mid_block) %>%
    mutate(midblock_cum=as.numeric(mid_block+chr_start))
    
    
df.tmp4 <- axisdf %>%
    left_join(df.tmp3, ., by=c("gene_chr"="gene_chr")) %>%
    arrange(gene_chr, mid_gene) %>%
    mutate(midgene_cum=as.numeric(mid_gene+chr_start.y))


cumtot=2105119857

cumtot=axisdf[axisdf$gene_chr==10,]$chr_end

minor_breaks=sort(c(0,axisdf$chr_start[-1],axisdf$chr_end))

cistrans=ggplot(df.tmp4,aes(x=midblock_cum,y=midgene_cum)) +
	geom_abline(slope=1, intercept=0,alpha=0.5) +
	scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center_shift,minor_breaks=minor_breaks,limits=c(0, cumtot)) +
    scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center_shift,minor_breaks=minor_breaks,limits=c(0, cumtot)) +
    geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
    geom_point(aes(color=time),size=0.5) + 
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




time="WD_0720"
subtmp=df.tmp4[df.tmp4$time==time,]


env1=subtmp[,c('Trait','SNP','CHR','block_start','block_end','mid_block','tot.x','midgene_cum')]
names(env1)=c('Trait.1','SNP.1','CHR','block_start','block_end','mid_block.1','tot.x.1','midgene_cum.1')
env1=as.data.table(env1)
env2=subtmp[,c('Trait','SNP','CHR','block_start','block_end','mid_block','tot.x','midgene_cum')]
names(env2)=c('Trait.2','SNP.2','CHR','block_start','block_end','mid_block.2','tot.x.2','midgene_cum.1')
env2=as.data.table(env2)
setkey(env2,CHR,block_start,block_end)
comparison4=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
comparison4=comparison4[comparison4$Trait.1!=comparison4$Trait.2,]

overlap=union(unique(comparison4$Trait.1),unique(comparison4$Trait.2))
subtmp$overlap=subtmp$Trait %in% overlap
cistrans=ggplot(subtmp,aes(x=midblock_cum,y=midgene_cum)) +
geom_abline(slope=1, intercept=0,alpha=0.5) +
scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center,minor_breaks=axisdf$tot,limits=c(0, cumtot)) +
geom_hline(yintercept=cumtot,colour="darkgrey") +
geom_vline(xintercept=cumtot,colour="darkgrey") +
annotate('rect', xmin=(comparison4$block_start+comparison4$tot.x.2), xmax=(comparison4$block_end+comparison4$tot.x.2), ymin=0, ymax=comparison4$midgene_cum.1, alpha=.5, fill='blue',alpha=0.3) +
geom_point(aes(color=overlap),size=1) + 
scale_color_manual(values = c("FALSE"="red", "TRUE"="blue")) + 
labs(x = "Position of eQTL variants",y="Position of eQTL Transcripts") +
theme_classic() +
theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_line(colour="black"))
 
pdf(sprintf('eqtl/images/%s_shadow_fdr.pdf',time))
print(cistrans)
dev.off()

