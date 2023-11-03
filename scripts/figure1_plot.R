#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')
library('ggnewscale')
library('RColorBrewer')


#eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
#eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
#eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
#eqtl=as.data.frame(eqtl2)

#all_founder_blocks=c()
#for(chr in 1:10){#
#  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
#  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
#}
#eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
#eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end
#
#
#
#trans=fread('eqtl/results/all_trans_fdr_SIs_FIXED.txt',data.table=F)
#names(trans)=c('Trait','time','CHR','X_ID','BP','block_start','block_end_bp','left_bound_snp','right_bound_snp','block_end')
#trans=trans[,c('Trait','X_ID','CHR','BP','time','block_start','block_end')]
#trans$class="trans"
#
#eqtl=eqtl[,c('Trait','X_ID','CHR','BP','time','block_start','block_end')]
#eqtl$class='cis'
#
#
#allhits=rbind(eqtl,trans)
##allhits=rbind(allhits,factoreqtl)
#names(allhits)[2]="SNP"
#fwrite(allhits,'eqtl/results/all_local_distal_eQTL_fdr_peaks.txt',row.names=F,quote=F,sep='\t')



allhits=fread('eqtl/results/all_local_distal_eQTL_fdr_peaks.txt',data.table=F)

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)


allhits$block_start=all_founder_blocks[match(allhits$SNP,all_founder_blocks$focal_snp),]$start
allhits$block_end=all_founder_blocks[match(allhits$SNP,all_founder_blocks$focal_snp),]$end

#fwrite(axisdf,'eqtl/data/chromosome_axis.txt',row.names=F,quote=F,sep='\t')
axisdf=fread('eqtl/data/chromosome_axis.txt',data.table=F)
axisdf$chr_start=as.numeric(axisdf$chr_start)
axisdf$chr_end=as.numeric(axisdf$chr_end)
axisdf$center_shift=as.numeric(axisdf$center_shift)

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
df.tmp4$time_f=factor(df.tmp4$time,levels=c("WD_0712","WD_0718","WD_0720","WD_0727"))

df.tmp4=df.tmp4[order(df.tmp4$midblock_cum),]
row.names(df.tmp4)=seq(1,nrow(df.tmp4))
df.tmp4$X_SNP=factor(df.tmp4$SNP,levels=c(unique(df.tmp4$SNP)))

snp_order=unique(df.tmp4$SNP)
#gene_snp=c()
#for(chr in 1:10){
#	testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
#	genes=unique(df.tmp4[df.tmp4$gene_chr==chr,]$Trait)
#	for(gene in genes){
#		snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
#		line=data.frame(chr=chr,gene=gene,snp=snp,stringsAsFactors=F)
#		gene_snp=rbind(gene_snp,line)
#	}
#}
#gene_snp=as.data.frame(gene_snp,stringsAsFactors=F)
#fwrite(gene_snp,'eqtl/data/gene_focal_snps.txt',row.names=F,quote=F,sep='\t')
gene_snp=fread('eqtl/data/gene_focal_snps.txt',data.table=F)
df.tmp4$gene_SNP=gene_snp[match(df.tmp4$Trait,gene_snp$gene),]$snp

df.tmp4=df.tmp4[order(df.tmp4$midgene_cum),]
row.names(df.tmp4)=seq(1,nrow(df.tmp4))

df.tmp4$Y_SNP=factor(df.tmp4$gene_SNP,levels=c(snp_order))

snps=unique(c(unique(df.tmp4$SNP),unique(df.tmp4$gene_SNP)))

subtmp=df.tmp4[df.tmp4$time=="WD_0712",]
mat=expand.grid(snps,snps)
ngenes=df.tmp4 %>% group_by(X_SNP,Y_SNP) %>% reframe(n=length(unique(Trait)))
ngenes=as.data.frame(ngenes)
#mat$snp_pair=paste0(mat$Var2,'_',mat$Var1)

ngenes$snp_pair=paste0(ngenes$X_SNP,'_',ngenes$Y_SNP)
#notin=mat$snp_pair[!(mat$snp_pair %in% ngenes$snp_pair)]
#mat$n=ngenes[match(c(mat$Var2,mat$Var1),c(ngenes$X_SNP,ngenes$Y_SNP)),]$n
#mat$n[is.na(mat$n),]=0

#nmat=mat[!(mat$snp_pair %in% ngenes$snp_pair),]
#nmat$n=0
#nmat=nmat[,c('Var2','Var1','n','snp_pair')]
#names(nmat)=c('X_SNP','Y_SNP','n','snp_pair')

#ngenes2=rbind(ngenes,nmat)


theme_set(theme_classic())
#theme_update(text=element_text(family="Times"))
theme_update(plot.caption = element_text(hjust = 0))
theme_update(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
theme_update(plot.title = element_text(size=20),axis.title=element_text(size=20))
theme_update(panel.background=element_blank())
theme_update(plot.caption=element_text(size=20))

#reds=c("#FEE5D9" ,"#FCAE91", "#FB6A4A", "#CB181D")
reds=brewer.pal(n=9,"Reds")
reds=reds[5:9]
ngenes$X_bin=as.numeric(ngenes$X_SNP)
ngenes$Y_bin=as.numeric(ngenes$Y_SNP)

ngenes$x_chr=all_founder_blocks[match(ngenes$X_SNP,all_founder_blocks$focal_snp),]$chr
ngenes$y_chr=all_founder_blocks[match(ngenes$Y_SNP,all_founder_blocks$focal_snp),]$chr

# Which SNP is the last one on the chromosome?
#for(i in 1:10){
#	print(tail(all_founder_blocks[all_founder_blocks$chr==i,],n=1))
#}

breaks=c('AX-91509113','AX-91553484','AX-91217538','AX-91221920','AX-91681652','AX-90594801','AX-91068340','AX-90636227','AX-91157044','AX-91834694')
xbreak=ngenes[match(breaks,ngenes$X_SNP),]$X_bin
#ybreak=ngenes[match(breaks,ngenes$Y_SNP),]$Y_bin

p1=ggplot(ngenes,aes(x=X_SNP,y=Y_SNP,fill=n)) + geom_raster() +
scale_fill_gradientn(colours=c(reds,'black')) +
#theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
#axis.ticks.x=element_blank(),axis.ticks.y=element_blank()) +
scale_x_discrete(label = axisdf$gene_chr,breaks=breaks) +
scale_y_discrete(label = axisdf$gene_chr,breaks=breaks) +
theme(panel.grid.major=element_line(colour="darkgrey"),panel.grid.minor=element_blank()) +
xlab('eQTL Variant') + ylab("eQTL Transcript")

#png('eqtl/images/WD_0712_heatmap.png')
#print(p1)
#dev.off()

pdf('eqtl/images/WD_0712_heatmap.pdf',width=12,height=12)
print(p1)
dev.off()


cistrans=ggplot(df.tmp4,aes(x=midblock_cum,y=midgene_cum)) +
	geom_abline(slope=1, intercept=0,alpha=0.5) +
	scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center_shift,minor_breaks=minor_breaks,limits=c(0, cumtot)) +
    scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center_shift,minor_breaks=minor_breaks,limits=c(0, cumtot)) +
    geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
    geom_point(aes(color=time_f),size=0.1) +
    scale_color_hue(name="Timepoint",labels=c("T12","T18","T20","T27")) + 
    labs(x = "Position of eQTL variants",y="Position of eQTL Transcripts") +
    theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_blank())
 
pdf('paper_figures/cis_trans.pdf')
print(cistrans)
dev.off()

output=df.tmp4[,c('Trait','SNP','CHR','block_start','block_end','gene_chr','gene_start','gene_end','class','time')]
fwrite(output,'eqtl/results/all_cis_trans_fdr_hits_FIXED.txt',row.names=F,quote=F,sep='\t')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")

for(time in times){
	subtmp=df.tmp4[df.tmp4$time==time,]
	cistrans=ggplot(subtmp,aes(x=midblock_cum,y=midgene_cum)) +
	geom_abline(slope=1, intercept=0,alpha=0.5) +
	scale_x_continuous(label = axisdf$gene_chr,breaks=axisdf$center_shift,minor_breaks=minor_breaks,limits=c(0, cumtot)) +
    scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center_shift,minor_breaks=minor_breaks,limits=c(0, cumtot)) +
    geom_hline(yintercept=cumtot,colour="darkgrey") +
    geom_vline(xintercept=cumtot,colour="darkgrey") +
 	geom_point(aes(color=class),size=0.1) + 
 	scale_color_hue(name="Type",labels=c("local",'distal'))
 	#scale_color_manual(values = c("cis" = "coral","trans" = "lightblue")) + 
    labs(x = "Position of eQTL variants",y="Position of eQTL Transcripts") +
    theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_blank())
 
	pdf(sprintf('paper_figures/%s_cis_trans_fdr.pdf',time))
	print(cistrans)
	dev.off()

}

#######




time="WD_0720"
subtmp=df.tmp4[df.tmp4$time==time,]
subtmp=subtmp[subtmp$CHR==8,]

env1=subtmp[,c('Trait','SNP','CHR','block_start','block_end','mid_block','tot.x','midgene_cum','mid_gene')]
names(env1)=c('Trait.1','SNP.1','CHR','block_start','block_end','mid_block.1','tot.x.1','midgene_cum.11','mid_gene.1')
env1=as.data.table(env1)
env2=subtmp[,c('Trait','SNP','CHR','block_start','block_end','mid_block','tot.x','midgene_cum','mid_gene')]
names(env2)=c('Trait.2','SNP.2','CHR','block_start','block_end','mid_block.2','tot.x.2','midgene_cum.2','mid_gene.2')
env2=as.data.table(env2)
setkey(env2,CHR,block_start,block_end)
comparison4=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
comparison4=comparison4[comparison4$Trait.1!=comparison4$Trait.2,]

overlap=union(unique(comparison4$Trait.1),unique(comparison4$Trait.2))
subtmp$overlap=subtmp$Trait %in% overlap

minor_breaks=sort(c(0,axisdf$chr_start[-1],axisdf$chr_end))

cistrans=ggplot(subtmp,aes(x=mid_block/1e6,y=midgene_cum)) +
geom_abline(slope=1, intercept=0,alpha=0.5) +
scale_y_continuous(label = axisdf$gene_chr,breaks=axisdf$center_shift,minor_breaks=minor_breaks,limits=c(0, cumtot)) +
annotate('rect', xmin=comparison4$block_start/1e6, xmax=comparison4$block_end/1e6, ymin=0, ymax=cumtot,color='black',fill='dodgerblue',alpha=0.25) +
geom_point(aes(color=overlap),size=1) + 
scale_color_manual(values = c("FALSE"="red", "TRUE"="blue")) + 
labs(x = "Position of eQTL variants",y="Position of eQTL Transcripts") +
theme(panel.grid.minor=element_line(colour="darkgrey"),panel.grid.major=element_line(colour="black"))
 
pdf(sprintf('eqtl/images/%s_shadow_fdr.pdf',time),width=10,height=8)
print(cistrans)
dev.off()

