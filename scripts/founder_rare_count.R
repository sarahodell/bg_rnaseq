#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('reshape2')
library('dplyr')
library('plyr')
# Tester rare alleles
tester_rares=fread('datasets/all_tester_rare_allele_probs.txt',data.table=F)

bimbam=c()
for(chr in 1:10){
	b=fread(sprintf('datasets/chr%s_biogemma_rare_allele_probs.txt',chr),data.table=F)
	bimbam=rbind(bimbam,b)
	rm(b)
}
bimbam=as.data.frame(bimbam,stringsAsFactors=F)

fwrite(bimbam,'datasets/all_biogemma_rare_allele_probs.txt',row.names=F,quote=F,sep='\t')

# How many rares does each founder have?

allfounders=c()
for(chr in 1:10){
	founder=fread(sprintf('datasets/hapmap/chr%s_founder_rare_alleles.txt',chr),data.table=F)
	names(founder)=c("SNP","CHR","POS","REF","ALT","A632_usa","A654_inra","B73_inra","B96","C103_inra","CO255_inra","D105_inra","DK63","EP1_inra","F492", "FV252_inra", "FV2_inra","MBS847","ND245","OH43_inra","VA85","W117_inra")
	#overlap=merge(rares,hapsnps,by='SNP')
	founders=c("MBS847","B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
	founder=founder[,c("SNP","CHR","POS","REF","ALT",founders)]
	founder[,founders]=apply(founder[,founders],MARGIN=2,function(x) ifelse(x=="0/0",0,ifelse(x=="1/1",1,ifelse(x=="2/2",2,NA))))
	# remove multiallelic sites?
	founder$ALT2=NA
	mult=which(grepl(',',founder$ALT))
	alts=founder$ALT[mult]
	alt1=sapply(seq(1,length(alts)),function(x) strsplit(alts[x],',')[[1]][[1]])
	alt2=sapply(seq(1,length(alts)),function(x) strsplit(alts[x],',')[[1]][[2]])
	founder[mult,]$ALT=alt1
	founder[mult,]$ALT2=alt2
	founder=founder[,c("SNP","CHR","POS","REF","ALT","ALT2",founders)]
	allfounders=rbind(allfounders,founder)
}
allfounders=as.data.frame(allfounders,stringsAsFactors=F)


totals=colSums(allfounders[,founders],na.rm=T)

#    MBS847   B73_inra   A632_usa CO255_inra FV252_inra  OH43_inra  A654_inra 
#    150832        535     108726     396490     247302     138968     277072 
#  FV2_inra  C103_inra   EP1_inra  D105_inra  W117_inra        B96       DK63 
#    436192     143049     457816     441348     285371     406276     186512 
#      F492      ND245       VA85 
#    338849     315699     230018 

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

# How many rare alleles does each founder have within a window relative to the window's size
get_rare_counts=function(rep){
	row1=all_founder_blocks[rep,]
	#gene=as.character(row1$Gene_ID)
	#id=as.character(row1$ID)
	#r=ranks[id,gene]
	chrom=row1$chr
	rstart=row1$start
	rend=row1$end
	tmp=allfounders[allfounders$CHR==chrom & allfounders$POS>=rstart & allfounders$POS<=rend,]
	tcounts=colSums(tmp[,founders],na.rm=T)
	
	line=as.data.frame(c(row1,tcounts),stringsAsFactors=F)
	return(line)
}

new_founder_blocks=sapply(seq(1,nrow(all_founder_blocks)),function(x) get_rare_counts(x))

tmp=t(new_founder_blocks)
tmp=as.data.frame(tmp,stringsAsFactors=F)

for(f in founders){
	tmp[,f]=as.numeric(tmp[,f])
}
tmp$chr=as.numeric(tmp$chr)
tmp$start=as.numeric(tmp$start)
tmp$end=as.numeric(tmp$end)
tmp$focal_snp=as.character(tmp$focal_snp)

fwrite(tmp,'datasets/founder_block_rare_counts.txt',row.names=F,quote=F,sep='\t')
# Start just with one chrom
plot_list=list()
count=1
for(chr in 1:10){
	tmp10=tmp[tmp$chr==chr,]
	tmp_df=reshape2::melt(tmp10,id=c('chr','start','end','focal_snp'))
	tmp_df$size=tmp_df$end-tmp_df$start
	tmp_df$adjusted_count=tmp_df$value/tmp_df$size
	mean_df=tmp_df %>% group_by(focal_snp) %>% summarize(mean_count=mean(adjusted_count))
	mean_df=as.data.frame(mean_df,stringsAsFactors=F)
	mean_df=merge(mean_df,all_founder_blocks,by='focal_snp')
	#ymax=max(tmp_df$adjusted_count)+0.0001

	p = ggplot(tmp_df,aes(x = start/1e6, y=adjusted_count)) +
	geom_segment(aes(xend = end/1e6, yend = adjusted_count,group=variable,color=variable),linewidth=2,alpha=0.5) 
	p = p + geom_segment(data=mean_df,aes(x=start/1e6, xend = end/1e6,y=mean_count, yend = mean_count),linewidth=1,color='black') +
	xlab(sprintf("Chr %s Position (Mb)",chr)) + ylim(0,0.005) +
	ylab("Rare Alleles Counts (per bp)")
	plot_list[[count]]=p
	count=count+1
}


pdf('images/chr10_rares_by_founder.pdf')
print(p)
dev.off()
 
pdf('images/rares_by_founder.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()


plot_list=list()
count=1
for(chr in 1:10){
	tmp10=tmp[tmp$chr==chr,]
	tmp10=tmp10[order(tmp10$start),]
	rownames(tmp10)=seq(1,nrow(tmp10))
	tmp10$snp_f=factor(tmp10$focal_snp,levels=tmp10$focal_snp)
	tmp_df=reshape2::melt(tmp10,id=c('chr','start','end','focal_snp','snp_f'))
	tmp_df$size=tmp_df$end-tmp_df$start
	tmp_df$adjusted_count=tmp_df$value/tmp_df$size
	tmp_df=tmp_df[tmp_df$size>1000,]

	mean_df=tmp_df %>% group_by(focal_snp) %>% summarize(mean_count=mean(adjusted_count))
	mean_df=as.data.frame(mean_df,stringsAsFactors=F)
	mean_df=merge(mean_df,all_founder_blocks,by='focal_snp')
	mean_df=mean_df[order(mean_df$start),]
	row.names(mean_df)=seq(1,nrow(mean_df))
	mean_df$snp_f=factor(mean_df$focal_snp,levels=c(mean_df$focal_snp))
	ymax=max(tmp_df$adjusted_count)+0.0001

	p = ggplot(tmp_df,aes(x = snp_f, y=adjusted_count)) +
	geom_line(aes(group=variable,color=variable),alpha=0.5) 
	p = p + geom_line(data=mean_df,aes(x=snp_f, y=mean_count),color='black') +
	xlab(sprintf("Chr %s Position (Mb)",chr)) + ylim(0,ymax) +
	ylab("Rare Alleles Counts (per bp)")
	plot_list[[count]]=p
	count=count+1
}

 
pdf('images/rares_by_founder.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()


plot_list=list()
count=1

allrares=c()
#chr=10
for(chr in 1:10){
	tmp10=allfounders[allfounders$CHR==chr,]
	minb=round_any(min(tmp10$POS), 1e5,f=floor) 
	maxb=round_any(max(tmp10$POS), 1e5,f=ceiling) 
	bbreaks=seq(minb,maxb,1e5)
	tmp10$bin=cut(tmp10$POS, breaks=bbreaks,labels=F)
	
	allbins=unique(tmp10$bin)
	rare_bins=c()
	#binstarts=c()
	#binends=c()
	for(bin in allbins){
		binstart=bbreaks[bin]
		binend=binstart+1e5
		#binstarts=c(binstarts,binstart)
		#binends=c(binends,binend)
		#row1=all_founder_blocks[rep,]
		#gene=as.character(row1$Gene_ID)
		#id=as.character(row1$ID)
		#r=ranks[id,gene]
		#chrom=row1$chr
		#rstart=row1$start
		#rend=row1$end
		tmp=tmp10[tmp10$bin==bin,]
		tcounts=colSums(tmp[,founders],na.rm=T)
		
		line=as.data.frame(tcounts,stringsAsFactors=F)
		line$bin=bin
		line$founder=names(tcounts)
		line$chr=chr
		line$binstart=binstart
		line$binend=binend
		rare_bins=rbind(rare_bins,line)
	}
	
	rare_bins=as.data.frame(rare_bins,stringsAsFactors=F)
	rare_bins=rare_bins[rare_bins$founder!="MBS847",]
	#rare_bins$perc=rare_bins$tcounts/rare_bins$total
	allrares=rbind(allrares,rare_bins)
	p = ggplot(rare_bins,aes(x = bin, y=tcounts)) +
	geom_line(aes(group=founder,color=founder),alpha=0.5) +
	#p = p + geom_line(data=mean_df,aes(x=snp_f, y=mean_count),color='black') +
	xlab(sprintf("Chr %s Position (Mb)",chr)) + #ylim(0,1) +
	ylab("Homozygous Rare Alleles Counts (per 100kb)")
	plot_list[[count]]=p
	count=count+1
}


fwrite(allrares,'datasets/founder_block_rare_bins_melt.txt',row.names=F,quote=F,sep='\t')


#allrares=fread('datasets/founder_block_rare_bins_melt.txt',data.table=F)
pdf('images/rares_by_founder_bins.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()

### Create plot

#shift=2.5e7
#
#tot_shift=c(0)
#center_shift=c(axisdf[axisdf$gene_chr==1,]$center)
#cums=c()
#for(i in 1:10){
#	seqs=length(1:i)-1
#	addon=sum(axisdf$chr_len[1:i])
#	addon=addon + (shift*seqs)
#	cums=c(cums,addon)
#}
#axisdf$len_shift=cums
#
#cstart=c(0)
#cend=c(306972432)
#cshift=c(153486216)
#for(i in 1:nrow(axisdf)){
#	if(i !=1){
#		row=axisdf[i,]
#		prev=axisdf[(i-1),]
#		prev_end=prev$len_shift
#		current_start=prev_end+shift
#		cstart=c(cstart,current_start)
#		current_end=row$len_shift
#		cend=c(cend,current_end)
#		center_shift=current_start+((current_end-current_start)/2)
#		cshift=c(cshift,center_shift)
#	}
#}
#axisdf$center_shift=cshift
#axisdf$chr_start=cstart
#axisdf$chr_end=cend
#
#fwrite(axisdf,'eqtl/data/chromosome_axis.txt',row.names=F,quote=F,sep='\t')

### just cis and trans plot first

tester_rares=allfounders[allfounders$MBS847==1,]
tester_rares=tester_rares[!is.na(tester_rares$MBS847),]
totals=colSums(tester_rares[,founders],na.rm=T)
# MBS847   B73_inra   A632_usa CO255_inra FV252_inra  OH43_inra  A654_inra 
#    150280        200       4924       8174       9275       5701      11007 
#  FV2_inra  C103_inra   EP1_inra  D105_inra  W117_inra        B96       DK63 
#      8233       5157       8264       9967       9093       6360      76942 
#      F492      ND245       VA85 
#      7130      10169       7460 

plot_list=list()
count=1

allhrares=c()
#chr=10
for(chr in 1:10){
	trares=tester_rares[tester_rares$CHR==chr,]
	minb=round_any(min(trares$POS), 1e5,f=floor) 
	maxb=round_any(max(trares$POS), 1e5,f=ceiling) 
	bbreaks=seq(minb,maxb,1e5)
	trares$bin=cut(trares$POS, breaks=bbreaks,labels=F)
	
	allbins=unique(trares$bin)
	hrare_bins=c()
	#binstarts=c()
	#binends=c()
	for(bin in allbins){
		binstart=bbreaks[bin]
		binend=binstart+1e5
		#binstarts=c(binstarts,binstart)
		#binends=c(binends,binend)
		#row1=all_founder_blocks[rep,]
		#gene=as.character(row1$Gene_ID)
		#id=as.character(row1$ID)
		#r=ranks[id,gene]
		#chrom=row1$chr
		#rstart=row1$start
		#rend=row1$end
		tmp=trares[trares$bin==bin,]
		tcounts=colSums(tmp[,founders],na.rm=T)
		
		line=as.data.frame(tcounts,stringsAsFactors=F)
		line$bin=bin
		line$founder=names(tcounts)
		line$total=tcounts['MBS847']
		line$chr=chr
		line$binstart=binstart
		line$binend=binend
		hrare_bins=rbind(hrare_bins,line)
	}
	
	hrare_bins=as.data.frame(hrare_bins,stringsAsFactors=F)
	hrare_bins=hrare_bins[hrare_bins$founder!="MBS847",]
	hrare_bins$perc=hrare_bins$tcounts/hrare_bins$total

	allhrares=rbind(allhrares,hrare_bins)
	
	p = ggplot(hrare_bins,aes(x = bin, y=tcounts)) +
	geom_line(aes(group=founder,color=founder),alpha=0.5) +
	#p = p + geom_line(data=mean_df,aes(x=snp_f, y=mean_count),color='black') +
	xlab(sprintf("Chr %s Position (Mb)",chr)) + #ylim(0,1) +
	ylab("% Homozygous Rare Alleles Counts (per 100kb)")
	plot_list[[count]]=p
	count=count+1
}

fwrite(allhrares,'datasets/founder_block_hrare_bins_melt.txt',row.names=F,quote=F,sep='\t')


pdf('images/hrares_by_founder.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()

pdf('images/chr10_hrares_by_founder.pdf')
print(p)
dev.off()