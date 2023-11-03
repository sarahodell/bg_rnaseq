#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')

# Disregulation smile plot

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10
### high GERP
totalrares=c()
for(time1 in times){
	for(chr in chroms){
		allrares=fread(sprintf('eqtl/results/%s_%s_5kb_rare_counts_high_GERP.txt',time1,chr),data.table=F)
		#allrares=allrares[allrares$max_f!="B73_inra",]
		allrares$time=time1
		allrares$chr=chr
		totalrares=rbind(totalrares,allrares)
	}
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)

totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
totalrares$gene_time_id=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$ID)

totalrares$total=rowSums(totalrares[,c('gerp4','gerp3','gerp2')])
gerp_rares=totalrares

# Get max_f and beta info from allrares files
totalrares=c()
for(time1 in times){
	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
	allrares=allrares[allrares$max_f!="B73_inra",]
	totalrares=rbind(totalrares,allrares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
totalrares$gene_time_id=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$ID)

# first I need to break up by gene_time_founder
totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

# Merge the two together, make nice, and 
gerp_rares=merge(gerp_rares,totalrares,by='gene_time_id')
#gerp_rares$gene_time_founder=paste0(gerp_rares$Gene_ID.x,'-',gerp_rares$time.x,'-',gerp_rares$max_f)

gerp_rares=gerp_rares[,c('gene_time_id','Gene_ID.x','ID.x','rank.x','gerp4','gerp3','gerp2','time.x','chr.x','gene_time.x','total','rare_count','max_f','beta','beta_rank','add_rank','gene_time_founder')]

names(gerp_rares)=c('gene_time_id','Gene_ID','ID','rank','gerp4','gerp3','gerp2','time','chr','gene_time','total','rare_count','max_f','beta','beta_rank','add_rank','gene_time_founder')

fwrite(gerp_rares,'eqtl/results/all_5kb_rare_counts_high_GERP.txt',row.names=F,quote=F,sep='\t')


gerp_rares=gerp_rares[!is.na(gerp_rares$max_f),]
gerp_rares=gerp_rares[gerp_rares$max_f!="",]
totf=gerp_rares %>% group_by(gene_time_founder) %>% summarize(Gene_ID=unique(Gene_ID.x),time=unique(time.x),chr=unique(chr.x),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(total),gene_time=unique(gene_time.x),max_f=unique(max_f))

beta_z=totf%>% group_by(gene_time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))

#beta_z=totalrares%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
# Extreme bins are >=3/<=-3 SD
bbreaks=c(min(beta_z$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_z$beta_z)+0.1)
beta_z$bin=cut(beta_z$beta_z, breaks=bbreaks, label=F)
beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)


p3=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/high_gerp_beta_z_bin_disregulation_smooth_smile.png')
print(p3)
dev.off()

p2=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(method = "lm",mapping=aes(weight=ngenes),formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))


png('eqtl/images/high_gerp_beta_z_bin_disregulation_smile.png')
print(p2)
dev.off()

quadratic_highgerp=lm(avg_rares~poly(bin,2),avg_z)
anova(quadratic_highgerp)
summary(quadratic_highgerp)
######### All genes

totalrares=c()
for(time1 in times){
	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
	allrares=allrares[allrares$max_f!="B73_inra",]
	totalrares=rbind(totalrares,allrares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
totalrares$gene_time_id=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$ID)

# first I need to break up by gene_time_founder
totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

#totalrares=totalrares[,-1]
#totalrares=rbind(totalrares,totallow)
#fwrite(totalrares,'eqtl/results/all_rare_counts_max_f_all_exp.txt',row.names=F,quote=F,sep='\t')

totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]
totf=totalrares %>% group_by(gene_time_founder) %>% summarize(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))

beta_z=totf%>% group_by(gene_time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))

#beta_z=totalrares%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
# Extreme bins are >=3/<=-3 SD
bbreaks=c(min(beta_z$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_z$beta_z)+0.1)
beta_z$bin=cut(beta_z$beta_z, breaks=bbreaks, label=F)
beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)
#avg_z=avg_z[!is.na(bin),]

p3=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/all_beta_z_bin_disregulation_smooth_smile.png')
print(p3)
dev.off()

p2=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(method = "lm",mapping=aes(weight=ngenes),formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))


png('eqtl/images/all_beta_z_bin_disregulation_smile.png')
print(p2)
dev.off()

quadratic_high=lm(avg_rares~poly(bin,2),avg_z)
anova(quadratic_high)
summary(quadratic_high)
# poly2 1.56310

###### Looking at lower expressed genes
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10

totalrares=c()
for(time1 in times){
	for(chr in chroms){
		allrares=fread(sprintf('eqtl/results/rare_counts_%s_c%s_max_f_lower_exp.txt',time1,chr),data.table=F)
		allrares=allrares[allrares$max_f!="B73_inra",]
		totalrares=rbind(totalrares,allrares)
	}
	
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
fwrite(totalrares,'eqtl/results/all_rare_counts_max_f_lower_exp.txt',row.names=F,quote=F,sep='\t')

# first I need to break up by gene_time_founder
totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]
totf=totalrares %>% group_by(gene_time_founder) %>% summarize(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))

beta_z=totf%>% group_by(gene_time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))

#beta_z=totalrares%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
# Extreme bins are >=3/<=-3 SD
bbreaks=c(min(beta_z$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_z$beta_z)+0.1)
beta_z$bin=cut(beta_z$beta_z, breaks=bbreaks, label=F)
beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)
#avg_z=avg_z[!is.na(bin),]

p3=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/all_beta_z_bin_disregulation_smooth_smile_lower_exp.png')
print(p3)
dev.off()

p2=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(method = "lm",mapping=aes(weight=ngenes),formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))


png('eqtl/images/all_beta_z_bin_disregulation_smile_lower_exp.png')
print(p2)
dev.off()


#avg_rares=avg_rares[!is.na(avg_rares$add_rank),]

# Build quadratic model and get estimate of smile
quadratic_low=lm(avg_rares~poly(bin,2),avg_z)
anova(quadratic_low)
# poly2 1.31064


# Grab FT genes
totalrares=fread('eqtl/results/all_rare_counts_max_f_all_exp.txt',data.table=F)
totf=totalrares %>% group_by(gene_time_founder) %>% summarize(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))

localcomp=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
ftgenes=ft_genelist$V4

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}
ft_genelist$focal_snp=sapply(seq(1,nrow(ft_genelist)),function(x) all_founder_blocks[all_founder_blocks$chr==ft_genelist$V1[x] & all_founder_blocks$start<=ft_genelist$V2[x] & all_founder_blocks$end>ft_genelist$V3[x],]$focal_snp)
ft_genelist$focal_snp=as.character(ft_genelist$focal_snp)
ft_genelist[ft_genelist$focal_snp=="character(0)",]$focal_snp=sapply(seq(1,nrow(ft_genelist[ft_genelist$focal_snp=="character(0)",])),function(x) tail(all_founder_blocks[all_founder_blocks$chr==ft_genelist[ft_genelist$focal_snp=="character(0)",]$V1[x] & all_founder_blocks$start<=ft_genelist[ft_genelist$focal_snp=="character(0)",]$V2[x],]$focal_snp,1))



fts=c("male_flowering_d6","female_flowering_d6")
cutoff=0.75
# local eQTL with high |r| values
localcomp=localcomp[localcomp$phenotype %in% fts,]
localcomp=localcomp[abs(localcomp$r)>=cutoff,]
ftgenes=c(ftgenes,unique(localcomp$Trait))
ftgenes=unique(ftgenes)
# 936 with local
#distalcomp=fread('QTT/QTL_trans_eQTL_interval_overlap.txt',data.table=F)
#distalcomp=distalcomp[distalcomp$phenotype %in% fts,]
#distalcomp=distalcomp[abs(distalcomp$r)>=cutoff,]

#ftgenes=c(ftgenes,unique(distalcomp$gene))
#ftgenes=unique(ftgenes)
# 1024 genes with distal

ftmarkers=unique(localcomp$SNP)
#ftmarkers=unique(ftmarkers) # 86 markers

ftmarkers=c(unique(ft_genelist$focal_snp),ftmarkers)
ftmarkers=unique(ftmarkers) # 757 (3959 markers left over)

#ftgenes2=c(unique(localcomp$Trait),unique(distalcomp$gene))
#ftgenes2=unique(ftgenes2)
# 153 genes from eQTL



# 18837 genes
allgenes=unique(totalrares$Gene_ID)
# To be a non-FT genes, it needs to not be in the FT genelist AND not have an eQTL
# with the same marker as an FT gene
nftgenes=allgenes[!(allgenes %in% ftgenes)] #5662 to start
# not localcomp
localcomp1=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
#distalcomp1=fread('QTT/QTL_trans_eQTL_interval_overlap.txt',data.table=F)

sub1=localcomp1[localcomp1$Trait %in% nftgenes,]
sub1=sub1[sub1$SNP %in% ftmarkers,]$Trait
nftgenes=nftgenes[!(nftgenes %in% sub1)]
# 16727 genes

# not in distalcomp
#sub2=distalcomp1[distalcomp1$gene %in% nftgenes,]
#sub2=sub2[sub2$SNP %in% ftmarkers,]$gene
#nftgenes=nftgenes[!(nftgenes %in% sub2)]
# 3851 genes

##### FT genes
totf=as.data.frame(totf,stringsAsFactors=F)
ft_rares=totf[totf$Gene_ID %in% ftgenes,]
ft_rares=ft_rares[!is.na(ft_rares$beta_rank),]
print(length(unique(ft_rares$Gene_ID)))
# 593 genes
### Beta Z-scores

beta_z=ft_rares%>% group_by(Gene_ID) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
# Extreme bins are >=3/<=-3 SD
bbreaks=c(min(beta_z$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_z$beta_z)+0.1)
beta_z$bin=cut(beta_z$beta_z, breaks=bbreaks, label=F)
beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)

avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)
#avg_z=avg_z[!is.na(bin),

p3=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme(axis.text.x=element_text(angle=-45)) + theme_classic()

png('eqtl/images/all_beta_z_bin_FT_disregulation_smooth_smile.png')
print(p3)
dev.off()

ft_df=data.frame(Gene_ID=ftgenes,stringsAsFactors=F)
fwrite(ft_df,'eqtl/data/FT_genelist.txt',row.names=F,quote=F,sep='\t')


#### Non-FT genes

nft_rares=totf[totf$Gene_ID %in% nftgenes,]
nft_rares=nft_rares[!is.na(nft_rares$beta_rank),]
print(length(unique(nft_rares$Gene_ID)))
#5225 genes
### Beta Z-scores

beta_z=nft_rares%>% group_by(Gene_ID) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
beta_z$gene_ind=paste0(beta_z$Gene_ID,'-',beta_z$ID)
# Extreme bins are >=3/<=-3 SD
bbreaks=c(min(beta_z$beta_z)-0.1,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,max(beta_z$beta_z)+0.1)
beta_z$bin=cut(beta_z$beta_z, breaks=bbreaks, label=F)
beta_z$gene_time=paste0(beta_z$Gene_ID,'-',beta_z$time)
beta_z$gene_founder=paste0(beta_z$Gene_ID,'-',beta_z$max_f)
avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(gene_founder)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)
nbins=max(avg_z$bin)
#avg_z=avg_z[!is.na(bin),

p3=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(mapping=aes(weight=ngenes),linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,nbins),breaks=c(seq(1,nbins)),labels=c(breaks=c("<-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]",">3.0)"))) +
theme_classic() + scale_size_continuous(name="# of Gene Alleles",breaks=c(500,1000,5000,10000,20000)) +
theme(axis.text.x=element_text(angle=-45))

png('eqtl/images/all_beta_z_bin_nonFT_disregulation_smooth_smile.png')
print(p3)
dev.off()

nft_df=data.frame(Gene_ID=nftgenes,stringsAsFactors=F)
fwrite(nft_df,'eqtl/data/Non_FT_genelist.txt',row.names=F,quote=F,sep='\t')



time1="WD_0727"

#allrares=c()
#for(chr in chroms){
#	rcount=fread(sprintf('eqtl/results/rare_counts_%s_c%s_max_f.txt',time1,chr),data.table=F)
#	allrares=rbind(allrares,rcount)
#}
#allrares=as.data.frame(allrares,stringsAsFactors=F)
#fwrite(allrares,sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),row.names=F,quote=F,sep='\t')



allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
allrares=allrares[allrares$max_f!="B73_inra",]

 #### x-axis 0-1 instead of 1-16
 
 maxrank=allrares %>% group_by(Gene_ID) %>% summarize(maxrank=max(beta_rank,na.rm=T))
 maxrank=as.data.frame(maxrank,stringsAsFactors=F)
 allrares$max_rank=maxrank[match(allrares$Gene_ID,maxrank$Gene_ID),]$maxrank
 
allrares$beta_rank0=1-(allrares$beta_rank-1)/allrares$max_rank
 
avg_rares0=allrares %>% group_by(beta_rank0) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(Gene_ID)))
avg_rares0=as.data.frame(avg_rares0,stringsAsFactors=F)
avg_rares0=avg_rares0[!is.na(avg_rares0$beta_rank0),]
#avg_rares0=avg_rares2[avg_rares2$ngenes>=200,]
quadratic2=lm(avg_rares~poly(beta_rank0,2),avg_rares0,weights=avg_rares0$ngenes)
anova(quadratic2)
 summary(quadratic2)
 
p2=ggplot(aes(x=beta_rank0,y=avg_rares),data=avg_rares0) + geom_point() + stat_smooth(method = "lm", mapping=aes(weight=ngenes),formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Rank (High to Low)") + ylab("Mean Rare Allele Count")
png(sprintf('eqtl/images/beta_rank0_disregulation_smile_%s.png',time1))
print(p2)
dev.off()

#WD_0712 0.68334
#WD_0718 1.10588
#WD_0720 1.09706
#WD_0727 1.27775 

# Z-score

beta_z=allrares%>% group_by(Gene_ID) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
beta_z$bin=cut(beta_z$beta_z, breaks=c(-5.6,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0), label=F)

avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(Gene_ID)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)

p2=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(method = "lm",mapping=aes(weight=ngenes),formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,20),breaks=c(seq(1,20)),labels=c(breaks=c("<-4.5]","(-4.5..-4.0]","(-4.0..-3.5]","(-3.5..-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]","(3.0..3.5]","(3.5..4.0]","(4.0..4.5]",">4.5)"))) +
theme(axis.text.x=element_text(angle=-45))

png(sprintf('eqtl/images/beta_z_bin_disregulation_smile_%s.png',time1))
print(p2)
dev.off()


p3=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,20),breaks=c(seq(1,20)),labels=c(breaks=c("<-4.5]","(-4.5..-4.0]","(-4.0..-3.5]","(-3.5..-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]","(3.0..3.5]","(3.5..4.0]","(4.0..4.5]",">4.5)"))) +
theme(axis.text.x=element_text(angle=-45))

png(sprintf('eqtl/images/beta_z_bin_disregulation_smooth_smile_%s.png',time1))
print(p3)
dev.off()


# Using additive rank
avg_rares=allrares %>% group_by(add_rank) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(Gene_ID)))
avg_rares=as.data.frame(avg_rares,stringsAsFactors=F)

avg_rares=avg_rares[!is.na(avg_rares$add_rank),]
quadratic=lm(avg_rares~poly(add_rank,2),avg_rares)
anova(quadratic)

p1=ggplot(aes(x=add_rank,y=avg_rares),data=avg_rares) + geom_point() + stat_smooth(method = "lm", formula = y ~ poly(x,2), linewidth = 1) +
xlab("Founder Effect Size Rank (High to Low)") + ylab("Mean Rare Allele Count")

png(sprintf('eqtl/images/add_rank_disregulation_smile_%s.png',time1))
print(p1)
dev.off()

# Using founder beta rank
avg_rares2=allrares %>% group_by(beta_rank) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(Gene_ID)))
avg_rares2=as.data.frame(avg_rares2,stringsAsFactors=F)
avg_rares2=avg_rares2[!is.na(avg_rares2$beta_rank),]
avg_rares2=avg_rares2[avg_rares2$ngenes>=200,]
quadratic2=lm(avg_rares~poly(beta_rank,2),avg_rares2)
anova(quadratic2)

p2=ggplot(aes(x=beta_rank,y=avg_rares),data=avg_rares2) + geom_point() + stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Rank (High to Low)") + ylab("Mean Rare Allele Count")
png(sprintf('eqtl/images/beta_rank_disregulation_smile_%s.png',time1))
print(p2)
dev.off()


# Separate out ft genes from non-ft genes
#eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
#eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
#eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
#eqtl=as.data.frame(eqtl2)

#alldf=fread('eqtl/data/founder_frequency.txt',data.table=F)
localcomp=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)

ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
ftgenes=ft_genelist$V4
fts=c("male_flowering_d6","female_flowering_d6")
ftgenes=c(ftgenes,unique(localcomp[localcomp$phenotype %in% fts,]$Trait))
ftgenes=unique(ftgenes)
# 2437 genes
allgenes=unique(allrares$Gene_ID)
nftgenes=allgenes[!(allgenes %in% ftgenes)]
# 16675 genes

# For local-eQTL related to FT, do more extreme founder effect sizes tend to be at higher frequency

#fteqtl=eqtl[eqtl$Trait %in% ftgenes,]
#fteqtl=merge(fteqtl,alldf,by.x='X_ID',by.y='snp')
#fsnps=unique(fteqtl$X_ID)
# 705 markers
# For local-eQTL, do more extreme founder effect sizes tend to be at lower frequency
#nft=eqtl[!(eqtl$Trait %in% ftgenes),]
# eQTL are not linked with FT genes
#nft=nft[!(nft$X_ID %in% fsnps),] #12945

#nft=merge(nft,alldf,by.x='X_ID',by.y='snp')

#nsnps = nft %>% group_by(X_ID) %>% summarize(n=length(unique(Trait)))
#### ALL GENES #####



 
 ######## FT rares  #########
ft_rares=allrares[allrares$Gene_ID %in% ftgenes,]
ft_rares=ft_rares[!is.na(ft_rares$beta_rank),]
print(length(unique(ft_rares$Gene_ID)))

### Beta Z-scores

beta_z=ft_rares%>% group_by(Gene_ID) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
beta_z$bin=cut(beta_z$beta_z, breaks=c(-5.6,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0), label=F)

avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(Gene_ID)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)

p2=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(method = "lm",mapping=aes(weight=ngenes),formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,20),breaks=c(seq(1,20)),labels=c(breaks=c("<-4.5]","(-4.5..-4.0]","(-4.0..-3.5]","(-3.5..-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]","(3.0..3.5]","(3.5..4.0]","(4.0..4.5]",">4.5)"))) +
theme(axis.text.x=element_text(angle=-45))

png(sprintf('eqtl/images/beta_z_bin_FT_disregulation_smile_%s.png',time1))
print(p2)
dev.off()


# Using additive rank
avg_ft_rares2=ft_rares %>% group_by(add_rank) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(Gene_ID)))
avg_ft_rares2=as.data.frame(avg_ft_rares2,stringsAsFactors=F)

avg_ft_rares2=avg_ft_rares2[!is.na(avg_ft_rares2$add_rank),]
ftquadratic2=lm(avg_rares~poly(add_rank,2),avg_ft_rares2)
anova(ftquadratic2)
ft_quad_res2=summary(ftquadratic2)
print(ft_quad_res2)
ft_quad_alpha2=ft_quad_res2$coefficients[3,1]
print(ft_quad_alpha2)

p1=ggplot(aes(x=add_rank,y=avg_rares),data=avg_ft_rares2) + geom_point() + stat_smooth(method = "lm", formula = y ~ poly(x,2), linewidth = 1) +
xlab("cis Additive Effect Rank (High to Low)") + ylab("Mean Rare Allele Count")

png(sprintf('eqtl/images/add_rank_FT_disregulation_smile_%s.png',time1))
print(p1)
dev.off()


# WD_0712 1.347262
# WD_0718 1.924343
# WD_0720 2.173846
# WD_0727 2.017173

# Using founder beta rank
avg_ft_rares=ft_rares %>% group_by(beta_rank) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(Gene_ID)))
avg_ft_rares=as.data.frame(avg_ft_rares,stringsAsFactors=F)
avg_ft_rares=avg_ft_rares[avg_ft_rares$ngenes>=200,]
ftquadratic=lm(avg_rares~poly(beta_rank,2),data=avg_ft_rares)
anova(ftquadratic)
ft_quad_res=summary(ftquadratic)
print(ft_quad_res)
ft_quad_alpha=ft_quad_res$coefficients[3,1]
print(ft_quad_alpha)
# WD_0712 0.4375026
# WD_0712 598 genes
# WD_0718 0.5571328
# WD_0718 595 genes
# WD_0720 0.5608608
# WD_0720 597 genes
# WD_0727 0.4426189
# WD_0727 595 genes


p2=ggplot(aes(x=beta_rank,y=avg_rares),data=avg_ft_rares) + geom_point() + stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Rank (High to Low)") + ylab("Mean Rare Allele Count")

png(sprintf('eqtl/images/beta_rank_FT_disregulation_smile_%s.png',time1))
print(p2)
dev.off()
 
####### Non-FT rares #######
nft_rares=allrares[!(allrares$Gene_ID %in% ftgenes),]
nft_rares=nft_rares[!is.na(nft_rares$beta_rank),]
print(length(unique(nft_rares$Gene_ID)))


### Beta Z-scores

beta_z=nft_rares%>% group_by(Gene_ID) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]
beta_z$bin=cut(beta_z$beta_z, breaks=c(-5.6,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0), label=F)

avg_z=beta_z %>% group_by(bin) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(Gene_ID)))
avg_z=as.data.frame(avg_z,stringsAsFactors=F)

p2=ggplot(aes(x=bin,y=avg_rares),data=avg_z) + geom_point(aes(size=ngenes)) + stat_smooth(method = "lm",mapping=aes(weight=ngenes),formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Z-Score (Low to High)") + ylab("Mean Rare Allele Count") +
scale_x_continuous(limits=c(1,20),breaks=c(seq(1,20)),labels=c(breaks=c("<-4.5]","(-4.5..-4.0]","(-4.0..-3.5]","(-3.5..-3.0]","(-3.0..-2.5]","(-2.5..-2.0]","(-2.0..-1.5]","(-1.5..-1.0]","(-1.0..-0.5]","(-0.5..0.0]","(0.0..0.5]","(0.5..1.0]","(1.0..1.5]","(1.5..2.0]","(2.0..2.5]","(2.5..3.0]","(3.0..3.5]","(3.5..4.0]","(4.0..4.5]",">4.5)"))) +
theme(axis.text.x=element_text(angle=-45))

png(sprintf('eqtl/images/beta_z_bin_nonFT_disregulation_smile_%s.png',time1))
print(p2)
dev.off()


# Using additive rank
avg_nft_rares2=nft_rares %>% group_by(add_rank) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(Gene_ID)))
avg_nft_rares2=as.data.frame(avg_nft_rares2,stringsAsFactors=F)

avg_nft_rares2=avg_nft_rares2[!is.na(avg_nft_rares2$add_rank),]
nftquadratic2=lm(avg_rares~poly(add_rank,2),avg_nft_rares2)
anova(nftquadratic2)
nft_quad_res2=summary(nftquadratic2)
print(nft_quad_res2)
nft_quad_alpha2=nft_quad_res2$coefficients[3,1]
print(nft_quad_alpha2)

p1=ggplot(aes(x=add_rank,y=avg_rares),data=avg_nft_rares2) + geom_point() + stat_smooth(method = "lm", formula = y ~ poly(x,2), linewidth = 1) +
xlab("cis Additive Effect Rank (High to Low)") + ylab("Mean Rare Allele Count")

png(sprintf('eqtl/images/add_rank_nonFT_disregulation_smile_%s.png',time1))
print(p1)
dev.off()
# WD_0712 0.6576233
# WD_0718 1.505524
# WD_0720 2.249892
# WD_0727 2.258254

# Using founder beta rank
avg_nft_rares=nft_rares %>% group_by(beta_rank) %>% summarize(avg_rares=mean(rare_count),sd_rares=sd(rare_count),n=length(rare_count),ngenes=length(unique(Gene_ID)))
avg_nft_rares=as.data.frame(avg_nft_rares,stringsAsFactors=F)
avg_nft_rares=avg_nft_rares[avg_nft_rares$ngenes>=200,]

nftquadratic=lm(avg_rares~poly(beta_rank,2),data=avg_nft_rares)
anova(nftquadratic)
nft_quad_res=summary(nftquadratic)
print(nft_quad_res)
nft_quad_alpha=nft_quad_res$coefficients[3,1]
print(nft_quad_alpha)
# WD_0712 0.2396853
# WD_0712 4384
# WD_0718 0.7219281
# WD_0718 4398 genes
# WD_0720 0.642901
# WD_0720 4396 genes
# WD_0727 0.7290002
# WD_0727 4392 genes

p2=ggplot(aes(x=beta_rank,y=avg_rares),data=avg_nft_rares) + geom_point() + stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1) +
xlab("Founder Effect Size Rank (High to Low)") + ylab("Mean Rare Allele Count")

png(sprintf('eqtl/images/beta_rank_nonFT_disregulation_smile_%s.png',time1))
print(p2)
dev.off()




# Is the relationship between extreme expression and allele frequency 
# different for flowering time genes than other genes?

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

#get_maxf=function(i){
#	row1=rcount[i,]
#	gene=row1$Gene_ID
#	id=row1$ID
#	snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps[1]
#	fprobs=unlist(lapply(X_list,function(x) x[id,snp]))
#	if(max(fprobs)>0.75){
#		max_f=names(which.max(fprobs))
#	}else{
#		max_f=NA
#	}
#	res=results[results$X_ID==snp & results$Trait==gene,]
#	betas=unlist(res[,founders])
#	wn=which(!is.na(betas))[1]
#	betas[-wn]=betas[-wn]+betas[wn]
#	betas=sort(betas,decreasing=TRUE)
#	beta_rank=match(max_f,names(betas))
#	return(list(max_f,betas[max_f],beta_rank))
#}


alldf=fread('eqtl/data/founder_frequency.txt',data.table=F)
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)
eqtl$gene_time_snp=paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)

#### 
fteqtl=eqtl[eqtl$Trait %in% ft_genelist$V4,]
# 530 genes, 1898 total eqtl
time1="WD_0720"

tmpdf=fread(sprintf('eqtl/results/eQTL_%s_freq_chi_data.txt',time1),data.table=F)
tmpdf$gene_time_snp=paste0(tmpdf$Trait,'-',tmpdf$time,'-',tmpdf$X_ID)
tmpdf=tmpdf[tmpdf$gene_time_snp %in% fteqtl$gene_time_snp,]


max_beta=function(row1){
	betas=unlist(row1[,founders])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	betas=betas[!is.na(betas)]
	w=which.max(abs(betas))
	return(betas[w]-mean(betas))
}

mbetas=sapply(seq(1,nrow(tmpdf)),function(x) max_beta(tmpdf[x,]))
tmpdf$max_betas=unlist(mbetas)
tmpdf$max_founder=names(mbetas)

p6=ggplot(tmpdf,aes(x=max_freq,y=max_betas)) + geom_point() + xlab("Max Founder Frequency") + 
ylab("Max Founder Effect Size") + geom_vline(xintercept=0.0625)

png(sprintf('eqtl/images/max_fbetas_freq_FT_% s.png',time1))
print(p6)
dev.off()


# Are more extreme founder alleles at higher frequencies?

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
df=c()
for(time in times){
	print(time)
	ngenes=0
	for(c in 1:10){
		d=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_weights_results_FIXED.txt',time,c))
		d$time=time
		pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
		d$CHR=c
		d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
		#d=d[,c('Trait','X_ID','p_value_ML','CHR','BP','time')]
		df=rbind(df,d)
		
		ngenes=ngenes+length(unique(d$Trait))
	}
	print(ngenes)
}
df=as.data.frame(df,stringsAsFactors=F)


chidf=c()
for(chr in 1:10){
	chi=fread(sprintf('../selection/founder_probs/bg%s_founder_chisq_results.txt',chr),data.table=F)
	chi$chr=chr
	chidf=rbind(chidf,chi)
}	

log_inverse=function(x){
  	return(2^x)
}

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

chidf$end=all_founder_blocks[match(chidf$marker,all_founder_blocks$focal_snp),]$end

####### Rank genes by highest and lowest expressed


time1="WD_0727"
exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
rownames(exp)=exp$V1
exp=exp[,-1]
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

### X-squared pvalue as a a test of segregation distortion


unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
avg_exp = apply(unlog,2,mean)
names(avg_exp)=names(exp)
avg_exp[avg_exp<1]=0
	# re-log the input data
avg_logexp=log2(avg_exp)
avg_logexp[is.infinite(avg_logexp)]=0

expdf=data.frame(Gene_ID=names(avg_logexp),avg_logexp=avg_logexp,stringsAsFactors=F)
expdf=expdf[order(expdf$avg_logexp),]
expdf$rank=seq(nrow(expdf),1)

expdf=merge(expdf,genetable,by="Gene_ID")


markers=c()
p_chis=c()


p1=ggplot(expdf,aes(x=rank,y=-log10(p_chi))) + geom_point() + xlab("Gene Expression Rank") +
ylab("Chi-Squared -log10(p-value)")

png(sprintf('eqtl/images/chi_sq_rank_%s.png',time1))
print(p1)
dev.off()

highrank=expdf[expdf$rank<=5000,]

p2=ggplot(highrank,aes(x=rank,y=-log10(p_chi))) + geom_point() + xlab("Gene Expression Rank") +
ylab("Chi-Squared -log10(p-value)")

png(sprintf('eqtl/images/chi_sq_rank5000_%s.png',time1))
print(p2)
dev.off()

#######  Grab the distance between largest B and avg expression

tmpdf=df[df$time==time1,]
tmpdf=merge(tmpdf,expdf,by.x='Trait',by.y='Gene_ID')


max_beta=function(row1){
	avg=row1$avg_logexp
	betas=unlist(row1[,founders])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	betas=betas[!is.na(betas)]
	diff=avg-betas
	w=which.max(abs(diff))
	return(diff[w])
}

mbetas=sapply(seq(1,nrow(tmpdf)),function(x) max_beta(tmpdf[x,]))
tmpdf$max_betas=unlist(mbetas)
tmpdf$max_founder=names(mbetas)

### Relationship between deviation of founder from mean and expression rank

p3=ggplot(tmpdf,aes(x=rank,y=max_betas)) + geom_point() + xlab("Gene Expression Rank") + 
ylab("Max Founder Deviation")

png(sprintf('eqtl/images/fbetas_rank5000_%s.png',time1))
print(p3)
dev.off()

### Relationship between deviation of founder from mean and chisq p-value

p4=ggplot(tmpdf,aes(x=-log10(p_chi),y=max_betas)) + geom_point() + xlab("Chi-Squared -log10(p-value)") + 
ylab("Max Founder Deviation")

png(sprintf('eqtl/images/fbetas_p_chi_%s.png',time1))
print(p4)
dev.off()


### If there is signfiicant segregation distortion, is the extreme founder higher or lower freq than expected?
max_freq=c()
for(i in 1:nrow(tmpdf)){
	row1=tmpdf[i,]
	snp=row1$X_ID
	f=row1$max_founder
	freq=alldf[alldf$snp==snp & alldf$variable==f,]$value
	max_freq=c(max_freq,freq)
}
tmpdf$max_freq=max_freq

p5=ggplot(tmpdf,aes(x=max_freq,y=max_betas)) + geom_point() + xlab("Max Founder Frequency") + 
ylab("Max Founder Deviation") + geom_vline(xintercept=0.0625)

png(sprintf('eqtl/images/fbetas_freq_%s.png',time1))
print(p5)
dev.off()

max_beta2=function(row1){
	avg=row1$avg_logexp
	snp=row1$X_ID
	freq=alldf[alldf$snp==snp,]
	mid=0.0625
	
	betas=unlist(row1[,founders])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	betas=betas[!is.na(betas)]
	diff=avg-betas
	freq=freq[freq$variable %in% names(betas),]
	freq$dist=freq$value-mid
	#w=which.max(abs(diff))
	w2=which.max(abs(freq$dist))
	wf=freq[w2,]$variable
	return(diff[wf])
}

mbetas2=sapply(seq(1,nrow(tmpdf)),function(x) max_beta2(tmpdf[x,]))
tmpdf$max_betas2=unlist(mbetas2)
tmpdf$max_founder2=names(mbetas2)

max_freq2=c()
for(i in 1:nrow(tmpdf)){
	row1=tmpdf[i,]
	snp=row1$X_ID
	f=row1$max_founder2
	freq=alldf[alldf$snp==snp & alldf$variable==f,]$value
	max_freq2=c(max_freq2,freq)
}
tmpdf$max_freq2=max_freq2

# What is the beta deviation for the most distorted founder?
p6=ggplot(tmpdf,aes(x=max_freq2,y=max_betas2)) + geom_point() + xlab("Max Founder Frequency") + 
ylab("Max Founder Deviation") + geom_vline(xintercept=0.0625)

png(sprintf('eqtl/images/fbetas_freq_deviant_%s.png',time1))
print(p6)
dev.off()

tmpdf2=tmpdf[tmpdf$rank<=5000,]
p7=ggplot(tmpdf2,aes(x=max_freq2,y=max_betas2)) + geom_point() + xlab("Max Founder Frequency") + 
ylab("Max Founder Deviation") + geom_vline(xintercept=0.0625)

png(sprintf('eqtl/images/fbetas_freq_deviant_%s_rank5000.png',time1))
print(p7)
dev.off()


fwrite(tmpdf,sprintf('eqtl/results/eQTL_%s_freq_chi_data.txt',time1),row.names=F,quote=F,sep='\t')

##### Looking specifically at fT genes

ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)

time1="WD_0712"
tmpdf=fread(sprintf('eqtl/results/eQTL_%s_freq_chi_data.txt',time1),data.table=F)
ftdf=tmpdf[tmpdf$Trait %in% ft_genelist$V4,]

p6=ggplot(ftdf,aes(x=max_freq2,y=max_betas2)) + geom_point() + xlab("Max Founder Frequency") + 
ylab("Max Founder Deviation") + geom_vline(xintercept=0.0625)

png(sprintf('eqtl/images/fbetas_freq_FT_%s.png',time1))
print(p6)
dev.off()
