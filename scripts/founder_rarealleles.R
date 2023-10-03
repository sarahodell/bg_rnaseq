#!/usr/bin/env Rscript
#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])

library('data.table')
library('dplyr')
library('stringr')
library('ggplot2')

#haprare=fread(sprintf('datasets/hapmap/chr%s_hapmap_biogemma_rare_alleles_ref.txt',chr),data.table=F)
#names(haprare)=c('ID','CHR','POS','MAJ','MIN')

#bgrare=fread(sprintf('datasets/hapmap/chr%s_founder_rare_alleles.txt',chr),data.table=F)
#names(bgrare)=c("SNP","CHR","POS","REF","ALT","A632_usa","A654_inra","B73_inra","B96","C103_inra","CO255_inra","D105_inra","DK63","EP1_inra","F492", "FV252_inra", "FV2_inra","MBS847","ND245","OH43_inra","VA85","W117_inra")

#raremerge=merge(haprare,bgrare,by='POS')
#print(dim(raremerge))
#print(dim(raremerge[raremerge$MAJ==raremerge$REF,]))


# ames filter for tropical
#tropical=fread('datasets/ames/tropical_vcf.txt',data.table=F,header=F)
#for(chr in 1:10){
#	ames=fread(sprintf('datasets/ames/AmesUSInbreds_AllZeaGBSv1.0_imputed_20130508_chr%s.hmp.txt.gz',chr),data.table=F)
#	tlist=intersect(names(ames),tropical$V1)
#	tropames=ames[,c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode",tlist)]
#	fwrite(tropames,sprintf('datasets/ames/AmesTropical_imputed_chr%s.hmp.txt',chr),row.names=F,quote=F,sep='\t')
#}

#chr="10"

#rares=fread(sprintf('datasets/hapmap/hapmap321_%s_biogemma_rare_alleles.txt',chr),data.table=F)

#names(rares)=c('CHR','POS','REF','ALT',',match')
#rares$SNP=paste0(rares$CHR,'-',rares$POS)

hapsnps=fread(sprintf('datasets/hapmap/chr%s_hapmap_biogemma_rare_alleles_ref.txt',chr),data.table=F)
names(hapsnps)=c("SNP","CHR","POS","MAJ","MIN")
hapsnps=hapsnps[,c("SNP","CHR","POS","MAJ","MIN")]

#hmpfreq=fread(sprintf('datasets/hapmap/hapmap321_%s.frq.counts',chr),data.table=F)
#hmpfreq$freq=hmpfreq$C1/(hmpfreq$C1+hmpfreq$C2)
#pos=sapply(seq(1,nrow(hmpfreq)),function(x) strsplit(hmpfreq$SNP[x],'-')[[1]][[2]])
#hmpfreq$POS=as.numeric(pos)
#hapsnps=merge(hapsnps,hmpfreq,by='SNP')

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

#founder[,founder]=apply(founder[,founder],MARGIN=1,function(x) ifelse(x==0,x$REF,ifelse))

overlap=merge(hapsnps,founder,by.x='POS',by.y='POS')

# remove reference minor alleles 
b73alts=unique(overlap[!is.na(overlap$B73_inra) & overlap$B73_inra!=0,]$B73_inra)
overlap=overlap[overlap$B73_inra!=1,]

# remove multi-allelic sites
overlap=overlap[is.na(overlap$ALT2),]

#overlap=overlap[overlap$ALT.x==overlap$ALT.y,]
# Is the reference allele minor or major allele in hapmap?
#refminor=c()
#for(i in 1:nrow(overlap)){
#	row1=overlap[i,]
#	minor=row1$A1
#	major=row1$A2
#	ref=row1$REF
#	alt1=row1$ALT
#	alt2=row1$ALT2
#	if(minor==ref){
#		refminor=c(refminor,'REF')
#	}else if(minor==alt1){
#		refminor=c(refminor,'ALT')
#	}else if(!is.na(alt2) & minor==alt2){
#		refminor=c(refminor,'ALT')
#	}else{
#		refminor=c(refminor,NA)
#	}
#}


#overlap$MINOR=refminor

#10-6787144

# Grab minor allele, grab founders with that allele

# Which founders have the minor allele?

#names(which(row1[,founders])==1)
#min_founder=list()
#for(i in 1:nrow(overlap)){#
#	row1=overlap[i,]
#	pos=paste0('S',row1$CHR.x,'_',row1$POS)
#	minf=founders[which(row1[,founders]==1)]
#	min_founder[[pos]]=minf
#}

#saveRDS(min_founder,sprintf('datasets/hapmap/minor_allele_founders_c%s.rds',chr))

min_founder=readRDS(sprintf('datasets/hapmap/minor_allele_founders_c%s.rds',chr))
bimbam=fread(sprintf('../genotypes/probabilities/allele_probs/bg%s_wgs_alleleprobs_bimbam.txt',chr),data.table=F)
K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

#Drop samples without phenotype data
names(bimbam)=c('marker','alt1','ref',rownames(K))
bimbam=bimbam[bimbam$marker %in% names(min_founder),]
bimbam$pos=sapply(seq(1,nrow(bimbam)),function(x) as.numeric(strsplit(bimbam$marker[x],'_')[[1]][[2]]))

fwrite(bimbam,sprintf('datasets/chr%s_biogemma_rare_allele_probs.txt',chr),row.names=F,quote=F,sep='\t')

#which lines are homozygous for a rare allele?

# Use # of homozygous rare alleles - correlation with fitness (grain yield or thousand kernel weight BLUPs)

# 

#grab top 5000 expressed genes, rank by individual expression
#time1="WD_0727"
#tmpdf=fread(sprintf('eqtl/results/eQTL_%s_freq_chi_data.txt',time1),data.table=F)
#tmpdf$gene_time_snp=paste0(tmpdf$Trait,'-',tmpdf$time,'-',tmpdf$X_ID)

#top5k=unique(tmpdf[tmpdf$rank<=5000,]$Trait)



# What are the frequencies of these rare alleles in the founders? In the MAGIC lines?

#ranks=fread(,sprintf('eqtl/results/%s_top5k_exp_ranks.txt',time1),data.table=F)
#rownames(ranks)=ranks$V1
#genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
#genetable=genetable[genetable$CHROM==chr,]
#chr5k=top5k[top5k %in% genetable$Gene_ID]

#ranks=ranks[,chr5k]

#min_founder=readRDS(sprintf('datasets/hapmap/minor_allele_founders_c%s.rds',chr))
#bimbam=fread(sprintf('../genotypes/probabilities/allele_probs/bg%s_wgs_alleleprobs.txt.gz',chr),data.table=F)
#bimbam=bimbam[bimbam$marker %in% names(min_founder),]
#bimbam$pos=sapply(seq(1,nrow(bimbam)),function(x) as.numeric(strsplit(bimbam$marker[x],'_')[[1]][[2]]))

#rcounts=as.data.frame(matrix(nrow=length(inds),ncol=length(chr5k)),stringsAsFactors=F)
#rownames(rcounts)=inds
#colnames(rcounts)=chr5k
#inds=rownames(ranks)
#rmelt=as.data.frame(expand.grid(chr5k,inds),stringsAsFactors=F)
#names(rmelt)=c("Gene_ID","ID")
#cutoff=0.75

#ranklist=c()
#countlist=c()
#for(i in 1:nrow(rmelt)){
#	row1=rmelt[i,]
#	gene=as.character(row1$Gene_ID)
#	id=as.character(row1$ID)
#	r=ranks[id,gene]
#	
#	rend=genetable[genetable$Gene_ID==gene,]$START
#	rstart=rend-5000
#	t=bimbam[bimbam$pos>=rstart & bimbam$pos<=rend,]
#	t=t[,c('marker','alt1','ref',id)]
#	rcount=sum(t[,id]>=cutoff)
#	ranklist=c(ranklist,r)
#	countlist=c(countlist,rcount)
#}

#rmelt$rank=ranklist
#rmelt$rare_count=countlist
#fwrite(rmelt,sprintf('eqtl/results/%s_%s_5kb_rare_counts.txt',time1,chr),row.names=T,quote=F,sep='\t')

#avg_rares=rmelt %>% group_by(rank) %>% summarize(avg=mean(rare_count),sd=sd(rare_count))
#
#p1=ggplot(avg_rares,aes(x=rank,y=avg)) + geom_point() + xlab("Expression Rank") + ylab("Average Rare Allele Count")
#
#png(sprintf('eqtl/images/%s_c%s_rare_alleles_by_exp.png',time1,chr))
#print(p1)
#dev.off()
##rare_count=function(gene,i){
##	rend=genetable[genetable$Gene_ID==gene,]$START
##	rstart=rend-5000
##	t=bimbam[bimbam$pos>=rstart & bimbam$pos<=rend,]
##	t=t[,c('marker','alt1','ref',i)]
##	rcount=sum(t[,i]>=0.75)
##	return(rcount)
##}
##inds=rownames(ranks)
##rcounts=sapply(chr5k,function(i) sapply(inds,function(j) rare_count(i,j)))
## Values are prob of alternate allele, # grab lines with high certainty
#
#
#eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
#eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
## Grab only the highest cis SNP
#eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
#eqtl=as.data.frame(eqtl2)
#eqtl$gene_time_snp=paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
#
#all_founder_blocks=c()
#for(chr in 1:10){#
#  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
#  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
#}
#eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
#eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end
#eqtl=eqtl[eqtl$CHR==chr,]
#
#genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
#eqtl=merge(eqtl,genetable,by.x='Trait',by.y='Gene_ID')
#rare_count=function(row1){
#	start=row1$START-100000
#	end=row1$END+100000
#	t=bimbam[bimbam$pos>=start & bimbam$pos<=end,]
#	return(nrow(t))
#}
#
#eqtl$n_rares=sapply(seq(1,nrow(eqtl)),function(x) rare_count(eqtl[x,]))
#
#tmpdf=fread('eqtl/results/eQTL_WD_0712_freq_chi_data.txt',data.table=F)
#tmpdf$gene_time_snp=paste0(tmpdf$Trait,'-',tmpdf$time,'-',tmpdf$X_ID)
#time1="WD_0712"
#teqtl=eqtl[eqtl$time==time1,]
#
#teqtl=merge(teqtl,tmpdf,by='gene_time_snp')
#
#bins=seq(0,17000,1000)
#names(bins)=seq(1,length(bins))
#
#teqtl=teqtl[order(teqtl$rank),]
#rownames(teqtl)=seq(1,nrow(teqtl))
#
#teqtl$bin=sapply(seq(1,nrow(teqtl)),function(x) (which(teqtl$rank[x]<bins)[1]-1))
#teqtl$bin=as.factor(teqtl$bin)
#
#p1=ggplot(teqtl,aes(x=bin,y=n_rares)) + geom_point() + xlab("Rank Gene Expression") +
#ylab("# of rares alleles")
#png(sprintf('eqtl/images/%s_c%s_rare_alleles_by_exp.png',time1,chr))
#print(p1)
#dev.off()
# Grabbing individuals genotypes at these sites
# 

# Are any of the eQTL SNPs rare variants?


#X_list=readRDS('../genotypes/probabilities/geno_probs/bg10_filtered_genotype_probs.rds')
#pmap=fread('../genotypes/qtl2/startfiles/Biogemma_pmap_c10.csv',data.table=F)
# Which individuals have the the rare allele? 


# 
#exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
#exp=exp[,c('V1',top5k)]

#inds=exp$V1

#find_rank=function(gene){
#	sub=exp[order(exp[,gene],decreasing=TRUE),c('V1',gene)]
#	r=match(inds,sub$V1)
#	return(r)
#}
#ranks=as.data.frame(unlist(sapply(top5k,function(x) find_rank(x))),stringsAsFactors=F)
#rownames(ranks)=inds
#colnames(ranks)=top5k

#fwrite(ranks,sprintf('eqtl/results/%s_top5k_exp_ranks.txt',time1),row.names=T,quote=F,sep='\t')
