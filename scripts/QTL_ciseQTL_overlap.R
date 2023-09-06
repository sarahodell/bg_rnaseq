library('abind')
library('data.table')
library('ggplot2')
library('dplyr')

##
time="WD_0712"
deg=fread(sprintf('limma_results/%s_top_MITE_DEGs.txt',time))
eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)


eqtl=eqtl[eqtl$time==time,]

eqtl= eqtl[eqtl$Trait %in% deg$V1,]

eqtl= eqtl[eqtl$Trait %in% deg$Gene_ID,]


all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}
eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end


qtl=fread('QTL/all_adjusted_QTL_support_intervals.txt',data.table=F)
qtl$block_start=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$start
qtl$block_end=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$end
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comp=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)


founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
plot_list=list()
count=1
#adj_chr=c(5,9)
cors=c()
pvals=c()
for(i in 1:nrow(comp)){
	row=comp[i,]
	pheno=row$phenotype
	env=row$environment
	chr=row$CHR
	gene=row$Trait
	qsnp=row$SNP
	#id=row$ID
	time=row$time
	#exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	gene=row$Trait
	esnp=row$X_ID
	#chr=row$CHR
	#if(chr %in% adj_chr){
	#	X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
	#}else{
	#	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	#
	#}
	#inds=rownames(X_list[[1]])
	#inter=intersect(exp$V1,inds)
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	#X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	#colnames(X) = founders
	#rownames(X) = inter
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	#w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	#results=results[[w]]
	results=results[results$X_ID==esnp & results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	pvals=c(pvals,p)
	
	df=data.frame(founder=founders,pheno=effect_size,eqts=betas,stringsAsFactors=F)
	df=df[complete.cases(df),]
	p1=ggplot(df,aes(x=eqts,y=pheno)) + geom_point(aes(color=founder)) +
	xlab("eQTL effect size (log2cpm)") + ylab("QTL effect size") +
	ggtitle(sprintf("%s %s and %s %s, r=%.2f",gene,time,pheno,env,r))
	plot_list[[count]]=p1
	count=count+1
	
}

comp$r=cors
comp$pvalue=pvals

pdf(sprintf('QTL/images/%s_MITE_DEG_effect_size_correlations.pdf',time))
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()



eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

# check that the genes with more than one SNP local eQTL that the effect sizes are correlated - of course they are, but check

#sub=fread('eqtl/results/all_cis_trans_fdr_hits.txt',data.table=F)
#factordf=fread('eqtl/results/all_factor_trans_eqtl_fdr_genes.txt',data.table=F)

eqtl %>% group_by(time) %>% count()
# A tibble: 4 × 2
# Groups:   time [4]
#  time        n
#  <chr>   <int>
#1 WD_0712 10674
#2 WD_0718 15081
#3 WD_0720 16381
#4 WD_0727 16114




gtime= eqtl %>% group_by(Trait) %>% summarize(length(unique(time)))

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}
eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end

# Just for qDTS3_2 aLL

qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
qtl=qtl[qtl$ID=="qDTS3_2",]
qtl=qtl[2,]
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
# 22 overlaps, 6 genes

#### Within the Founder LD r2=0.5 cutoff interval
qtl=fread('QTL/all_adjusted_QTL_peaks.txt',data.table=F)
qtl$block_start=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$start
qtl$block_end=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$end
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,leftmost,alt_rightmost)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)
# Overlap within every QTL ID

comparison2 %>% group_by(ID) %>% summarize(length(unique(Trait)))
# A tibble: 23 × 2
#   ID      `length(unique(Trait))`
#   <chr>                     <int>
# 1 qDT33_2                     130
# 2 qDTA3_1                      14
# 3 qDTA3_2                     163
# 4 qDTA7                         9
# 5 qDTA8                       435
# 6 qDTA9                        62
# 7 qDTS3_1                       1
# 8 qDTS3_2                     228
# 9 qDTS8                       320
#10 qDTS9                        57


env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,block_start,block_end)
comparison1=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
# 18 of 23 QTL have cis-eQTL in the same LD window
# A tibble: 18 × 2
#   ID      `length(unique(Trait))`
#   <chr>                     <int>
# 1 qDT33_2                       4
# 2 qDTA3_1                       3
# 3 qDTA3_2                      54
# 4 qDTA8                        84
# 5 qDTA9                        11
# 6 qDTS3_1                       1
# 7 qDTS3_2                       9
# 8 qDTS8                        67
# 9 qDTS9                         8
#10 qHGM1                         2
#11 qHGM3_1                       1
#12 qHGM3_2                      15
#13 qHGM7                         2
#14 qTKW2                         4
#15 qTKW7_1                       9
#16 qTKW7_2                       2
#17 qTPH6                         1
#18 qTPH8                         2

# Just in EXP_STPAUL_2017_WD

 inenv %>% group_by(ID) %>% summarize(length(unique(Trait)))
# A tibble: 4 × 2
#  ID      `length(unique(Trait))`
#  <chr>                     <int>
#1 qDTA3_2                      11
#2 qDTA8                         8
#3 qDTA9                         7
#4 qDTS8                        23


### DO any of these have correlated effect sizes?
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
plot_list=list()
count=1
adj_chr=c(5,9)
cors=c()
pvals=c()
for(i in 1:nrow(inenv)){
	row=inenv[i,]
	pheno=row$phenotype
	env=row$environment
	chr=row$CHR
	gene=row$Trait
	snp=row$SNP
	time=row$time
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	gene=row$Trait
	snp=row$SNP
	chr=row$CHR
	if(chr %in% adj_chr){
		X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
	}else{
		X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))

	}
	inds=rownames(X_list[[1]])
	inter=intersect(exp$V1,inds)
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	#w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	#results=results[[w]]
	results=results[results$X_ID==snp& results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	pvals=c(pvals,p)
	
	df=data.frame(founder=founders,pheno=effect_size,eqts=betas,stringsAsFactors=F)
	if(p<=0.1){
		p1=ggplot(df,aes(x=eqts,y=pheno)) + geom_point(aes(color=founder))
		plot_list[[count]]=p1
		count=count+1
	}
	
}
inenv$r=cors
inenv$pvalue=pvals

ids=unique(inenv$ID)
for(i in ids){
	sub=inenv[inenv$ID==i,]
	loc=which.max(abs(sub$r))
	hval=sub[loc,]$i.value
	hr=sub[loc,]$r
	h1=ggplot(sub,aes(x=i.value)) + geom_histogram() + geom_vline(xintercept) +
	ggtitle(sprintf('eQTL qvalues within %s',i)) + theme(label=sprintf('%s r=%.2f'))
}


pdf('QTT/QTL_eQTL_effect_size_correlations.png')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()
##### Within the same LD block

##### cis-eQTL overlap with QTL in EXP_STPAUL_2017
interval=c()
highest=c()
# Overlap of cis eQTL and St.Paul recomb blocks peaks
# in the QTL interval
qtl2=fread('QTL/adjusted/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,leftmost,alt_rightmost)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)
comparison2$phenotype="male_flowering_d6"
interval=comparison2

# 766 instances of overlap, 281 genes
comparison2 %>% group_by(SNP) %>% count()
# A tibble: 4 × 2
# Groups:   SNP [4]
#  SNP               n
#  <chr>         <int>
#1 AX-91100761     176 # qDTA8_1
#2 AX-91145110      60	qDTA9
#3 AX-91771656     502	qDTA8_2
#4 PZE-103093413    28	qDTA3_2

# for the highest QTL peak
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
comparison2$phenotype="male_flowering_d6"
#78 instances of overlap of 28 genes
highest=comparison2
comparison2 %>% group_by(SNP) %>% count()

# A tibble: 4 × 2
# Groups:   SNP [4]
#  SNP               n
#  <chr>         <int>
#1 AX-91100761    13 # qDTA8_1
#2 AX-91145110    25 qDTA9
#3 AX-91771656    12	qDTA8_2
#4 PZE-103093413  28	qDTA3_2

#DTS
qtl2=fread('QTL/adjusted/female_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,leftmost,alt_rightmost)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)
comparison2$phenotype="female_flowering_d6"
# 244 instances of overlap, 92 genes
comparison2 %>% group_by(SNP) %>% count()
# A tibble: 2 × 2
# Groups:   SNP [2]
#  SNP             n
#  <chr>       <int>
#1 AX-91772415   141 qDTA8_2
#2 PZA01038.1    103 qDTA8_1

interval=rbind(interval,comparison2)


env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,block_start,block_end)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
comparison2$phenotype="female_flowering_d6"
# 52 overlaps of 23 genes
highest=rbind(highest,comparison2)

comparison2 %>% group_by(SNP) %>% count()
# A tibble: 2 × 2
# Groups:   SNP [2]
#  SNP             n
#  <chr>       <int>
#1 AX-91772415     3 qDTA8_2
#2 PZA01038.1     49 qDTA8_1

#tkw_15
qtl2=fread('QTL/adjusted/tkw_15_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',data.table=F)
qtl2$block_start=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$start
qtl2$block_end=all_founder_blocks[match(qtl2$SNP,all_founder_blocks$focal_snp),]$end
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl2)
setkey(env2,CHR,leftmost,alt_rightmost)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','leftmost','alt_rightmost'),nomatch=NULL)
# 9 overlaps with 4 genes
comparison2$phenotype="tkw_15"

interval=rbind(interval,comparison2)
highest=rbind(highest,comparison2)


fwrite(highest,'QTT/cis_eQTL_STPAUL_QTL_peak_overlaps.txt',row.names=F,quote=F,sep='\t')
fwrite(interval,'QTT/cis_eQTL_STPAUL_QTL_overlaps.txt',row.names=F,quote=F,sep='\t')


# Other phenotypes

phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:9]
thresh=0.10
all_qtl=c()
for(p in phenos){
	for(e in envs){
		hitfile=sprintf('QTL/adjusted/%s_%s_QTL_scan_%.2f_peaks.txt',p,e,thresh)
		if(file.exists(hitfile)){
			df=fread(hitfile,data.table=F)
			df$phenotype=p
			df$environment=e
			all_qtl=rbind(all_qtl,df)
		}
	}
}
all_qtl=all_qtl[,c("phenotype","environment","CHR","BP","SNP","leftmost","alt_rightmost","leftmost_SNP","rightmost_SNP","rightmost","value")]
#fwrite(all_qtl,'QTL/all_adjusted_QTL_peaks.txt',row.names=F,quote=F,sep='\t')

# Pleiotropic QTL
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

all_qtl=fread('QTL/all_adjusted_QTL_peaks.txt',data.table=F)

####### qDTS3_2 and qHGM3_2 ####### 
chr="3"
pheno1="female_flowering_d6"
env1="ALL"
pheno2="harvest_grain_moisture"
env2="ALL"
qtl1="qDTS3_2"
qtl2="qHGM3_2"

snp1="PZE-103084819"
snp2="PZE-103084819"
#r=0.6841095

####### #qDTA7 and qTKW7_1 ####### 
chr="7"
pheno1="male_flowering_d6"
env1="GRANEROS_2015_OPT"
pheno2="tkw_15"
env2="BLOIS_2017_OPT"
qtl1="qDTA7"
qtl2="qTKW7_1"

snp1="AX-91715233"
#snp2="SYN16002"
#r=0.1599592
snp2="AX-91714901"
#r=0.5285355

####### #qHGM7 and qTKW7_2 in GRANEROS ####### 
chr="7"
pheno1="harvest_grain_moisture"
env1="GRANEROS_2015_OPT"
pheno2="tkw_15"
env2="GRANEROS_2015_OPT"

snp1="AX-91053569"
snp2="AX-91053569"
qtl1="qHGM7"
qtl2="qTKW7_2"
#r=0.2688263

##########
res1=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno1,env1),data.table=F)
peaks1=fread(sprintf('QTL/adjusted/%s_%s_QTL_scan_0.10_peaks.txt',pheno1,env1),data.table=F)
peaks1=peaks1[peaks1$CHR==chr,]


res2=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno2,env2),data.table=F)
peaks2=fread(sprintf('QTL/adjusted/%s_%s_QTL_scan_0.10_peaks.txt',pheno2,env2),data.table=F)
peaks2=peaks2[peaks2$CHR==chr,]



res1=res1[res1$X_ID==snp1,]
beta1=unlist(res1[,founders])
beta1[-1]=beta1[-1]+beta1[1]
res2=res2[res2$X_ID==snp2,]
beta2=unlist(res2[,founders])
beta2[-1]=beta2[-1]+beta2[1]

r=cor(beta1,beta2,use="complete.obs")

df=data.frame(founder=founders,beta1=beta1,beta2=beta2)
p1=ggplot(df,aes(x=beta1,y=beta2,color=founder)) + geom_point() +
xlab(sprintf('%s %s effect size',env1,qtl1)) + ylab(sprintf('%s %s effect size',env2,qtl2)) + ggtitle(sprintf("Founder Effect Size Correlation r = %.2f",r))

png(sprintf('QTL/images/%s_%s_effect_sizes.png',qtl1,qtl2))
print(p1)
dev.off()
####


qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
fqtl=qtl[qtl$Method=="Founder_probs",]
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

fqtl$block_start=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$start
fqtl$block_end=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$end


#Look 5kb upstream or downstream of gene
#genetable$window_START=genetable$START-5000
#genetable$window_END=genetable$END+5000

# Does it make more sense to use 10kb region around the significant SNP
#ciseqtl$window_START=ciseqtl$BP-5000
#ciseqtl$window_END=ciseqtl$BP+5000
# Or use the founder LD block? I think this one is more accurate

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
phenos=unique(full_comparison$phenotype)
cors=c()
ps=c()
#for(i in 1:nrow(full_comparison)){
for(i in 1:nrow(all_comparison)){
	row=all_comparison[i,]
	#row=full_comparison[i,]

	time=row$time
	pheno=row$phenotype
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	gene=row$Trait
	snp=row$SNP
	chr=row$CHR
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])
	inter=intersect(exp$V1,inds)
	effect_sizes=fread(sprintf('QTL/Biogemma_chr%s_%s_x_EXP_STPAUL_2017_WD_vst_founderprobs.txt',chr,pheno),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	#w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	#results=results[[w]]
	results=results[results$X_ID==snp& results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	ps=c(ps,p)
	
	df=data.frame(founder=founders,pheno=effect_size,eqts=betas,stringsAsFactors=F)
	if(p<=0.1){
		p1=ggplot(df,aes(x=eqts,y=pheno)) + geom_point(aes(color=founder))
		png(sprintf('QTT/%s_%s_%s_%s_effect_sizes.png',time,gene,snp,pheno))
		print(p1)
		dev.off()
	}
	
	#subpheno=phenotypes[phenotypes$ID %in% inter,c('ID','female_flowering_d6','male_flowering_d6','tkw_15')]
	#fdata$dta=subpheno[match(fdata$ID,subpheno$ID),]$male_flowering_d6
	#fdata$dts=subpheno[match(fdata$ID,subpheno$ID),]$female_flowering_d6

}

all_comparison$cor=cors
all_comparison$pvalue=ps

all_comparison=all_comparison[order(all_comparison$pvalue),]
rownames(all_comparison)=seq(1,nrow(all_comparison))
all_comparison$p_adjust=p.adjust(all_comparison$pvalue,method="fdr")

png('QTT/overlap_qqplot.png')
print(qqman::qq(all_comparison$pvalue))
dev.off()

fwrite(all_comparison,'QTT/cis_eQTL_STPAUL_QTL_overlaps.txt',row.names=F,quote=F,sep='\t')
# most are related to DTA
# Two most correlated are Zm00001d011459 in WD0727 (r=-0.778)- a AAA ATPase related to hypersensitivity
# and Zm00001d011225 (WD_0720), a giberrellin sensitivity gene (r=0.742)

# Zm00001d010804 (WD_0727), Proteosome maturation factor UMP1
# What is a good null expectation?



# What about QTL from other environments?
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)


phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

# First, DTA and trans-eQTL peaks
phenos=unique(envcomp$Phenotype)
cors=c()
ps=c()
#for(i in 1:nrow(full_comparison)){
for(i in 1:nrow(envcomp)){
	row=envcomp[i,]
	#row=full_comparison[i,]
	env=row$Environment
	time=row$time
	pheno=row$Phenotype
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	gene=row$Trait
	snp=row$highest_SNP
	chr=row$CHR
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	inds=rownames(X_list[[1]])
	inter=intersect(exp$V1,inds)
	effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%s_%s_x_%s_vst_founderprobs.rds',chr,pheno,env))
	#effect_sizes=fread(sprintf('../GridLMM/QTL/Biogemma_chr%s_%s_x_EXP_STPAUL_2017_WD_vst_founderprobs.txt',chr,pheno),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
	effect_size=unlist(effect_size[,6:21])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))
	colnames(X) = founders
	rownames(X) = inter
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	#w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene)))
	#results=results[[w]]
	results=results[results$X_ID==snp& results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	cors=c(cors,r)
	ps=c(ps,p)
	
	df=data.frame(founder=founders,pheno=effect_size,eqts=betas,stringsAsFactors=F)
	if(p<=0.05){
		p1=ggplot(df,aes(x=eqts,y=pheno)) + geom_point(aes(color=founder))
		png(sprintf('QTT/%s_%s_%s_%s_%s_effect_sizes.png',time,gene,snp,pheno,env))
		print(p1)
		dev.off()
	}
	
	#subpheno=phenotypes[phenotypes$ID %in% inter,c('ID','female_flowering_d6','male_flowering_d6','tkw_15')]
	#fdata$dta=subpheno[match(fdata$ID,subpheno$ID),]$male_flowering_d6
	#fdata$dts=subpheno[match(fdata$ID,subpheno$ID),]$female_flowering_d6

}

envcomp$cor=cors
envcomp$pvalue=ps

all_comparison=all_comparison[order(all_comparison$pvalue),]
rownames(all_comparison)=seq(1,nrow(all_comparison))
all_comparison$p_adjust=p.adjust(all_comparison$pvalue,method="fdr")

png('QTT/overlap_qqplot.png')
print(qqman::qq(all_comparison$pvalue))
dev.off()

fwrite(all_comparison,'QTT/cis_eQTL_STPAUL_QTL_overlaps.txt',row.names=F,quote=F,sep='\t')

#full_comparison$cor=cors
#full_comparison$pvalue=ps
#fwrite(full_comparison,'QTT/cis_eQTL_STPAUL_QTL_peak_overlaps.txt',row.names=F,quote=F,sep='\t')


#full_comparison$p_adjust=p.adjust(full_comparison$pvalue,method="fdr")


# Overlap of eQTL variants with eQTL genes 
env1=sub[,c('Trait','CHR','block_start','block_end','SNP','class')]
names(env1)=c('Trait.variant','CHR','block_start','block_end','SNP.variant','class.variant')
env1=as.data.table(env1)
env2=sub[,c('Trait','SNP','class','gene_chr','gene_start','gene_end')]
names(env2)=c('Trait.transcript','SNP.transcript','class.transcript','gene_chr','gene_start','gene_end')
env2=as.data.table(env2)
setkey(env2,gene_chr,gene_start,gene_end)
comparison3=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('gene_chr','gene_start','gene_end'),nomatch=NULL)

comparison3$gene1_gene2=paste0(comparison3$Trait.transcript,'-',comparison3$Trait.variant)
length(unique(comparison3$gene1_gene2))
# 90 (45) instances of eQTL variants overlapping with eQTL genes for other trans-eQTL

# Overlap of eQTL variants with each other
env1=sub[,c('Trait','CHR','block_start','block_end','SNP','class','time')]
names(env1)=c('Trait.1','CHR','block_start','block_end','SNP.1','class.1','time.1')
env1=as.data.table(env1)
env2=sub[,c('Trait','CHR','block_start','block_end','SNP','class','time')]
names(env2)=c('Trait.2','CHR','block_start','block_end','SNP.2','class.2','time.2')
env2=as.data.table(env2)
setkey(env2,CHR,block_start,block_end)
comparison4=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)
comparison4=comparison4[comparison4$Trait.1!=comparison4$Trait.2,]

geneinter4=intersect(comparison4$Trait.1,comparison4$Trait.2)
length(unique(comparison4$gene1_gene2))
# 418 (209) instances of eQTL variants controlling multiple genes in trans
# 904 (452) of the shared SNPs are in the same timepoint
# 324 (162) genes share a trans-eQTL variant in the same time point

#Across all timepoints
snp_avg=comparison4 %>% group_by(SNP.1) %>% summarize(avg=length(unique(gene1_gene2)))
summary(snp_avg$avg/2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   1.000   1.000   3.116   3.000  28.000	

# The WD_0712 cis-eQTL variant on CHR 1 overlaps with a WD_0727 trans-eQTL variant for Zm00001d047592
# This gene also has two large trans-eQTL on chr 5 in WD_0727 (this is the big inter-LD region)

# how often is a SNP an eQTL in multiple timepoints?
snptime=sub %>% group_by(SNP) %>% summarize(ntimes=length(unique(time)))
summary(snptime$ntimes)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   1.000   1.000   1.065   1.000   3.000
# There were 60 eQTL variants that were eQTL in more than one timepoint
# There was one SNP that was an eQTL in three different timepoints (AX-90719118)
# It controlled three different genes in WD_0712, WD_0718, and WD_0727


# How many genes were different SNPs associated with
genesnp= sub %>% group_by(SNP,time) %>% summarize(ngenes=length(unique(Trait)))
summary(genesnp$ngenes)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.00    1.00    1.00    1.23    1.00    7.00 
genesnp[genesnp$ngenes==7,]
#  SNP         time    ngenes
#  <chr>       <chr>    <int>
#1 AX-90610207 WD_0727      7
#2 AX-90794853 WD_0727      7
#3 AX-90965986 WD_0727      7
#4 AX-91673845 WD_0727      7

# 173 SNPs associated with more than one gene, (111 with 2, 14 with 3, 12 with 4,two with 5, 9 with 6, 3 with 7, 2 with 8)


# Were any genes trans-eQTL across timepoints?
genetime=sub %>% group_by(Trait) %>% summarize(ntimes=length(unique(time)))
summary(genetime$ntimes)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   1.000   1.000   1.009   1.000   2.000
#genetime[genetime$ntimes==2,]
# A tibble: 3 × 2
#  Trait          ntimes
#  <chr>           <int>
#1 Zm00001d049691      2
#2 Zm00001d051459      2
#3 Zm00001d051801      2


## Factors
factorgenes=fread('eqtl/results/all_factor_trans_eqtl_fdr_genes.txt',data.table=F)
egenes=unique(sub$Trait)
sum(egenes %in% factorgenes$Gene)
infactors=egenes[egenes %in% factorgenes$Gene]
# There are 77 trans-eQTL genes that load onto factors with factor trans-eQTL
subf=factorgenes[factorgenes$Gene %in% egenes,]
subf %>% group_by(Trait) %>% count()
# A tibble: 3 × 2
# Groups:   Trait [3]
#  Trait        n
#  <chr>    <int>
#1 Factor14   235
#2 Factor16     2
#3 Factor2     50


### overlap of eQTL and factor eQTL
env1=sub[,c('Trait','CHR','block_start','block_end','SNP','class')]
names(env1)=c('Trait.variant','CHR','block_start','block_end','SNP.variant','class.variant')
env1=as.data.table(env1)
env2=factoreqtl[,c('Trait','CHR','SNP','class','time','block_start','block_end')]
names(env2)=c('Trait.factor','CHR','SNP.factor','class.factor','time.factor','block_start','block_end')
env2=as.data.table(env2)
setkey(env2,CHR,block_start,block_end)
comparison5=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','block_start','block_end'),nomatch=NULL)

factorg=comparison5 %>% group_by(Trait.factor) %>% summarize(ngenes=length(unique(Trait.variant)))
#  Trait.factor ngenes
#  <chr>         <int>
#1 Factor14          9
#2 Factor2           1

# In factor 2, "Zm00001d010440" is overlapping with the factor trans-eQTL (prop_var = 0.1194212)

f14g=c("Zm00001d027804","Zm00001d034716", "Zm00001d049577", "Zm00001d050770",
 "Zm00001d013614", "Zm00001d015242", "Zm00001d021130" ,"Zm00001d021854",
  "Zm00001d009256")
f14props=c(0.436,0.671,0.862,0.756,0.503,0,0.879,0.698,0.453) # all but "Zm00001d015242" are loaded on Factor 14

# Overlap of cis and trans eQTL genes with factor genes
#env1=sub[,c('Trait','gene_chr','gene_start','gene_end','SNP','class','time')]
#names(env1)=c('Trait.1','gene_chr','gene_start','gene_end','SNP.1','class.1','time1')
#env1=as.data.table(env1)
#env2=factorgenes[,c('Trait','Gene','gene_chr','gene_start','gene_end','SNP','class','time')]
#names(env2)=c('Trait.factor','Gene','gene_chr','gene_start','gene_end','SNP.factor','class.factor','time.factor')
#env2=as.data.table(env2)
#setkey(env2,gene_chr,gene_start,gene_end)
#comparison5=foverlaps(env1,env2,by.x=c('gene_chr','gene_start','gene_end'),by.y=c('gene_chr','gene_start','gene_end'),nomatch=NULL)


# Overlap of cis and trans variants with factor genes
#env1=eqtl[,c('Trait','CHR','block_start','block_end','SNP','class','time')]
#names(env1)=c('Trait.1','CHR','block_start','block_end','SNP.1','class.1','time1')
#env1=as.data.table(env1)
#env2=factorgenes[,c('Trait','Gene','gene_chr','gene_start','gene_end','SNP','class','time')]
#names(env2)=c('Trait.factor','Gene','gene_chr','gene_start','gene_end','SNP.factor','class.factor','time.factor')
#env2=as.data.table(env2)
#setkey(env2,gene_chr,gene_start,gene_end)
#comparison6=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('gene_chr','gene_start','gene_end'),nomatch=NULL)



# Overlap of factor eQTL variants with QTL support intervals
qtl10=fread('../GridLMM/Biogemma_10p_QTL.csv',data.table=F)
env1=factoreqtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(qtl10)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

# Factor 2 trans-eQTL on chr 8 overlaps with qTPH8 in Graneros 2015 OPt




#overlap of gene location +1 10kb? Look at papers. what do they use?
#ciseqtl=ciseqtl[ciseqtl$time=="WD_0718",]

ciseqtl_genetable=genetable[genetable$Gene_ID %in% ciseqtl$Gene,]
#ciseqtl_genetable=ciseqtl_genetable[grepl('T001',ciseqtl_genetable$TXNAME),]
#i=47920
#ciseqtl_genetable=rbind(ciseqtl_genetable,genetable[i,])
#ciseqtl_genetable$TXCHROM=as.integer(ciseqtl_genetable$TXCHROM)
rownames(ciseqtl_genetable)=seq(1,nrow(ciseqtl_genetable))
env1=ciseqtl_genetable
env1=as.data.table(env1)
env2=as.data.table(qtl)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison2=foverlaps(env1,env2,by.x=c('CHROM','window_START','window_END'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

#overlap of founder recombination blocks - does this make sense?
ciseqtl$window_START=all_founder_blocks[match(ciseqtl$SNP,all_founder_blocks$focal_snp),]$start
ciseqtl$window_END=all_founder_blocks[match(ciseqtl$SNP,all_founder_blocks$focal_snp),]$end


env1=ciseqtl
#env1$BP_end=env1$BP
env1=as.data.table(env1)
env2=as.data.table(qtl)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','window_START','window_END'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

#overlap of gene location +1 10kb? Look at papers. what do they use?
ciseqtl=ciseqtl[ciseqtl$time=="WD_0712",]

ciseqtl_genetable=genetable[genetable$Gene_ID %in% ciseqtl$Gene,]
#ciseqtl_genetable=ciseqtl_genetable[grepl('T001',ciseqtl_genetable$TXNAME),]
#i=47920
#ciseqtl_genetable=rbind(ciseqtl_genetable,genetable[i,])
#ciseqtl_genetable$TXCHROM=as.integer(ciseqtl_genetable$TXCHROM)
rownames(ciseqtl_genetable)=seq(1,nrow(ciseqtl_genetable))
env1=ciseqtl_genetable
env1=as.data.table(env1)
env2=as.data.table(qtl)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison2=foverlaps(env1,env2,by.x=c('CHROM','window_START','window_END'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)



## Correlation of effect sizes
# Grab founder effect sizes from QTL_F -
# Grab founder effect sizes from cis-eQTL test from GridLMM


# Change betas based on time point of eQTL
#eqtl_betas=c()
#for(chr in 1:10){
#  eqtl_beta=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_vst_results.txt',time,chr),data.table=F)
#  eqtl_betas=rbind(eqtl_betas,eqtl_beta)
#}
#eqtl_betas$snp_gene=paste0(eqtl_betas$X_ID,'_',eqtl_betas$Trait)
pvalue=c()
es_cor=c()

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder

comp_list=list()
count=1
for(q in 1:nrow(comparison)){
  row=comparison[q,]
  id=row$pheno_env_id
  chr=row$CHR
  pheno=row$Phenotype
  env=row$Environment
  gene=row$Gene
  #time=row$time
  time=row$time
  eqtl_betas=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_fkeep_results.txt',time,chr),data.table=F)
  eqtl_betas$snp_gene=paste0(eqtl_betas$X_ID,'_',eqtl_betas$Trait)
#  effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_vst_founderprobs.rds',chr,pheno,env))

  effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',chr,pheno,env))
  effect_sizes=unlist(unname(effect_sizes[effect_sizes$X_ID==row$highest_SNP,6:21]))
  effect_sizes[-1]=effect_sizes[1] + effect_sizes[-1]
  #test=which(unlist(unname(lapply(effect_sizes,function(x) x$i==id))))
  #es=effect_sizes[[test]]$values
  tmp=eqtl_betas[eqtl_betas$X_ID==row$SNP & eqtl_betas$Trait==row$Gene,c(6,10:24)]
  eqs=unlist(unname(tmp))
  intercept=min(which(!is.na(eqs)))
  df=data.frame(founder=founders,exp_effects=eqs,pheno_effects=effect_sizes,stringsAsFactors=F)
  df=df[order(df$exp_effects),]
  df$founder_f=factor(df$founder,levels=c(df$founder))
  #eqs[-intercept]=eqs[intercept]+eqs[-intercept]
  t=cor.test(eqs,effect_sizes)
  pvalue=c(pvalue,t$p.value)
  es_cor=c(es_cor,t$estimate)
  p1=ggplot(df,aes(x=exp_effects,y=pheno_effects)) + geom_point(aes(color=founder_f)) + xlab("Expression effect sizes") +
  scale_color_manual(values=colorcodes[levels(df$founder_f),]$hex_color,labels=levels(df$founder_f))+
  ylab("Phenotype effect sizes") + ggtitle(sprintf('r=%.2f Founder %s and %s effects',t$estimate,gene,pheno))
  comp_list[[count]]=p1
  count=count+1

  drop=c("VA85","A632_usa")
  dloc=c(2,16)
  subf=founders[-dloc]
  subdf= df[df$founder %in% subf,]
  newt=cor.test(subdf$exp_effects,subdf$pheno_effects)
  p1=ggplot(subdf,aes(x=exp_effects,y=pheno_effects)) + geom_point(aes(color=founder_f)) + xlab("Expression effect sizes") +
  scale_color_manual(values=colorcodes[levels(df$founder_f),]$hex_color,labels=levels(df$founder_f))+
  ylab("Phenotype effect sizes") + ggtitle(sprintf('r=%.2f Founder %s and %s effects, dropped',newt$estimate,gene,pheno))
  comp_list[[count]]=p1
  count=count+1
}

drop=c("VA85","A632_usa")
dloc=c(2,16)
subf=founders[-dloc]
subdf= df[df$founder %in% subf,]


pdf(sprintf('images/%s_eQTL_QTL_effect_sizes.pdf',time))
for(i in 1:length(comp_list)){
  print(comp_list[[i]])
}
dev.off()


comparion=as.data.frame(comparison,stringsAsFactors=F)
comparison$es_cor=es_cor
comparison$cortest_pvalue=pvalue

fwrite(comparison,'eqtl/results/all_eQTL_QTL_overlap.txt',row.names=F,quote=F,sep='\t')


comparison=fread('eqtl/results/all_eQTL_QTL_overlap.txt',data.table=F)
#hasf=c('S_F_H','F_and_H','F_only','S_and_F')
# Null expectation
infosnps=unique(comparison$highest_SNP) # SNPs for QTL
eqtlsnps=unique(comparison$SNP) # SNPs for eQTL



eqtl_betas=c()
for(chr in 1:10){
  eqtl_beta=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_fkeep_results.txt',time,chr),data.table=F)
  eqtl_betas=rbind(eqtl_betas,eqtl_beta)
}
eqtl_betas$snp_gene=paste0(eqtl_betas$X_ID,'_',eqtl_betas$Trait)

allcorsdf=c()
testeqtls=ciseqtl[!(ciseqtl$SNP %in% eqtlsnps),]
for(q in 1:length(infosnps)){
  #row=comparison[q,]
  infosnp=infosnps[q]
  row=comparison[comparison$highest_SNP==infosnp,]

  #id=row$pheno_env_id
  dontcomp=unique(comparison[comparison$highest_SNP==infosnp,]$SNP)
  testesnps=testeqtls[!(testeqtls$SNP %in% dontcomp),]
  testesnps$snp_gene=paste0(testesnps$SNP,'_',testesnps$Gene)
  tmp1=eqtl_betas[match(testesnps$snp_gene,eqtl_betas$snp_gene),]
  tmp1=tmp1[which(unlist(unname(sapply(seq(1,nrow(tmp1)),function(x) sum(is.na(tmp1[x,]))<16)))),]
  tmp=tmp1[,c(6,10:24)]
  rownames(tmp)=seq(1,nrow(tmp))
  #for(t in 1:nrow(tmp)){
  #  crow=tmp[t,]
    #intercept=min(which(!is.na(crow)))
  #  eqs=unlist(unname(crow))
    #eqs[-intercept]=eqs[intercept]+eqs[-intercept]
  #  tmp[t,]=eqs
  #}

  for(r in 1:nrow(row)){
    chr=row$CHR[r]
    pheno=row$Phenotype[r]
    env=row$Environment[r]
    effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_vst_founderprobs.rds',chr,pheno,env))
    effect_sizes=unlist(unname(effect_sizes[effect_sizes$X_ID==infosnp,6:21]))
    effect_sizes[-1]=effect_sizes[1] + effect_sizes[-1]

    allcors=apply(tmp,MARGIN=1,function(x) cor.test(unlist(unname(x)),effect_sizes)$estimate)
    allps=apply(tmp,MARGIN=1,function(x) cor.test(unlist(unname(x)),effect_sizes)$p.value)

    nullp=quantile(allps,0.05)
    obvp=comparison[comparison$SNP %in% dontcomp & comparison$highest_SNP==infosnp,]$cortest_pvalue
    if(sum(obvp < nullp)>0){
      print(infosnp)
      print(env)
    }
    nullcor=quantile(abs(allcors),0.95)
    obvcor=comparison[comparison$SNP %in% dontcomp & comparison$highest_SNP==infosnp,]$es_cor
    if(sum(abs(obvcor) > nullcor)>0){
      print(infosnp)
      print(env)
    }
    png(sprintf('images/null_correlation_dist_%s_%s_x_%s.png',infosnp,pheno,env))
    hist(abs(allcors))
    dev.off()
    x=which.max(abs(allcors))
    allcorsdf=rbind(allcorsdf,c(infosnp,pheno,env,unlist(unname(tmp1[x,])),allcors[x],allps[x]))
  }
}

allcorsdf=as.data.frame(allcorsdf,stringsAsFactors=F)
names(allcorsdf)=c('qtl_snp','pheno','env',names(tmp1),'es_cor','cor_pvalue')
allcorsdf$cor_pvalue=as.numeric(allcorsdf$cor_pvalue)
  #test=which(unlist(unname(lapply(effect_sizes,function(x) x$fsnp==infosnp))))
  #for now, just use effect sizes from the first QTL. Maybe check to see if
  # the effect sizes across the environemnt_phenotypes is highly correlated
  #if(length(test)>1){
  #  test=test[1]
  #}
  #test=which(unlist(unname(lapply(effect_sizes,function(x) !(x$fsnp %in% dontcomp) & x$label %in% hasf))))

  #es=effect_sizes[[test]]$values



#None of the cis-eQTL that overlap with QTL have correlated effect sizes that are
# more correlated with QTL effect sizes than an other eQTL that don't colocalize

# I should double check my work here though...
# What are these very high correlations? Could these be due to LD? cis-eQTL effecting the trait in trans?
# Need to look at factors from MegaLMM
# Based off of allele frequency as well? Or does that not matter?


#ZmRap2.7 effect size
pmap=fread('../genotypes/qtl2/startfiles/Biogemma_pmap_c8.csv',data.table=F)

rap27="Zm00001d010987"
#genetable[genetable$Gene_ID==rap27,]
#             Gene_ID CHROM     START       END
#32330 Zm00001d010987     8 136009216 136012084
chr=8
env="ALL"
pheno="male_flowering_d6"
qdta8_snp="AX-91102763"
#pmap[pmap$marker==qdta8_snp,]
#           marker chr       pos
#28709 AX-91102763   8 135165296
effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_vst_founderprobs.rds',chr,pheno,env))
betas=unlist(unname(effect_sizes[effect_sizes$X_ID==qdta8_snp,6:21]))
betas[-1]=betas[1]+betas[-1]
names(betas)=founders

eqtl_snp="AX-91102912"
#marker chr       pos
#28918 AX-91102912   8 135736942
eqtl_betas=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_results2.txt',time,chr),data.table=F)
eqtl_betas$snp_gene=paste0(eqtl_betas$X_ID,'_',eqtl_betas$Trait)
eqtl_rap27=eqtl_betas[eqtl_betas$Trait==rap27,]

eqtl_beta=unlist(unname(eqtl_betas[eqtl_betas$Trait==rap27,6:21]))
#eqtl_beta[-1]=eqtl_beta[1]+eqtl_beta[-1]
names(eqtl_beta)=founders
cor(eqtl_beta,betas)


# How correlated with phenotype effect sizes are genes within the support interval?
times=c("WD_0712","WD_0718",'WD_0720','WD_0727')

qtl_genes=fread('metadata/QTL_support_interval_genes.txt',data.table=F)
qtl_genes$pheno_env=paste0(qtl_genes$Phenotype,'-',qtl_genes$Environment)
pheno_envs=unique(qtl_genes$pheno_env)
for(pe in pheno_envs){
  df=qtl_genes[qtl_genes$pheno_env==pe,]
  chroms=unique(df$CHROM)
  pheno=strsplit(pe,'-')[[1]][1]
  env=strsplit(pe,'-')[[1]][2]
  for(chr in chroms){
    effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_vst_founderprobs.rds',chr,pheno,env))
    tmp=df[df$CHROM==chr,]
    qs=unique(tmp$ID)
    for(q in qs){
      subdf=tmp[tmp$ID==q,]
      infosnp=unique(subdf[subdf$ID==q,]$highest_SNP)
      effect_size=unlist(unname(effect_sizes[effect_sizes$X_ID==infosnp,6:21]))
      effect_size[-1]=effect_size[1] + effect_size[-1]

      for(time in times){
        eqtl_betas=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_fkeep_results.txt',time,chr),data.table=F)
        eqtl_betas$snp_gene=paste0(eqtl_betas$X_ID,'_',eqtl_betas$Trait)
        #eqtl_rap27=eqtl_betas[eqtl_betas$Trait==rap27,]

        eqtl_beta=eqtl_betas[eqtl_betas$Trait %in% subdf$Gene_ID,c(1,6,10:24)]
        rownames(eqtl_beta)=eqtl_beta$Trait
        eqtl_beta=eqtl_beta[,-1]
        #tmp1=eqtl_betas[match(testesnps$snp_gene,eqtl_betas$snp_gene),]
        #tmp1=tmp1[which(unlist(unname(sapply(seq(1,nrow(tmp1)),function(x) sum(is.na(tmp1[x,]))<16)))),]
        #tmp=tmp1[,c(6,10:24)]
        #rownames(tmp)=seq(1,nrow(tmp))
        results=apply(eqtl_beta,MARGIN=1,function(x) cor.test(effect_size,x))
        pvalues=unlist(lapply(results,function(x) x$p.value))
        cors=unlist(lapply(results,function(x) x$estimate))
        correlations=data.frame(Gene_ID=rownames(eqtl_beta),pheno_env=pe,cor=cors,pvalue=pvalues,stringsAsFactors=F)
        coln=paste0(time,'-correlation')
        coln2=paste0(time,'-pvalue')
        names(correlations)=c('Gene_ID','pheno_env',coln,coln2)
        #qtl_genes[(match(correlations$Gene,qtl_genes$Gene_ID) & qtl_genes$pheno_env==pe),coln]=correlations[,coln]

        # How do I update qtl_genes for correlations and fill in as I go?
        for(x in 1:nrow(eqtl_beta)){
          qtl_genes[(qtl_genes$pheno_env==pe & qtl_genes$Gene_ID==rownames(eqtl_beta)[x]),coln]=correlations[correlations$Gene_ID==rownames(eqtl_beta)[x],coln]
          qtl_genes[(qtl_genes$pheno_env==pe & qtl_genes$Gene_ID==rownames(eqtl_beta)[x]),coln2]=correlations[correlations$Gene_ID==rownames(eqtl_beta)[x],coln2]

        }
      }
    }
  }
}
# Include cor.test - correct for multiple testing
# This is looking at additive effects only - should do with full expression as well
#(after correcting for PCs?)
ntests=sum(!is.na(qtl_genes[,c(26,28,30,32)]))
#23894
threshold=0.05/ntests

sig=qtl_genes[sum(qtl_genes[,c(27,29,31,33)]<=threshold) > 0,]
which.main
fwrite(qtl_genes,'eqtl/QTL_GeneExp_correlations.txt',row.names=F,quote=F,sep='\t')
