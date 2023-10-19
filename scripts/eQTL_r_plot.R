#!/usr/bin/env Rscript

library('ggplot2')
library('data.table')
library('dplyr')

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}
eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end

qtl=fread('QTL/all_adjusted_QTL_SIs.txt',data.table=F)

#qtl=fread('QTL/all_adjusted_QTL_peaks_trimmed.txt',data.table=F)
qtl$block_start=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$start
qtl$block_end=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$end
env1=eqtl
env1=as.data.table(env1)
env2=as.data.table(qtl)
setkey(env2,CHR,left_bound_bp,alt_right_bound_bp)
comparison2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)
fwrite(comparison2,'QTT/QTL_cis_eQTL_overlap.txt',row.names=F,quote=F,sep='\t')

comparison2=fread('QTT/QTL_cis_eQTL_overlap.txt',data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
plot_list=list()
count=1
#adj_chr=c(5,9)
cors=c()
pvals=c()
for(i in 1:nrow(comparison2)){
	row=comparison2[i,]
	pheno=row$phenotype
	env=row$environment
	chr=row$CHR
	gene=row$Trait
	snp=row$SNP
	id=row$ID
	time=row$time
	#exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	gene=row$Trait
	snp=row$SNP
	esnp=row$X_ID
	chr=row$CHR
	#if(chr %in% adj_chr){
	#	X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
	#}else{
	#	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	#
	#}
	#inds=rownames(X_list[[1]])
	#inter=intersect(exp$V1,inds)
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==snp,]
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
	
	#df=data.frame(founder=founders,pheno=effect_size,eqts=betas,stringsAsFactors=F)
	#if(abs(r)>0.6){
	#	p1=ggplot(df,aes(x=eqts,y=pheno)) + geom_point(aes(color=founder)) +
	#	xlab("eQTL effect size (log2cpm)") + ylab("QTL effect size") +
	#	ggtitle(sprintf("%s %s and %s %s, r=%.2f",gene,time,id,env,r))
	#	plot_list[[count]]=p1
	#	count=count+1
	#}
	
}

comparison2$r=cors
comparison2$pvalue=pvals
fwrite(comparison2,'QTT/QTL_cis_eQTL_interval_overlap.txt',row.names=F,quote=F,sep='\t')


pdf('QTL/images/QTL_eQTL_effect_size_correlations.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()



# For each QTL, plot r by log10p, and r by var explained
prop_var=c()
for(i in 1:nrow(comparison2)){
	row=comparison2[i,]
	time=row$time
	chr=row$CHR
	pv=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_prop_var_FIXED.txt',time,chr),data.table=F)
	gene=row$Trait
	esnp=row$X_ID
	p=pv[pv$snp==esnp & pv$gene==gene,]$prop_var
	prop_var=c(prop_var,p)
}
comparison2$prop_var=prop_var
comparison2$qtl_prop_var=qtl[match(comparison2$pheno_env_ID,qtl$pheno_env_ID),]$prop_var

fwrite(comparison2,'QTT/QTL_cis_eQTL_interval_overlap.txt',row.names=F,quote=F,sep='\t')


# How does correlation change over time?

# For each QTL, grab the top 50 most significant eQTL - what is the max abs(r) for those
comparison2=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
comparison2$pheno_env_ID=paste0(comparison2$phenotype,'-',comparison2$environment,'-',comparison2$ID)
genecount=comparison2 %>% group_by(pheno_env_ID) %>% count()
pheno_env_ids=unique(comparison2$pheno_env_ID)

get_cor=function(row1,d2){
	row2=qtl[d2,]
	pheno=row2$phenotype
	env=row2$environment
	chr1=row1$CHR
	chr2=row2$CHR
	gene=row1$Trait
	esnp=row1$X_ID
	qsnp=row2$SNP
	id=row2$ID
	time=row1$time
	gene=row1$Trait
	
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr2,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr1),data.table=F)
	results=results[results$X_ID==esnp & results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	return(test$estimate)
}


random_ten=function(pei){
	d2=which(qtl$pheno_env_ID==pei)
	ndraws=genecount[genecount$pheno_env_ID==pei,]$n
	dropgenes=unique(comparison2[comparison2$pheno_env_ID==pei,]$gene_time_SNP)
	subcis=eqtl[!(eqtl$gene_time_snp %in% dropgenes),]
	drawl=sample(seq(1,nrow(subcis)),ln)
	subeqtl=subcis[drawl,]
	n=min(10,ndraws)
	top10=subeqtl %>% slice_max(value, n = n)
	top10=as.data.frame(top10,stringsAsFactors=F)
	all_r=sapply(seq(1,n),function(x) get_cor(top10[x,],d2))
	highest_r=max(abs(all_r))
	return(highest_r)
}

qtl$pheno_env_ID=paste0(qtl$phenotype,'-',qtl$environment,'-',qtl$ID)

all_perms=c()
for(pei in pheno_env_ids){
	bootstrap=sapply(seq(1,300),function(x) random_ten(pei))
	perms=data.frame(pei=pei,rep=seq(1,300),max_r=bootstrap,stringsAsFactors=F)
	all_perms=rbind(all_perms,perms)
}
names(all_perms)=c('pei','rep','perm_max_r')
fwrite(all_perms,'QTT/top10_eqtl_permutations.txt',row.names=F,quote=F,sep='\t')

all_perms=fread('QTT/top10_eqtl_permutations.txt',data.table=F)
n=10
top_r=c()
for(pei in pheno_env_ids){
	subdf=comparison2[comparison2$pheno_env_ID==pei,]
	subdf=subdf[order(subdf$value),]
	rownames(subdf)=seq(1,nrow(subdf))
	subdf$rank=seq(nrow(subdf),1)
	max_loc=which.max(abs(subdf$r))
	max_r=max(abs(subdf$r))
	top10=subdf %>% slice_max(value, n = n)
	top10=as.data.frame(top10)
	#top10_r=top10_r[order(top10_r$i.value),]
	#top10_r$rank=seq(1,10)
	top10_r_loc=unlist(which.max(abs(top10$r)))
	top_rank=top10[top10_r_loc,]$rank
	top10_r=abs(top10[top10_r_loc,]$r)
	newline=data.frame(pheno_env_ID=pei,max_r=max_r,top10_r=top10_r,max_rank=top_rank,stringsAsFactors=F)
	top_r=rbind(top_r,newline)
}

all_perms$max_r=top_r[match(all_perms$pei,top_r$pheno_env_ID),]$max_r




all_perms=all_perms[order(all_perms$max_r),]
rownames(all_perms)=seq(1,nrow(all_perms))
all_perms$pei_f=factor(all_perms$pei,levels=c(unique(all_perms$pei)))

top_r=top_r[order(top_r$max_r),]
rownames(top_r)=seq(1,nrow(top_r))
top_r$pei_f=factor(top_r$pheno_env_ID,levels=c(unique(top_r)$pheno_env_ID))

max_r_f=sort(unique(top_r$max_r2))
top_r$max_r2=top_r$max_r**2
top_r$max_r_f=factor(top_r$max_r2,levels=c(max_r_f))

all_perms$max_r2=all_perms$max_r**2
all_perms$max_r_f=factor(all_perms$max_r2,levels=c(max_r_f))


all_perms$max_r100=round(all_perms$max_r2*100)

top_r$max_r100=round(top_r$max_r2*100)

p2=ggplot(aes(x=max_r,y=top10_r),data=top_r) + geom_point() + xlab("Highest correlation eQTL (|r|)") + 
ylab("Highest correlation gene of 10 most significant eQTL (|r|)") + geom_abline(slope=1) +
theme_classic()

png('paper_figures/local_eQTL_r_by_sig.png')
print(p2)
dev.off()


p2=ggplot(all_perms,aes(x=max_r100,y=perm_max_r**2,group=max_r100))  + 
scale_x_continuous(limits=c(0,100),breaks=c(seq(0,100,10)),labels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9",'1.0')) +
scale_y_continuous(limits=c(0,1),breaks=c(seq(0,1,0.1)),labels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9",'1.0')) +
geom_boxplot(position =position_dodge(0),alpha=0.2) +
geom_point(data=top_r,aes(x=max_r100,y=top10_r**2),color="red",size=2) +
xlab("Highest correlation eQTL (r^2)") + 
ylab("Highest correlation of 10 most significant eQTL (r^2)") + geom_abline(slope=0.01) +
ggtitle("Correlation of QTL effect sizes with overlapping eQTL") +
theme_classic()

png('paper_figures/local_eQTL_r_by_sig.png',width=800,height=800)
print(p2)
dev.off()


top_r$cand=top_r$pheno_env_ID %in% cand$pei
# highlight candidate genes
p2=ggplot(all_perms,aes(x=max_r100,y=perm_max_r**2,group=max_r100))  + 
scale_x_continuous(limits=c(30,100),breaks=c(seq(30,100,10)),labels=c("0.3","0.4","0.5","0.6","0.7","0.8","0.9",'1.0')) +
scale_y_continuous(limits=c(0,1),breaks=c(seq(0,1,0.1)),labels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9",'1.0')) +
geom_boxplot(position =position_dodge(0),alpha=0.2) +
geom_point(data=top_r,aes(x=max_r100,y=top10_r**2,color=cand),size=2) +
scale_color_manual(values = c("TRUE" = "red","FALSE"="black")) +
xlab(bquote('Highest correlation local-eQTL'(r^2))) + 
ylab(bquote('Highest correlation of 10 most significant local-eQTL' (r^2))) + geom_abline(slope=0.01) +
ggtitle("Correlation of QTL effect sizes with overlapping eQTL") +
theme_classic() + guides(color="none")

png('paper_figures/local_eQTL_r_by_sig.png',width=800,height=800)
print(p2)
dev.off()

# is this the plot I want?
# 

# generally, how many are higher or lower than expected by chance
all_perms$top10_r=top_r[match(all_perms$pei,top_r$pheno_env_ID),]$top10_r

#lowerbound=all_perms %>% group_by(pei) %>% summarize(cutoff=quantile(perm_max_r,0.025))
upperbound=all_perms %>% group_by(pei) %>% summarize(cutoff=quantile(perm_max_r,0.95))

sig=c()
#lowers=c()
highers=c()
for(i in 1:nrow(top_r)){
	row=top_r[i,]
	pei=row$pheno_env_ID
	top10_r=row$top10_r
	#low=lowerbound[lowerbound$pei==pei,]$cutoff
	high=upperbound[upperbound$pei==pei,]$cutoff
	#test=(top10_r < low | top10_r > high)
	test=top10_r > high
	sig=c(sig,test)
	#lowers=c(lowers,low)
	highers=c(highers,high)
}
top_r$top10_sig=sig
#top_r$perm_lower=lowers
top_r$perm_higher=highers
# How does this compare to random?
top_r[top_r$top10_sig==TRUE,]
#                                    pheno_env_ID     max_r   top10_r max_rank
#19   male_flowering_d6-GRANEROS_2015_OPT-qDTA3_2 0.7700711 0.7029747        2
#23 female_flowering_d6-GRANEROS_2015_OPT-qDTS3_2 0.7785476 0.7237072        2
#30      female_flowering_d6-BLOIS_2014_OPT-qDTS8 0.7962589 0.7962589       10

#                                           pei_f    max_r2 max_r_f max_r100
#19   male_flowering_d6-GRANEROS_2015_OPT-qDTA3_2 0.5930095    <NA>       59
#23 female_flowering_d6-GRANEROS_2015_OPT-qDTS3_2 0.6061364    <NA>       61
#30      female_flowering_d6-BLOIS_2014_OPT-qDTS8 0.6340283    <NA>       63

#   top10_sig perm_lower perm_higher
#19      TRUE  0.2837239   0.6561011
#23      TRUE  0.2815416   0.6528896
#30      TRUE  0.3379347   0.7076527


# these three QTL, the highest r 
pei="female_flowering_d6-GRANEROS_2015_OPT-qDTS3_2"
subdf=comparison2[comparison2$pheno_env_ID==pei,]
top10=subdf %>% slice_max(abs(r), n = n)

unique(top10$Trait)
#[1]"Zm00001d042307" "Zm00001d042291" "Zm00001d000636" "Zm00001d042263"
#[5] "Zm00001d042315" "Zm00001d042306" "Zm00001d042146" "Zm00001d042192"
#[9] "Zm00001d042166"


pei="female_flowering_d6-BLOIS_2014_OPT-qDTS8"
subdf=comparison2[comparison2$pheno_env_ID==pei,]
top10=subdf %>% slice_max(abs(r), n = n)

unique(top10$Trait)
#[1] "Zm00001d010946" "Zm00001d010975" "Zm00001d010804" "Zm00001d010974"
#[5] "Zm00001d010892" "Zm00001d010982" "Zm00001d010758"

pei="male_flowering_d6-GRANEROS_2015_OPT-qDTA3_2"
subdf=comparison2[comparison2$pheno_env_ID==pei,]
top10=subdf %>% slice_max(abs(r), n = n)

unique(top10$Trait)
#[1] "Zm00001d042307" "Zm00001d042291" "Zm00001d042315" "Zm00001d000636"
#[5] "Zm00001d042306" "Zm00001d042214" "Zm00001d042263" "Zm00001d042228"

# Grab the max r for all genes within the interval
# Plot out the points

plot_list=list()
count=1
for(pei in pheno_env_ids){
	subdf=comparison2[comparison2$pheno_env_ID==pei,]
	qtlp=unique(subdf$value)
	p1=ggplot(aes(x=abs(r),y=i.value),data=subdf) + geom_point() + xlab("eQTL correlation with QTL effect size (|r|)") + 
	ylab("log10(q-value) of eQTL") + ggtitle(sprintf('%s overlapping eQTL (log10(pvalue)= %.2f)',pei,qtlp))
	plot_list[[count]]=p1
	count=count+1
}

pdf('QTT/images/eQTL_r_by_qvalue.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()


# What is the relationship between QTL prop_var (effect size) and max eQTL r?
max_r= comparison2 %>% group_by(pheno_env_ID) %>% slice(which.max(abs(r)))
ft=c("male_flowering_d6","female_flowering_d6")
max_r$ft=max_r$phenotype %in% ft
m1=lm(abs(r) ~ qtl_prop_var,max_r)


m2=lm(abs(r) ~ as.factor(ft)+qtl_prop_var,max_r)
p3=ggplot(aes(x=abs(r),y=qtl_prop_var),data=max_r) + geom_point(aes(color=ft)) + xlab('Max eQTL correlation |r|') +
ylab("Proportion of phenotypic variance explained by QTL")
png('QTT/images/eQTL_r_by_qtl_propvar.png')
print(p3)
dev.off()

p3=ggplot(aes(x=abs(r),y=qtl_prop_var),data=max_r) + geom_point(aes(color=ID)) + xlab('Max eQTL correlation |r|') +
ylab("Proportion of phenotypic variance explained by QTL")
png('QTT/images/eQTL_r_by_ID_propvar.png')
print(p3)
dev.off()

# Make file with the samples used in each timepoint
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	norm=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
	inter=intersect(norm$V1,names(K)[-1])
	inter=data.frame(id=inter,stringsAsFactors=F)
	fwrite(inter,sprintf('eqtl/data/%s_samples_FIXED.txt',time))
}

# top eQTL expression by founder
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

plot_list=list()
n=10
top10=eqtl %>% slice_max(value, n = n)
count=1
adj_chr=c(5,9)
for(i in 1:nrow(top10)){
	gene=top10[i,]$Trait
	time=top10[i,]$time
	samplefile=fread(sprintf('eqtl/data/%s_samples_FIXED.txt',time),data.table=F)
	inter=samplefile$id
	snp=top10[i,]$X_ID
	chr=top10[i,]$CHR
	if(chr %in% adj_chr){
		X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
	}else{
		X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
	}
	X_list=lapply(X_list,function(x) x[inter,])
	ex=data.frame(ID=norm$V1,exp=unlist(norm[,gene]),stringsAsFactors=F)
	ex=ex[ex$ID %in% inter,]
	ex=ex[order(ex$exp),]
	rownames(ex)=seq(1,nrow(ex))
	ex$ID_f=factor(ex$ID,levels=c(unique(ex$ID)))
	X = do.call(cbind,lapply(X_list,function(x) x[,snp]))
    colnames(X) = founders
    rownames(X) = dimnames(X_list[[1]])[[1]]
    X=X[inter,]
    founder_id=apply(X,MARGIN=1,function(x) names(x[which.max(x)]))
    ex$founder=founder_id[match(ex$ID,names(founder_id))]
    p=ggplot(aes(x=ID_f,y=exp),data=ex) + 
    geom_point(aes(color=founder)) +
     xlab('Sample') +
     ylab('Expression (log2CPM)') + geom_hline(yintercept=1) +
     ggtitle(sprintf('%s %s Expression by Sample',gene,time))
    plot_list[[count]]=p
    count=count+1
}

pdf('eqtl/cis/images/founder_eqtl_by_ind_FIXED.pdf')
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()


#### For eQTL within QTL, how does correlation plot with location
comparison2=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
comparison2$pheno_env_ID=paste0(comparison2$phenotype,'-',comparison2$environment,'-',comparison2$ID)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
compmerge=merge(comparison2,genetable,by.x='Trait',by.y="Gene_ID")
compmerge$midgene=compmerge$START + ((compmerge$END-compmerge$START)/2)
pheno_env_ids=unique(compmerge$pheno_env_ID)

plot_list=list()
count=1
for(pei in pheno_env_ids){
	subdf=compmerge[compmerge$pheno_env_ID==pei,]
	#maxr=max(ab)
	#times=unique(subdf$time)
	p1=ggplot(subdf, aes(xmin = START/1e6, xmax = END/1e6, ymin = r-0.02, ymax = r+0.02)) + 
  	geom_rect(stat = "identity") + facet_grid(rows = vars(time)) +
#	p1=ggplot(subdf,aes(x=midgene/1e6,y=r)) + geom_bar(aes(x=STARt/1e6,)) + xlab("Physical position (bp)") +
	ylab("Correlation of eQTL") + ggtitle(sprintf("%s eQTL correlation",pei)) +
	ylim(-1,1) + geom_hline(yintercept=0,color="black")
	plot_list[[count]]=p1
	count=count+1

}

pdf('QTT/images/eQTL_cor_by_pos.pdf')
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()



########### trans-eQLT STPUAL corrleatio

qtl=fread('QTL/all_adjusted_QTL_all_methods.txt',data.table=F)
qtl$pheno_env_ID=paste0(qtl$phenotype,'-',qtl$environment,'-',qtl$ID)
qtl=qtl[qtl$method=="Founder_probs",]
subcomp=fread('QTT/QTL_trans_eQTL_EXPSTPAUL_overlap.txt',data.table=F)
trans=fread('eqtl/results/all_trans_fdr_SIs_FIXED.txt',data.table=F)

prop_var=c()
for(i in 1:nrow(subcomp)){
	row=subcomp[i,]
	time=row$time
	chr=row$CHR
	pv=fread(sprintf('eqtl/trans/results/eQTL_%s_c%s_weights_prop_var_FIXED.txt',time,chr),data.table=F)
	gene=row$gene
	esnp=row$SNP
	p=pv[pv$snp==esnp & pv$gene==gene,]$prop_var
	prop_var=c(prop_var,p)
}
subcomp$prop_var=prop_var
subcomp$qtl_prop_var=qtl[match(subcomp$pheno_env_ID,qtl$pheno_env_ID),]$prop_var


genecount=subcomp %>% group_by(pheno_env_ID) %>% count()
pheno_env_ids=unique(subcomp$pheno_env_ID)

get_cor=function(row1,d2){
	row2=qtl[d2,]
	pheno=row2$phenotype
	env=row2$environment
	chr1=row1$CHR
	chr2=row2$CHR
	gene=row1$gene
	esnp=row1$SNP
	qsnp=row2$SNP
	id=row2$ID
	time=row1$time
	
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr2,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	
	results=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_filtered_FIXED.txt',time,chr1),data.table=F)
	results=results[results$X_ID==esnp & results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	return(test$estimate)
}


random_ten=function(pei){
	d2=which(qtl$pheno_env_ID==pei)
	ndraws=genecount[genecount$pheno_env_ID==pei,]$n
	draw=sample(seq(1,nrow(trans)),ndraws)
	subeqtl=trans[draw,]
	n=min(10,ndraws)
	top10=subeqtl %>% slice_max(value, n = n)
	top10=as.data.frame(top10,stringsAsFactors=F)
	all_r=sapply(seq(1,n),function(x) get_cor(top10[x,],d2))
	highest_r=max(abs(all_r))
	return(highest_r)
}

qtl$pheno_env_ID=paste0(qtl$phenotype,'-',qtl$environment,'-',qtl$ID)

all_perms=c()
for(pei in pheno_env_ids){
	bootstrap=sapply(seq(1,300),function(x) random_ten(pei))
	perms=data.frame(pei=pei,rep=seq(1,300),max_r=bootstrap,stringsAsFactors=F)
	all_perms=rbind(all_perms,perms)
}
fwrite(all_perms,'QTT/top10_eqtl_permutations.txt',row.names=F,quote=F,sep='\t')
names(all_perms)=c('pei','rep','perm_max_r')

all_perms=fread('QTT/top10_eqtl_permutations.txt',data.table=F)
n=10
top_r=c()


pheno_env_ids=unique(subcomp$pheno_env_ID)

for(pei in pheno_env_ids){
	subdf=subcomp[subcomp$pheno_env_ID==pei,]
	subdf=subdf[order(subdf$value),]
	rownames(subdf)=seq(1,nrow(subdf))
	subdf$rank=seq(nrow(subdf),1)
	max_loc=which.max(abs(subdf$r))
	max_r=max(abs(subdf$r))
	top10=subdf %>% slice_max(value, n = n)
	top10=as.data.frame(top10)
	#top10_r=top10_r[order(top10_r$i.value),]
	#top10_r$rank=seq(1,10)
	top10_r_loc=unlist(which.max(abs(top10$r)))
	top_rank=top10[top10_r_loc,]$rank
	top10_r=abs(top10[top10_r_loc,]$r)
	newline=data.frame(pheno_env_ID=pei,max_r=max_r,top10_r=top10_r,max_rank=top_rank,stringsAsFactors=F)
	top_r=rbind(top_r,newline)
}

all_perms$max_r=top_r[match(all_perms$pei,top_r$pheno_env_ID),]$max_r



p2=ggplot(aes(x=max_r,y=top10_r),data=top_r)  + xlab("Highest correlation eQTL (|r|)") + 
ylab("Highest correlation gene of 10 most significant eQTL (|r|)") + geom_abline(slope=1)

png('QTT/images/eQTL_r_by_sig.png')
print(p2)
dev.off()


plot_list=list()
count=1
for(pei in pheno_env_ids){
	subdf=subcomp[subcomp$pheno_env_ID==pei,]
	qtlp=unique(subdf$value)
	p1=ggplot(aes(x=abs(r),y=value),data=subdf) + geom_point() + xlab("eQTL correlation with QTL effect size (|r|)") + 
	ylab("log10(q-value) of eQTL") + ggtitle(sprintf('%s overlapping eQTL (log10(pvalue)= %.2f)',pei,qtlp))
	plot_list[[count]]=p1
	count=count+1
}

pdf('QTT/images/trans_eQTL_r_by_qvalue.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()

max_r= subcomp %>% group_by(pheno_env_ID) %>% slice(which.max(abs(r)))
max_r=as.data.frame(max_r)

p3=ggplot(aes(x=abs(r),y=qtl_prop_var),data=max_r) + geom_point(aes(color=prop_var)) + xlab('Max eQTL correlation |r|') +
ylab("Proportion of phenotypic variance explained by QTL")
png('QTT/images/trans_eQTL_r_by_ID_propvar.png')
print(p3)
dev.off()

names(all_perms)=c('pei','rep','perm_max_r')

all_perms=fread('QTT/top10_trans_eqtl_permutations.txt',data.table=F)
n=10
top_r=c()


pheno_env_ids=unique(subcomp$pheno_env_ID)

for(pei in pheno_env_ids){
	subdf=subcomp2[subcomp2$pheno_env_ID==pei,]
	subdf=subdf[order(subdf$value),]
	rownames(subdf)=seq(1,nrow(subdf))
	subdf$rank=seq(nrow(subdf),1)
	max_loc=which.max(abs(subdf$r))
	max_r=max(abs(subdf$r))
	top10=subdf %>% slice_max(value, n = n)
	top10=as.data.frame(top10)
	#top10_r=top10_r[order(top10_r$i.value),]
	#top10_r$rank=seq(1,10)
	top10_r_loc=unlist(which.max(abs(top10$r)))
	top_rank=top10[top10_r_loc,]$rank
	top10_r=abs(top10[top10_r_loc,]$r)
	newline=data.frame(pheno_env_ID=pei,max_r=max_r,top10_r=top10_r,max_rank=top_rank,stringsAsFactors=F)
	top_r=rbind(top_r,newline)
}

all_perms$max_r=top_r[match(all_perms$pei,top_r$pheno_env_ID),]$max_r



p2=ggplot(aes(x=max_r,y=top10_r),data=top_r)  + xlab("Highest correlation eQTL (|r|)") + 
ylab("Highest correlation gene of 10 most significant eQTL (|r|)") + geom_abline(slope=1)

png('QTT/images/trans_eQTL_r_by_sig.png')
print(p2)
dev.off()


plot_list=list()
count=1
for(pei in pheno_env_ids){
	subdf=subcomp2[subcomp2$pheno_env_ID==pei,]
	qtlp=unique(subdf$value)
	p1=ggplot(aes(x=abs(r),y=value),data=subdf) + geom_point() + xlab("eQTL correlation with QTL effect size (|r|)") + 
	ylab("log10(q-value) of eQTL") + ggtitle(sprintf('%s overlapping eQTL (log10(pvalue)= %.2f)',pei,qtlp))
	plot_list[[count]]=p1
	count=count+1
}

pdf('QTT/images/trans_eQTL_r_by_qvalue.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()

max_r= subcomp2 %>% group_by(pheno_env_ID) %>% slice(which.max(abs(r)))
max_r=as.data.frame(max_r)

p3=ggplot(aes(x=abs(r),y=qtl_prop_var),data=max_r) + geom_point(aes(color=prop_var)) + xlab('Max eQTL correlation |r|') +
ylab("Proportion of phenotypic variance explained by QTL")
png('QTT/images/trans_eQTL_r_by_ID_propvar.png')
print(p3)
dev.off()

ft=c("male_flowering_d6","female_flowering_d6")
max_r$ft=ifelse(max_r$phenotype %in% ft,TRUE,FALSE)
p3=ggplot(aes(x=abs(r),y=qtl_prop_var),data=max_r) + geom_point(aes(color=ft)) + xlab('Max eQTL correlation |r|') +
ylab("Proportion of phenotypic variance explained by QTL")
png('QTT/images/trans_eQTL_r_by_FT_propvar.png')
print(p3)
dev.off()

fwrite(max_r,'QTT/trans_max_r_QTL.txt',row.names=F,quote=F,sep='\t')