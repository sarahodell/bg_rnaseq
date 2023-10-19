#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('stringr')
library('ggrepel')
library('DescTools')
library('dplyr')

all_perms=fread('QTT/z_permutations.txt',data.table=F)

truediff=fread('QTT/local_distal_max_z.txt',data.table=F)


pot=truediff[truediff$zdiff>0,]
#pot=truediff
pheno_env_ids=pot$pei

high=c()
low=c()
sig=c()
for(i in 1:nrow(truediff)){
	row1=truediff[i,]
	pei=row1$pei
	subperm=all_perms[all_perms$pei==pei,]
	quant=quantile(subperm$difference,0.95)
	quant2=quantile(subperm$difference,0.05)
	space=row1$zdiff-quant
	perm_high=row1$local_z-quant
	perm_low=row1$local_z-quant2
	
	# distance of true z_diff from higher perm z_diff
	
	high=c(high,perm_high)
	low=c(low,perm_low)
	# Figure out how to get slope of permutations 95% CI 
	true=row1$zdiff
	if(true>=quant){
		sig=c(sig,TRUE)
	}else{
		sig=c(sig,FALSE)
	}
}

truediff$sig=sig
truediff$perm_high=high
truediff$perm_low=low


p1=ggplot(truediff,aes(x=distal_z,y=local_z)) + geom_point(aes(color=sig)) +
	xlab("Max |z| of distal-eQTL") + ylab("Max |z| of local-eQTL") +
	geom_abline(slope=1) + 
	geom_segment(data=truediff,aes(x=perm_high,xend = perm_low)) +
	theme_classic() + 
	#geom_text_repel(data=subset(truediff, sig==TRUE),aes(x=local_z,y=distal_z,label=pei)) +
	scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) + 
	ggtitle("Most Correlated eQTL per QTL") + guides(color="none")


png('paper_figures/eQTL_z_values.png')
print(p1)
dev.off()

cand=fread('QTT/sig_candidate_genes.txt',data.table=F)
comparison2=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
genecount=comparison2 %>% group_by(pheno_env_ID) %>% count()
pheno_env_ids=unique(comparison2$pheno_env_ID)
max_r_l=fread('QTT/local_eQTL_candidates.txt',data.table=F)
max_r_l$qtl_prop_var=max_r_l$prop_var

eqtl_prop_vars= comparison2 %>% group_by(gene_time_SNP) %>% summarize(prop_var=unique(prop_var))
max_r_l$eqtl_prop_var=eqtl_prop_vars[match(max_r_l$gene_time_SNP,eqtl_prop_vars$gene_time_SNP),]$prop_var

fwrite(max_r_l,'QTT/local_eQTL_candidates.txt',row.names=F,quote=F,sep='\t')

p2=ggplot(aes(x=max_r,y=top10_r),data=top_r) + geom_point() + xlab("Highest correlation eQTL (|r|)") + 
ylab("Highest correlation gene of 10 most significant eQTL (|r|)") + geom_abline(slope=1) +
theme_classic()

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


top_r$top10_z=FisherZ(abs(top_r$top10_r))
top_r$max_z=FisherZ(abs(top_r$max_r))

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

top_r$sig=top_r$pheno_env_ID %in% cand$pei
top_r$ngenes=genecount[match(top_r$pheno_env_ID,genecount$pheno_env_ID),]$n

p2=ggplot(aes(x=max_z,y=top10_z),data=top_r) + geom_point(data=subset(top_r,sig==FALSE),aes(size=ngenes),color='black') +
geom_point(data=subset(top_r,sig==TRUE),aes(size=ngenes),color='red') +
 xlab("Highest correlation eQTL (|z|)") + 
ylab("Highest correlation gene of 10 most significant eQTL (|z|)") + geom_abline(slope=1) +
theme_classic() + scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) + guides(color="none") +
labs(size= "# Genes")

png('paper_figures/local_eQTL_r_by_sig.png')
print(p2)
dev.off()


ft=c("male_flowering_d6","female_flowering_d6")
max_r_l$ft=max_r_l$phenotype %in% ft
max_r_l$max_z=FisherZ(abs(max_r_l$r))

p3=ggplot(aes(x=max_z,y=qtl_prop_var),data=max_r_l) + #geom_point(aes(color=ft)) +
geom_point(aes(color=sig)) + scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) +
theme_classic() +
 xlab('Max eQTL correlation |z|') +
ylab("Proportion of phenotypic variance explained by QTL")  #+ labs(color="Flowering Time QTL")
png('paper_figures/eQTL_z_by_qtl_propvar.png')
print(p3)
dev.off()

#cand=pot[pot$sig==TRUE,]
fwrite(cand,'QTT/sig_candidate_genes.txt',row.names=F,quote=F,sep='\t')

plot_list=list()
count=1
#### Look at correaltion of bv with FT
for(i in 1:nrow(cand)){
	row1=cand[i,]
	chr=row1$CHR
	pheno_env_ID=cand$pei
	#split1=strsplit(pheno_env_ID,'_')
	pheno=row1$phenotype
	env=row1$environment
	id=row1$ID
	time1=row1$time
	gene=row1$Trait
	exp=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_bv_FIXED.txt',time1,chr),data.table=F)
	rownames(exp)=exp$V1
	exp=exp[,-1]
	phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
	phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
	rownames(phenotypes)=phenotypes$Genotype_code
	
	inter=intersect(rownames(phenotypes),rownames(exp))
	phenotypes=phenotypes[inter,]
	exp=exp[inter,]
	df=phenotypes[,c('Genotype_code',pheno)]
	df$exp_bv=exp[match(df$Genotype_code,rownames(exp)),gene]
	beta_r=row1$r
	names(df)=c('ID','pheno','exp_bv')
	bv_r=cor(df$pheno,df$exp_bv,use="complete.obs")

	p1=ggplot(df,aes(x=exp_bv,y=pheno)) + geom_point() + xlab("Additive Local Effect on Expression") +
	ylab('Phenotype Value') + ggtitle(sprintf('%s %s',pheno_env_ID,gene),subtitle=sprintf("beta r=%.2f, bv r = %.2f",beta_r,bv_r)) 
	
	plot_list[[count]]=p1
	count=count+1
	
}

pdf('QTT/candidate_bv_pheno.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()

#topcand=pot[pot$sig==TRUE,]

#                                            pei  local_z  distal_z     zdiff
#4  male_flowering_d6-EXP_STPAUL_2017_WD-qDTA3_2 1.256848 1.0084874 0.2483607
#12     male_flowering_d6-STPAUL_2017_WD-qDTA3_2 1.239797 0.9018092 0.3379882
#14     male_flowering_d6-BLOIS_2017_OPT-qDTA3_2 1.431885 1.0378326 0.3940528
#30       male_flowering_d6-STPAUL_2017_WD-qDTA8 1.552356 1.1160330 0.4363226
#35   male_flowering_d6-EXP_STPAUL_2017_WD-qDTA8 1.267919 1.0293634 0.2385558
#    sig
#4  TRUE
#12 TRUE
#14 TRUE
#30 TRUE
#35 TRUE

max_r_l=fread('QTT/local_eQTL_candidates.txt',data.table=F)
pei=cand$pei
max_r_l=max_r_l[max_r_l$pheno_env_ID %in% pei,]

localcomp=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
cgenes=max_r_l$Trait

for(gene in cgenes){
	id=max_r_l[max_r_l$Trait==gene,]$ID
	print(gene)
	print(localcomp[localcomp$Trait==gene & localcomp$ID==id ,c('environment','time','r')])
}
#"Zm00001d011123" "Zm00001d011294" had basially no corrleation with flowering time BLOIS_2014_OPT
# For all other environments, the WD_0712 correlations were also high
# Plot out correlation of candidate genes

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
plot_list=list()
count=1

for(i in 1:nrow(max_r_l)){
	row=max_r_l[i,]
	pheno=row$phenotype
	env=row$environment
	chr=row$CHR
	gene=row$Trait
	qsnp=row$SNP
	id=row$ID
	pei=row$pheno_env_ID
	gts=row$gene_time_SNP
	time=row$time
	gene=row$Trait
	esnp=row$X_ID
	
	effect_sizes=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_unscaled_founderprobs.txt',chr,pheno,env),data.table=F)
	effect_size=effect_sizes[effect_sizes$X_ID==qsnp,]
	effect_size=unlist(effect_size[,c(6:21)])
	wn=which(!is.na(effect_size))[1]
	effect_size[-wn]=effect_size[-wn]+effect_size[wn]
	
	results=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time,chr),data.table=F)
	results=results[results$X_ID==esnp & results$Trait==gene,]
	betas=unlist(results[,c(6,10:24)])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[-wn]+betas[wn]
	
	test=cor.test(effect_size,betas,use="complete.obs")
	r=test$estimate
	p=test$p.value
	
	df=data.frame(founder=founders,qtl_es=effect_size,eqtl_es=betas,stringsAsFactors=F)
	
	p1=ggplot(df,aes(x=eqtl_es,y=qtl_es)) + geom_point(aes(color=founder)) +
	xlab("eQTL effect size (log2cpm)") + ylab("FT effect size (gdd)") +
	ggtitle(sprintf("%s and %s,",gts,pei),subtitle=sprintf('r=%.2f',r)) +
	theme(plot.title = element_text(size = 9))
	plot_list[[count]]=p1
	count=count+1
	
}

pdf('QTT/images/candidate_gene_effect_size_correlations.pdf')
for(i in 1:length(plot_list)){
	print(plot_list[[i]])
}
dev.off()

prop_var=fread('MegaLMM/MegaLMM_WD_0712_prop_variance_FIXED.txt',data.table=F)
# Which factor are these genes most loaded on?
rownames(prop_var)=prop_var$V1
prop_var=prop_var[,-1]

apply(prop_var,MARGIN=2,function(x) sum(cgenes %in% rownames(prop_var)[which(x>0.1)]))
# Factor 4 and 5 have 3 of the 5 genes loaded on them
# In WD, Factor 3, Factor 8 are enriched for FT genes (from the list)
# which 3
# 2 and 12 have factoreqtl

# Factor 4
f4=prop_var[prop_var$Factor4>0.1,c('Factor4'),drop=F]
f4[cgenes,]
#"Zm00001d042291" 0.167
#"Zm00001d041900" 0.361
#"Zm00001d042306" 0.1109

# Factor 5
f5=prop_var[prop_var$Factor5>0.1,c('Factor5'),drop=F]
#"Zm00001d042291" 0.3430424
#"Zm00001d011294" 0.2088616
#"Zm00001d042306" 0.1805666

# Factor 6 
#Zm00001d042291 0.143712571
#Zm00001d042306 0.489615495

# Factor 10
#Zm00001d041900 0.1502910533

# factor 13 
#Zm00001d011123 0.456330076


# Looking at DE

prop_var=fread('MegaLMM/MegaLMM_residauls_WD_0712_prop_variance_FIXED.txt',data.table=F)
rownames(prop_var)=prop_var$V1
prop_var=prop_var[,-1]

apply(prop_var,MARGIN=2,function(x) sum(cgenes %in% rownames(prop_var)[which(x>0.1)]))

# factor eqtl for Factor2, (chr 2), Factor 9, (chr 5) Factor 23 (chr 4,6,7)

# Factor 1
#Zm00001d042291 0.161617299
#Zm00001d041900 0.207356319
#Zm00001d042306 0.235207523

# Factor 2
#Zm00001d042291 0.10637893

# Factor 4
#Zm00001d042291 0.1231471134
#Zm00001d041900 0.1590809554

# Factor 6
#Zm00001d011294 0.1065955176

# Factor 7
#Zm00001d042306 0.224203031

# Factor 9
#Zm00001d042291 0.2176312512

# Factor 23
#Zm00001d042291 0.1061234

# Factor 34
#Zm00001d041900 0.1411403

# Factor 35 
#Zm00001d011294 0.1408263328
#Zm00001d011123 0.5836578358

# Factor 42 
#Zm00001d041900 0.301303567
#Zm00001d011294 0.107176388

# Factor 44
#Zm00001d011123	0.2596519888

# Are the F-values of these factor-eQTL 
