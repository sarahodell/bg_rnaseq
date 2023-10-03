#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')

all_perms=fread('QTT/%s_z_permutations.txt',data.table=F)

truediff=fread('QTT/local_distal_max_z.txt',data.table=F)


pot=truediff[truediff$zdiff>0,]
#pot=truediff
pheno_env_ids=pot$pei

sig=c()
for(i in 1:nrow(pot)){
	row1=pot[i,]
	pei=row1$pei
	subperm=all_perms[all_perms$pei==pei,]
	quant=quantile(subperm$difference,0.95)
	true=pot[pot$pei==pei,]$zdiff
	if(true>=quant){
		sig=c(sig,TRUE)
	}else{
		sig=c(sig,FALSE)
	}
}

pot$sig=sig
cand=pot[pot$sig==TRUE,]
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
