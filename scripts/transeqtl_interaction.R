#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('ggplot2')
library('lme4')
library('emmeans')

# trans eQTL locations
#hits=fread('eqtl/results/all_eQTL_hits.txt',data.table=F)
#hits$genesnp=paste0(hits$Trait,'_',hits$SNP)

eqtl=fread('eqtl/results/cis_trans_eQTL_hits.txt',data.table=F)
eqtl=eqtl[eqtl$class=="trans",]
#eqtl$genesnp=paste0(eqtl$Trait,'_',eqtl$SNP)
#eqtl$time=hits[match(eqtl$genesnp,hits$genesnp),]$time
#fwrite(eqtl,'eqtl/results/cis_trans_eQTL_hits.txt',row.names=F,quote=F,sep='\t')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

sig=c()

times=c('WD_0712','WD_0718','WD_0720','WD_0727')
for(time in times){
	print(time)
	samples=fread(sprintf('eqtl/data/%s_samples.txt',time),data.table=F,header=F)
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
	rownames(exp)=exp$V1
	subeqtl=eqtl[eqtl$time==time,]
	for(i in 1:nrow(subeqtl)){
		row=subeqtl[i,]
		snp=row$SNP
		snpchr=row$CHR
		gene=row$Trait
		genechr=row$gene_chr
		start=row$gene_start
		snp_X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',snpchr))
		snp_X_list=lapply(snp_X_list,function(x) x[samples$V1,])
		snpX = do.call(cbind,lapply(snp_X_list,function(x) x[,snp]))
    	frep1=apply(snpX,MARGIN=2,function(x) round(sum(x[x>0.75])))
    	fkeep1=founders[frep1>3]
    	certain1=apply(snpX,MARGIN=1,function(x) sum(x>0.75)>0)
    	xexp=exp[,c('V1',gene)]
    	names(xexp)=c('ind','gene')
    	gene_X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',genechr))
		gene_X_list=lapply(gene_X_list,function(x) x[samples$V1,])
		pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',genechr),data.table=F)
		pmap=pmap[pmap$marker %in% colnames(gene_X_list[[1]]),]
		pmap$end=pmap$pos+1
		lower=pmap[pmap$pos<=start,]
		leftindex=which.min(start-lower$pos)
		leftmarker=lower[leftindex,]
		geneX = do.call(cbind,lapply(gene_X_list,function(x) x[,leftmarker$marker]))
		frep2=apply(geneX,MARGIN=2,function(x) round(sum(x[x>0.75])))
    	fkeep2=founders[frep2>3]
    	certain2=apply(geneX,MARGIN=1,function(x) sum(x>0.75)>0)
    	comb=intersect(names(certain1[certain1==T]),names(certain2[certain2==T]))
    	snpX=snpX[comb,]
    	geneX=geneX[comb,]
		xexp=xexp[comb,]
		xexp$snpfounder=unlist(unname(apply(snpX,MARGIN=1,function(x) colnames(snpX)[which.max(x)])))
		xexp$genefounder=unlist(unname(apply(geneX,MARGIN=1,function(x) colnames(geneX)[which.max(x)])))
		m1=lm(gene~snpfounder*genefounder,xexp)
		print(gene)
		print(snp)
		print(anova(m1))  
	}
}

sig=c()

sig=c("Zm00001d039421_AX-90780424","Zm00001d022126_AX-91451546","Zm00001d047592_AX-90965571","Zm00001d047592_AX-90972801")


interactions=eqtl[eqtl$genesnp %in% sig,]

models=vector("list",length=4)
for(i in 1:nrow(interactions)){
		row=interactions[i,]
		time=row$time
		samples=fread(sprintf('eqtl/data/%s_samples.txt',time),data.table=F,header=F)
		exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
		rownames(exp)=exp$V1
		snp=row$SNP
		snpchr=row$CHR
		gene=row$Trait
		genechr=row$gene_chr
		start=row$gene_start
		snp_X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',snpchr))
		snp_X_list=lapply(snp_X_list,function(x) x[samples$V1,])
		snpX = do.call(cbind,lapply(snp_X_list,function(x) x[,snp]))
    	frep1=apply(snpX,MARGIN=2,function(x) round(sum(x[x>0.75])))
    	fkeep1=founders[frep1>3]
    	certain1=apply(snpX,MARGIN=1,function(x) sum(x>0.75)>0)
    	xexp=exp[,c('V1',gene)]
    	names(xexp)=c('ind','gene')
    	gene_X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',genechr))
		gene_X_list=lapply(gene_X_list,function(x) x[samples$V1,])
		pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',genechr),data.table=F)
		pmap=pmap[pmap$marker %in% colnames(gene_X_list[[1]]),]
		pmap$end=pmap$pos+1
		lower=pmap[pmap$pos<=start,]
		leftindex=which.min(start-lower$pos)
		leftmarker=lower[leftindex,]
		geneX = do.call(cbind,lapply(gene_X_list,function(x) x[,leftmarker$marker]))
		frep2=apply(geneX,MARGIN=2,function(x) round(sum(x[x>0.75])))
    	fkeep2=founders[frep2>3]
    	certain2=apply(geneX,MARGIN=1,function(x) sum(x>0.75)>0)
    	comb=intersect(names(certain1[certain1==T]),names(certain2[certain2==T]))
    	snpX=snpX[comb,]
    	geneX=geneX[comb,]
		xexp=xexp[comb,]
		xexp$snpfounder=unlist(unname(apply(snpX,MARGIN=1,function(x) colnames(snpX)[which.max(x)])))
		xexp$genefounder=unlist(unname(apply(geneX,MARGIN=1,function(x) colnames(geneX)[which.max(x)])))
		m1=lm(gene~snpfounder*genefounder,xexp)
		print(gene)
		print(snp)
		print(anova(m1))  
		models[[i]]$model=m1
		models[[i]]$data=xexp
		models[[i]]$gene=gene
		models[[i]]$snp=snp
		
		p1=emmip(m1,genefounder~snpfounder,CIs = T)

		pdf(sprintf('eqtl/trans/images/%s_%s_%s_interaction.pdf',gene,snp,time))
		print(p1)
		dev.off()
}

model1=models[[1]]$model

means_by_snp = emmeans(model1,specs = 'snpfounder',by = 'genefounder')
snp_founder_effects = contrast(means_by_snp,
                           method = 'pairwise',
                           name = 'snp_founder_effect')
print(summary(snp_founder_effects,infer = c(T,T)))
regrouped_snp_effects = update(snp_founder_effects,by = 'snp_founder_effect',adjust = 'none')

interaction_effects = contrast(regrouped_snp_effects,
                                          'pairwise',
                                          name='gene_founder_effects_on_snp_effects')
                                          
print(summary(interaction_effects,infer = c(T,T),level = 1-0.05/16))

plot(interaction_effects)

emmip(model1,genefounder~snpfounder,CIs = T)
#incidence_mat = is_crossed(snpfounder~genefounder,model1)
#incidence_mat
#crossprod(incidence_mat)

ggplot(RCBD_data,aes(x=Variety,y=Yield)) + 
  geom_jitter(aes(color = Chamber),width = 0.1) + 
  facet_wrap(~Soil,labeller = label_both)

# FV252_inra - OH43_inra
# OH43_inra - VA85 


print(regrouped_variety_effects)



