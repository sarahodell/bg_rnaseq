#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('abind')

# How many genes are in QTL support intervals?
find_qtts=function(x,e,p){
	xexp=exp[,x,drop=F]
	pheno=phenotypes[phenotypes$Loc.Year.Treat==e,c('Genotype_code',p)]
	pheno$exp=xexp[match(pheno$Genotype_code,rownames(xexp)),x]
	test=cor.test(pheno[,p],pheno$exp,use="complete.obs")
	return(test)
}



qtl=fread('QTL/male_flowering_d6_ALL_QTL_scan_hits.txt',data.table=F)
all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}


# Just looking at DTA in EXP_STPAUL
qtl$block_end=all_founder_blocks[match(qtl$SNP,all_founder_blocks$focal_snp),]$end
env2=qtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env2=as.data.table(env2)
env1=as.data.table(genetable)
#env2$end=env2$end-1
setkey(env2,CHR,BP,block_end)
comparison=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('CHR','BP','block_end'),nomatch=NULL,type="within")

time="WD_0712"
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
rownames(exp)=exp$V1
exp=exp[,-1]
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)

metadata=metadata[metadata$experiment==time,]

geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
not_genes=geneh2s[geneh2s$h2==0 ,]$gene

genes=kept_genes
genes=names(exp)

e="ALL"
p="male_flowering_d6"
qgenes=intersect(genes,unique(comparison$Gene_ID))
if(length(qgenes)!=0){
	correlations=sapply(qgenes,function(x) find_qtts(x,e,p))
	tline=data.frame(pheno=p,env=e,gene=qgenes,r=unlist(correlations['estimate',]),pvalue=unlist(correlations['p.value',]),stringsAsFactors=F)
}

tline=as.data.frame(tline,stringsAsFactors=F)
tline=tline[order(tline$pvalue),]
tline$padjust=p.adjust(tline$pvalue,method='BH')


#n=length(unique(comparison$Gene_ID))
quant5=c()
for(i in 1:1000){
	n=length(qgenes)
	draw=sample(not_genes,n)
	correlations=sapply(draw,function(x) find_qtts(x,e,p))
	line=data.frame(pheno=p,env=e,gene=draw,r=unlist(correlations['estimate',]),pvalue=unlist(correlations['p.value',]),stringsAsFactors=F)
	quant5=c(quant5,min(line$pvalue))
}
cutoff=quantile(quant5,0.05)
print(cutoff)
#0.002895952
tline[tline$pvalue<=cutoff ,]



#exp=exp[,kept_genes]


#tline$padjust=p.adjust(tline$pvalue,method='BH')
#fwrite(all_qtts,sprintf('MegaLMM/QTT_QTL_overlap_%s.txt',time),row.names=F,quote=F,sep='\t')

notgenes=genes[!(genes %in% qgenes)]

quant5=c()
for(i in 1:1000){
	n=length(qgenes)
	draw=sample(notgenes,n)
	correlations=sapply(draw,function(x) find_qtts(x,e,p))
	line=data.frame(pheno=p,env=e,gene=draw,r=unlist(correlations['estimate',]),pvalue=unlist(correlations['p.value',]),stringsAsFactors=F)
	quant5=c(quant5,min(line$pvalue))
}
cutoff=quantile(quant5,0.05)
print(cutoff)

tline[tline$pvalue<=cutoff ,]


#### All QTL
qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

fqtl=qtl[qtl$Method=="Founder_probs",]


env2=fqtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env2=as.data.table(env2)
env1=as.data.table(genetable)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL,type="within")
gene_counts=comparison %>% group_by(pheno_env_id) %>% count()

comparison=as.data.frame(comparison,stringsAsFactors=F)
fwrite(comparison,'metadata/QTL_support_interval_genes.txt',row.names=F,quote=F,sep='\t')

find_qtts=function(x,e,p){
	xexp=exp[,x,drop=F]
	
	pheno=phenotypes[phenotypes$Loc.Year.Treat==e,c('Genotype_code',p)]
	pheno$exp=xexp[match(pheno$Genotype_code,rownames(xexp)),x]
	test=cor.test(pheno[,p],pheno$exp,use="complete.obs")
	return(test)
}

time="WD_0712"

phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)


exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
rownames(exp)=exp$V1
exp=exp[,-1]
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)

metadata=metadata[metadata$experiment==time,]


geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
not_genes=geneh2s[geneh2s$h2==0 ,]$gene

genes=kept_genes




#0.002895952
tline[tline$pvalue<=cutoff ,]


envs=unique(comparison$Environment)
phenos=unique(comparison$Phenotype)

all_qtts=c()

for(e in envs){
	for(p in phenos){
		print(p)
		print(e)
		subdf=comparison[comparison$Environment==e & comparison$Phenotype==p,]
		qgenes=intersect(kept_genes,unique(subdf$Gene_ID))
		if(length(qgenes)!=0){
			quant5=c()
			for(i in 1:1000){
				n=length(qgenes)
				draw=sample(not_genes,n)
				correlations=sapply(draw,function(x) find_qtts(x,e,p))
				line=data.frame(pheno=p,env=e,gene=draw,r=unlist(correlations['estimate',]),pvalue=unlist(correlations['p.value',]),stringsAsFactors=F)
				quant5=c(quant5,min(line$pvalue))
			}
			cutoff=quantile(quant5,0.05)
			print(cutoff)
			correlations=sapply(qgenes,function(x) find_qtts(x,e,p))
			tline=data.frame(pheno=p,env=e,gene=qgenes,r=unlist(correlations['estimate',]),pvalue=unlist(correlations['p.value',]),stringsAsFactors=F)
			sig=tline[tline$pvalue<=cutoff,]
			if(nrow(sig)!=0){
				print(sig)
				all_qtts=rbind(all_qtts,sig)
			}
			
		}
	}
	
}

all_qtts=as.data.frame(all_qtts,stringsAsFactors=F)
all_qtts=all_qtts[order(all_qtts$pvalue),]
all_qtts$padjust=p.adjust(all_qtts$pvalue,method='BH')
fwrite(all_qtts,sprintf('MegaLMM/QTT_QTL_overlap_%s.txt',time),row.names=F,quote=F,sep='\t')

qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
fqtl=qtl[qtl$Method=="Founder_probs",]
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
df=fread('eqtl/results/all_eQTL_hits.txt',data.table=F)
df=df[df$class=="factor",]
all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

df$block_start=all_founder_blocks[match(df$SNP,all_founder_blocks$focal_snp),]$start
df$block_end=all_founder_blocks[match(df$SNP,all_founder_blocks$focal_snp),]$end

qtl$block_start=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$start

qtl=fread('../GridLMM/Biogemma_10p_QTL.csv',data.table=F)
df=fread('eqtl/results/all_eQTL_hits.txt',data.table=F)

df$block_start=all_founder_blocks[match(df$SNP,all_founder_blocks$focal_snp),]$start
df$block_end=all_founder_blocks[match(df$SNP,all_founder_blocks$focal_snp),]$end
#qtl$block_start=qtl
#qtl$block_end=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$end

#overlap of SNP 5kb upstream or downstream of SNP
env1=df
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(qtl)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

#F values corr with BLUPs
time="WD_0727"
fvalues=fread(sprintf('MegaLMM/MegaLMM_%s_all_F_means.txt',time),data.table=F)
rownames(fvalues)=fvalues$V1
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:9]
#envs=c('ALL','EXP_STPAUL_2017_WD','GRANEROS_2015_OPT')

all_qtts=c()
for(e in envs){
	for(p in phenos){
		print(p)
		print(e)
		subdf=phenotypes[phenotypes$Loc.Year.Treat==e,c('Genotype_code',p)]
		rownames(subdf)=subdf$Genotype_code
		factors=names(fvalues)[-1]
		inds=intersect(subdf$Genotype_code,fvalues$V1)
		subdf=subdf[inds,]
		fvalues=fvalues[inds,]
		for(factor in factors){
			subdf[,factor]=fvalues[,factor]
			test=cor.test(subdf[,p],subdf[,factor],use="complete.obs")
			#correlations=sapply(qgenes,function(x) find_qtts(x,e,p))
			tline=data.frame(pheno=p,env=e,factor=factor,r=unlist(test['estimate']),pvalue=unlist(test['p.value']),stringsAsFactors=F)			
			all_qtts=rbind(all_qtts,tline)
		}
	}
}

all_qtts=all_qtts[order(all_qtts$pvalue),]
fwrite(all_qtts,sprintf('MegaLMM/F_QTTs_%s.txt',time),row.names=F,quote=F,sep='\t')


		qgenes=intersect(kept_genes,unique(subdf$Gene_ID))
		if(length(qgenes)!=0){
			quant5=c()
			for(i in 1:1000){
				n=length(qgenes)
				draw=sample(not_genes,n)
				correlations=sapply(draw,function(x) find_qtts(x,e,p))
				line=data.frame(pheno=p,env=e,gene=draw,r=unlist(correlations['estimate',]),pvalue=unlist(correlations['p.value',]),stringsAsFactors=F)
				quant5=c(quant5,min(line$pvalue))
			}
			cutoff=quantile(quant5,0.05)
			print(cutoff)
			correlations=sapply(qgenes,function(x) find_qtts(x,e,p))
			tline=data.frame(pheno=p,env=e,gene=qgenes,r=unlist(correlations['estimate',]),pvalue=unlist(correlations['p.value',]),stringsAsFactors=F)
			sig=tline[tline$pvalue<=cutoff,]
			if(nrow(sig)!=0){
				print(sig)
				all_qtts=rbind(all_qtts,sig)
			}
			
		}
	}
	
}


# WD_0712 2010 tests max 0.328
# WD_0718 2483 max 0.255 r
# WD_0720 2298 max -0.22 r
# WD_0727 1376 max -0.222

qtts %>% group_by(env_pheno) %>% count()



rap27="Zm00001d010987"
qtts[qtts$Gene==rap27,]
#env_pheno                                    n
#<chr>                                    <int>
#1 BLOIS_2014_OPT-female_flowering_d6           8
#2 BLOIS_2014_OPT-male_flowering_d6             4
#3 BLOIS_2017_OPT-female_flowering_d6          20
#4 BLOIS_2017_OPT-male_flowering_d6           152
#5 BLOIS_2017_OPT-tkw_15                       76
#6 GRANEROS_2015_OPT-female_flowering_d6      124
#7 GRANEROS_2015_OPT-harvest_grain_moisture    30
#8 GRANEROS_2015_OPT-male_flowering_d6         56
#9 GRANEROS_2015_OPT-tkw_15                    13
#10 NERAC_2016_WD-female_flowering_d6            2
#11 NERAC_2016_WD-male_flowering_d6            109
#12 STPAUL_2017_WD-male_flowering_d6            77
#13 SZEGED_2017_OPT-female_flowering_d6         44
#14 SZEGED_2017_OPT-male_flowering_d6            9
