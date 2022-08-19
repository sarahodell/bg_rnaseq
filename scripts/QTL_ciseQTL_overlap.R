library('abind')
library('data.table')

time="WD_0712"

ciseqtl=fread(sprintf('eqtl/results/%s_cis_eQTL_vst_hits.txt',time),data.table=F)
qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
qtl=qtl[qtl$Method=="Founder_probs",]
genetable=fread('eqtl/data/Zea_mays.B73_v4_generanges.txt',data.table=F)

eqtl_betas=c()
for(chr in 1:10){
  eqtl_beta=fread(sprintf('eqtl/cis/results/eQTL_%s_c%.0f_vst_results.txt',time,chr),data.table=F)
  eqtl_betas=rbind(eqtl_betas,eqtl_beta)
}
eqtl_betas$snp_gene=paste0(eqtl_betas$X_ID,'_',eqtl_betas$Trait)



#all_founder_blocks=c()
#for(chr in 1:10){#
#  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
#  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
#}

#overlap of SNP

env1=ciseqtl
env1$BP_end=env1$BP
env1=as.data.table(env1)
env2=as.data.table(qtl)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','BP','BP_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

#overlap of gene location +1 10kb? Look at papers. what do they use?

ciseqtl_genetable=genetable[genetable$Gene_ID %in% ciseqtl$Gene,]
ciseqtl_genetable=ciseqtl_genetable[grepl('T001',ciseqtl_genetable$TXNAME),]
i=47920
ciseqtl_genetable=rbind(ciseqtl_genetable,genetable[i,])
ciseqtl_genetable$TXCHROM=as.integer(ciseqtl_genetable$TXCHROM)
rownames(ciseqtl_genetable)=seq(1,nrow(ciseqtl_genetable))
env1=ciseqtl_genetable
env1=as.data.table(env1)
env2=as.data.table(qtl)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison2=foverlaps(env1,env2,by.x=c('TXCHROM','TXSTART','TXEND'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

#overlap of founder recombination blocks - does this make sense?

## Correlation of effect sizes
# Grab founder effect sizes from QTL_F -
# Grab founder effect sizes from cis-eQTL test from GridLMM
pvalue=c()
es_cor=c()

for(q in 1:nrow(comparison)){
  row=comparison[q,]
  id=row$pheno_env_id
  chr=row$CHR
  pheno=row$Phenotype
  env=row$Environment
  effect_sizes=readRDS(sprintf('../GridLMM/GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_vst_founderprobs.rds',chr,pheno,env))
  effect_sizes=unlist(unname(effect_sizes[effect_sizes$X_ID==row$highest_SNP,6:21]))
  effect_sizes[-1]=effect_sizes[1] + effect_sizes[-1]
  #test=which(unlist(unname(lapply(effect_sizes,function(x) x$i==id))))
  #es=effect_sizes[[test]]$values
  tmp=eqtl_betas[eqtl_betas$X_ID==row$SNP & eqtl_betas$Trait==row$Gene,c(6,10:24)]
  eqs=unlist(unname(tmp))
  intercept=min(which(!is.na(eqs)))
  eqs[-intercept]=eqs[intercept]+eqs[-intercept]
  t=cor.test(eqs,effect_sizes)
  pvalue=c(pvalue,t$p.value)
  es_cor=c(es_cor,t$estimate)
}

comparion=as.data.frame(comparison,stringsAsFactors=F)
comparison$es_cor=es_cor
comparison$cortest_pvalue=pvalue

fwrite(comparison,sprintf('eqtl/results/%s_eQTL_QTL_overlap.txt',time),row.names=F,quote=F,sep='\t')


comparison=fread(sprintf('eqtl/results/%s_eQTL_QTL_overlap.txt',time),data.table=F)
#hasf=c('S_F_H','F_and_H','F_only','S_and_F')
# Null expectation
infosnps=unique(comparison$highest_SNP)
eqtlsnps=unique(comparison$SNP)

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
  tmp=tmp1[,c(6,10:24)]
  rownames(tmp)=seq(1,nrow(tmp))
  for(t in 1:nrow(tmp)){
    crow=tmp[t,]
    intercept=min(which(!is.na(crow)))
    eqs=unlist(unname(crow))
    eqs[-intercept]=eqs[intercept]+eqs[-intercept]
    tmp[t,]=eqs
  }

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
