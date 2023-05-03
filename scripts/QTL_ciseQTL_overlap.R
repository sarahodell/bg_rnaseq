library('abind')
library('data.table')
library('ggplot2')
library('dplyr')

time="WD_0718"
#ciseqtl=fread('eqtl/results/all_cis_eQTL_vst_hits.txt',data.table=F)
ciseqtl=fread('eqtl/results/WD_0712_cis_eQTL_fkeep_hits.txt',data.table=F)

times=c("WD_0712","WD_0720","WD_0727")
all_cis_eqtl=c()
for(t in times){
  ciseqtl=fread(sprintf('eqtl/results/%s_cis_eQTL_fkeep_hits.txt',t),data.table=F)
  ciseqtl$time=t
  all_cis_eqtl=rbind(all_cis_eqtl,ciseqtl)
}
fwrite(all_cis_eqtl,'eqtl/results/all_cis_eQTL_fkeep_hits2.txt',row.names=F,quote=F,sep='\t')

ciseqtl=fread('eqtl/results/all_cis_eQTL_fkeep_hits2.txt',data.table=F)

#ciseqtl=all_cis_eqtl
#ciseqtl=fread(sprintf('eqtl/results/%s_cis_eQTL_vst_hits.txt',time),data.table=F)
qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
fqtl=qtl[qtl$Method=="Founder_probs",]
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)


#Look 5kb upstream or downstream of gene
genetable$window_START=genetable$START-5000
genetable$window_END=genetable$END+5000

# Does it make more sense to use 10kb region around the significant SNP
ciseqtl$window_START=ciseqtl$BP-5000
ciseqtl$window_END=ciseqtl$BP+5000
# Or use the founder LD block? I think this one is more accurate
all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

ciseqtl$block_start=all_founder_blocks[match(ciseqtl$SNP,all_founder_blocks$focal_snp),]$start
ciseqtl$block_end=all_founder_blocks[match(ciseqtl$SNP,all_founder_blocks$focal_snp),]$end

#overlap of SNP 5kb upstream or downstream of SNP
env1=ciseqtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(fqtl)
#env2$end=env2$end-1
setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
comparison=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL)

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
