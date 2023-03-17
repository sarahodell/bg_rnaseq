#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('abind')

# How many genes are in QTL support intervals?

qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

fqtl=qtl[qtl$Method=="Founder_probs",]
fids=unique(fqtl$ID)
notin=qtl[!(qtl$ID %in% fids),]
fqtl=rbind(fqtl,notin)
fwrite(fqtl,'metadata/unique_QTL.csv',row.names=F,quote=F,sep=',')

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

#test_qtl=qtl[qtl$Chromosome==3,]
#test_genes=genetable[genetable$CHROM==3,]
#env2=qtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
#env2=as.data.table(env2)
#env1=as.data.table(test_genes)
##env2$end=env2$end-1
#setkey(env2,left_bound_bp,alt_right_bound_bp)
#comparison=foverlaps(env1,env2,by.x=c('START','END'),by.y=c('left_bound_bp','alt_right_bound_bp'),nomatch=NULL,type="within")
qtts=c()

# How many QTTs are in the support intervals of relevant QTL?
factor_groups=readRDS('MegaLMM/pheno_MegaLMM_WD_0718_factor_groups.rds')
qtl=qtl[qtl$Environment!="ALL",]
qtl$env_pheno=paste0(qtl$Environment,'-',qtl$Phenotype)
env_phenos=unique(paste0(qtl$Environment,'-',qtl$Phenotype))
pheno_df=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_sig_factors.txt',time),data.table=F)
for(pe in env_phenos){
  factors=which(unlist(unname(lapply(factor_groups,function(x) pe %in% x$phenotype))))
  print(pe)
  print(factors)
  for(f in factors){
    genelist=factor_groups[[f]]$genes
    subgenes=genetable[genetable$Gene_ID %in% genelist,]
    subqtl=qtl[qtl$env_pheno==pe,]

    env2=subqtl
    env2=as.data.table(env2)
    env1=as.data.table(subgenes)
    #env2$end=env2$end-1
    setkey(env2,Chromosome,left_bound_bp,alt_right_bound_bp)
    comparison=foverlaps(env1,env2,by.x=c('CHROM','START','END'),by.y=c('Chromosome','left_bound_bp','alt_right_bound_bp'),nomatch=NULL,type="within")
    comparison$factor=f
    #gene_counts=comparison %>% group_by(pheno_env_id) %>% count()
    qtts=rbind(qtts,comparison)

  }
}

qtts=as.data.frame(qtts,stringsAsFactors=F)
fwrite(qtts,'MegaLMM/QTT_QTL_overlap_WD_0718.txt',row.names=F,quote=F,sep='\t')

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
