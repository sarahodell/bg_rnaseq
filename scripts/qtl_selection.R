#!/usr/bin/env Rscript

library('data.table')
library('dplyr')


ciseqtl=fread('eqtl/results/all_cis_eQTL_fkeep_hits.txt',data.table=F)
qtl=fread('../GridLMM/Biogemma_QTL.csv',data.table=F)
fqtl=qtl[qtl$Method=="Founder_probs",]

chi_segs=c()
for(chr in 1:10){
  rep=fread(sprintf('../selection/founder_probs/bg%.0f_chi_peak_seqments.bed',chr))
  chi_segs=rbind(chi_segs,rep)
}
chi_segs=as.data.frame(chi_segs,stringsAsFactors=F)
names(chi_segs)=c('CHR','START','END','markers')

#Overlap of chi_seg with qtl
fqtl$block_start=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$start
fqtl$block_end=all_founder_blocks[match(fqtl$highest_SNP,all_founder_blocks$focal_snp),]$end



env1=fqtl
env1=as.data.table(env1)
env2=as.data.table(chi_segs)
#env2$end=env2$end-1
setkey(env2,CHR,START,END)
chi_seg_qtl=foverlaps(env1,env2,by.x=c('Chromosome','block_start','block_end'),by.y=c('CHR','START','END'),nomatch=NULL)
length(unique(chi_seg_qtl$pheno_env_id))
# 19 (out of 36)    ONLY 6 with founder blocks of highest SNP, rather than support intervals
length(unique(chi_seg_qtl$ID))
# 8 (out of 14)

#Overlap of chi_seg with cis-eqtl
all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

ciseqtl$block_start=all_founder_blocks[match(ciseqtl$SNP,all_founder_blocks$focal_snp),]$start
ciseqtl$block_end=all_founder_blocks[match(ciseqtl$SNP,all_founder_blocks$focal_snp),]$end



env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(chi_segs)
#env2$end=env2$end-1
setkey(env2,CHR,START,END)
chi_seg_eqtl=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR','START','END'),nomatch=NULL)
length(unique(chi_seg_eqtl$Gene))
# 15 (out of 40)


# What about high inter-chromosomal LD regions?
ld=fread('../stats/ld_decay/Biogemma_DHLines_rsquared_all_chroms_r2_0.9.ld',data.table=F)
ld=ld[ld$CHR_A != ld$CHR_B,]

ld2=fread('../stats/ld_decay/circos/ld_bundled_links_filtered.txt',data.table=F)
names(ld2)=c('CHR_A','START_A','END_A','CHR_B','START_B','END_B','SIZE')
ld2$CHR_A=sapply(seq(1,nrow(ld2)),function(x) as.numeric(strsplit(ld2$CHR_A[x],'chr')[[1]][2]))
ld2$CHR_B=sapply(seq(1,nrow(ld2)),function(x) as.numeric(strsplit(ld2$CHR_B[x],'chr')[[1]][2]))

#QTL
env1=fqtl
env1=as.data.table(env1)
env2=as.data.table(ld2)
#env2$end=env2$end-1
setkey(env2,CHR_A,START_A,END_A)
ld_qtl1=foverlaps(env1,env2,by.x=c('Chromosome','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR_A','START_A','END_A'),nomatch=NULL)
setkey(env2,CHR_B,START_B,END_B)
ld_qtl2=foverlaps(env1,env2,by.x=c('Chromosome','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR_B','START_B','END_B'),nomatch=NULL)
length(unique(ld_qtl1$pheno_env_id))
# 1 male_flowering_d6_GRANEROS_2015_OPT_qDTA7, qDTA7
length(unique(ld_qtl2$pheno_env_id))

#only one overlap with founder block aroundhighest SNP
# 14 with QTL support intervals qDTS8, qDTA8, qDTA7, qTKW7_1,qHGM3_1,qHGM3_2

# cis-eQTL
env1=ciseqtl
env1=as.data.table(env1)
env2=as.data.table(ld2)
#env2$end=env2$end-1
setkey(env2,CHR_A,START_A,END_A)
ld_eqtl1=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR_A','START_A','END_A'),nomatch=NULL)
setkey(env2,CHR_B,START_B,END_B)
ld_eqtl2=foverlaps(env1,env2,by.x=c('CHR','block_start','block_end'),by.y=c('CHR_B','START_B','END_B'),nomatch=NULL)
length(unique(ld_eqtl1$Gene))
# 12 total genes - none of them on chr7
length(unique(ld_eqtl2$Gene))
