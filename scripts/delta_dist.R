#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('stringr')

truediff=fread('QTT/local_distal_max_z.txt',data.table=F)

truediff$pheno=sapply(seq(1,nrow(truediff)),function(x) strsplit(truediff$pei[x],'-')[[1]][1])
truediff$env=sapply(seq(1,nrow(truediff)),function(x) strsplit(truediff$pei[x],'-')[[1]][2])

theme_set(theme_classic())
#theme_update(text=element_text(family="Times"))
theme_update(plot.caption = element_text(hjust = 0))
theme_update(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))
theme_update(plot.title = element_text(size=20),axis.title=element_text(size=18))
theme_update(panel.background=element_blank())
theme_update(plot.caption=element_text(size=20))

p1=ggplot(aes(x=zdiff),data=truediff) + geom_histogram()


genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
distal=fread('QTT/distal_eQTL_candidates.txt',data.table=F)
# tkw_15-GRANEROS_2015_OPT-qTKW7_2  Zm00001d052810 T12  r=0.9632258  4 201667424 201671532
# also has a distal-eQTL overlapping qDTA8, qDTS8, qHGM7

# harvest_grain_moisture-GRANEROS_2015_OPT-qHGM7 Zm00001d047181 T12 r=0.9409329 9 121092224 121099700 
# also has a distal-eQTL overlapping qHGM2, qTKW7_2

# male_flowering_d6-BLOIS_2017_OPT-qDTA3_1 Zm00001d015469 r=0.8342525 5 91973542 91980501 (dys2)
# also has distal eQTL overlapping qDTA9, qDTS3_1, and qHGM3_1 

# worth noting that zdiff was negative for 35/54 QTL?

ld=fread('../stats/ld_decay/circos/ld_bundled_links_filtered.txt',data.table=F)
names(ld)=c('CHR_A','START_A','END_A','CHR_B','START_B','END_B','size')
chra=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(ld[x,]$CHR_A,'chr')[[1]][2]))
ld$CHR_A=chra
chrb=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(ld[x,]$CHR_B,'chr')[[1]][2]))
ld$CHR_B=chrb

## Are regions with highly correlated distal eQTL in high interchromosomal LD regions?
# Try for first three, then generally distal eQTL table? Not for these 3

peis=c("tkw_15-GRANEROS_2015_OPT-qTKW7_2","harvest_grain_moisture-GRANEROS_2015_OPT-qHGM7","male_flowering_d6-BLOIS_2017_OPT-qDTA3_1")

subdiff=truediff[truediff$pei %in% peis,]
subdiff=merge(subdiff,distal,by.x='pei',by.y='pheno_env_ID')
subdiff=merge(subdiff,genetable,by.x='gene',by.y='Gene_ID')

# Location of QTL in Ld region

env1=subdiff
env1=as.data.table(env1)
env2=as.data.table(ld)
setkey(env2,CHR_A,START_A,END_A)
compA=foverlaps(env1,env2,by.x=c('CHR','left_bound_bp','alt_right_bound_bp'),by.y=c('CHR_A','START_A','END_A'),nomatch=NULL)

compA=as.data.frame(compA,stringsAsFactors=F)
env1=compA
env1=as.data.table(env1)
env2=as.data.table(subdiff)
setkey(env2,CHROM,START,END)
compA2=foverlaps(env1,env2,by.x=c('CHR_B','START_B','END_B'),by.y=c('CHROM','START','END'),nomatch=NULL)
