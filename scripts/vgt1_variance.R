#!/usr/bin/env Rscript

library('data.table')

miteprob=fread('../GridLMM/mite_probabilities.txt')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra",
           "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra",
           "W117_inra","B96","DK63","F492","ND245","VA85")
has_mite=list("B73_inra"=F,"A632_usa"=T,"CO255_inra"=T,"FV252_inra"=T,"OH43_inra"=F,
           "A654_inra"=T,"FV2_inra"=T,"C103_inra"=T,"EP1_inra"=T,"D105_inra"=T,
           "W117_inra"=T,"B96"=F,"DK63"=T,"F492"=T,"ND245"=T,"VA85"=F)

founder_freq=(12*1 + 4*0)/16
founder_var=var(c(rep(1,12),rep(0,4)))
#0.75
magic_freq=sum(miteprob$final)/nrow(miteprob)
#0.77
results=fread('eqtl/cis/results/eQTL_WD_0712_c8_results2.txt',data.table=F)
rap27="Zm00001d010987"
results=results[results$Trait==rap27,]
chr="8"
founder_probs = readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
snp="AX-91102912"
X = do.call(cbind,lapply(founder_probs,function(x) x[,snp]))
time="WD_0712"
phenotype=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
K = fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
inter=intersect(rownames(K),phenotype$ID)
X=X[inter,]
miteprob=miteprob[miteprob$ID %in% inter,]
magic_freq=sum(miteprob$final)/nrow(miteprob)
#0.7341837 in subset of 79 individuals for timepoint WD_0712
betas=unlist(unname(results[1,c(6,10:24)]))
est= X %*% betas

obs_var=var(est)
#0.217

var(phenotype[phenotype$ID %in% inter,rap27])
# 3.866943

#2.7 fold increase in expression for those without the MITE (Salvi et al. 2007)
miteprob$effect=ifelse(miteprob$final==0,log2(2.7),0)
mite_var=var(miteprob$effect)
#0.3932132


# What is the variance of a cis eQTL we do detect

cis=fread('eqtl/results/WD_0718_cis_eQTL_hits.txt',data.table=F)
# What is a cis-eQTL that just passed significance
test=cis[which.min(cis$value),]
chr=as.character(test$CHR)
gene=test$Gene
snp=test$SNP
founder_probs = readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
X = do.call(cbind,lapply(founder_probs,function(x) x[,snp]))
results=fread('eqtl/cis/results/eQTL_WD_0712_c8_results2.txt',data.table=F)
results=results[results$Trait==gene,]
betas=unlist(unname(results[1,c(6,10:24)]))
est= X %*% betas
test_var=var(est)
#0.51

var(phenotype[phenotype$ID %in% inter,gene])
# 0.8884748

# Is the window around rap27 including 70kb upstream?
chr="8"
rap27="Zm00001d010987"
snp="AX-91102912"

testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%s.rds',chr))
founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
testsnp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==rap27)))]]$focal_snps
pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
pos=pmap[pmap$marker==snp,]$pos
founder_blocks[founder_blocks$start<=pos & founder_blocks$end > pos,]

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
genetable=genetable[genetable$CHROM==chr,]

rap27pos=genetable[genetable$Gene_ID==rap27,]
miteloc=rap27pos$START - 70000
# test snp is 20kb upstream of the MITE (within the block tested, so that's not the issue)
