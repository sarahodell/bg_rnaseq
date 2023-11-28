#!/usr/bin/env Rscript

library('ggplot2')
library('data.table')
library('dplyr')

times=c('WD_0712',"WD_0718","WD_0720","WD_0727")
allgenes=c()
for(time in times){
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	genes=names(exp)[-1]
	allgenes=c(allgenes,genes)
}


trans=fread('eqtl/results/all_trans_fdr_SIs_FIXED.txt',data.table=F)

time1="WD_0727"
tmp=trans[trans$time==time1,]

hotspot=trans[trans$time=="WD_0727" & trans$SNP=="AX-91102858",]
prop_var=fread('MegaLMM/MegaLMM_WD_0727_prop_variance_FIXED.txt',data.table=F)
rownames(prop_var)=prop_var$V1
tgenes=unique(hotspot$gene)
prop_var=prop_var[tgenes,]
apply(prop_var[,-1], MARGIN=2,function(x) sum(x>0.1))

f4=prop_var[prop_var$Factor4>=0.1,c('V1','Factor4')]
local=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)
sub1=local[local$Trait %in% tgenes,]
#Factor2  Factor3  Factor4  Factor5  Factor6  Factor7  Factor8  Factor9 
#       8       21       42        4        6        1        1       16 
#Factor10 Factor11 Factor12 Factor13 
#      13       38       13        1 

# WD_0727 factor 4 go terms  - 223 floral development!


prop_var=fread('MegaLMM/MegaLMM_residuals_WD_0727_prop_variance_FIXED.txt',data.table=F)
rownames(prop_var)=prop_var$V1
tgenes=unique(hotspot$gene)
prop_var=prop_var[tgenes,]
apply(prop_var[,-1], MARGIN=2,function(x) sum(x>0.1))
#Factor1  Factor2  Factor3  Factor4  Factor5  Factor6  Factor7  Factor8 
#       4       15       24       15       15        5        7        7 
# Factor9 Factor10 Factor11 Factor12 Factor13 Factor14 Factor15 Factor16 
#       8       28       17        0        0        3        4        0 
#Factor17 Factor18 Factor19 
#       2        0        0 

bysnp=trans %>% group_by(time,SNP) %>% reframe(ngenes=length(unique(gene)))

# 18252 time-SNPs
summary(bysnp$ngenes)
# Just WD_0712
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   3.000   5.000   6.291   8.000  45.000

# ALL times
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   3.000   6.000   6.669   9.000  62.000

hotspots=bysnp[bysnp$ngenes>1,]
# 4699 SNPs, 17036 time-SNPs

fwrite(hotspots,'paper_figures/data/trans_hotspot_SNP_counts.txt',row.names=F,quote=F,sep='\t')

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

#8 AX-91102858 135530058 135656887  - 126829bp window

sub=trans[trans$SNP=="AX-91102858" & trans$time=="WD_0727",]
cand=fread('QTT/sig_candidate_genes.txt',data.table=F)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

sub=merge(sub,genetable,by.x="gene",by.y="Gene_ID")


# GO IDs of these genes


library('goseq')
library('rols')
time1="WD_0727"
genetable$LENGTH=genetable$END-genetable$START
genetable=genetable[order(genetable$CHROM,genetable$START),]
rownames(genetable)=seq(1,nrow(genetable))

annotation=fread('GO/maize.B73.AGPv4.aggregate.withancestors.csv',data.table=F)
annotation=merge(annotation,genetable,by.x="Gene",by.y="Gene_ID")
annotation=annotation[order(annotation$CHROM,annotation$START),]

log_inverse=function(x){
  	return(2^x)
}

genes=unique(annotation$Gene)
	
exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
rownames(exp)=exp$V1
exp=exp[,-1]


unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
avg_exp = apply(unlog,2,mean)
names(avg_exp)=names(exp)
avg_exp[avg_exp<1]=0
# re-log the input data
avg_logexp=log2(avg_exp)
avg_logexp[is.infinite(avg_logexp)]=0
	
#all_genes=intersect(unique(genetable$Gene_ID),unique(annotation$Gene))


pdf(sprintf('images/%s_nullp_plots.pdf',time1))

genes=names(avg_logexp)

test=unique(sub$gene)
deg=ifelse(genes %in% test,1,0)
names(deg)=genes
inter=intersect(test,genes[which(genes %in% test)])
test=inter
go=nullp(DEgenes=deg,bias.data=avg_logexp,plot.fit=T)
GO.samp=goseq(go,gene2cat=annotation[,c('Gene','GO')],use_genes_without_cat=TRUE,test.cats=c("GO:BP"))
GO.samp=as.data.frame(GO.samp,stringsAsFactors=F)
GO.samp$time=time1

dev.off()



#fwrite(all_go,sprintf('eqtl/results/chr8_trans_hotspot_GO_terms.txt',cutoff),row.names=F,quote=F,sep='\t')

GO.samp$adjusted_p=p.adjust(GO.samp$over_represented_pvalue,method="fdr")
enriched_GO=GO.samp[GO.samp$adjusted_p<=0.05,]
fwrite(enriched_GO,'eqtl/results/WD_0727_chr8_trans_hotspot_GO_terms.txt',row.names=F,quote=F,sep='\t')

# All chloroplast & thylakoid/plastid localized, no BP or MF terms

# Where are they located?
# Are these genes just in interchromosomal LD with this region?

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
eqtl$gene_time_SNP = paste0(eqtl$Trait,'-',eqtl$time,'-',eqtl$X_ID)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

eqtl[eqtl$X_ID=="AX-91102858",]
# Zm00001d010974 Zm00001d010975 Zm00001d010976
# Which of these has more genes with correlated effect sizes with trans-eQTL?
