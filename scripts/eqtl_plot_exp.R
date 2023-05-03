#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')

time="WD_0712"
ciseqtl=fread('eqtl/results/%s_cis_eQTL_fkeep_hits.txt',time),data.table=F)
norm=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
rownames(norm)=norm$V1
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

colorcodes=fread('../GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder

plot_list=list()
count=1
for(i in 1:nrow(ciseqtl)){
  gene=ciseqtl[i,]$Gene
  if(gene %in% names(norm)){
    chr=genetable[genetable$Gene_ID==gene,]$CHROM
    X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',chr))
    testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%.0f.rds',chr))
    inter=intersect(norm$V1,dimnames(X_list[[1]])[[1]])
    norm=norm[inter,]

    snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==gene)))]]$focal_snps
    X = do.call(cbind,lapply(X_list,function(x) x[inter,snp]))

    colnames(X) = founders
    rownames(X) = inter
    founder=unlist(unname(apply(X,MARGIN=1,function(x) colnames(X)[which.max(x)])))

    #print("True")
    ex=data.frame(ID=norm$V1,exp=norm[,gene],founder=founder,stringsAsFactors=F)
    ex=ex[ex$ID!="",]
    ex=ex[order(ex$exp),]
    rownames(ex)=seq(1,nrow(ex))
    ex$ID_f=factor(ex$ID,levels=c(unique(ex$ID)))
    ex$founder_f=factor(ex$founder,levels=c(unique(ex$founder)))
    #subex=subset(ex, ID %in% drop_ind)
    p=ggplot(aes(x=ID_f,y=exp),data=ex) + geom_point() +
    geom_point(aes(color=founder_f)) +
    scale_color_manual(values=colorcodes[levels(ex$founder_f),]$hex_color,labels=levels(ex$founder_f))+
     xlab('Sample') +
     ylab('Expression (log2CPM)') + geom_hline(yintercept=1)
    plot_list[[count]]=p
    count=count+1
  }
}

pdf('images/founder_eqtl_by_ind.pdf')
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()
