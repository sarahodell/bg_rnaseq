#!/usr/bin/env Rscript

library('data.table')
library('dplyr')

base=c(7,7,8,6,7,6,7,7,6,8)
qtl=fread(sprintf('eqtl/results/%s_cis_eQTL_fkeep_hits.txt',time),data.table=F)
eqtl_haps=c()
for(i in 1:nrow(qtl)){
  chr=qtl[i,]$CHR
  snp=qtl[i,]$SNP
  pos=qtl[i,]$BP
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  gmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_gmap_c%.0f.csv',chr),data.table=F)
  for(h in seq(base[chr],16)){
    hfile=sprintf('../genotypes/probabilities/haplotype_probs/RefinedIBD_600K/dropped/bg%.0f_haplogroup%.0f_dropped_markers.rds',chr,h)
    if(file.exists(hfile)){
      g=readRDS(hfile)
      all_linked=unlist(lapply(g,function(x) x$linked))
      markers=unlist(lapply(g,function(x) x$marker))
      if(snp %in% all_linked){
        print(qtl[i,]$Gene)
        print("in")
        print(h)
        eqtl_haps=c(eqtl_haps,h)
      }else if(snp %in% markers){
        print(qtl[i,]$Gene)
        print("in")
        print(h)
        eqtl_haps=c(eqtl_haps,h)
      }else{
        print("not in")
        print(h)
      }
    }else{print("no file")}
	}
}


qtl$haps=eqtl_haps
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

for(i in 1:nrow(qtl)){
  chr=qtl[i,]$CHR
  snp=qtl[i,]$SNP
  pos=qtl[i,]$BP
  geneinfo=genetable[genetable$Gene_ID==qtl[i,]$Gene,]
  genestart=geneinfo$START
  geneend=geneinfo$END
  ibd=fread(sprintf('../ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  within=which(ibd$start<=pos & ibd$end>pos)
  row=ibd[within,]

  within2=which(ibd$start<=genestart & ibd$end>genestart)
  within3=which(ibd$start<=geneend & ibd$end>geneend)
  row=ibd[unique(c(within,within2,within3)),]
  print(qtl[i,]$Gene)
  print(row)
}

#Gene
#Outlier parent
#SNP IBD
#Gene IBD (Start and End)

#"Zm00001d031961"
# FV2_inra
#    chrom     start       end B73_inra A632_usa CO255_inra FV252_inra OH43_inra
#567     1 207691974 208050077        1        2          3          4         5
#568     1 208050077 208357564        1        2          3          4         5
#    A654_inra FV2_inra C103_inra EP1_inra D105_inra W117_inra B96 DK63 F492
#567         6        3         7        3         3         8   9    2   10
#568         6        7         8        3         3         9  10    2   11
#    ND245 VA85 n_haps
#567    11    7     11
#568    12    8     12

#[1] "Zm00001d032971"
# DK63
#    chrom     start       end B73_inra A632_usa CO255_inra FV252_inra OH43_inra
#704     1 244574866 244857001        1        2          3          2         4
#705     1 244857001 245057363        1        2          3          4         5
#    A654_inra FV2_inra C103_inra EP1_inra D105_inra W117_inra B96 DK63 F492
#704         5        6         7        3         8         9   2   10   11
#705         6        7         8        3         9        10   4   11   12
#    ND245 VA85 n_haps
#704    12    4     12
#705    13    5     13


#[1] "Zm00001d007378"
# F492
#    chrom     start       end B73_inra A632_usa CO255_inra FV252_inra OH43_inra
#794     2 229257688 230365709        1        2          3          4         5
#    A654_inra FV2_inra C103_inra EP1_inra D105_inra W117_inra B96 DK63 F492
#794         6        7         8        9        10        11  12   13   14
#    ND245 VA85 n_haps
#794    15   16     16


#[1] "Zm00001d051632"
# D105_inra
#    chrom     start       end B73_inra A632_usa CO255_inra FV252_inra OH43_inra
#357     4 164769860 164826383        1        1          2          3         4
#359     4 164932123 165362670        1        1          2          3         4
#    A654_inra FV2_inra C103_inra EP1_inra D105_inra W117_inra B96 DK63 F492
#357         5        6         7        2         8         9  10   11   12
#359         3        5         6        2         5         7   8    9   10
#    ND245 VA85 n_haps
#357    13   14     14
#359    11   12     12

#[1] "Zm00001d037902"
# W117_inra
#    chrom     start       end B73_inra A632_usa CO255_inra FV252_inra OH43_inra
#247     6 140470476 141555276        1        2          3          2         2
#    A654_inra FV2_inra C103_inra EP1_inra D105_inra W117_inra B96 DK63 F492
#247         4        5         5        6         6         7   8    2    9
#    ND245 VA85 n_haps
#247     2    5      9


#[1] "Zm00001d019592"
# None are different -all the ones that change are dropped
#    chrom    start      end B73_inra A632_usa CO255_inra FV252_inra OH43_inra
#164     7 41806410 44682779        1        1          2          3         4
#165     7 44682779 49829061        1        1          2          3         4
#    A654_inra FV2_inra C103_inra EP1_inra D105_inra W117_inra B96 DK63 F492
#164         5        6         7        8         9        10  11    4   12
#165         5        6         7        8         9        10  11   12   13
#    ND245 VA85 n_haps
#164    13   14     14
#165    14   15     15


#[1] "Zm00001d020799"
# A632
#    chrom     start       end B73_inra A632_usa CO255_inra FV252_inra OH43_inra
#257     7 132498742 133378157        1        2          3          4         5
#    A654_inra FV2_inra C103_inra EP1_inra D105_inra W117_inra B96 DK63 F492
#257         4        6         7        8         9         3  10    4   11
#    ND245 VA85 n_haps
#257    12    5     12

# How often is this the case for other tested genes?
time="WD_0712"
phenotypes=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
phenotypes=phenotypes[,c('V1',kept_genes)]
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

genes=intersect(unique(genetable$Gene_ID),kept_genes)
#phenotypes=phenotypes[,c('V1',genes)]

one_block=c()
nhaps=c()
for(i in genes){
  geneinfo=genetable[genetable$Gene_ID==i,]

  chr=geneinfo$CHR
  genestart=geneinfo$START
  geneend=geneinfo$END
  testsnps=readRDS(sprintf('eqtl/data/gene_focal_snps_c%.0f.rds',chr))

  snp=testsnps[[which(unlist(lapply(testsnps,function(x) x$gene==i)))]]$focal_snps
  if(length(snp)>1){
    snp=snp[1]
  }
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos

  ibd=fread(sprintf('../ibd_segments/refinedibd/600K/bg%.0f_refined_ibd_blocks.txt',chr),data.table=F)
  within=which(ibd$start<=pos & ibd$end>pos)
  row=ibd[within,]

  within2=which(ibd$start<=genestart & ibd$end>genestart)
  within3=which(ibd$start<=geneend & ibd$end>geneend)
  row=ibd[unique(c(within,within2,within3)),]
  if(nrow(row)==1){
    one_block=c(one_block,TRUE)
    nhaps=c(nhaps,row$n_haps)
  }
  else{
    one_block=c(one_block,FALSE)
    nhaps=c(nhaps,row$n_haps[1])
  }
  #print(qtl[i,]$Gene)
  #print(row)
}

data=data.frame(Gene=genes,sameIBD=one_block,n_haps=nhaps,stringsAsFactors=F)

# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
