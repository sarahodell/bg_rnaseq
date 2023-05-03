#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
#library('biomaRt')
#library('AnnotationHub')
#library('AnnotationForge')
#library('org.Zm.eg.db')
library('goseq')
library('data.table')
#library('clusterProfiler')
library('rols')
#library('TissueEnrich')
library('dplyr')

#time="WD_0727"

#factor_groups=readRDS(sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_factor_groups.rds',time))


genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
genetable$LENGTH=genetable$END-genetable$START
genetable=genetable[order(genetable$CHROM,genetable$START),]
rownames(genetable)=seq(1,nrow(genetable))

annotation=fread('annotations/GO/maize.B73.AGPv4.aggregate.withancestors.csv',data.table=F)
annotation=merge(annotation,genetable,by.x="Gene",by.y="Gene_ID")
annotation=annotation[order(annotation$CHROM,annotation$START),]

#inter1=intersect(annotation$Gene,genetable$Gene_ID)
#annotation=annotation[annotation$Gene %in% inter1,]
#rownames(genetable)=genetable$Gene_ID
#genetable=genetable[inter1,]
#rownames(annotation)=seq(1,nrow(annotation))
#rownames(genetable)=seq(1,nrow(genetable))

#pheno_df=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_sig_factors.txt',time),data.table=F)
#pheno_df=fread(sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_sig_factors.txt',time),data.table=F)

#'MegaLMM/pheno_MegaLMM_residuals_%s_sig_factors.txt',time)
#testv5=glist[match(test,glist$v4_Gene_ID),]$v5_Gene_ID
#kept_test=test[which(!is.na(testv5))]
#testv5=testv5[!is.na(testv5)]
pheno_factors=c('Factor14')
genes=unique(annotation$Gene)
genelength=genetable[match(genes,genetable$Gene_ID),]$LENGTH
names(genelength)=genes

exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
rownames(exp)=exp$ID
exp=exp[,-1]

geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
exp=exp[,c(kept_genes)]

log_inverse=function(x){
  return(2^x)
}
unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
avg_exp = apply(unlog,2,mean)
names(avg_exp)=names(exp)
avg_exp[avg_exp<1]=0
# re-log the input data
avg_logexp=log2(avg_exp)
avg_logexp[is.infinite(avg_logexp)]=0

#avg_exp2=apply(exp,2,mean)

#genes=genes[1:10]
#go_list = sapply(genes,function(x) NULL)
#go_list = sapply(seq(1,length(genes)),function(x) go_list[[x]]=annotation[annotation$Gene==names(go_list)[x],]$GO)
#for(i in 1:length(genes)){
#  go_list[[i]]=annotation[annotation$Gene==genes[i],]$GO
#  names(go_list[[i]])=genes[i]
#}
#for(i in 1:length(go_list)){
#  names(go_list[[i]])=rep(genes[i],length(go_list[[i]]))
#}
#saveRDS(go_list,'annotations/GO/maize_annotations_list.rds')
#go_list=readRDS('annotations/GO/maize_annotations_list.rds')
####GOSeq


##### TissueEnrich ######

#expressionData<-fread('annotations/expression_atlas/MATURE_LEAF_TISSUE_leaf_8.fpkm_tracking',data.table=F)
#rownames(expressionData)=expressionData$V1
#expressionData=expressionData[,1:2]
#expressionData$V2= log2(expressionData$V2)
#expressionData[is.infinite(expressionData$V2),]$V2=0
#names(expressionData)=c('Gene','Mature_Leaf_8')
#seed=fread('annotations/expression_atlas/B73_12DAP_Whole_seed.fpkm_tracking',data.table=F)
#rownames(seed)=seed$V1
#seed=seed[,1:2]
#names(seed)=c('Gene','DAP12_Whole_seed')


#root=fread('annotations/expression_atlas/CrownRoot_Node4_V7.fpkm_tracking',data.table=F)
#rownames(root)=seed$root
#root=root[,1:2]
#seed$V2=log2(seed$V2)
#seed[is.infinite(seed$V2),]$V2=0
#names(root)=c('Gene','CrownRoot_Node4_V7')

#expressionData=merge(expressionData,seed,by="Gene",all=T)
#expressionData=merge(expressionData,root,by="Gene",all=T)
#rownames(expressionData)=expressionData$Gene
#expressionData[is.na(expressionData$Mature_Leaf_8),]$Mature_Leaf_8=0
#expressionData=expressionData[grepl('Zm',expressionData$Gene),]
#expressionData=expressionData[,-1]

#se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = rownames(expressionData),colData = names(expressionData))
#output<-teGeneRetrieval(se,maxNumberOfTissues=2)
#output_df=as.data.frame(assay(output),stringsAsFactors=F)
#leaf_enriched=output_df[output_df$Tissue=="Mature_Leaf_8" & output_df$Group=="Tissue-Enriched",]
#fwrite(leaf_enriched,'annotations/Mature_Leaf_8_enriched_genes.txt',row.names=F,quote=F,sep='\t')

#f=pheno_factors[1]
#indices=which(unlist(unname(lapply(factor_groups,function(x) x$factor==f))))
#test=factor_groups[[indices]]$genes
#test=test[test %in% rownames(expressionData)]
#gs<-GeneSet(geneIds=test)
#output2<-teEnrichmentCustom(gs,output,multiHypoCorrection=F)
#Plotting the P-Values
#enrichmentOutput<-setNames(data.frame(assay(output2[[1]]),row.names = rowData(output2[[1]])[,1]),colData(output2[[1]])[,1])
#enrichmentOutput$Tissue<-row.names(enrichmentOutput)

prop_var=fread('MegaLMM/MegaLMM_WD_0727_prop_variance.txt',data.table=F)
# How much variation is enough? 
subvar=prop_var[,c('V1','Factor2')]
cutoff=0.9
test=subvar[subvar$Factor2>=cutoff,]$V1

testlist=data.frame(genes=test,stringsAsFactors=F)
fwrite(testlist,sprintf('WD_0727_Factor2_%.1f_gene_list.txt',cutoff),row.names=F,col.names=F,quote=F,sep='\t')



inter2=length(intersect(names(avg_exp),annotation$Gene))
#genelength=genetable[match(genes,genetable$Gene_ID),]$LENGTH
#names(genelength)=genes

fulllist=data.frame(genes=genes,stringsAsFactors=F)
fwrite(fulllist,'WD_0727_full_gene_list.txt',row.names=F,col.names=F,quote=F,sep='\t')


ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
# Enrichment of flowering time genes?
cutoff=0.1
test=subvar[subvar$Factor2>=cutoff,]$V1
gtable=genetable[genetable$Gene_ID %in% test,]


find_nearest_snp=function(row){
    index=which.min(abs(row$START-pmap[pmap$chr==row$CHROM,]$pos))
    return(pmap[index,]$marker)
}

nearest_snps=sapply(seq(1,nrow(gtable)),function(x) find_nearest_snp(gtable[x,]))
gtable$SNP=nearest_snps

snplist=gtable[,'SNP',drop=F]
fwrite(snplist,'WD_0727_Factor2_snplist.txt',row.names=F,col.names=F,sep='\t',quote=F)

true_nft=length(intersect(ft_genelist$V4,test))

ngenes=length(test)
allgenes=length(genes)
nfts=c()
for(i in 1:1000){
	draw=sample(seq(1,allgenes),ngenes)
	nft=intersect(ft_genelist$V4,genes[draw])
	nfts=c(nfts,length(nft))
}

# 0.1 cutoff 60 genes
#summary(nfts)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  39.00   53.75   58.00   57.36   61.00   74.00 
# quantile(nfts,0.95)
#95% 
# 67 

# 0.5 cutoff 5 genes
summary(nfts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   7.000   9.000   9.166  11.000  23.000 
#quantile(nfts,0.05)
#5% 
# 4

######### GOSeq ########
pheno_factors=c('Factor16')
factor_groups=readRDS(sprintf('MegaLMM/MegaLMM_%s_factor_groups.rds',time))
#pheno_factors=unique(pheno_df$factor)
genes=names(avg_exp)

nullp_list=vector("list",length(pheno_factors))



all_go=c()
enriched_go=c()
for(f in pheno_factors){
  indices=which(unlist(unname(lapply(factor_groups,function(x) x$factor==f))))
  test=factor_groups[[indices]]$genes
  deg=ifelse(genes %in% test,1,0)
  names(deg)=genes
  inter=intersect(test,genes[which(genes %in% test)])
  test=inter
  #genelength2=genelength[genes]
  go=nullp(DEgenes=deg,bias.data=avg_exp,plot.fit=T)
  GO.samp=goseq(go,gene2cat=annotation[,c('Gene','GO')],use_genes_without_cat=TRUE,test.cats=c("GO:BP"))
  #GO.samp=goseq(go,gene2cat=annotation[,c('Gene','GO')],method="Sampling",repcnt=10000,use_genes_without_cat=TRUE)
  GO.samp=as.data.frame(GO.samp,stringsAsFactors=F)
  GO.samp$factor=f
  #enriched.GO=GO.samp[p.adjust(GO.samp$over_represented_pvalue,method="BH")<.05,]
  all_go=rbind(all_go,GO.samp)
  GO.samp$adjusted_p=p.adjust(GO.samp$over_represented_pvalue,method="fdr")
  enriched_go=rbind(enriched_go,GO.samp[GO.samp$adjusted_p<(0.05/length(pheno_factors)),])
}
all_go=as.data.frame(all_go,stringsAsFactors=F)
fwrite(all_go,sprintf('MegaLMM/GO/MegaLMM_%s_GOSeq_all.txt',time),row.names=F,quote=F,sep='\t')
#all_go=all_go[all_go$ontology!="CC",]
#all_go=all_go[!is.na(all_go$ontology),]
#rownames(all_go)=seq(1,nrow(all_go))
#enriched.GO=all_go[all_go$adjusted_p<(0.05/length(pheno_factors)),]


fwrite(enriched_go,sprintf('MegaLMM/GO/MegaLMM_%s_GOSeq_enriched.txt',time),row.names=F,quote=F,sep='\t')

pdf('images/WD_0727_nullp_plots.pdf', onefile=TRUE)
for (my.plot in my.plots) {
    replayPlot(my.plot)
}
graphics.off()
#### Generally leaf-tissue enriched GO terms
#deg=ifelse(genes %in% leaf_enriched$Gene,1,0)
#names(deg)=genes
#inter=intersect(test,genes[which(genes %in% test)])
#test=inter
#go=nullp(DEgenes=deg,bias.data=genelength)
#GO.samp=goseq(go,gene2cat=annotation[,c('Gene','GO')],method="Hypergeometric")
#leaf_enriched.GO=GO.samp[p.adjust(GO.samp$over_represented_pvalue,method="BH")<.05,]

#enriched.GO$leaf_enriched = enriched.GO$category %in% leaf_enriched.GO$category




#f=pheno_factors[1]
#indices=which(unlist(unname(lapply(factor_groups,function(x) x$factor==f))))
#test=factor_groups[[indices]]$genes

#deg=ifelse(genes %in% test,1,0)
#names(deg)=genes
#inter=intersect(test,genes[which(genes %in% test)])
#test=inter
#genelength=genetable[match(genes,genetable$Gene_ID),]$LENGTH
#names(genelength)=genes
#go=nullp(DEgenes=deg,bias.data=genelength)
#GO.samp=goseq(go,gene2cat=annotation[,c('Gene','GO')],method="Sampling",repcnt=1000)

#enriched.GO=GO.samp[p.adjust(GO.samp$over_represented_pvalue,method="BH")<.05,]

####
#all_go=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_GOSeq_enriched.txt',time),data.table=F)


pheno_loc=c()
for(f in pheno_factors){
  indices=which(unlist(unname(lapply(factor_groups,function(x) x$factor==f))))
  pheno_loc=c(pheno_loc,indices)
}
for(i in pheno_loc){
  print(i)
  test=factor_groups[[i]]$genes
  factor=factor_groups[[i]]$factor
  #results=enricher(test,TERM2GENE=annotation)
  factor_groups[[i]]$go=all_go[all_go$factor==factor,]
}



#Grab ontology descriptions
ol=rols::Ontologies()
go=ol[["po"]]
gotrms=rols::terms(go)

library('GO.db')
#for(go in enriched.GO[1:10,]$category){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n")
#}

missing_go=c()
for(i in pheno_loc){
  sub=as.data.frame(factor_groups[[i]]$go,stringsAsFactors=F)
  if(nrow(sub)!=0){
    res=unlist(unname(sub$category))
    labels=c()
    for(r in res){
      a=Term(r)
      if(!is.null(a)){
        label=a
      }else{
        label=NA
        missing_go=c(missing_go,res)
      }
      labels=c(labels,label)
    }
    #a=gotrms[c(res)]
    sub$label=labels
    factor_groups[[i]]$go=sub
  }
}


saveRDS(factor_groups,sprintf('MegaLMM/MegaLMM_%s_factor_groups.rds',time))
#saveRDS(factor_groups,sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_factor_groups.rds',time))

all_go_df=c()
for(i in pheno_loc){
  if(nrow(sub)!=0){
    factor=factor_groups[[i]]$factor
    sub=as.data.frame(factor_groups[[i]]$go,stringsAsFactors=F)
    #sub=sub[sub$over_represented_pvalue<=0.05,]

    sub$factor=factor
    #all_go_df=rbind(all_go_df,sub)
    sub=sub[,c('category','over_represented_pvalue')]
    fwrite(sub,sprintf('MegaLMM/GO/MegaLMM_%s_%s_GO_terms.txt',time,factor),row.names=F,quote=F,sep='\t')
  }
}


#enriched.GO=GO[p.adjust(GO$over_represented_pvalue,method="BH")<.05,]

#fwrite(all_go_df,sprintf('MegaLMM/pheno_MegaLMM_%s_GO_terms.txt',time),row.names=F,quote=F,sep='\t')

#group_by(A, B) %>%
#             filter(value == max(value)) %>%
#             arrange(A,B,C)
#maxfactor=pheno_df %>% group_by(phenotype) %>% filter(prop_var==max(prop_var)) %>% arrange(factor,phenotype,prop_var)
#saveRDS(factor_groups,sprintf('eqtl/results/%s_factor_groupings.rds',time))




#listMarts(host="plants.ensembl.org")
#ensembl <- useEnsembl(host="plants.ensembl.org",biomart = "plants_mart",dataset="zmays_eg_gene")

#go <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id", "start_position", "end_position","go_id","name_1006"), mart=ensembl)

#listAttributes()
#attributes=c("entrezgene_id","ensembl_gene_id","chromosome_name","start_position",
#"end_position","go_id","name_1006","definition_1006")
#gomap <- getBM(attributes=attributes,
#               filters="ensembl_gene_id",
#               values=testv5,
#               mart=ensembl)
#buildGOmap(gomap)

#gomap=gomap[!is.na(gomap$entrezgene_id),]
#kept_testv5=testv5[which(!is.na(go))]

#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#datasets <- listDatasets(ensembl)
#searchDatasets(mart = ensembl, pattern = "zea")

#all_go_mappings = fread('maize.B73.AGPv4.aggregate.gaf',data.table=F)
#gomap = buildGOmap(data.frame(GO = all_go_mappings$term_accession,Gene = all_go_mappings$db_object_id))



#ah = AnnotationHub()
#dm <- query(ah, c("Zea mays"))
#db=dm[["AH72174"]]
#file.copy(AnnotationHub::cache(dm["AH72174"]), ./"org.Zm.eg.sqlite")
#seed <- new("AnnDbPkgSeed", Package = "org.Zm.eg.db", Version = "0.0.1",Author = "Sarah Odell", Maintainer = "Sarah Odell <sgodell@ucdavis.edu>", PkgTemplate = "NOSCHEMA.DB", AnnObjPrefix = "org.Zm.eg", organism = "Zea mays", species = "Zea mays", biocViews = "annotation", manufacturerUrl = "none", manufacturer = "none", chipName = "none")
#makeAnnDbPkg(seed, "org.Zm.eg.sqlite")

#install.packages("org.Zm.eg.db/", type = "source", repos = NULL)


#annotations_orgDb <-AnnotationDbi::select(db, # database
#                                     keys = ,  # data to use for retrieval
#                                     columns = c("SYMBOL", "ENTREZID","GENENAME"), # information to retreive for given data
#                                     keytype = "ENSEMBL")

#testv5_entrez = gomap[match(testv5,gomap$ensembl_gene_id),]$entrezgene_id
#testv5_entrez=testv5_entrez[!is.na(testv5_entrez)]
#GOe <- enrichGO(gene=testv5_entrez,OrgDb='org.Zm.eg.db',keyType="ENTREZID",ont = "BP",pvalueCutoff = 0.05, qvalue = 0.1, readable = TRUE)

# make a plot
#p1 <- plot(GOe, type = "bar", order = TRUE, showCategory = 15)
#print(p1)

# write the results to a file
#write.table(summary(GOe), file = "out_GOenrichment.txt", sep = "\t")
