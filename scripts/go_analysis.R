#!/usr/bin/env Rscript
library('biomaRt')
library('AnnotationHub')
library('AnnotationForge')
library('org.Zm.eg.db')

library('data.table')
library('clusterProfiler')
library('rols')

time="WD_0720"

#glist=fread('Zea_mays_v4_to_v5_gene_list.txt',data.table=F)
factor_groups=readRDS(sprintf('eqtl/results/%s_factor_groupings.rds',time))
annotation=fread('annotations/GO/maize.B73.AGPv4.aggregate.withancestors.csv',data.table=F)

#testv5=glist[match(test,glist$v4_Gene_ID),]$v5_Gene_ID
#kept_test=test[which(!is.na(testv5))]
#testv5=testv5[!is.na(testv5)]
for(i in 1:length(factor_groups)){
  print(i)
  test=factor_groups[[i]]$genes
  results=enricher(test,TERM2GENE=annotation)
  factor_groups[[i]]$go=results
}
#Grab ontology descriptions
ol=Ontologies
go=ol[["go"]]
gotrms=terms(go)

for(i in 1:length(factor_groups)){
  sub=as.data.frame(factor_groups[[i]]$go,stringsAsFactors=F)
  if(nrow(sub)!=0){
    res=unlist(unname(sub$ID))
    a=gotrms[c(res)]
    labels=unlist(unname(termLabel(a)))
    sub$label=labels
    factor_groups[[i]]$go=sub
  }
}


saveRDS(factor_groups,sprintf('eqtl/results/%s_factor_groupings.rds',time))


#listMarts(host="plants.ensembl.org")
ensembl <- useEnsembl(host="plants.ensembl.org",biomart = "plants_mart",dataset="zmays_eg_gene")

#go <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id", "start_position", "end_position","go_id","name_1006"), mart=ensembl)

#listAttributes()
attributes=c("entrezgene_id","ensembl_gene_id","chromosome_name","start_position",
"end_position","go_id","name_1006","definition_1006")
gomap <- getBM(attributes=attributes,
               filters="ensembl_gene_id",
               values=testv5,
               mart=ensembl)
buildGOmap(gomap)

gomap=gomap[!is.na(gomap$entrezgene_id),]
#kept_testv5=testv5[which(!is.na(go))]

#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#datasets <- listDatasets(ensembl)
#searchDatasets(mart = ensembl, pattern = "zea")

all_go_mappings = fread('maize.B73.AGPv4.aggregate.gaf',data.table=F)
gomap = buildGOmap(data.frame(GO = all_go_mappings$term_accession,Gene = all_go_mappings$db_object_id))



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

testv5_entrez = gomap[match(testv5,gomap$ensembl_gene_id),]$entrezgene_id
testv5_entrez=testv5_entrez[!is.na(testv5_entrez)]
GOe <- enrichGO(gene=testv5_entrez,OrgDb='org.Zm.eg.db',keyType="ENTREZID",ont = "BP",pvalueCutoff = 0.05, qvalue = 0.1, readable = TRUE)

# make a plot
p1 <- plot(GOe, type = "bar", order = TRUE, showCategory = 15)
print(p1)

# write the results to a file
write.table(summary(GOe), file = "out_GOenrichment.txt", sep = "\t")
