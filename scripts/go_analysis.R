#!/usr/bin/env Rscript
library('goseq')
library('data.table')
library('rols')
library('dplyr')


genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
genetable$LENGTH=genetable$END-genetable$START
genetable=genetable[order(genetable$CHROM,genetable$START),]
rownames(genetable)=seq(1,nrow(genetable))

annotation=fread('GO/maize.B73.AGPv4.aggregate.withancestors.csv',data.table=F)
annotation=merge(annotation,genetable,by.x="Gene",by.y="Gene_ID")
annotation=annotation[order(annotation$CHROM,annotation$START),]

log_inverse=function(x){
  	return(2^x)
}

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
all_go=c()
for(time in times){
	data=fread(sprintf('MegaLMM/MegaLMM_%s_factor_correlations_FIXED.txt',time),data.table=F)
	count=c()
	for(i in 1:nrow(data)){
  		count=c(count,sum(abs(data[i,c(4:6)])>=0.9,na.rm=T))
	}
	strong=which(count==3)
	sdata=data[strong,]

	pheno_factors=sdata$r1_factor_1
	genes=unique(annotation$Gene)
	
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	rownames(exp)=exp$V1
	exp=exp[,-1]

	#geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
	#kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
	#exp=exp[,c(kept_genes)]
	unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
	avg_exp = apply(unlog,2,mean)
	names(avg_exp)=names(exp)
	avg_exp[avg_exp<1]=0
	# re-log the input data
	avg_logexp=log2(avg_exp)
	avg_logexp[is.infinite(avg_logexp)]=0
	
	all_genes=intersect(unique(genetable$Gene_ID),unique(annotation$Gene))
	#tp_genes=kept_genes

	genelength=genetable[match(all_genes,genetable$Gene_ID),]$LENGTH
	names(genelength)=all_genes

	#deg=ifelse(all_genes %in% tp_genes,1,0)
	#names(deg)=all_genes
	#go=nullp(DEgenes=deg,bias.data=genelength,plot.fit=F)
	#GO.full=goseq(go,gene2cat=annotation[,c('Gene','GO')],use_genes_without_cat=TRUE,test.cats=c("GO:BP"))
	#GO.full=as.data.frame(GO.full,stringsAsFactors=F)
	#GO.full$adjusted_p=p.adjust(GO.full$over_represented_pvalue,method="fdr")
	#enriched_GO.full=GO.full[GO.full$adjusted_p<0.05,]
	

	#full_cats=unique(enriched_GO.full$category)
	######### GOSeq ########

	
	pdf(sprintf('images/%s_nullp_plots.pdf',time))
	prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance_FIXED.txt',time),data.table=F)
	# How much variation is enough? 
	subvar=prop_var[,c('V1',pheno_factors)]
	cutoff=0.1

	nullp_list=vector("list",length(pheno_factors))
	genes=names(avg_exp)
	for(f in pheno_factors){
  		test=subvar[subvar[,f]>=cutoff,]$V1
  		if(length(test)>0){
  			deg=ifelse(genes %in% test,1,0)
  			names(deg)=genes
  			inter=intersect(test,genes[which(genes %in% test)])
  			test=inter
  			go=nullp(DEgenes=deg,bias.data=avg_logexp,plot.fit=T)
  			GO.samp=goseq(go,gene2cat=annotation[,c('Gene','GO')],use_genes_without_cat=TRUE,test.cats=c("GO:BP"))
  			GO.samp=as.data.frame(GO.samp,stringsAsFactors=F)
  			GO.samp$factor=f
  			GO.samp$time=time
  			#GO.samp$in_full=GO.samp$category %in% full_cats
  			all_go=rbind(all_go,GO.samp)
  		}else{
  			print("No genes loaded")
  			print(time)
  			print(f)
  		}	
  	}
  	dev.off()
}


all_go=as.data.frame(all_go,stringsAsFactors=F)
fwrite(all_go,sprintf('GO/MegaLMM_%.2f_GOSeq_all_FIXED.txt',cutoff),row.names=F,quote=F,sep='\t')

all_go$adjusted_p=p.adjust(all_go$over_represented_pvalue,method="fdr")
enriched_GO=all_go[all_go$adjusted_p<0.05,]
fwrite(enriched_GO,sprintf('GO/MegaLMM_%.2f_GOSeq_enriched_FIXED.txt',cutoff),row.names=F,quote=F,sep='\t')

########### Residuals ##############

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
all_go=c()
for(time in times){
	data=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_factor_correlations_FIXED.txt',time),data.table=F)
	count=c()
	for(i in 1:nrow(data)){
  		count=c(count,sum(abs(data[i,c(4:6)])>=0.9,na.rm=T))
	}
	strong=which(count==3)
	sdata=data[strong,]

	pheno_factors=sdata$r1_factor_1
	genes=unique(annotation$Gene)
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	rownames(exp)=exp$V1
	exp=exp[,-1]

	unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
	avg_exp = apply(unlog,2,mean)
	names(avg_exp)=names(exp)
	avg_exp[avg_exp<1]=0
	# re-log the input data
	avg_logexp=log2(avg_exp)
	avg_logexp[is.infinite(avg_logexp)]=0
	
	all_genes=intersect(unique(genetable$Gene_ID),unique(annotation$Gene))
	#tp_genes=kept_genes

	genelength=genetable[match(all_genes,genetable$Gene_ID),]$LENGTH
	names(genelength)=all_genes
	
	pdf(sprintf('images/%s_residuals_nullp_plots.pdf',time))
	prop_var=fread(sprintf('MegaLMM/MegaLMM_residuals_%s_prop_variance_FIXED.txt',time),data.table=F)
	# How much variation is enough? 
	subvar=prop_var[,c('V1',pheno_factors)]
	cutoff=0.1

	nullp_list=vector("list",length(pheno_factors))
	genes=names(avg_exp)
	for(f in pheno_factors){
  		test=subvar[subvar[,f]>=cutoff,]$V1
  		if(length(test)>0){
  			deg=ifelse(genes %in% test,1,0)
  			names(deg)=genes
  			inter=intersect(test,genes[which(genes %in% test)])
  			test=inter
  			go=nullp(DEgenes=deg,bias.data=avg_logexp,plot.fit=T)
  			GO.samp=goseq(go,gene2cat=annotation[,c('Gene','GO')],use_genes_without_cat=TRUE,test.cats=c("GO:BP"))
  			GO.samp=as.data.frame(GO.samp,stringsAsFactors=F)
  			GO.samp$factor=f
  			GO.samp$time=time
  			all_go=rbind(all_go,GO.samp)
  		}else{
  			print("No genes loaded")
  			print(time)
  			print(f)
  		}	
  	}
  	dev.off()
}


all_go=as.data.frame(all_go,stringsAsFactors=F)
fwrite(all_go,sprintf('GO/MegaLMM_residuals_%.2f_GOSeq_all_FIXED.txt',cutoff),row.names=F,quote=F,sep='\t')

all_go$adjusted_p=p.adjust(all_go$over_represented_pvalue,method="fdr")
enriched_GO=all_go[all_go$adjusted_p<0.05,]
fwrite(enriched_GO,sprintf('GO/MegaLMM_residuals_%.2f_GOSeq_enriched_FIXED.txt',cutoff),row.names=F,quote=F,sep='\t')









#testlist=data.frame(genes=test,stringsAsFactors=F)
#fwrite(testlist,sprintf('WD_0727_Factor2_%.1f_gene_list.txt',cutoff),row.names=F,col.names=F,quote=F,sep='\t')

#inter2=length(intersect(names(avg_exp),annotation$Gene))
#fulllist=data.frame(genes=genes,stringsAsFactors=F)
#fwrite(fulllist,'WD_0727_full_gene_list.txt',row.names=F,col.names=F,quote=F,sep='\t')

## Are any of the ~5000 genes we are comparing GO enriched for anything compared to the rest 
# of the genes in the genome?


enriched_GO %>% group_by(time,factor) %>% count()

#WD_0727

# A tibble: 14 × 2
# Groups:   factor [14]
#   factor       n
#   <chr>    <int>
# 1 Factor1     17
# 2 Factor13    21
# 3 Factor15     6
# 4 Factor19    36
# 5 Factor2      1
# 6 Factor21     9
# 7 Factor23    82
# 8 Factor3    118
# 9 Factor4      3
#10 Factor5     10
#11 Factor6      2
#12 Factor7     28
#13 Factor8    141
#14 Factor9      1


all_go$under_adjusted_p=p.adjust(all_go$under_represented_pvalue,method="fdr")
depleted_GO=all_go[all_go$under_adjusted_p<0.05,]











########### MITE-No MITE DEG GO

all_go=c()
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	degs=fread(sprintf('limma_results/%s_top_MITE_late_DEGs.txt',time),data.table=F)
	degs$time=time
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
	rownames(exp)=exp$V1
	exp=exp[,-1]
	
	unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
	avg_exp = apply(unlog,2,mean)
	names(avg_exp)=names(exp)
	avg_exp[avg_exp<1]=0
	# re-log the input data
	avg_logexp=log2(avg_exp)
	avg_logexp[is.infinite(avg_logexp)]=0
	
	all_genes=intersect(unique(genetable$Gene_ID),unique(annotation$Gene))
	#tp_genes=kept_genes

	genelength=genetable[match(all_genes,genetable$Gene_ID),]$LENGTH
	names(genelength)=all_genes

	
	pdf(sprintf('images/%s_MITE_nullp_plots.pdf',time))
	genes=names(avg_exp)
	
  	test=unique(degs$Gene_ID)
  	if(length(test)>1){
  		deg=ifelse(genes %in% test,1,0)
  		names(deg)=genes
  		inter=intersect(test,genes[which(genes %in% test)])
  		test=inter
  		go=nullp(DEgenes=deg,bias.data=avg_logexp,plot.fit=T)
  		GO.samp=goseq(go,gene2cat=annotation[,c('Gene','GO')],use_genes_without_cat=TRUE,test.cats=c("GO:BP"))
  		GO.samp=as.data.frame(GO.samp,stringsAsFactors=F)
  		GO.samp$time=time
  		#GO.samp$in_full=GO.samp$category %in% full_cats
  		all_go=rbind(all_go,GO.samp)
  	}
  	dev.off()
}


#all_go=as.data.frame(all_go,stringsAsFactors=F)
#fwrite(all_go,sprintf('GO/MegaLMM_%.2f_GOSeq_all_FIXED.txt',cutoff),row.names=F,quote=F,sep='\t')

all_go$adjusted_p=p.adjust(all_go$over_represented_pvalue,method="fdr")
enriched_GO=all_go[all_go$adjusted_p<0.05,]
fwrite(enriched_GO,sprintf('GO/MITE_%.2f_GOSeq_enriched_FIXED.txt',cutoff),row.names=F,quote=F,sep='\t')


alldegs=c()
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
for(time in times){
	#degs=fread(sprintf('limma_results/%s_top_MITE_DEGs.txt',time),data.table=F)
	degs=fread(sprintf('limma_results/%s_top_MITE_late_DEGs.txt',time),data.table=F)

	degs$time=time
	alldegs=rbind(alldegs,degs)
}

test=unique(alldegs$Gene_ID)
deg=ifelse(genes %in% test,1,0)
names(deg)=genes
inter=intersect(test,genes[which(genes %in% test)])
test=inter
genelength=genetable[match(genes,genetable$Gene_ID),]$LENGTH
names(genelength)=genes
go=nullp(DEgenes=deg,bias.data=genelength,plot.fit=T)
GO.samp=goseq(go,gene2cat=annotation[,c('Gene','GO')],use_genes_without_cat=TRUE,test.cats=c("GO:BP"))
GO.samp=as.data.frame(GO.samp,stringsAsFactors=F)

enriched_GO=GO.samp[GO.samp$over_represented_pvalue<=0.05,]
fwrite(enriched_GO,'GO/Early_Late_MITE_GOSeq_enriched_FIXED.txt',row.names=F,quote=F,sep='\t')
#fwrite(enriched_GO,'GO/MITE_NOMITE_GOSeq_enriched_FIXED.txt',row.names=F,quote=F,sep='\t')

#GO.samp$time=time
  		#GO.samp$in_full=GO.samp$category %in% full_cats
#all_go=rbind(all_go,GO.samp)
#all_go=all_go[all_go$ontology!="CC",]
#all_go=all_go[!is.na(all_go$ontology),]
#rownames(all_go)=seq(1,nrow(all_go))


#pdf('images/WD_0727_nullp_plots.pdf', onefile=TRUE)
#for (my.plot in my.plots) {
#    replayPlot(my.plot)
#}
#graphics.off()
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


#pheno_loc=c()
#for(f in pheno_factors){
#  indices=which(unlist(unname(lapply(factor_groups,function(x) x$factor==f))))
#  pheno_loc=c(pheno_loc,indices)
#}
#for(i in pheno_loc){
#  print(i)
#  test=factor_groups[[i]]$genes
#  factor=factor_groups[[i]]$factor
#  #results=enricher(test,TERM2GENE=annotation)
#  factor_groups[[i]]$go=all_go[all_go$factor==factor,]
#}



#Grab ontology descriptions
#ol=rols::Ontologies()
#go=ol[["po"]]
#gotrms=rols::terms(go)

#library('GO.db')
#for(go in enriched.GO[1:10,]$category){
#  print(GOTERM[[go]])
#  cat("--------------------------------------\n")
#}

#missing_go=c()
#res=unlist(unname(enriched_GO$category))
#labels=c()
#for(i in 1:nrow(enriched_GO)){#
#	row=enriched_GO[i,]
#	r=row$category
#	a=row$term
#	if(is.na(a)){
#		missing_go=c(missing_go,r)
#		a=Term(r)
#	}
#	labels=c(labels,a)
#}
#enriched_GO$term=labels


#saveRDS(factor_groups,sprintf('MegaLMM/MegaLMM_%s_factor_groups.rds',time))
#saveRDS(factor_groups,sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_factor_groups.rds',time))

#all_go_df=c()
#for(i in pheno_loc){
#  if(nrow(sub)!=0){
#    factor=factor_groups[[i]]$factor
#    sub=as.data.frame(factor_groups[[i]]$go,stringsAsFactors=F)
    #sub=sub[sub$over_represented_pvalue<=0.05,]

#    sub$factor=factor
    #all_go_df=rbind(all_go_df,sub)
#    sub=sub[,c('category','over_represented_pvalue')]
#    fwrite(sub,sprintf('MegaLMM/GO/MegaLMM_%s_%s_GO_terms.txt',time,factor),row.names=F,quote=F,sep='\t')
#  }
#}


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
