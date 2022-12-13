#!/usr/bin/env Rscript

library('data.table')

times=c("WD_0718","WD_0720","WD_0727","WD_0712")
metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)

for(time in times){
  resids=fread(sprintf('eqtl/results/cis_eQTL_%s_all_vst_residuals.txt',time),data.table=F)
  norm=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts.txt',time),data.table=F)
  samples=colnames(norm)[-1]
  genos=metadata[match(samples,metadata$sample_name),]$dh_genotype
  rownames(norm)=norm$V1
  norm=norm[,-1]
  names(norm)=c(genos)
  norm=as.data.frame(t(norm),stringsAsFactors=F)
  genes=names(norm)
  norm$ID=rownames(norm)
  norm=norm[,c('ID',genes)]
  fwrite(norm,sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),row.names=F,quote=F,sep='\t')
}


phenos=c("female_flowering_d6","male_flowering_d6","total_plant_height","harvest_grain_moisture",
"grain_yield_15","tkw_15",'asi')
envs=c("BLOIS_2014_OPT","BLOIS_2017_OPT","GRANEROS_2015_OPT","NERAC_2016_WD",
"STPAUL_2017_WD","SZEGED_2017_OPT")

phenotypes=fread('../GridLMM/phenotypes_asi.csv',data.table=F)
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)


for(time in times){
  norm=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)

  phenotypes2=phenotypes[phenotypes$Genotype_code %in% norm$ID,]
  phenotypes2=phenotypes2[,c('Genotype_code',"Loc.Year.Treat",phenos)]
  newpheno=data.frame(ID=norm$ID,stringsAsFactors=F)
  for(e in envs){
    df=phenotypes2[phenotypes2$Loc.Year.Treat==e,]
    names(df)=c('ID','Loc.Year.Treat',paste0(e,'-',phenos))
    df=df[match(newpheno$ID,df$ID),]
    rownames(df)=seq(1,nrow(df))
    newpheno=cbind(newpheno,df[3:9])
  }
  norm=cbind(norm,newpheno[,-1])
  fwrite(norm,sprintf('eqtl/results/%s_vst_counts_x_phenotypes.txt',time),row.names=F,quote=F,sep='\t')

}
