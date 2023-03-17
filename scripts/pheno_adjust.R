#!/usr/bin/env Rscript

library('xlsx')
library('data.table')
library('dplyr')

data=read.xlsx('phenotypes/Template_STPAUL_BALANCE_WD_GS17_transfert_Sarah.xlsx',sheetName='BALANCE_TRAITS_calculated')
names(data)=data[1,]
data1=fread('phenotypes/phenotypes_asi.csv',data.table=F)

data$Loc.Year.Treat=interaction(data$LOCATION,data$YEAR,data$TREATMENT,sep='_')

full_data=rbind(data1,data[,c('Loc.Year.Treat','LOCATION','')])

lct="STPAUL_2017_WD"
subdata1=data1[data1$Loc.Year.Treat==lct,]
subdata=data[data$Loc.Year.Treat==lct,]
subdata=subdata[,c('Loc.Year.Treat','LOCATION','YEAR','GENOTYPE','TREATMENT','sowing_date','female_flowering_d6','male_flowering_d6','total_plant_height','harvest_grain_moisture','grain_yield_15','tkw_15','asi_d6')]

metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
subdata$ID=metadata[match(subdata$GENOTYPE,metadata$genotype),]$dh_genotype
subdata1$Genotype_code=gsub('-','.',subdata1$Genotype_code)

subdata$female_flowering_d6=as.numeric(as.character(subdata$female_flowering_d6))
subdata$male_flowering_d6=as.numeric(as.character(subdata$male_flowering_d6))
subdata$total_plant_height=as.numeric(as.character(subdata$total_plant_height))
subdata$harvest_grain_moisture=as.numeric(as.character(subdata$harvest_grain_moisture))
subdata$grain_yield_15=as.numeric(as.character(subdata$grain_yield_15))
subdata$tkw_15=as.numeric(as.character(subdata$tkw_15))
subdata$asi_d6=as.numeric(as.character(subdata$asi_d6))

subdata=subdata[!is.na(subdata$ID),]
subdata=subdata[subdata$ID!="",]
rownames(subdata)=subdata$ID

inter=intersect(subdata$ID,subdata1$Genotype_code)

data_means=subdata %>% group_by(ID) %>% summarize(female_flowering_d6=mean(female_flowering_d6),male_flowering_d6=mean(male_flowering_d6),total_plant_height=mean(total_plant_height),harvest_grain_moisture=mean(harvest_grain_moisture),grain_yield_15=mean(grain_yield_15),tkw_15=mean(tkw_15),asi=mean(asi_d6))
data_means=as.data.frame(data_means,stringsAsFactors=F)
row.names(data_means)=data_means$ID
data_means=data_means[inter,]


data_means1=subdata1 %>% group_by(Genotype_code) %>% summarize(female_flowering_d6=mean(female_flowering_d6),male_flowering_d6=mean(male_flowering_d6),total_plant_height=mean(total_plant_height),harvest_grain_moisture=mean(harvest_grain_moisture),grain_yield_15=mean(grain_yield_15),tkw_15=mean(tkw_15),asi=mean(asi))
data_means1=as.data.frame(data_means1,stringsAsFactors=F)
row.names(data_means1)=data_means1$Genotype_code
data_means1=data_means1[inter,]

original=fread


#female_flowering_days
cor(data_means$female_flowering_d6,data_means1$female_flowering_d6,use="complete.obs")
#[1] 0.9237243

#male_flowering_days
# 0.9221134

# total plant height
cor(data_means$total_plant_height,data_means1$total_plant_height,use="complete.obs")
# 0.7877814

# harvest_grain_moisture
cor(data_means$harvest_grain_moisture,data_means1$harvest_grain_moisture,use="complete.obs")
# 0.9528051

# grain_yield_15
cor(data_means$grain_yield_15,data_means1$grain_yield_15,use="complete.obs")
# 0.9129779

# tkw_15
cor(data_means$tkw_15,data_means1$tkw_15,use="complete.obs")
# 0.9913168

# asi
cor(data_means$asi,data_means1$asi,use="complete.obs")
# 0.8625361

fwrite(data_means,'phenotypes/EXP_STPAUL_2017_WD_phenotypes.csv',row.names=F,quote=F,sep='\t')

data_means=fread('phenotypes/EXP_STPAUL_2017_WD_phenotypes.csv',data.table=F)
new_time="EXP_STPAUL_2017_WD"
names(data_means)=c('ID',paste0(new_time,'-',names(data_means)[-1]))
#exp=fread('eqtl/normalized/WD_0712_voom_normalized_gene_counts_formatted.txt',data.table=F)

times=c('WD_0712','WD_0718','WD_0720','WD_0727')
for(time in times){
  exp=fread(sprintf('eqtl/results/%s_vst_counts_x_phenotypes.txt',time),data.table=F)
  new_exp=merge(exp,data_means,by.x='ID',by.y='ID',all.x=TRUE)
  fwrite(new_exp,sprintf('eqtl/results/%s_vst_counts_x_phenotypes.txt',time),row.names=F,quote=F,sep='\t')
}


K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

phenos=c('female_flowering_d6','male_flowering_d6','total_plant_height','harvest_grain_moisture','grain_yield_15','tkw_15','asi')
phenotypes=fread('phenotypes/phenotypes_asi.csv',data.table=F)

data_means=fread('phenotypes/EXP_STPAUL_2017_WD_phenotypes.csv',data.table=F)
data_means$Loc.Year.Treat="EXP_STPAUL_2017_WD"

phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)

inter=intersect(unique(phenotypes$Genotype_code),rownames(K))
blups=data.frame(ID=inter,stringsAsFactors=F)
for(pheno in phenos){
  phenodf=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
  phenodf=phenodf[phenodf$Genotype_code %in% rownames(K),]
  names(phenodf)=c('ID','Loc.Year.Treat',pheno)
  phenodf=rbind(phenodf,data_means[,c('ID','Loc.Year.Treat',pheno)])
  data=data.frame(ID=phenodf$ID,ID2=phenodf$ID,Loc.Year.Treat=phenodf$Loc.Year.Treat,y=phenodf[,c(pheno)],stringsAsFactors=F)

  m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
  data_blup = as.data.frame(ranef(m1)$ID2,stringsAsFactors=F)
  data_blup$ID = rownames(data_blup)
  data_blup$y=data_blup$`(Intercept)`
  data_blup=data_blup[,c('ID','y')]
  data_blup$pheno=pheno
  blups[,paste0('BLUP','-',pheno)]=data_blup[match(inter,data_blup$ID),]$y
}
fwrite(blups,'phenotypes/phenotypes_BLUPs.csv',row.names=F,quote=F,sep='\t')


times=c('WD_0712','WD_0718','WD_0720','WD_0727')
for(time in times){
  exp=fread(sprintf('eqtl/results/%s_vst_counts_x_phenotypes.txt',time),data.table=F)
  new_exp=merge(exp,blups,by.x='ID',by.y='ID',all.x=TRUE)
  fwrite(new_exp,sprintf('eqtl/results/%s_vst_counts_x_phenotypes.txt',time),row.names=F,quote=F,sep='\t')
}

#vst <- function(x){
#  return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
#}

#vst_vect=unlist(data %>% group_by(Loc.Year.Treat) %>% group_map(~vst(.x$y)))
#vst_data=data[,c('ID','ID2','Loc.Year.Treat')]
#vst_data$y=vst_vect
