#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
factor=as.character(args[[2]])
run=as.character(args[[3]])
cores=as.numeric(args[[4]])

library('data.table')
library('MegaLMM')
library('GridLMM')
library('dplyr')
library('parallel')
library('MASS')
library('stringr')

library('ggplot2')
library('reshape2')
library('tibble')
library('tidyr')
#library('cowplot')
# Libraries ====
library('readr')
library('ggrepel')
library('RColorBrewer')
library('qqman')



#time="WD_0712"
#factor=paste0("Factor",f)
#n="1000"
#run="1"
run_id=sprintf('MegaLMM/MegaLMM_%s_%s',time,run)

#MegaLMM_state=readRDS(sprintf('%s/MegaLMM_state_base.rds',run_id))
#MegaLMM_state$current_state=readRDS(sprintf('%s/current_state.rds',run_id))
#MegaLMM_state$Posterior=reload_Posterior(MegaLMM_state,'F')
#F_mean=get_posterior_mean(MegaLMM_state,'F')


#exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
#rownames(exp)=exp$V1
#K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
#rownames(K)=K[,1]
#rownames(K)=gsub("-",".",rownames(K))
#K=as.matrix(K[,-1])
#colnames(K)=rownames(K)

#inter=intersect(rownames(exp),rownames(K))
#n_k=nrow(exp)-1

#F_mean=as.data.frame(F_mean,stringsAsFactors=F)
#colnames(F_mean)=paste0('Factor',seq(1,n_k))
#rownames(F_mean)=inter

# drop factors that didn't load
#F_mean=F_mean[,1:50]

#fwrite(F_mean,paste0(run_id,'/F_means.txt'),row.names=T,quote=F,sep='\t')

#### GridLMM

phenotype=fread(sprintf('%s/F_means.txt',run_id),data.table=F)

metadata=fread('metadata/BG_completed_sample_list.txt',data.table=F)
metadata=metadata[metadata$experiment==time,]
plate=metadata[match(phenotype$V1,metadata$dh_genotype),]$plate

drop=list('Factor19'="EB.09S.H.00429",'Factor38'="EB.09S.H.00424",
'Factor48'="EB.09S.H.00494",'Factor49'="EB.09S.H.00402")
#factor="Factor7"

for(chr in 1:10){
  #data=data.frame(ID=phenotype$V1,ID2=phenotype$V1,y=phenotype[,factor],stringsAsFactors=F)
  #data=data[!is.na(data$y),]
  data=data.frame(ID=phenotype$V1,ID2=phenotype$V1,plate=plate,y=phenotype[,factor],stringsAsFactors=F)
  data=data[!is.na(data$y),]
  data$plate=as.factor(data$plate)
  K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%.0f.txt',chr),data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)

  X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%.0f_filtered_genotype_probs.rds',chr))
  if(factor %in% names(drop)){
    drop_ind=drop[factor]
    X_list=lapply(X_list,function(x) x[rownames(x)!=drop_ind,])
  }

  inds=rownames(X_list[[1]])
  inter=intersect(inds,data$ID)

  null_model = GridLMM_ML(y~1+plate+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

  nmarkers=dim(X_list[[1]])[2]
  frep2=sapply(seq(1,nmarkers),function(i) lapply(X_list,function(j) sum(j[,i]>0.75)))
  founders=names(X_list)
  fkeep=apply(frep2,MARGIN=2,function(x) x>2)
  markers=dimnames(X_list[[1]])[[2]]
  colnames(fkeep)=markers
  colnames(frep2)=markers
  fgroups=unique(colSums(fkeep))

  all_gwas=data.frame(matrix(ncol=27,nrow=0))
  names(all_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','plate',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')


  for(g in fgroups){
    subm=colnames(fkeep[,colSums(fkeep)==g])
    subfkeep=fkeep[,subm]
    X_list_sub=lapply(X_list,function(x) x[inter,subm])
    if(g==16){
      h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
      names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
      h2_start
      V_setup=null_model$setup
      Y=as.matrix(data$y)
      X_cov=null_model$lmod$X
      X_list_null=NULL

      gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
      gwas$Trait=factor
      names(gwas)[c(6,8:22)]=founders
      names(gwas)[7]='plate'
      gwas=gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','plate',founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
      all_gwas=rbind(all_gwas,gwas)
    }else{
      pattern=apply(subfkeep,MARGIN=2,function(x) str_flatten(c(unlist(founders[x])),'-'))
  #pattern=
     fdf=data.frame(marker=subm,fpattern=pattern,stringsAsFactors=F)
      fpatterns=unique(fdf$fpattern)
      for(i in fpatterns){
        subm2=fdf[fdf$fpattern==i,]$marker
        subf=subfkeep[,subm2,drop=F]
      #m=marker[i]
        fk=founders[subf[,1]]
        nfk=founders[!subf[,1]]
        X_list_sub2=X_list_sub[ - which(names(X_list_sub) %in% nfk)]
        X_list_sub2=lapply(X_list_sub2,function(x) x[,subm2,drop=F])
        h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
        names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
        h2_start
        V_setup=null_model$setup
        Y=as.matrix(data$y)
        X_cov=null_model$lmod$X

        X_list_null=NULL

        gwas=run_GridLMM_GWAS(Y,X_cov,X_list_sub2[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
        gwas$Trait=factor
        if(!("B73_inra" %in% fk)){
          end=8+length(fk)-2
          #betas=gwas[,c(6,8:end)]
          #betas=unlist(unname(betas))
          #betas[-1]=betas[1]+betas[-1]
          #names(gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',fkeep,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
          #names(betas)=fkeep
          new_gwas=gwas[,1:5]
          #new_gwas=cbind(new_gwas,betas[1,"B73_inra"])
          ncol1=data.frame(matrix(ncol=1,nrow=nrow(gwas)))
          names(ncol1)='B73_inra'
          new_gwas=cbind(new_gwas,ncol1)
          new_gwas=cbind(new_gwas,gwas[,7])
          new_gwas=cbind(new_gwas,gwas[,c(6,8:end)])
          names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra','plate',fk)
          #fdrop=founders[!(founders %in% fk)]
          fdrop=nfk[nfk!="B73_inra"]
          nacol=data.frame(matrix(ncol=length(fdrop),nrow=nrow(gwas)))
          names(nacol)=fdrop
          new_gwas=cbind(new_gwas,nacol)
          new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
          new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',"B73_inra","plate",founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
          #new_X=cbind(X,data[,c('PC1','PC2','PC3')])
          #new_X=new_X[,c(fkeep[1],'PC1','PC2','PC3',fkeep[-1])]
        }else{
          end=8+length(fk)-2
          #betas=gwas[,c(6,8:end)]
          #betas=unlist(unname(betas))
          #betas[-1]=betas[1]+betas[-1]
          #names(gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML',fkeep,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')
          #names(betas)=fk
          new_gwas=gwas[,1:5]
          new_gwas=cbind(new_gwas,gwas[,6])
          #names(new_gwas)[6]="B73_inra"
          new_gwas=cbind(new_gwas,gwas[,7])
          #names(new_gwas)[7]="plate"
          new_gwas=cbind(new_gwas,gwas[,8:end])
          names(new_gwas)=c('Trait','X_ID','s2','ML_logLik','ID.ML','B73_inra',"plate",fk[-1])
          #fdrop=founders[!(founders %in% fk)]
          nacol=data.frame(matrix(ncol=length(nfk),nrow=nrow(gwas)))
          names(nacol)=nfk
          new_gwas=cbind(new_gwas,nacol)
          new_gwas=cbind(new_gwas,gwas[,(end+1):ncol(gwas)])
          new_gwas=new_gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',"B73_inra","plate",founders[-1],'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
          #new_X=cbind(X,data[,c('PC1','PC2','PC3')])
          #new_X=new_X[,c('B73_inra','PC1','PC2','PC3',fkeep[-1])]
        }

        #names(gwas)[6:(6+length(fk)-1)]=fk
        #gwas[,nfk]=NA
        #gwas=gwas[,c('Trait','X_ID','s2','ML_logLik','ID.ML',founders,'n_steps','Df_X','ML_Reduced_logLik','Reduced_Df_X','p_value_ML')]
        #all_gwas=rbind(all_gwas,gwas)
        all_gwas=rbind(all_gwas,new_gwas)

      }
    }
  }
  all_gwas=as.data.frame(all_gwas,stringsAsFactors=F)
  fwrite(all_gwas,sprintf('%s/test/test_c%.0f_%s_results.txt',run_id,chr,factor),row.names=F,quote=F,sep='\t')
}


#### Plot Manhattan plot of results
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


mypalette = c("#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FF61C9")
greypalette=gray.colors(5)

df=c()
for(c in 1:10){
  #d=fread(sprintf('eqtl/trans/results/%s_c%s_pheno_residuals_factor_trans_eQTL.txt',time,c))
  d=fread(sprintf('%s/test/test_c%.0f_%s_results.txt',run_id,c,factor))
  pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
  d$CHR=c
  d$BP=pmap[match(d$X_ID,pmap$marker),]$pos
  df=rbind(df,d)
}
df=df[!is.na(df$p_value_ML),]
df=df[!is.infinite(df$ML_logLik),]

threshold=-log10(0.05/15)
print(threshold)

p_adjusted=p.adjust(df$p_value_ML,method='fdr')
df$value=-log10(p_adjusted)


df=df[,c('Trait','X_ID','p_value_ML','CHR','BP','value')]
names(df)=c('Factor','SNP','p_value_ML','CHR','BP','value')

subdf=df[,c('Factor','CHR','BP','SNP','value')]
subdf=subdf[order(subdf$CHR,subdf$BP),]

gg.manhattan2 <- function(df, threshold, col, ylims,bounds){
  df.tmp <- df %>%

    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=as.numeric(BP+tot)) %>%
    gather(key, value, -BP,-SNP,-CHR,-BPcum,-tot,-Factor)
  df.tmp$sig=df.tmp$value > threshold
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  df.tmp=df.tmp[!is.na(df.tmp$value),]
  rownames(df.tmp)=seq(1,nrow(df.tmp))
  ggplot(df.tmp, aes(x=BPcum, y=value)) +
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22 )) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    geom_hline(yintercept = threshold,linetype="dashed") +
    geom_point(data=subset(df.tmp, sig==T), color="coral2", size=5) +
    theme_classic() +
    theme(
      text = element_text(size=30),
      plot.title = element_text(hjust = 0.5),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + guides(color="none")
}

ymax=max(round(max(subdf$value))+1,threshold+1)

title=sprintf("eQTL for %s at timepoint %s",factor,time)
a2<-gg.manhattan2(subdf,threshold,
             col=greypalette,
             ylims=c(0,ymax)) + labs(caption = title)

sigs=subdf[subdf$value>=threshold,]
#print(dim(sigs))

if(dim(sigs)[1]!=0){
  #png(sprintf('eqtl/trans/images/%s_pheno_residuals_trans_eQTL_manhattan_%s.png',factor,time),width=2000,height=1500)
  png(sprintf('%s/images/test_%s_manhattan.png',run_id,factor),width=2000,height=1500)
  print(a2)
  dev.off()

  png(sprintf('%s/images/test_%s_qq_plot.png',run_id,factor))
  #png(sprintf('eqtl/trans/images/%s_%s_pheno_residuals_trans_qqplot.png',time,factor))
  qqman::qq(df$p_value_ML)
  dev.off()
  #fwrite(sigs2,sprintf('eqtl/results/%s_pheno_trans_%s_eQTL_fkeep_hits.txt',factor,time),row.names=F,quote=F,sep='\t')
  #fwrite(sigs2,sprintf('eqtl/results/%s_pheno_residuals_trans_%s_eQTL_hits.txt',factor,time),row.names=F,quote=F,sep='\t')
}
