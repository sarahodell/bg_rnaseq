#!/usr/bin/env Rscript

library('generics',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('MegaLMM',lib='/home/sodell/R/x86_64-conda-linux-gnu-library/4.2')
library('data.table')

time="WD_0720"
#n_ind=79
#n_ind=246
#n_k=100
#n_genes=5000
#n_k=83 #WD_0712
#n_k=143 #WD_0718
#n_k=220 #WD_0720
#n_k=194 #WD_0727

n_ks=list("WD_0712"=83,"WD_0718"=143,"WD_0720"=220,"WD_0727"=194)
#K=K[inter,inter]

#'pheno_MegaLMM_residuals_%s','pheno_MegaLMM_%s',
run_ids = c(sprintf('MegaLMM_%s',time))#,sprintf('MegaLMM_%s_residuals',time))
for(id in run_ids){
  if(grepl('residuals',id)){
    exp=fread(sprintf('eqtl/results/cis_eQTL_%s_all_vst_residuals.txt',time),data.table=F)
  }else{
    exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
  }
  rownames(exp)=exp$V1
  K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)
  inter=intersect(rownames(exp),rownames(K))
  n_k=n_ks[[time]]
  for(r in 1:3){
    run_id=paste0('MegaLMM/',id,'_',r,'_FIXED')
    #Reload model and current state
    MegaLMM_state=readRDS(sprintf('%s/MegaLMM_state_base.rds',run_id))
    MegaLMM_state$current_state=readRDS(sprintf('%s/current_state.rds',run_id))
    #Reload posterior
    MegaLMM_state$Posterior=reload_Posterior(MegaLMM_state,params=c("Lambda","F"))
    Lambda_mean=get_posterior_mean(MegaLMM_state$Posterior$Lambda)
    Lambda_mean=as.data.frame(Lambda_mean,stringsAsFactors=F)

    rownames(Lambda_mean)=paste0('Factor',seq(1,n_k))
    fwrite(Lambda_mean,paste0(run_id,'/Lambda_means.txt'),row.names=T,quote=F,sep='\t')

    F_mean=get_posterior_mean(MegaLMM_state$Posterior$F)
    F_mean=as.data.frame(F_mean,stringsAsFactors=F)
    colnames(F_mean)=paste0('Factor',seq(1,n_k))
    rownames(F_mean)=inter
    fwrite(F_mean,paste0(run_id,'/F_means.txt'),row.names=T,quote=F,sep='\t')
  }
}


#Get F_h2


#run_id=sprintf("MegaLMM/pheno_MegaLMM_%s_3",time)
#MegaLMM_state=readRDS(paste0(run_id,'/','MegaLMM_state_base.rds'))
#MegaLMM_state$current_state=readRDS(paste0(run_id,'/','current_state.rds'))
#MegaLMM_state$Posterior=reload_Posterior(MegaLMM_state,params=c("F_h2"))
#h2_mean=get_posterior_mean(MegaLMM_state$Posterior$F_h2)

#png(paste0(run_id,'/',sprintf('%s_posterior_F_h2_boxplot.png',time)))
#print(boxplot(MegaLMM_state$Posterior$F_h2[,1,]))
#dev.off()


# Average across runs

# Find correlated factors across runs
# MegaLMM gene counts
#run 1


# Average Lambda and F across mcmc runs
#WD_0712
r1_r=55
r2_r=57
r3_r=59

#WD_0718
r1_r=25
r2_r=24
r3_r=27

#WD_0720
r1_r=21
r2_r=18
r3_r=20

#WD_0727
r1_r=17
r2_r=16
r3_r=16

#f_all_means=c()

# Average F values across runs
#run 1
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f_FIXED',time,1)
r1=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
r1_probs=paste0('Factor',c(r1_r:n_k))
r1=r1[,!names(r1) %in% (r1_probs)]


#run 2
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f_FIXED',time,2)
r2=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
r2_probs=paste0('Factor',c(r2_r:n_k))
r2=r2[,!names(r2) %in% (r2_probs)]

#run 3
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f_FIXED',time,3)
r3=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
r3_probs=paste0('Factor',c(r3_r:n_k))
r3=r3[,!names(r3) %in% (r3_probs)]

# Average across runs
#nr=max(c(ncol(r1),ncol(r2),ncol(r3)))


data=data.frame(matrix(ncol = 5, nrow = ncol(r1)))
names(data)=c('r1_factor_1','r2_factor_1','r3_factor_1','r1_r2_cor','r2_r3_cor')
data$r1_factor_1=colnames(r1)
#data$r2_factor_2=colnames(r2)
for(i in 1:ncol(r1)){
  max_factor=NA
  max_cor=0
  for(j in 1:ncol(r2)){
    cur_cor=cor(r1[,i],r2[,j],use="complete.obs")
    if(abs(cur_cor)>abs(max_cor)){
      max_cor=cur_cor
      max_factor=colnames(r2)[j]
    }
  }
  data$r1_r2_cor[i]=max_cor
  data$r2_factor_1[i]=max_factor
}

for(i in 1:ncol(r1)){
  max_factor=NA
  max_cor=0
  for(j in 1:ncol(r3)){
    cur_cor=cor(r1[,i],r3[,j],use="complete.obs")
    if(abs(cur_cor)>abs(max_cor)){
      max_cor=cur_cor
      max_factor=colnames(r3)[j]
    }
  }
  data$r1_r3_cor[i]=max_cor
  data$r3_factor_1[i]=max_factor
}

r2_1_r3_cor=c()
for(i in 1:nrow(data)){
  r=cor(r2[,data$r2_factor_1[i]],r3[,data$r3_factor_1[i]])
  r2_1_r3_cor=c(r2_1_r3_cor,r)
}
data$r2_r3_cor=r2_1_r3_cor

count=c()
for(i in 1:nrow(data)){
  count=c(count,sum(abs(data[i,c(4:6)])>=0.90,na.rm=T))
}

fwrite(data,sprintf('MegaLMM/MegaLMM_%s_factor_correlations_FIXED.txt',time),row.names=F,quote=F,sep='\t')


data=fread(sprintf('MegaLMM/MegaLMM_%s_factor_correlations_FIXED.txt',time),data.table=F)
count=c()
for(i in 1:nrow(data)){
  count=c(count,sum(abs(data[i,c(4:6)])>=0.90,na.rm=T))
}
strong=which(count==3)
sdata=data[strong,]

factors=sdata$r1_factor_1
df=data.frame(time=time,chr=rep(seq(1,10),length(factors)),stringsAsFactors=F)
f=c()
for(i in factors){f=c(f,rep(i,10))}
df$factor=f
fwrite(df,sprintf('eqtl/trans/%s_chrom_factor_FIXED.txt',time),row.names=F,col.names=F,quote=F,sep=',')
fwrite(as.data.frame(factors),sprintf('eqtl/trans/%s_factors_FIXED.txt',time),row.names=F,quote=F,col.names=F,sep=',')


f_all_means=data.frame(matrix(ncol = nrow(sdata), nrow = nrow(r1)),stringsAsFactors=F)
rownames(f_all_means)=rownames(r1)
names(f_all_means)=sdata$r1_factor_1

for(i in 1:nrow(sdata)){
  row=sdata[i,]
  isneg=row[,c(4:6)]<0
  mult=c(1,1,1)
  if(sum(isneg==c(T,T,F))==3){
    mult=c(-1,1,1)
  }else if(sum(isneg==c(F,T,T))==3){
    mult=c(1,1,-1)
  }else if(sum(isneg==c(T,F,T))==3){
    mult=c(1,-1,1)
  }
  r1f=row$r1_factor_1
  r2f=row$r2_factor_1
  r3f=row$r3_factor_1
  newdf=cbind(r1[,r1f],r2[,r2f],r3[,r3f])
  newdf=sapply(1:ncol(newdf),function(x) newdf[,x] * mult[x] )
  fmean=as.data.frame(apply(newdf,MARGIN=1,mean),stringsAsFactors=F)
  f_all_means[,i]=fmean
}

fwrite(f_all_means,sprintf('MegaLMM/MegaLMM_%s_all_F_means_FIXED.txt',time),row.names=T,quote=F,sep='\t')

#Lambda

run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f_FIXED',time,1)

r1=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
r1_probs=paste0('Factor',c(r1_r:n_k))
r1=r1[!rownames(r1) %in% (r1_probs),]
r1=t(as.matrix(r1))

#run 2
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f_FIXED',time,2)
r2=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
r2_probs=paste0('Factor',c(r2_r:n_k))
r2=r2[!rownames(r2) %in% (r2_probs),]
r2=t(as.matrix(r2))

#run 3
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f_FIXED',time,3)
r3=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
r3_probs=paste0('Factor',c(r3_r:n_k))
r3=r3[!rownames(r3) %in% (r3_probs),]
r3=t(as.matrix(r3))



lambda_all_means=data.frame(matrix(ncol = nrow(sdata), nrow = nrow(r1)),stringsAsFactors=F)
rownames(lambda_all_means)=rownames(r1)
names(lambda_all_means)=sdata$r1_factor_1

for(i in 1:nrow(sdata)){
  row=sdata[i,]
  isneg=row[,c(4:6)]<0
  mult=c(1,1,1)
  if(sum(isneg==c(T,T,F))==3){
    mult=c(-1,1,1)
  }else if(sum(isneg==c(F,T,T))==3){
    mult=c(1,1,-1)
  }else if(sum(isneg==c(T,F,T))==3){
    mult=c(1,-1,1)
  }
  #print(mult)
  #print(row)
  #mult=isneg*-1
  r1f=row$r1_factor_1
  r2f=row$r2_factor_1
  r3f=row$r3_factor_1
  newdf=cbind(r1[,r1f],r2[,r2f],r3[,r3f])
  newdf=sapply(1:ncol(newdf),function(x) newdf[,x] * mult[x] )
  #newdf=newdf*mult
  lmean=as.data.frame(apply(newdf,MARGIN=1,mean),stringsAsFactors=F)
  lambda_all_means[,i]=lmean
}

fwrite(lambda_all_means,sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means_FIXED.txt',time),row.names=T,quote=F,sep='\t')


######## MegaLMM using residuals ###########
#############################################

time="WD_0712"

#'pheno_MegaLMM_residuals_%s','pheno_MegaLMM_%s',
run_ids = c(sprintf('MegaLMM_%s_residuals',time))#,sprintf('MegaLMM_%s_residuals',time))
for(id in run_ids){
  #if(grepl('residuals',id)){
  exp=fread(sprintf('eqtl/results/cis_eQTL_%s_all_residuals_FIXED.txt',time),data.table=F)
  #}else{
  #  exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time),data.table=F)
  #}
  rownames(exp)=exp$V1
  K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)
  inter=intersect(rownames(exp),rownames(K))
  n_k=60
  for(r in 1:3){
    run_id=paste0('MegaLMM/',id,'_',r,'_FIXED')
    #Reload model and current state
    MegaLMM_state=readRDS(sprintf('%s/MegaLMM_state_base.rds',run_id))
    MegaLMM_state$current_state=readRDS(sprintf('%s/current_state.rds',run_id))
    #Reload posterior
    MegaLMM_state$Posterior=reload_Posterior(MegaLMM_state,params=c("Lambda","F"))
    Lambda_mean=get_posterior_mean(MegaLMM_state$Posterior$Lambda)
    Lambda_mean=as.data.frame(Lambda_mean,stringsAsFactors=F)

    rownames(Lambda_mean)=paste0('Factor',seq(1,n_k))
    fwrite(Lambda_mean,paste0(run_id,'/Lambda_means.txt'),row.names=T,quote=F,sep='\t')

    F_mean=get_posterior_mean(MegaLMM_state$Posterior$F)
    F_mean=as.data.frame(F_mean,stringsAsFactors=F)
    colnames(F_mean)=paste0('Factor',seq(1,n_k))
    rownames(F_mean)=inter
    fwrite(F_mean,paste0(run_id,'/F_means.txt'),row.names=T,quote=F,sep='\t')
  }
}


#WD_0712
r1_r=53
r2_r=54
r3_r=49

#WD_0718
r1_r=39
r2_r=35
r3_r=33

#WD_0720
r1_r=21
r2_r=50
r3_r=21

#WD_0727
r1_r=22
r2_r=21
r3_r=23

#f_all_means=c()

# Average F values across runs
#run 1
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f_FIXED',time,1)
r1=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
r1_probs=paste0('Factor',c(r1_r:n_k))
r1=r1[,!names(r1) %in% (r1_probs)]


#run 2
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f_FIXED',time,2)
r2=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
r2_probs=paste0('Factor',c(r2_r:n_k))
r2=r2[,!names(r2) %in% (r2_probs)]

#run 3
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f_FIXED',time,3)
r3=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
r3_probs=paste0('Factor',c(r3_r:n_k))
r3=r3[,!names(r3) %in% (r3_probs)]

# Average across runs
#nr=max(c(ncol(r1),ncol(r2),ncol(r3)))


data=data.frame(matrix(ncol = 5, nrow = ncol(r1)))
names(data)=c('r1_factor_1','r2_factor_1','r3_factor_1','r1_r2_cor','r2_r3_cor')
data$r1_factor_1=colnames(r1)
#data$r2_factor_2=colnames(r2)
for(i in 1:ncol(r1)){
  max_factor=NA
  max_cor=0
  for(j in 1:ncol(r2)){
    cur_cor=cor(r1[,i],r2[,j],use="complete.obs")
    if(abs(cur_cor)>abs(max_cor)){
      max_cor=cur_cor
      max_factor=colnames(r2)[j]
    }
  }
  data$r1_r2_cor[i]=max_cor
  data$r2_factor_1[i]=max_factor
}

for(i in 1:ncol(r1)){
  max_factor=NA
  max_cor=0
  for(j in 1:ncol(r3)){
    cur_cor=cor(r1[,i],r3[,j],use="complete.obs")
    if(abs(cur_cor)>abs(max_cor)){
      max_cor=cur_cor
      max_factor=colnames(r3)[j]
    }
  }
  data$r1_r3_cor[i]=max_cor
  data$r3_factor_1[i]=max_factor
}

r2_1_r3_cor=c()
for(i in 1:nrow(data)){
  r=cor(r2[,data$r2_factor_1[i]],r3[,data$r3_factor_1[i]])
  r2_1_r3_cor=c(r2_1_r3_cor,r)
}
data$r2_r3_cor=r2_1_r3_cor

count=c()
for(i in 1:nrow(data)){
  count=c(count,sum(abs(data[i,c(4:6)])>=0.90,na.rm=T))
}

fwrite(data,sprintf('MegaLMM/MegaLMM_%s_residuals_factor_correlations_FIXED.txt',time),row.names=F,quote=F,sep='\t')


data=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_factor_correlations_FIXED.txt',time),data.table=F)
count=c()
for(i in 1:nrow(data)){
  count=c(count,sum(abs(data[i,c(4:6)])>=0.90,na.rm=T))
}
strong=which(count==3)
sdata=data[strong,]

factors=sdata$r1_factor_1
df=data.frame(time=time,chr=rep(seq(1,10),length(factors)),stringsAsFactors=F)
f=c()
for(i in factors){f=c(f,rep(i,10))}
df$factor=f
fwrite(df,sprintf('eqtl/trans/%s_residuals_chrom_factor_FIXED.txt',time),row.names=F,col.names=F,quote=F,sep=',')
fwrite(as.data.frame(factors),sprintf('eqtl/trans/%s_residuals_factors_FIXED.txt',time),row.names=F,quote=F,col.names=F,sep=',')


f_all_means=data.frame(matrix(ncol = nrow(sdata), nrow = nrow(r1)),stringsAsFactors=F)
rownames(f_all_means)=rownames(r1)
names(f_all_means)=sdata$r1_factor_1

for(i in 1:nrow(sdata)){
  row=sdata[i,]
  isneg=row[,c(4:6)]<0
  mult=c(1,1,1)
  if(sum(isneg==c(T,T,F))==3){
    mult=c(-1,1,1)
  }else if(sum(isneg==c(F,T,T))==3){
    mult=c(1,1,-1)
  }else if(sum(isneg==c(T,F,T))==3){
    mult=c(1,-1,1)
  }
  r1f=row$r1_factor_1
  r2f=row$r2_factor_1
  r3f=row$r3_factor_1
  newdf=cbind(r1[,r1f],r2[,r2f],r3[,r3f])
  newdf=sapply(1:ncol(newdf),function(x) newdf[,x] * mult[x] )
  fmean=as.data.frame(apply(newdf,MARGIN=1,mean),stringsAsFactors=F)
  f_all_means[,i]=fmean
}

fwrite(f_all_means,sprintf('MegaLMM/MegaLMM_%s_residuals_all_F_means_FIXED.txt',time),row.names=T,quote=F,sep='\t')

#Lambda

run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f_FIXED',time,1)

r1=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
r1_probs=paste0('Factor',c(r1_r:n_k))
r1=r1[!rownames(r1) %in% (r1_probs),]
r1=t(as.matrix(r1))

#run 2
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f_FIXED',time,2)
r2=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
r2_probs=paste0('Factor',c(r2_r:n_k))
r2=r2[!rownames(r2) %in% (r2_probs),]
r2=t(as.matrix(r2))

#run 3
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f_FIXED',time,3)
r3=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
r3_probs=paste0('Factor',c(r3_r:n_k))
r3=r3[!rownames(r3) %in% (r3_probs),]
r3=t(as.matrix(r3))



lambda_all_means=data.frame(matrix(ncol = nrow(sdata), nrow = nrow(r1)),stringsAsFactors=F)
rownames(lambda_all_means)=rownames(r1)
names(lambda_all_means)=sdata$r1_factor_1

for(i in 1:nrow(sdata)){
  row=sdata[i,]
  isneg=row[,c(4:6)]<0
  mult=c(1,1,1)
  if(sum(isneg==c(T,T,F))==3){
    mult=c(-1,1,1)
  }else if(sum(isneg==c(F,T,T))==3){
    mult=c(1,1,-1)
  }else if(sum(isneg==c(T,F,T))==3){
    mult=c(1,-1,1)
  }
  #print(mult)
  #print(row)
  #mult=isneg*-1
  r1f=row$r1_factor_1
  r2f=row$r2_factor_1
  r3f=row$r3_factor_1
  newdf=cbind(r1[,r1f],r2[,r2f],r3[,r3f])
  newdf=sapply(1:ncol(newdf),function(x) newdf[,x] * mult[x] )
  #newdf=newdf*mult
  lmean=as.data.frame(apply(newdf,MARGIN=1,mean),stringsAsFactors=F)
  lambda_all_means[,i]=lmean
}

fwrite(lambda_all_means,sprintf('MegaLMM/MegaLMM_%s_residuals_all_Lambda_means_FIXED.txt',time),row.names=T,quote=F,sep='\t')



#### F_h2

# What are the F_h2 estimates of these factors?
time="WD_0712"

data=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_factor_correlations_FIXED.txt',time),data.table=F)
count=c()
for(i in 1:nrow(data)){
  count=c(count,sum(abs(data[i,c(4:6)])>=0.90,na.rm=T))
}
strong=which(count==3)
sdata=data[strong,]

#'pheno_MegaLMM_residuals_%s','pheno_MegaLMM_%s',
id = c(sprintf('MegaLMM_%s_residuals',time))#,sprintf('MegaLMM_%s_residuals',time))
n_k=60
for(r in 1:3){
	run_id=paste0('MegaLMM/',id,'_',r,'_FIXED')
    #Reload model and current state
    MegaLMM_state=readRDS(sprintf('%s/MegaLMM_state_base.rds',run_id))
    MegaLMM_state$current_state=readRDS(sprintf('%s/current_state.rds',run_id))
    #Reload posterior
    MegaLMM_state$Posterior=reload_Posterior(MegaLMM_state,params=c("F_h2"))
    Fh2_mean=get_posterior_mean(MegaLMM_state$Posterior$F_h2)
    Fh2_mean=as.data.frame(t(Fh2_mean),stringsAsFactors=F)

    rownames(Fh2_mean)=paste0('Factor',seq(1,n_k))
    fwrite(Fh2_mean,paste0(run_id,'/Fh2_means.txt'),row.names=T,quote=F,sep='\t')
}


#WD_0712
r1_r=53
r2_r=54
r3_r=49

#WD_0718
r1_r=39
r2_r=35
r3_r=33

#WD_0720
r1_r=21
r2_r=50
r3_r=21

#WD_0727
r1_r=22
r2_r=21
r3_r=23

run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f_FIXED',time,1)

r1=fread(sprintf('%s/Fh2_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1,drop=F]
r1_probs=paste0('Factor',c(r1_r:n_k))
r1=r1[!rownames(r1) %in% (r1_probs),,drop=F]
#r1=t(as.matrix(r1))

#run 2
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f_FIXED',time,2)
r2=fread(sprintf('%s/Fh2_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1,drop=F]
r2_probs=paste0('Factor',c(r2_r:n_k))
r2=r2[!rownames(r2) %in% (r2_probs),,drop=F]
#r2=t(as.matrix(r2))

#run 3
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f_FIXED',time,3)
r3=fread(sprintf('%s/Fh2_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1,drop=F]
r3_probs=paste0('Factor',c(r3_r:n_k))
r3=r3[!rownames(r3) %in% (r3_probs),,drop=F]
#r3=t(as.matrix(r3))



Fh2_all_means=data.frame(matrix(nrow = nrow(sdata), ncol=1),stringsAsFactors=F)
names(Fh2_all_means)="F_h2"
rownames(Fh2_all_means)=sdata$r1_factor_1

for(i in 1:nrow(sdata)){
  row=sdata[i,]
  isneg=row[,c(4:6)]<0
  mult=c(1,1,1)
  #if(sum(isneg==c(T,T,F))==3){
  #  mult=c(-1,1,1)
  #}else if(sum(isneg==c(F,T,T))==3){
  #  mult=c(1,1,-1)
  #}else if(sum(isneg==c(T,F,T))==3){
  #  mult=c(1,-1,1)
  #}
  #print(mult)
  #print(row)
  #mult=isneg*-1
  r1f=row$r1_factor_1
  r2f=row$r2_factor_1
  r3f=row$r3_factor_1
  newdf=cbind(r1[r1f,],r2[r2f,],r3[r3f,])
  #newdf=sapply(1:ncol(newdf),function(x) newdf[,x] * mult[x] )
  #newdf=newdf*mult
  lmean=as.data.frame(apply(newdf,MARGIN=1,mean),stringsAsFactors=F)
  Fh2_all_means[i,]=lmean
}

fwrite(Fh2_all_means,sprintf('MegaLMM/MegaLMM_%s_residuals_all_Fh2_means_FIXED.txt',time),row.names=T,quote=F,sep='\t')

