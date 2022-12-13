#!/usr/bin/env Rscript

library('generics',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.1')
library('MegaLMM')
library('data.table')

time="WD_0712"
#n_ind=79
n_ind=246
n_k=100
n_genes=5000


#K=K[inter,inter]

#'pheno_MegaLMM_residuals_%s','pheno_MegaLMM_%s',
run_ids = c(sprintf('MegaLMM_%s',time),sprintf('MegaLMM_%s_residuals',time))
for(id in run_ids){
  if(grepl('residuals',id)){
    exp=fread(sprintf('eqtl/results/cis_eQTL_%s_all_vst_residuals.txt',time),data.table=F)
  }else{
    exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
  }
  rownames(exp)=exp$ID
  K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)
  inter=intersect(rownames(exp),rownames(K))
  n_k=nrow(exp)
  for(r in 1:3){
    run_id=paste0('MegaLMM/',id,'_',r)
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

# Average across runs

# Find correlated factors across runs
# MegaLMM gene counts
#run 1
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f',time,1)

r1=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
r1_probs=paste0('Factor',c(96))
r1=r1[!rownames(r1) %in% (r1_probs),]
r1=t(as.matrix(r1))

#run 2
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f',time,2)
r2=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
r2_probs=paste0('Factor',c(95,96))
r2=r2[!rownames(r2) %in% (r2_probs),]
r2=t(as.matrix(r2))

#run 3
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f',time,3)
r3=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
r3_probs=paste0('Factor',c(96))
r3=r3[!rownames(r3) %in% (r3_probs),]
r3=t(as.matrix(r3))



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
  count=c(count,sum(abs(data[i,c(4:6)])>=0.75,na.rm=T))
}

fwrite(data,sprintf('MegaLMM/MegaLMM_%s_factor_correlations.txt',time),row.names=F,quote=F,sep='\t')

# Average Lambda and F across mcmc runs

strong=which(count==3)
sdata=data[strong,]

#f_all_means=c()
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

fwrite(lambda_all_means,sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means.txt',time),row.names=T,quote=F,sep='\t')

# Average F values across runs
#run 1
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f',time,1)
r1=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
r1_probs=paste0('Factor',c(96))
r1=r1[,!names(r1) %in% (r1_probs)]


#run 2
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f',time,2)
r2=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
r2_probs=paste0('Factor',c(95,96))
r2=r2[,!names(r2) %in% (r2_probs)]

#run 3
run_id=sprintf('MegaLMM/MegaLMM_%s_%.0f',time,3)
r3=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
r3_probs=paste0('Factor',c(96))
r3=r3[,!names(r3) %in% (r3_probs)]

# Average across runs

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

fwrite(f_all_means,sprintf('MegaLMM/MegaLMM_%s_all_F_means.txt',time),row.names=T,quote=F,sep='\t')

######## MegaLMM using residuals ###########
#############################################

run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f',time,1)

r1=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
#r1_probs=paste0('Factor',c(96))
#r1=r1[!rownames(r1) %in% (r1_probs),]
r1=t(as.matrix(r1))

#run 2
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f',time,2)
r2=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
#r2_probs=paste0('Factor',c(95,96))
#r2=r2[!rownames(r2) %in% (r2_probs),]
r2=t(as.matrix(r2))

#run 3
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f',time,3)
r3=fread(sprintf('%s/Lambda_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
#r3_probs=paste0('Factor',c(96))
#r3=r3[!rownames(r3) %in% (r3_probs),]
r3=t(as.matrix(r3))



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
  count=c(count,sum(abs(data[i,c(4:6)])>=0.75,na.rm=T))
}

fwrite(data,sprintf('MegaLMM/MegaLMM_%s_residuals_factor_correlations.txt',time),row.names=F,quote=F,sep='\t')

# Average Lambda and F across mcmc runs
#count=apply(data,MARGIN=1,function(x) sum(x[c(6,8,9)]>=0.75,na.rm=T))

strong=which(count==3)
sdata=data[strong,]

#f_all_means=c()
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

fwrite(lambda_all_means,sprintf('MegaLMM/MegaLMM_%s_residuals_all_Lambda_means.txt',time),row.names=T,quote=F,sep='\t')

# Average F values across runs
#run 1
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f',time,1)
r1=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
#r1_probs=paste0('Factor',c(96))
#r1=r1[,!names(r1) %in% (r1_probs)]


#run 2
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f',time,2)
r2=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
#r2_probs=paste0('Factor',c(95,96))
#r2=r2[,!names(r2) %in% (r2_probs)]

#run 3
run_id=sprintf('MegaLMM/MegaLMM_%s_residuals_%.0f',time,3)
r3=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
#r3_probs=paste0('Factor',c(96))
#r3=r3[,!names(r3) %in% (r3_probs)]

# Average across runs

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

fwrite(f_all_means,sprintf('MegaLMM/MegaLMM_%s_residuals_all_F_means.txt',time),row.names=T,quote=F,sep='\t')






#### With phenotypes included ###############
#############################################

#
run_ids = c(sprintf('pheno_MegaLMM_residuals_%s',time),sprintf('pheno_MegaLMM_%s',time))
for(id in run_ids){
  if(grepl('residuals',id)){
    exp=fread(sprintf('eqtl/results/cis_eQTL_%s_all_vst_residuals.txt',time),data.table=F)
  }else{
    exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
  }
  rownames(exp)=exp$ID
  K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)
  inter=intersect(rownames(exp),rownames(K))
  n_k=nrow(exp)
  for(r in 1:3){
    run_id=paste0('MegaLMM/',id,'_',r)
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



##### With residuals########################
############################################
# run 1
r1=fread(sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_1/Lambda_means.txt',time),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
r1=t(as.matrix(r1))

# run 2
r2=fread(sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_2/Lambda_means.txt',time),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
r2=t(as.matrix(r2))

# run 3

r3=fread(sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_3/Lambda_means.txt',time),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
r3=t(as.matrix(r3))



# Correlations between factors
data=data.frame(matrix(ncol = 6, nrow = ncol(r1)))
names(data)=c('r1_factor_1','r2_factor_1','r3_factor_1','r1_r2_cor','r2_r3_cor','r1_r3_cor')
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
  count=c(count,sum(abs(data[i,c(4:6)])>=0.75,na.rm=T))
}

strong=which(count==3)

fwrite(data,sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_factor_correlations.txt',time),row.names=F,quote=F,sep='\t')

sdata=data[strong,]

# Average Lambda across runs

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
  r1f=row$r1_factor_1
  r2f=row$r2_factor_1
  r3f=row$r3_factor_1
  newdf=cbind(r1[,r1f],r2[,r2f],r3[,r3f])
  newdf=sapply(1:ncol(newdf),function(x) newdf[,x] * mult[x] )
  #newdf=newdf*mult
  lmean=as.data.frame(apply(newdf,MARGIN=1,mean),stringsAsFactors=F)
  lambda_all_means[,i]=lmean
}

phenotypes=rownames(lambda_all_means)[!grepl('Zm',rownames(lambda_all_means))]

fwrite(lambda_all_means,sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_all_Lambda_means.txt',time),row.names=T,quote=F,sep='\t')


# Average F values across runs
#run 1
run_id=sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_%.0f',time,1)
r1=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
#r1_probs=paste0('Factor',c(96))
#r1=r1[,!names(r1) %in% (r1_probs)]


#run 2
run_id=sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_%.0f',time,2)
r2=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
#r2_probs=paste0('Factor',c(95,96))
#r2=r2[,!names(r2) %in% (r2_probs)]

#run 3
run_id=sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_%.0f',time,3)
r3=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
#r3_probs=paste0('Factor',c(93,96))
#r3=r3[,!names(r3) %in% (r3_probs)]

# Average across runs

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

fwrite(f_all_means,sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_all_F_means.txt',time),row.names=T,quote=F,sep='\t')


# Identify Factor membership

#factors=names(lambda_all_means)
#factor_groups=vector("list",length=nrow(lambda_all_means))
#for(i in 1:length(factors)){#
#  factor_groups[[i]]$factor=factors[i]
#  factor_groups[[i]]$genes=c(NA)
#}

#for(f in 1:nrow(lambda_all_means)){
#  subl=lambda_all_means[f,,drop=F]
#  gene=rownames(subl)
#  var_exp=apply(subl,MARGIN=1,function(x) x**2)
#  tot_var=sum(var_exp)
#  prop_var=var_exp/tot_var
#  fkeep=names(subl[,which(prop_var>=0.2),drop=F])
#  for(k in fkeep){
#    x=unlist(unname(lapply(factor_groups2,function(x) which(x$factor==k))))
#    factor_groups[[x]]$genes=c(factor_groups[[x]]$genes,gene)
#
#  }
#}
#saveRDS(factor_groups)
#pheno_lambdas=lambda_all_means[phenotypes,]
#pfactors=which(unlist(unname(lapply(factor_groups,function(x) x$genes %in% phenotypes))))


##### With raw gene counts########################
############################################
# run 1
r1=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_1/Lambda_means.txt',time),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
r1=t(as.matrix(r1))

# run 2
r2=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_2/Lambda_means.txt',time),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
r2=t(as.matrix(r2))

# run 3

r3=fread(sprintf('MegaLMM/pheno_MegaLMM_%s_3/Lambda_means.txt',time),data.table=F)
rownames(r3)=r3$V1
r3_probs=paste0('Factor',c(93,96))
r3=r3[!rownames(r3) %in% (r3_probs),]
r3=r3[,-1]
r3=t(as.matrix(r3))

# Correlations between factors
data=data.frame(matrix(ncol = 6, nrow = ncol(r1)))
names(data)=c('r1_factor_1','r2_factor_1','r3_factor_1','r1_r2_cor','r2_r3_cor','r1_r3_cor')
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
  count=c(count,sum(abs(data[i,c(4:6)])>=0.75,na.rm=T))
}

strong=which(count==3)

fwrite(data,sprintf('MegaLMM/pheno_MegaLMM_%s_factor_correlations.txt',time),row.names=F,quote=F,sep='\t')

sdata=data[strong,]

# Average Lambda across runs

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
  r1f=row$r1_factor_1
  r2f=row$r2_factor_1
  r3f=row$r3_factor_1
  newdf=cbind(r1[,r1f],r2[,r2f],r3[,r3f])
  newdf=sapply(1:ncol(newdf),function(x) newdf[,x] * mult[x] )
  #newdf=newdf*mult
  lmean=as.data.frame(apply(newdf,MARGIN=1,mean),stringsAsFactors=F)
  lambda_all_means[,i]=lmean
}

phenotypes=rownames(lambda_all_means)[grepl('_',rownames(lambda_all_means))]

fwrite(lambda_all_means,sprintf('MegaLMM/pheno_MegaLMM_%s_all_Lambda_means.txt',time),row.names=T,quote=F,sep='\t')

# Average F values across runs
#run 1
run_id=sprintf('MegaLMM/pheno_MegaLMM_%s_%.0f',time,1)
r1=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]
#r1_probs=paste0('Factor',c(96))
#r1=r1[,!names(r1) %in% (r1_probs)]


#run 2
run_id=sprintf('MegaLMM/pheno_MegaLMM_%s_%.0f',time,2)
r2=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]
#r2_probs=paste0('Factor',c(95,96))
#r2=r2[,!names(r2) %in% (r2_probs)]

#run 3
run_id=sprintf('MegaLMM/pheno_MegaLMM_%s_%.0f',time,3)
r3=fread(sprintf('%s/F_means.txt',run_id),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]
r3_probs=paste0('Factor',c(93,96))
r3=r3[,!names(r3) %in% (r3_probs)]

# Average across runs

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

fwrite(f_all_means,sprintf('MegaLMM/pheno_MegaLMM_%s_all_F_means.txt',time),row.names=T,quote=F,sep='\t')



######## Overlap between pheno & non-pheno residuals

pheno=fread('MegaLMM/pheno_MegaLMM_residuals_WD_0712_all_Lambda_means.txt',data.table=F)
rownames(pheno)=pheno$V1
pheno=pheno[,-1]
resid=fread('MegaLMM/MegaLMM_WD_0712_residuals_all_Lambda_means.txt',data.table=F)
rownames(resid)=resid$V1
resid=resid[,-1]
inter=intersect(rownames(pheno),rownames(resid))
pheno=pheno[inter,]
resid=resid[inter,]
# Correlations between factors
data=data.frame(matrix(ncol = 3, nrow = ncol(pheno)))
names(data)=c('r1_factor_1','r2_factor_1','r1_r2_cor')
data$r1_factor_1=colnames(pheno)
#data$r2_factor_2=colnames(r2)
for(i in 1:ncol(pheno)){
  max_factor=NA
  max_cor=0
  for(j in 1:ncol(resid)){
    cur_cor=cor(pheno[,i],resid[,j],use="complete.obs")
    if(abs(cur_cor)>abs(max_cor)){
      max_cor=cur_cor
      max_factor=colnames(resid)[j]
    }
  }
  data$r1_r2_cor[i]=max_cor
  data$r2_factor_1[i]=max_factor
}

sdata=data[abs(data$r1_r2_cor)>0.75,]
sdata[sdata$r1_factor_1=="Factor3",]
sdata[sdata$r1_factor_1=="Factor10",]



fwrite(data,sprintf('MegaLMM/pheno_x_genecounts_MegaLMM_residuals_%s_factor_correlations.txt',time),row.names=F,quote=F,sep='\t')

##### Overlap between pheno residuals & pheno gene counts

pheno=fread('MegaLMM/pheno_MegaLMM_WD_0712_all_Lambda_means.txt',data.table=F)
rownames(pheno)=pheno$V1
pheno=pheno[,-1]
resid=fread('MegaLMM/pheno_MegaLMM_residuals_WD_0712_all_Lambda_means.txt',data.table=F)
rownames(resid)=resid$V1
resid=resid[,-1]
inter=intersect(rownames(pheno),rownames(resid))
pheno=pheno[inter,]
resid=resid[inter,]
# Correlations between factors
data=data.frame(matrix(ncol = 3, nrow = ncol(pheno)))
names(data)=c('r1_factor_1','r2_factor_1','r1_r2_cor')
data$r1_factor_1=colnames(pheno)
#data$r2_factor_2=colnames(r2)
for(i in 1:ncol(pheno)){
  max_factor=NA
  max_cor=0
  for(j in 1:ncol(resid)){
    cur_cor=cor(pheno[,i],resid[,j],use="complete.obs")
    if(abs(cur_cor)>abs(max_cor)){
      max_cor=cur_cor
      max_factor=colnames(resid)[j]
    }
  }
  data$r1_r2_cor[i]=max_cor
  data$r2_factor_1[i]=max_factor
}

sdata=data[abs(data$r1_r2_cor)>0.75,]
vsdata=data[abs(data$r1_r2_cor)>0.9,]
sdata[sdata$r1_factor_1=="Factor3",]
sdata[sdata$r1_factor_1=="Factor10",]



fwrite(data,sprintf('MegaLMM/pheno_resid_x_pheno_genecounts_MegaLMM_%s_factor_correlations.txt',time),row.names=F,quote=F,sep='\t')
