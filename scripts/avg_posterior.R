#!/usr/bin/env Rscript

library('generics',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.1')
library('MegaLMM')
library('data.table')

time="WD_0712"
n_ind=79
n_k=100

for(r in 1:3){

  run_id=sprintf('vst_MegaLMM_%s_%.0f',time,r)

  exp=fread(sprintf('eqtl/results/cis_eQTL_%s_all_vst_residuals.txt',time),data.table=F)
  key=exp[,c('ID'),drop=F]
  rownames(exp)=exp$ID
  Y = exp[,-c(1,2)]
  data = key

  K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)

  inter=intersect(data$ID,rownames(K))
  K=K[inter,inter]
  run_parameters = MegaLMM_control(
    max_NA_groups = 3,
    scale_Y = FALSE,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
    h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
    h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
    burn = 10000,  # number of burn in samples before saving posterior samples
    thin = 40,
    K = 100 # number of factors
  )

  #Set the prior hyperparameters of the BSFG model
  priors = MegaLMM_priors(
    tot_Y_var = list(V = 0.5,   nu = 5),      # Prior variance of trait residuals after accounting for fixed effects and factors
    tot_F_var = list(V = 18/20, nu = 20),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
    Lambda_prior = list(
      sampler = sample_Lambda_prec_horseshoe, # function that implements the horseshoe-based Lambda prior described in Runcie et al 2020. See code to see requirements for this function.
      prop_0 = 0.1,    # prior guess at the number of non-zero loadings in the first and most important factor
      delta = list(shape = 3, scale = 1),    # parameters of the gamma distribution giving the expected change in proportion of non-zero loadings in each consecutive factor
      delta_iterations_factor = 100   # parameter that affects mixing of the MCMC sampler. This value is generally fine.
    ),
    h2_priors_resids_fun = function(h2s,n)  1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
    h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
  )


  MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                                ~(1|ID),  # RHS of base model for factors and residuals. Fixed effects defined here only apply to the factor residuals.
                                data = data,         # the data.frame with information for constructing the model matrices
                                relmat = list(ID = K), # covariance matrices for the random effects. If not provided, assume uncorrelated
                                run_parameters=run_parameters,
                                run_ID = run_id
  )

  #Reload posterior
  MegaLMM_state$Posterior=reload_Posterior(MegaLMM_state,params=c("Lambda","F"))
  Lambda_mean=get_posterior_mean(MegaLMM_state$Posterior$Lambda)
  Lambda_mean=as.data.frame(Lambda_mean,stringsAsFactors=F)

  rownames(Lambda_mean)=paste0('Factor',seq(1,100))
  fwrite(Lambda_mean,paste0(run_id,'/Lambda_means.txt'),row.names=T,quote=F,sep='\t')

  F_mean=get_posterior_mean(MegaLMM_state$Posterior$F)
  F_mean=as.data.frame(F_mean,stringsAsFactors=F)
  colnames(F_mean)=paste0('Factor',seq(1,100))
  rownames(F_mean)=inter
  fwrite(F_mean,paste0(run_id,'/F_means.txt'),row.names=T,quote=F,sep='\t')
}


#n_genes=29478
# Lambda average across samples
#for(r in 1:3){
#  all_l=array(numeric(),c(250,n_k,n_genes))
#  for(i in seq(25,250,25)){
#    l=readRDS(sprintf('MegaLMM_%s_%.0f/Posterior/Lambda_%.0f.rds',time,r,i))
#    h=i-24
    #print(h:i)
    #print("next")
#    all_l[h:i,,]=l
#  }

  # calculate mean and variance of factor loadings across mcmc samples

#  l_means=sapply(seq(1:n_k),function(i) sapply(seq(1:n_genes),function(j) mean(all_l[,i,j])))
#  l_var=sapply(seq(1:n_k),function(i) sapply(seq(1:n_genes),function(j) var(all_l[,i,j])))


#  data=fread(sprintf('eqtl/results/cis_eQTL_%s_all_residuals.txt',time),data.table=F)
#  data=colnames(data)[-1]

#  l_means=as.data.frame(l_means,stringsAsFactors=F)
#  l_var=as.data.frame(l_var,stringsAsFactors=F)
#  rownames(l_means)=data
#  colnames(l_means)=paste0('Factor',seq(1,100))
#  rownames(l_var)=data
#  colnames(l_var)=paste0('Factor',seq(1,100))


#  fwrite(l_means,sprintf('MegaLMM_%s_%.0f/Lambda_means.txt',time,r),row.names=T,quote=F,sep='\t')
#  fwrite(l_var,sprintf('MegaLMM_%s_%.0f/Lambda_var.txt',time,r),row.names=T,quote=F,sep='\t')
#}


r1=fread(sprintf('vst_MegaLMM_%s_1/Lambda_means.txt',time),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]

#r1_probs=paste0('Factor',c(76,77,78,79,80,83,92,93,99))
#vst
r1_probs=paste0('Factor',c(37,38,39,40,42,43,44,45,46,47,48,49,50,51,52,
53,54,55,56,57,58,59,60,63,64,65,67,69,70,71,72,74,75,76,77,78,79,80,81,82,83,
84,85,86))

r1=r1[!rownames(r1) %in% (r1_probs),]
#rownames(r1)=seq(1,nrow(r1))
r1=t(as.matrix(r1))

#Do this with Lambda! Not U_F


r2=fread(sprintf('vst_MegaLMM_%s_2/Lambda_means.txt',time),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]

#r2_probs=paste0('Factor',c(67,68,72,76,77,81,82,89,99))
#vst
r2_probs=paste0('Factor',c(35,36,37,38,39,40,41,42,43,44,45,46,47,48,51,52,53,54,
55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,80,81,82,83,84,
85,86,87,88,89,90,91))


r2=r2[!rownames(r2) %in% (r2_probs),]
#rownames(r2)=seq(1,nrow(r2))
r2=t(as.matrix(r2))

r3=fread(sprintf('vst_MegaLMM_%s_3/Lambda_means.txt',time),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]

#r3_probs=paste0('Factor',c(69,77,81,85,86,97,88,89,90,91,99))

#vst
r3_probs=paste0('Factor',c(26,27,29,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,
48,49,51,52,53,54,55,56,57,60,61,62,64,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
80,81,82,83,84,85,90))

r3=r3[!rownames(r3) %in% (r3_probs),]
#rownames(r3)=seq(1,nrow(r3))
r3=t(as.matrix(r3))

#colnames(r1)=paste0('Factor',colnames(r1))
#colnames(r2)=paste0('Factor',colnames(r2))
#colnames(r3)=paste0('Factor',colnames(r3))


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

#for(i in 1:ncol(r2)){
  #r2_names=data$r2_factor_1
#  max_factor=NA
#  max_cor=0
#  for(j in 1:ncol(r3)){
#    cur_cor=cor(r2[,i],r3[,j],use="complete.obs")
#   if(abs(cur_cor)>abs(max_cor)){
#      max_cor=cur_cor
#      max_factor=colnames(r3)[j]
#    }
#  }
#  data$r2_r3_cor[i]=max_cor
#  data$r3_factor_2[i]=max_factor
#}

r2_1_r3_cor=c()
for(i in 1:nrow(data)){
  r=cor(r2[,data$r2_factor_1[i]],r3[,data$r3_factor_1[i]])
  r2_1_r3_cor=c(r2_1_r3_cor,r)
}
data$r2_r3_cor=r2_1_r3_cor

#r1_1_r3_cor=c()
#for(i in 1:nrow(data)){
#  r=cor(r1[,data$r1_factor_1[i]],r3[,data$r3_factor_2[i]])
#  r1_1_r3_cor=c(r1_1_r3_cor,r)
#}
#data$r1_1_r3_cor=r1_1_r3_cor

count=c()
for(i in 1:nrow(data)){
  count=c(count,sum(abs(data[i,c(4:6)])>=0.75,na.rm=T))
}
#count=apply(data,MARGIN=1,function(x) sum(abs(x[c(6,8,9)])>=0.75,na.rm=T))

strong=which(count==3)
data[strong,]
med=which(count>=2)
data[med,]


#data[,c('r2_factor_2','r2_r3_cor')]=data[match(data$r2_factor_2,data$r2_factor1),c('r2_factor_2','r2_r3_cor')]

fwrite(data,'vst_MegaLMM_WD_0712_factor_correlations.txt',row.names=F,quote=F,sep='\t')

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

fwrite(lambda_all_means,'vst_WD_0712_all_Lambda_means.txt',row.names=T,quote=F,sep='\t')

# Average F values across runs

r1=fread(sprintf('vst_MegaLMM_%s_1/F_means.txt',time),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]

#r1_probs=paste0('Factor',c(76,77,78,79,80,83,92,93,99))

#vst
r1_probs=paste0('Factor',c(37,38,39,40,42,43,44,45,46,47,48,49,50,51,52,
53,54,55,56,57,58,59,60,63,64,65,67,69,70,71,72,74,75,76,77,78,79,80,81,82,83,
84,85,86))
r1=r1[,!names(r1) %in% (r1_probs)]
#rownames(r1)=seq(1,nrow(r1))
#r1=t(as.matrix(r1))

r2=fread(sprintf('vst_MegaLMM_%s_2/F_means.txt',time),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]

#r2_probs=paste0('Factor',c(67,68,72,76,77,81,82,89,99))

#vst
r2_probs=paste0('Factor',c(35,36,37,38,39,40,41,42,43,44,45,46,47,48,51,52,53,54,
55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,80,81,82,83,84,
85,86,87,88,89,90,91))
r2=r2[,!names(r2) %in% (r2_probs)]
#rownames(r2)=seq(1,nrow(r2))
#r2=t(as.matrix(r2))

r3=fread(sprintf('vst_MegaLMM_%s_3/F_means.txt',time),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]

#r3_probs=paste0('Factor',c(69,77,81,85,86,97,88,89,90,91,99))

#vst
r3_probs=paste0('Factor',c(26,27,29,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,
48,49,51,52,53,54,55,56,57,60,61,62,64,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
80,81,82,83,84,85,90))
r3=r3[,!names(r3) %in% (r3_probs)]

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
  #print(mult)
  #print(row)
  #mult=isneg*-1
  r1f=row$r1_factor_1
  r2f=row$r2_factor_1
  r3f=row$r3_factor_1
  newdf=cbind(r1[,r1f],r2[,r2f],r3[,r3f])
  newdf=sapply(1:ncol(newdf),function(x) newdf[,x] * mult[x] )
  #newdf=newdf*mult
  fmean=as.data.frame(apply(newdf,MARGIN=1,mean),stringsAsFactors=F)
  f_all_means[,i]=fmean
}

fwrite(f_all_means,'vst_WD_0712_all_F_means.txt',row.names=T,quote=F,sep='\t')




#l=readRDS('pheno_MegaLMM_residuals_WD_0712_1/Posterior/Lambda_250.rds')]
#n_k=100
#n_genes=dim(l)[3]
#l_means=sapply(seq(1:n_k),function(i) sapply(seq(1:n_genes),function(j) mean(l[,i,j])))

#data=fread('eqtl/results/WD_0712_residuals_x_phenotypes.txt',data.table=F)
#data=colnames(data)[-1]

#l_means=as.data.frame(l_means,stringsAsFactors=F)
#colnames(l_means)=paste0('Factor',seq(1,100))
#l_means$Gene=data
#l_means=l_means[,c('Gene',paste0('Factor',seq(1,100)))]

#phenotypes=unique(data[!grepl('Zm',data)])

# With phenotypes included
for(r in 1:3){

  run_id=sprintf('pheno_MegaLMM_vst_residuals_%s_%.0f',time,r)

  exp=fread(sprintf('eqtl/results/%s_vst_residuals_x_phenotypes.txt',time),data.table=F)
  key=exp[,c('ID'),drop=F]
  rownames(exp)=exp$ID
  Y = exp[,-c(1,2)]
  data = key

  K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
  rownames(K)=K[,1]
  rownames(K)=gsub("-",".",rownames(K))
  K=as.matrix(K[,-1])
  colnames(K)=rownames(K)

  inter=intersect(data$ID,rownames(K))
  K=K[inter,inter]
  run_parameters = MegaLMM_control(
    max_NA_groups = 3,
    scale_Y = FALSE,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
    h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
    h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
    burn = 10000,  # number of burn in samples before saving posterior samples
    thin = 40,
    K = 100 # number of factors
  )

  #Set the prior hyperparameters of the BSFG model
  priors = MegaLMM_priors(
    tot_Y_var = list(V = 0.5,   nu = 5),      # Prior variance of trait residuals after accounting for fixed effects and factors
    tot_F_var = list(V = 18/20, nu = 20),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
    Lambda_prior = list(
      sampler = sample_Lambda_prec_horseshoe, # function that implements the horseshoe-based Lambda prior described in Runcie et al 2020. See code to see requirements for this function.
      prop_0 = 0.1,    # prior guess at the number of non-zero loadings in the first and most important factor
      delta = list(shape = 3, scale = 1),    # parameters of the gamma distribution giving the expected change in proportion of non-zero loadings in each consecutive factor
      delta_iterations_factor = 100   # parameter that affects mixing of the MCMC sampler. This value is generally fine.
    ),
    h2_priors_resids_fun = function(h2s,n)  1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
    h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
  )


  MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                                ~(1|ID),  # RHS of base model for factors and residuals. Fixed effects defined here only apply to the factor residuals.
                                data = data,         # the data.frame with information for constructing the model matrices
                                relmat = list(ID = K), # covariance matrices for the random effects. If not provided, assume uncorrelated
                                run_parameters=run_parameters,
                                run_ID = run_id
  )

  #Reload posterior
  MegaLMM_state$Posterior=reload_Posterior(MegaLMM_state,params=c("Lambda","F"))
  Lambda_mean=get_posterior_mean(MegaLMM_state$Posterior$Lambda)
  Lambda_mean=as.data.frame(Lambda_mean,stringsAsFactors=F)

  rownames(Lambda_mean)=paste0('Factor',seq(1,100))
  fwrite(Lambda_mean,paste0(run_id,'/Lambda_means.txt'),row.names=T,quote=F,sep='\t')

  F_mean=get_posterior_mean(MegaLMM_state$Posterior$F)
  F_mean=as.data.frame(F_mean,stringsAsFactors=F)
  colnames(F_mean)=paste0('Factor',seq(1,100))
  rownames(F_mean)=inter
  fwrite(F_mean,paste0(run_id,'/F_means.txt'),row.names=T,quote=F,sep='\t')
}





r1=fread(sprintf('pheno_MegaLMM_vst_residuals_%s_1/Lambda_means.txt',time),data.table=F)
rownames(r1)=r1$V1
r1=r1[,-1]

r1_probs=paste0('Factor',c(28,30,31,32,33,34,35,36,37,39,40,41,42,43,46,47,48,49,
50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,66,67,68,69,70,71,72,73,74,75,76,77,
78,79,80,81,82,83,84,85,86,87,88,89,90,94,98))

r1=r1[!rownames(r1) %in% (r1_probs),]
#rownames(r1)=seq(1,nrow(r1))
r1=t(as.matrix(r1))
#r1_probs=c()
#r1=r1[!rownames(r1) %in% (r1_probs),]
#rownames(r1)=seq(1,nrow(r1))
#r1=t(as.matrix(r1))

#Do this with Lambda! Not U_F


r2=fread(sprintf('pheno_MegaLMM_vst_residuals_%s_2/Lambda_means.txt',time),data.table=F)
rownames(r2)=r2$V1
r2=r2[,-1]

r2_probs=paste0('Factor',c(31,33,34,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,
  52,53,54,55,56,57,58,59,60,61,62,63,64,66,67,69,70,71,72,73,74,75,76,78,80,81,
  82,83,84,85,86,87,88,89,90,91,92))
r2=r2[!rownames(r2) %in% (r2_probs),]
#rownames(r2)=seq(1,nrow(r2))
r2=t(as.matrix(r2))

r3=fread(sprintf('pheno_MegaLMM_vst_residuals_%s_3/Lambda_means.txt',time),data.table=F)
rownames(r3)=r3$V1
r3=r3[,-1]

r3_probs=paste0('Factor',c(40,41,47,48,49,50,51,52,53,55,56,57,58,60,61,62,63,64,
  65,66,67,68,70,71,72,73,74,75,76,77,79,80,81,82,84,85,86,87))
r3=r3[!rownames(r3) %in% (r3_probs),]
#rownames(r3)=seq(1,nrow(r3))
r3=t(as.matrix(r3))

#colnames(r1)=paste0('Factor',colnames(r1))
#colnames(r2)=paste0('Factor',colnames(r2))
#colnames(r3)=paste0('Factor',colnames(r3))


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

#for(i in 1:ncol(r2)){
  #r2_names=data$r2_factor_1
#  max_factor=NA
#  max_cor=0
#  for(j in 1:ncol(r3)){
#    cur_cor=cor(r2[,i],r3[,j],use="complete.obs")
#    if(abs(cur_cor)>abs(max_cor)){
#      max_cor=cur_cor
#      max_factor=colnames(r3)[j]
#    }
#  }
#  data$r2_r3_cor[i]=max_cor
#  data$r3_factor_2[i]=max_factor
#}

r2_1_r3_cor=c()
for(i in 1:nrow(data)){
  r=cor(r2[,data$r2_factor_1[i]],r3[,data$r3_factor_1[i]])
  r2_1_r3_cor=c(r2_1_r3_cor,r)
}
data$r2_r3_cor=r2_1_r3_cor

#r1_1_r3_cor=c()
#for(i in 1:nrow(data)){
#  r=cor(r1[,data$r1_factor_1[i]],r3[,data$r3_factor_2[i]])
#  r1_1_r3_cor=c(r1_1_r3_cor,r)
#}
#data$r1_1_r3_cor=r1_1_r3_cor

count=c()
for(i in 1:nrow(data)){
  count=c(count,sum(abs(data[i,c(4:6)])>=0.75,na.rm=T))
}
#count=apply(data,MARGIN=1,function(x) sum(abs(x[c(6,8,9)])>=0.75,na.rm=T))

strong=which(count==3)
data[strong,]
med=which(count>=2)
data[med,]


#data[,c('r2_factor_2','r2_r3_cor')]=data[match(data$r2_factor_2,data$r2_factor1),c('r2_factor_2','r2_r3_cor')]

fwrite(data,'pheno_MegaLMM_vst_residuals_WD_0712_factor_correlations.txt',row.names=F,quote=F,sep='\t')

sdata=data[strong,]
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

phenotypes=rownames(lambda_all_means)[!grepl('Zm',rownames(lambda_all_means))]

fwrite(lambda_all_means,'pheno_MegaLMM_vst_residuals_WD_0712_all_Lambda_means.txt',row.names=T,quote=F,sep='\t')


factors=names(lambda_all_means)
factor_groups=vector("list",length=nrow(lambda_all_means))
for(i in 1:length(factors)){
  factor_groups[[i]]$factor=factors[i]
  factor_groups[[i]]$genes=c(NA)
}
#v=apply(lambda_all_means,MARGIN=1,function(x) x**2)
#v2=apply(v,MARGIN=2,function(x) which(x>0.2))

for(f in 1:nrow(lambda_all_means)){
  #index=f-1
  subl=lambda_all_means[f,,drop=F]
  gene=rownames(subl)
  var_exp=apply(subl,MARGIN=1,function(x) x**2)
  tot_var=sum(var_exp)
  prop_var=var_exp/tot_var
  fkeep=names(subl[,which(prop_var>=0.2),drop=F])
  for(k in fkeep){
    x=unlist(unname(lapply(factor_groups2,function(x) which(x$factor==k))))
    factor_groups[[x]]$genes=c(factor_groups[[x]]$genes,gene)

  }
  #print(length(fkeep))
}

# total variance
#lambda[,i]^2 + 1/tot_Eta_prec[i]

#factor_groups=vector("list",length=nrow(sdata))
#for(f in 1:nrow(sdata)){
  #index=f-1
#  subl=lambda_all_means[,f,drop=F]
#  var_exp=subl[,1]**2
#  tot_var=sum(var_exp)
#  prop_var=var_exp/tot_var
#  fgenes=rownames(subl[(subl[,1]**2)>=0.2,,drop=F])
#  factor_groups[[f]]=list(factor=paste0('Factor',f),genes=fgenes)
#  print(length(fgenes))
#}

pheno_lambdas=lambda_all_means[phenotypes,]

pfactors=which(unlist(unname(lapply(factor_groups,function(x) x$genes %in% phenotypes))))
