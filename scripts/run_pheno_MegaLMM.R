#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
run=as.character(args[[2]])
cores=as.numeric(args[[3]])

#Running MegaLMM
#installed in R/4.1.0
#devtools::install_github('deruncie/MegaLMM',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.1')
library('generics',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.1')
library('rlang',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.1')
library('vctrs',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.1')
library('glue',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.1')
library('tibble',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.1')
library('tidyselect',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.1')
library('pillar',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.1')

library('MegaLMM')
library('data.table')
#options(warn=2)
#Read in data
run_id=sprintf('MegaLMM/pheno_MegaLMM_residuals_%s_%s',time,run)

phenos=c("female_flowering_d6","male_flowering_d6","total_plant_height","harvest_grain_moisture",
"grain_yield_15","tkw_15",'asi')
envs=c("BLOIS_2014_OPT","BLOIS_2017_OPT","GRANEROS_2015_OPT","NERAC_2016_WD",
"STPAUL_2017_WD","SZEGED_2017_OPT")

exp=fread(sprintf('eqtl/results/%s_vst_residuals_x_phenotypes.txt',time),data.table=F)

iterations=50000
n_iter = 1000 # how many samples to collect at once?
runs=50
#burn_drop=0.5
burn=0
burnin=30
k=100
#Increase thinning?
thin=40



start_time <- Sys.time()

#exp=fread(sprintf('eqtl/results/cis_eQTL_%s_all_residuals.txt',time),data.table=F)
# full gene counts
#exp=fread(sprintf('eqtl/%s_voom_normalized_gene_counts.txt',time),data.table=F)
# need to t() and add genotype codes for this to work
#meta=fread('metadata/BG_completed_sample_list.txt',data.table=F)
#snames=colnames(exp)[-1]
#gnames=meta[match(snames,meta$sample_name),]$dh_genotype
#znames=exp$V1
#colnames(exp)=c('ID',gnames)
#rownames(exp)=exp$ID
#exp=exp[,-1]
#exp=as.matrix(exp)
#exp=t(exp)
#exp=as.data.frame(exp,stringsAsFactors=F)
#exp$ID=rownames(exp)
#exp=exp[,c('ID',znames)]


#phenotypes=fread('../GridLMM/phenotypes_asi.csv',data.table=F)
#phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
#phenotypes2=phenotypes[phenotypes$Genotype_code %in% exp$ID,]
#phenotypes2=phenotypes2[,c('Genotype_code',"Loc.Year.Treat",phenos)]
#newpheno=data.frame(ID=exp$ID,stringsAsFactors=F)
#for(e in envs){
#  df=phenotypes2[phenotypes2$Loc.Year.Treat==e,]
#  names(df)=c('ID','Loc.Year.Treat',paste0(e,'-',phenos))
#  df=df[match(newpheno$ID,df$ID),]
#  rownames(df)=seq(1,nrow(df))
#  newpheno=cbind(newpheno,df[3:9])
#}
#newpheno_scaled=as.data.frame(sapply(seq(2,ncol(newpheno)),function(x) (newpheno[,x]-mean(newpheno[,x],na.rm=T))/sd(newpheno[,x],na.rm=T)),stringsAsFactors=F)
#names(newpheno_scaled)=names(newpheno)[-1]

#exp=cbind(exp,newpheno_scaled)
#residuals
#fwrite(exp,sprintf('eqtl/results/%s_residuals_x_phenotypes.txt',time),row.names=F,quote=F,sep='\t')
#gene counts
#fwrite(exp,sprintf('eqtl/results/%s_genecounts_x_phenotypes.txt',time),row.names=F,quote=F,sep='\t')

#testing
#exp=exp[,1:100]
key=exp[,c('ID'),drop=F]
rownames(exp)=exp$ID
#roots=read.table('/home/sodell/projects/BSFG/RootTraitMatrix.txt',sep='\t',header=TRUE)
#key=read.table('/home/sodell/projects/BSFG/RootTraitPlantKey.txt',sep='\t',header=TRUE)

Y = exp[,-1]
#Y=as.matrix(Y)
Y =as.data.frame(apply(Y,MARGIN=2,function(x)  (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
data = key

K=fread('../GridLMM/K_matrices/K_matrix_full.txt',data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

inter=intersect(data$ID,rownames(K))
K=K[inter,inter]
#K = setup$K covariance matrix?

#cf_DFinf2NA = function(x){
#  for(i in 1:ncol(x)){
#    x[,i][is.infinite(x[,i])]=NA
#  }
#  return(x)
#}

#Y=cf_DFinf2NA(Y)

run_parameters = MegaLMM_control(
  max_NA_groups = 3,
  scale_Y = FALSE,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
  h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
  h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
  burn = burn,  # number of burn in samples before saving posterior samples
  thin = thin,
  K = k # number of factors
)

#Set the prior hyperparameters of the BSFG model
priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.5,   nu = 5),      # Prior variance of trait residuals after accounting for fixed effects and factors
  tot_F_var = list(V = 1, nu = 100000),
  #tot_F_var = list(V = 18/20, nu = 20),
  Lambda_prior = list(
      sampler = sample_Lambda_prec_ARD,
      Lambda_df = 3,
      delta_1 = list(shape = 20, rate = 1/2),
      delta_2 = list(shape = 3, rate = 1)
    ),
    # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
  #Lambda_prior = list(
  #  sampler = sample_Lambda_prec_horseshoe, # function that implements the horseshoe-based Lambda prior described in Runcie et al 2020. See code to see requirements for this function.
  #  prop_0 = 0.1,    # prior guess at the number of non-zero loadings in the first and most important factor
    #delta = list(shape = 3, scale = 1),    # parameters of the gamma distribution giving the expected change in proportion of non-zero loadings in each consecutive factor
  #  delta = list(shape = 30, scale = 1/10),
  #  delta_iterations_factor = 100   # parameter that affects mixing of the MCMC sampler. This value is generally fine.
  #),
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

maps = make_Missing_data_map(MegaLMM_state)
MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map)
MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)  # apply the priors
saveRDS(MegaLMM_state,sprintf('%s/MegaLMM_state_base.rds',run_id))

MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state) # initialize the model
MegaLMM_state = initialize_MegaLMM(MegaLMM_state) # run the initial calculations
MegaLMM_state = clear_Posterior(MegaLMM_state) # prepare the output directories


MegaLMM_state$Posterior$posteriorSample_params = c('Lambda',
                                                   'F','U_F',
                                                   'F_h2','resid_h2','tot_Eta_prec',
                                                   'delta')
MegaLMM_state = clear_Posterior(MegaLMM_state) # prepare the output directories


# The following code is optional, but tries to guess for your system how many CPUs to use for fastest processing
#n_threads = optimize_n_threads(MegaLMM_state,seq(1,RcppParallel::defaultNumThreads(),by=1),times=2)
#set_MegaLMM_nthreads(n_threads$optim)

#set_MegaLMM_nthreads(cores)

set_omp_nthreads(RcppParallel::defaultNumThreads()/2)
set_MegaLMM_nthreads(RcppParallel::defaultNumThreads()/2)
# now do sampling is smallish chunks

#Burnin
for(i  in 1:burnin) {
  if(i %% 2==0){
    print(sprintf('Run %d',i))
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_samples iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
    MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
    print(MegaLMM_state) # print status of current chain
    plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
    MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
    #MegaLMM_state$current_state$F = scale(MegaLMM_state$current_state$F)
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    print(MegaLMM_state$run_parameters$burn)
  }else{
    print(sprintf('Run %d',i))
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_samples iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
    MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
    print(MegaLMM_state) # print status of current chain
    plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
    #MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
    MegaLMM_state$current_state$F = scale(MegaLMM_state$current_state$F)
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    print(MegaLMM_state$run_parameters$burn)
  }
}
#Start sampling from posterior
for(i  in (burnin+1):runs) {
  print(sprintf('Run %d',i))
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_samples iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
  MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
  print(MegaLMM_state) # print status of current chain
  plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
}

end_time <- Sys.time()

# MegaLMM_state$Posterior = reload_Posterior(MegaLMM_state)
# U_hat = get_posterior_mean(MegaLMM_state,U_R + U_F %*% Lambda)
# all parameter names in Posterior
MegaLMM_state$Posterior$posteriorSample_params
MegaLMM_state$Posterior$posteriorMean_params  # these ones only have the posterior mean saved, not individual posterior samples
print("F_var")
print(apply(MegaLMM_state$current_state$F,2,var))
print("delta")
print(MegaLMM_state$current_state$delta)
print("Lambda_tau2")
print(MegaLMM_state$current_state$Lambda_tau2)

print("Cumprod(delta[1:20])")
print(cumprod(MegaLMM_state$current_state$delta[1:20]))
# instead, load only a specific parameter
# Lambda = load_posterior_param(MegaLMM_state,'Lambda')
# boxplots are good ways to visualize Posterior distributions on sets of related parameters

#MegaLMM_state$Posterior$F_h2 = load_posterior_param(MegaLMM_state,'F_h2')
#png(paste0(run_id,'/',sprintf('%s_posterior_F_h2_boxplot.png',time)))
#print(boxplot(MegaLMM_state$Posterior$F_h2[,1,]))
#dev.off()

# Run info file
run_time=end_time-start_time


line0=paste0(run_id,'\n')
line1="Full gene set and phenotypes"
line2=paste0(time,' Timepoint\n')
line3=paste0(cores,' CPU')
line4="60G"
line5=paste0(run_time,' minute run time\n')
line6="Parameters\n"
line7=paste0("burn: ",burn)
line8=paste0("iterations: ",iterations)
line9=paste0("n_iter: ",n_iter)
line10=paste0("runs: ",runs,'\n')
line11=paste0("thin: ",thin)
line12=paste0("K: ",k)


fileConn<-file(paste0(run_id,"/about.txt"))
writeLines(c(line0,line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,line12), fileConn)
close(fileConn)
