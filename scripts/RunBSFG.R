#Running BSFG on Busch Root Data
library('BSFG')

#Read in data
roots=read.table('/home/sodell/projects/BSFG/RootTraitMatrix.txt',sep='\t',header=TRUE)
key=read.table('/home/sodell/projects/BSFG/RootTraitPlantKey.txt',sep='\t',header=TRUE)

Y = roots[,-1]
data = key
#K = setup$K covariance matrix?

run_parameters = BSFG_control(
  sampler='fast_BSFG',
  scale_Y = TRUE,
  simulation = FALSE,
  h2_divisions = 2,
  h2_step_size = NULL,
  burn = 100,
  k_init = 15
)

#Set the prior hyperparameters of the BSFG model
priors = BSFG_priors(
  fixed_var=list(V=1,nu=3),
  tot_Y_var = list(V = 0.5,   nu = 5),
  tot_F_var = list(V = 18/20, nu = 20),
  delta_1   = list(shape = 2.1,  rate = 1/20),
  delta_2   = list(shape = 3, rate = 1),
  Lambda_df=3,
  B_df=3,
  B_F_df=3,
  h2_priors_resids_fun = function(h2s,n) 1,
  h2_priors_factors_fun = function(h2s,n) 1
)

#Construct the model
BSFG_state = BSFG_init(Y,
                       model=~1 + (1|ACC_ID),
                       factor_model_fixed=NULL, 
                       data = data, # the data.frame with information for constructing the model matrices
                       #K_mats = list(animal = K), # covariance matrices for the random effects. If not provided, assume uncorrelated
                       run_parameters=run_parameters,
                       priors=priors #,
                       #setup = setup  # only if running simulated data. Stores simulated values for comparison
)


n_samples = 100;  # how many samples to collect at once?
for(i  in 1:70) {
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,grainSize=1)  # run MCMC chain n_samples iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
  # set of commands to run during burn-in period to help chain converge
  if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
    BSFG_state = reorder_factors(BSFG_state) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
    BSFG_state = clear_Posterior(BSFG_state)
    # BSFG_state$current_state = update_k(BSFG_state) # use to drop insignificant factors
    BSFG_state$run_parameters$burn = max(BSFG_state$run_parameters$burn,BSFG_state$current_state$nrun+100) # if you made changes, set a new burn-in period
    print(BSFG_state$run_parameters$burn)
  }
  BSFG_state = save_posterior_chunk(BSFG_state)  # save any accumulated posterior samples in the database to release memory
  print(BSFG_state) # print status of current chain
  plot(BSFG_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
}