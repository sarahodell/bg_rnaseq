#devtools::install_github("deruncie/SparseFactorMixedModel",subdir="BSFG",ref="one_general_model")
devtools::install_github('deruncie/MegaLMM')
library('MegaLMM')

seed = 1  # for reproducibility
nSire = 50 # simulation design is a half-sib design: each of nSire fathers has nRep children, each with a different female
nRep = 10
nTraits = 100
nFixedEffects = 2 # fixed effects are effects on the factors
nFactors = 10
factor_h2s = c(rep(0,nFactors/2),rep(0.3,nFactors/2))  # factor_h2 is the % of variance in each factor trait that is explained by additive genetic variation.
# This is after accounting for the fixed effects on the factors
Va = 2 # residual genetic variance in each of the observed traits after accounting for the factors
Ve = 2 # residual microenvironmental variance in each of the observed traits after accounting for the factors
Vb = 2 # magnitude of the fixed effects (just factors)
new_halfSib_simulation('Sim_FE_1', nSire=nSire,nRep=nRep,p=nTraits, b=nFixedEffects, factor_h2s= factor_h2s,Va = Va, Ve = Ve,Vb = Vb)
set.seed(seed)
load('setup.RData')

Y = setup$Y
data = setup$data
K = setup$K

#Set the parameters of the BSFG model
run_parameters = BSFG_control(
  #sampler = 'fast_BSFG',  #
  #sampler = 'general_BSFG',
  scale_Y = FALSE,
  simulation = TRUE,
  h2_divisions = 20,
  h2_step_size = NULL,
  burn = 00,
  k_init = 15
)

#Set the prior hyperparameters of the BSFG model
priors = BSFG_priors(
  tot_Y_var = list(V = 0.5,   nu = 5),
  tot_F_var = list(V = 18/20, nu = 20),
  #delta_1   = list(shape = 2.1,  rate = 1/20),
  #delta_2   = list(shape = 3, rate = 1),
  Lambda_prior = list(sampler =
                        sample_Lambda_prec_ARD, Lambda_df = 3,delta_1   = list(shape = 2.1,  rate = 1/20),
                      delta_2   = list(shape = 3, rate = 1)),
  B_prior = list(sampler = sample_B_prec_ARD, global = list(V = 1, nu = 3),
                 global_F = list(V = 1, nu = 3), B_df = 3, B_F_df = 3),
  h2_priors_resids_fun = function(h2s,n) 1,
  h2_priors_factors_fun = function(h2s,n) 1
)


#BSFG_priors(tot_Y_var = list(V = 0.5, nu = 3), tot_F_var = list(V = 18/20,
#nu = 20), h2_priors_resids_fun = function(h2s, n) 1,
#h2_priors_factors_fun = function(h2s, n) 1, Lambda_prior = list(sampler =
#sample_Lambda_prec_ARD, Lambda_df = 3, delta_1 = list(shape = 2, rate = 1),
#delta_2 = list(shape = 3, rate = 1), delta_iteractions_factor = 100),
#B_prior = list(sampler = sample_B_prec_ARD, global = list(V = 1, nu = 3),
#global_F = list(V = 1, nu = 3), B_df = 3, B_F_df = 3),
#QTL_prior = list(sampler = sample_QTL_prec_horseshoe, separate_QTL_shrinkage
# = T, cauchy_iteractions_factor = 10))

#Construct the model
BSFG_state = BSFG_init(Y,
                       model=~Fixed1+(1|animal), # This model has one fixed and one random term.
                       #factor_model_fixed = ~0, # we could specify a different fixed effect model for the factors
                       data = data, # the data.frame with information for constructing the model matrices
                       K_mats = list(animal = K), # covariance matrices for the random effects. If not provided, assume uncorrelated
                       run_parameters=run_parameters,
                       priors=priors,
                       setup = setup  # only if running simulated data. Stores simulated values for comparison
)

#Run MCMC

#A MCMC chain is a way of fitting parameters from a Bayesian model. The chain is a sequence of draws from
#the Posterior distribution of each parameter. BSFG uses a Gibbs sampler, which means that we iterate
#through all of the model's parameters and for each parameter draw a new value from its posterior holding
#all other parameters constant. This works, but the individual draws are not independent draws from
#the joint posterior of all parameters. Therefore, to get independent draws, we have to collect many
#posterior samples, and then save only a portion of them.

#Also, it may take many iterations for the MCMC chain to converge to the central part of the distribution
#This is the burnin period, and we do not stor these values (they are not useful). At the end,
#we have a collection of samples from the posterior distribution of each parameter. We can use these
#samples to characterize each distribution: mean, SD, histogram, HDinterval, etc.
#The way BFSG works is that you ask the program to collect a small number of samples and then can
#assess how the chain is performing (how well it is mixing, is it converged?)
#And then either save the samples as posterior samples, or declare them as burn-in and discard them.
#Then, you can ask for more samples, and repeat until you have enough.

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


dim(BSFG_state$Posterior$Lambda)
#Work with the posterior samples
