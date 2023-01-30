#!/usr/bin/env Rscript

library(data.table)
library(JointGWAS)
library(Matrix)
library(foreach)
library(parallel)
library(doParallel)
​
run = commandArgs(t=T)
dataset = as.character(run[1])
​
# set ncores. I use the RcppParallel function below, but you can set manually.
ncores = RcppParallel::defaultNumThreads()-1
​
​
#Read in genotype data
# Need to format them as PLINK?
#Y

#Read in phenotype data (F values from MegaLMM)
​
​
​
SampleIDs = data[, 1]
# print(Y)
​
#dimnames(G_hat) = NULL
#dimnames(R_hat) = NULL
​
# check that colnames(Y) == colnames(G_hat) == colnames(R_hat) == rownames(G_hat) == rownames(R_hat)
​
# run GWAS for each chromosome
​
#TC_info = fread('/group/runciegrp2/Projects/SeeD_GWAS_Gates/prepped_data/TestCross_passport.csv',data.table=F)
​
for(chr in 1:5) {
  print(chr)

  # load genotypes

  mat = fread(sprintf('%s/chr%02d_DS.csv',genotype_dir,chr),data.table=F)
  rownames(mat) = mat[,1]
  mat = as.matrix(mat[,-1])
  mat = mat[TC_info$DNA_Name[match(SampleIDs,TC_info$SampleID)],]
  # rownames(mat) = SampleIDs
  # for permutation purposes, shuffle rownames of markers
  rownames(mat) = sample(SampleIDs)
​
  print('finish loading genotypes')

  # discard genotypes with load maf or too much missing data (imputed)
  mean_notNA = colMeans(apply(mat,2,function(x) x %in% 0:2))
  maf = .5 - abs(.5-colMeans(mat)/2)
  select_cols = maf > 0.01 & mean_notNA>0.25
  mat = mat[,select_cols]
​
  print('finish filtering genotypes')

  # form map info for genotypes
  map = fread(sprintf('%s/All_SeeD_2.7_chr%d_nofilter.unimputed.012.pos',genotype_dir,chr),data.table=F)
  map$snp = sprintf('S%d_%d',map$V1,map$V2)
  map = map[match(colnames(mat),map$snp),]
  colnames(map) = c('Chr','V2_pos','V4_pos','snp')
  map$maf = maf[select_cols]
  map$mean_imputed = 1-mean_notNA[select_cols]
​
  print('finish forming map')

  # load K matrix
  K = fread(sprintf('/group/runciegrp2/Projects/SeeD_GWAS_Gates/Genetic_data/allLines/LOO_Ks/K_chr_%02d.csv',chr),data.table=F)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  K = K[SampleIDs,SampleIDs]

  print('finish loading K matrix')
​
  # prep covariance matrices
  sK = svd(K)
  n = nrow(Y)
  sK2 = simultaneous_diagonalize(K,diag(1,n))
  sGR = simultaneous_diagonalize(G_hat,R_hat)
​
  print('finish calcing cov matrices')

  # As a test, to check that everything is working,
  # subset to just the first 10 markers:
  # mat = mat[,1:10]
  # map = map[1:10,]
​
if(dataset == 'PC') {
  results = EMMAX_ANOVA_matrix(cbind(TemperaturePC, PrecipitationPC, TemperatureVarPC, PrecipitationVarPC)~X,
                            Y,mat,'SampleID',svd_matrices = list(sK,sGR),
                            mc.cores = ncores,verbose=T)
} else if (dataset == 'clim' | dataset == 'clim_perm') {
  results = EMMAX_ANOVA_matrix(cbind(altitude, meanTemp, meanTempWarmQ, meanTempColdQ, annualPrecipitation, precipWet, precipDry,
                                precipSeasonal, precipWetQ, precipDryQ, precipWarmQ, precipColdQ, meanDiurnalRange, Isothermality,
                                tempSeasonality, maxTempWarmestMonth, minTempColdestMonth, TempRangeQ, meanTempWetQ, meanTempDryQ)~X,
                            Y,mat,'SampleID',svd_matrices = list(sK,sGR),
                            mc.cores = ncores,verbose=T)
}
​
  write.csv(cbind(map,results$anova),file = sprintf('%s/eGWAS_results_chr_%02d.csv',output_dir,chr),row.names=F)
  write.csv(results$beta_hats,file = sprintf('%s/eGWAS_beta_hats_chr_%02d.csv',output_dir,chr),row.names=F)
  write.csv(results$SEs,file = sprintf('%s/eGWAS_SEs_chr_%02d.csv',output_dir,chr),row.names=F)
}
