#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J glmnet
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --array=0-3
#SBATCH --ntasks=4
#SBATCH --mem=8G

module load R

phenotypes=("tkw_15" "tkw_15" "grain_yield_15" "grain_yield_15")
environments=('ALL' "EXP_STPAUL_2017_WD" 'ALL' "EXP_STPAUL_2017_WD")

pheno=${phenotypes[$SLURM_ARRAY_TASK_ID]}
env=${environments[$SLURM_ARRAY_TASK_ID]}
echo $pheno
echo $env

echo "Z-score of effect sizes"
Rscript scripts/z_glmnet.R $pheno $env
echo "Total expression"
Rscript scripts/exp_glmnet.R $pheno $env

#tkw_15
#ALL
#Z-score of effect sizes
#[1] "WD_0712 4 Fold Validation"
#[1] "WD_0718 7 Fold Validation"
#[1] "WD_0720 11 Fold Validation"
#[1] "WD_0727 9 Fold Validation"
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712       0.0793         0.997        0.330 ***
#2 WD_0718      -0.0314         0.316        0.0202
#3 WD_0720       0.0457         0.958        0.348  ***
#4 WD_0727      -0.0747         0.749        0.145 
#              [,1]
#WD_0712  0.2368214
#WD_0718 -0.2233250
#WD_0720  0.6471975
#WD_0727  0.2043001
#Total expression
#[1] "WD_0712 4 Fold Validation"
#[1] "WD_0718 7 Fold Validation"
#[1] "WD_0720 11 Fold Validation"
#[1] "WD_0727 9 Fold Validation"
#               [,1]
#WD_0712 -0.22656757
#WD_0718 -0.16281887
#WD_0720  0.16590178
#WD_0727  0.03373725
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.0350         0.619      -0.0540 
#2 WD_0718      -0.102          0.377       0.00741
#3 WD_0720      -0.0802         0.160      -0.0888 
#4 WD_0727      -0.0260         0.181       0.0653 
#
#tkw_15
#EXP_STPAUL_2017_WD
#Z-score of effect sizes
#[1] "WD_0712 3 Fold Validation"
#[1] "WD_0718 6 Fold Validation"
#[1] "WD_0720 9 Fold Validation"
#[1] "WD_0727 8 Fold Validation"
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      0.0260          0.838        0.355  ***
#2 WD_0718     -0.0524          0.320        0.0636
#3 WD_0720      0.0955          0.949        0.335 ***
#4 WD_0727     -0.00185         0.581        0.156 
#               [,1]
#WD_0712  0.29914309
#WD_0718 -0.04822599
#WD_0720  0.84522164
#WD_0727  0.56540592
#Total expression
#[1] "WD_0712 3 Fold Validation"
#[1] "WD_0718 6 Fold Validation"
#[1] "WD_0720 9 Fold Validation"
#[1] "WD_0727 8 Fold Validation"
#              [,1]
#WD_0712 -0.2730968
#WD_0718 -0.1514908
#WD_0720 -0.1274734
#WD_0727 -0.2262693
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.136          0.601      -0.160  
#2 WD_0718      -0.120          0.476      -0.00662
#3 WD_0720      -0.0399         0.247       0.0913 
#4 WD_0727      -0.0230         0.193       0.0634 
#
#grain_yield_15
#ALL
#Z-score of effect sizes
#[1] "WD_0712 4 Fold Validation"
#[1] "WD_0718 7 Fold Validation"
#[1] "WD_0720 11 Fold Validation"
#[1] "WD_0727 9 Fold Validation"
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.168          0.859        0.0433
#2 WD_0718      -0.0445         0.368        0.0882
#3 WD_0720      -0.0844         0.522        0.143 
#4 WD_0727      -0.0725         0.846        0.236 
#              [,1]
#WD_0712 0.06828986
#WD_0718 0.22834440
#WD_0720 0.60634839
#WD_0727 0.56597587
#Total expression
#[1] "WD_0712 4 Fold Validation"
#[1] "WD_0718 7 Fold Validation"
#[1] "WD_0720 11 Fold Validation"
#[1] "WD_0727 9 Fold Validation"
#               [,1]
#WD_0712 -0.08746021
#WD_0718 -0.28128119
#WD_0720 -0.60664562
#WD_0727 -0.01247630
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.0213         0.551        0.0528
#2 WD_0718      -0.110          0.324        0.0529
#3 WD_0720      -0.0988         0.161       -0.128 
#4 WD_0727      -0.104          0.511        0.109
##
#grain_yield_15
#EXP_STPAUL_2017_WD
#Z-score of effect sizes
#[1] "WD_0712 3 Fold Validation"
#[1] "WD_0718 6 Fold Validation"
#[1] "WD_0720 9 Fold Validation"
#[1] "WD_0727 8 Fold Validation"
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.0836         0.691         0.100
#2 WD_0718      -0.217          0.617         0.164
#3 WD_0720      -0.0486         0.341         0.179
#4 WD_0727      -0.0834         0.683         0.205
#               [,1]
#WD_0712 -0.04070231
#WD_0718 -0.11517746
#WD_0720  0.39952883
#WD_0727  0.31904738
#Total expression
#[1] "WD_0712 3 Fold Validation"
#[1] "WD_0718 6 Fold Validation"
#[1] "WD_0720 9 Fold Validation"
#[1] "WD_0727 8 Fold Validation"
#              [,1]
#WD_0712 -0.1868796
#WD_0718 -0.5040925
#WD_0720 -0.2981355
#WD_0727  0.1359260
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.0521         0.685      -0.0178 
#2 WD_0718      -0.108          0.341       0.00832
#3 WD_0720      -0.0334         0.341       0.0716 
#4 WD_0727      -0.0938         0.230       0.0489 
