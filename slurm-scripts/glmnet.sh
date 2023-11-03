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
#1 WD_0712      0.137           0.997         0.520 ****
#2 WD_0718     -0.0785          0.275         0.106
#3 WD_0720     -0.00345         0.925         0.330
#4 WD_0727     -0.0486          0.735         0.172
#              [,1]
#WD_0712  0.4295083
#WD_0718 -0.1931142
#WD_0720  0.5228916
#WD_0727  0.2561534
#Total expression
#[1] "WD_0712 4 Fold Validation"
#[1] "WD_0718 7 Fold Validation"
#[1] "WD_0720 11 Fold Validation"
#[1] "WD_0727 9 Fold Validation"
#               [,1]
#WD_0712 -0.10442449
#WD_0718 -0.05372081
#WD_0720  0.02910795
#WD_0727 -0.19673955
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.223          0.462       -0.115 
#2 WD_0718      -0.0408         0.234        0.0343
#3 WD_0720      -0.0548         0.162       -0.0409
#4 WD_0727      -0.0408         0.179       -0.0375
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
#1 WD_0712      -0.113          0.773        0.311 
#2 WD_0718      -0.166          0.457        0.0318
#3 WD_0720      -0.0341         0.911        0.270 
#4 WD_0727      -0.0908         0.553        0.225 
#                [,1]
#WD_0712 -0.252740597
#WD_0718  0.006215678
#WD_0720  0.638558681
#WD_0727  0.313408926
#Total expression
#[1] "WD_0712 3 Fold Validation"
#[1] "WD_0718 6 Fold Validation"
#[1] "WD_0720 9 Fold Validation"
#[1] "WD_0727 8 Fold Validation"
#               [,1]
#WD_0712 -0.04886811
#WD_0718 -0.37867951
#WD_0720  0.12500611
#WD_0727 -0.09278258
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.187          0.902       -0.0757
#2 WD_0718      -0.225          0.477       -0.0371
#3 WD_0720      -0.0255         0.223        0.0638
#4 WD_0727      -0.114          0.234       -0.0503
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
#1 WD_0712      -0.0693         0.693        0.145 
#2 WD_0718      -0.120          0.446        0.0485
#3 WD_0720      -0.0651         0.504        0.125 
#4 WD_0727      -0.0223         0.874        0.270 
#              [,1]
#WD_0712  0.1459184
#WD_0718 -0.4307078
#WD_0720  0.5870080
#WD_0727  0.5557064
#Total expression
#[1] "WD_0712 4 Fold Validation"
#[1] "WD_0718 7 Fold Validation"
#[1] "WD_0720 11 Fold Validation"
#[1] "WD_0727 9 Fold Validation"
#              [,1]
#WD_0712 -0.2254243
#WD_0718 -0.3099085
#WD_0720 -0.4086227
#WD_0727  0.3124282
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.261          0.637       -0.0733
#2 WD_0718      -0.0811         0.424        0.111 
#3 WD_0720      -0.0592         0.160       -0.0588
#4 WD_0727      -0.225          0.506        0.0849
#
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
#1 WD_0712      -0.0340         0.682         0.263
#2 WD_0718      -0.0249         0.428         0.203
#3 WD_0720      -0.117          0.337         0.134
#4 WD_0727      -0.0315         0.769         0.148
#               [,1]
#WD_0712  0.19733904
#WD_0718 -0.06780383
#WD_0720  0.54327443
#WD_0727  0.19187122
#Total expression
#[1] "WD_0712 3 Fold Validation"
#[1] "WD_0718 6 Fold Validation"
#[1] "WD_0720 9 Fold Validation"
#[1] "WD_0727 8 Fold Validation"
#               [,1]
#WD_0712 -0.08468207
#WD_0718 -0.34562232
#WD_0720 -0.62541671
#WD_0727  0.47717416
## A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.0826         0.742       -0.0143
#2 WD_0718      -0.0267         0.325       -0.0149
#3 WD_0720      -0.0225         0.303        0.131 
#4 WD_0727      -0.0797         0.431        0.0198