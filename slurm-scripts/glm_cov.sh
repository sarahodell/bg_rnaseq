#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J glm_fp
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/%A_%a.error
#SBATCH -t 48:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=4
#SBATCH --mem=4G

module load R

#pheno="male_flowering_d6"
#env="EXP_STPAUL_2017_WD"
#chr=$SLURM_ARRAY_TASK_ID
#pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" QTL/pheno_env_chr.txt | cut -f2 -d,)"
#env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" QTL/pheno_env_chr.txt | cut -f3 -d,)"
#chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" QTL/pheno_env_chr.txt | cut -f1 -d,)"
#pheno="male_flowering_d6"
#env="EXP_STPAUL_2017_WD"
#thresh=0.10
#chr=$SLURM_ARRAY_TASK_ID

pheno="male_flowering_d6"
env="ALL"
chr=$SLURM_ARRAY_TASK_ID
id='qDTA3_2'
echo $pheno
echo $env
#echo $thresh
echo $chr
echo $id

#Rscript scripts/GridLMM_fp_updated.R $pheno $env $chr 1
Rscript scripts/GridLMM_covariate.R $pheno $env $chr $id 4

#Rscript scripts/GridLMM_fp_pheno_x_env.R $pheno $env $chr 1
#if [ $env == "ALL" ]
#then
#  echo "ALL environments"
#  Rscript GridLMM_run_founders.R $pheno $chr 4
#else
#  echo $env
#  Rscript GridLMM_pheno_x_env.R $pheno $env $chr 4
#fi
