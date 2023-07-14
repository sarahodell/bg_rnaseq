#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression/
#SBATCH -J glm_perm
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=1-560%50
#SBATCH --ntasks=8
#SBATCH --mem=8G

module load R

pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" QTL/pheno_env_chr.txt | cut -f2 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" QTL/pheno_env_chr.txt | cut -f3 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" QTL/pheno_env_chr.txt | cut -f1 -d,)"

echo $pheno
echo $env
echo $chr

#Rscript scripts/transeqtl_permute.R $time $factor $chr 4 1000

Rscript scripts/GridLMM_randomized_pheno_x_env.R $pheno $env $chr 8 1000

#echo "Finished array task $SLURM_ARRAY_TASK_ID" >> report.txt
