#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J qtl_plot
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-56%20
#SBATCH --ntasks=1
#SBATCH --mem 7G

module load R

pheno="$(sed "${SLURM_ARRAY_TASK_ID}q;d" QTL/pheno_env.txt | cut -f1 -d,)"
env="$(sed "${SLURM_ARRAY_TASK_ID}q;d" QTL/pheno_env.txt | cut -f2 -d,)"

echo $pheno
echo $env

Rscript scripts/manhattan_plots.R $pheno $env
