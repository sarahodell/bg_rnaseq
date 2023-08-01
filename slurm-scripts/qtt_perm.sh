#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J qtt_perm
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 10:00:00
#SBATCH --array=1-181%20
#SBATCH --ntasks=16
#SBATCH --mem 16G

module load R

#Rscript scripts/ciseqtl_weights_update.R $time $chr 1
Rscript scripts/qtt_perm.R $SLURM_ARRAY_TASK_ID 16