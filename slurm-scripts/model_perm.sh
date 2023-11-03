#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J perm
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-100%20
#SBATCH --ntasks=4
#SBATCH --mem 16G

module load R

Rscript scripts/model_perm.R $SLURM_ARRAY_TASK_ID 4