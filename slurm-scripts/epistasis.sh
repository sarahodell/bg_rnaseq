#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J epistasis
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 12:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=8
#SBATCH --mem 16G

module load R

chr=$SLURM_ARRAY_TASK_ID
echo $chr
ENV=ALL
echo $ENV
Rscript ../scripts/epistasis_scan.R $chr $ENV 7