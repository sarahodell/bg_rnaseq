#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J fdr
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=22
#SBATCH --mem 128G

module load R

#Rscript scripts/all_trans_weighted_plots.R
Rscript scripts/trans_fixed.R
