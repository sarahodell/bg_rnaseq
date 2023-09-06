#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J eqtl
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 36:00:00
#SBATCH --ntasks=8
#SBATCH --mem 16G

module load R

Rscript scripts/trans_eqtl_bounds.R

#Rscript scripts/ciseqtl_weights.R $time 5 4
#Rscript scripts/ciseqtl_weights.R $time 9 4
