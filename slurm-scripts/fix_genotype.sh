#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression/metadata
#SBATCH -J fix_genos
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=8
#SBATCH --mem 16G

module load R


Rscript fix_RNA_genotypes.R 
