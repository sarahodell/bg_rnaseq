#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J ld
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 96:00:00
#SBATCH --ntasks=22
#SBATCH --mem 64G

module load R

#c1="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/data/chr_combos.txt | cut -f1 -d,)"
#c2="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/data/chr_combos.txt | cut -f2 -d,)"

#Rscript scripts/filter_geno.R $c1 $c2 21

Rscript scripts/combine_ld.R