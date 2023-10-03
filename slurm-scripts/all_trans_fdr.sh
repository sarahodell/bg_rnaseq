#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J perm
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=21
#SBATCH --mem 168G

module load R


#pei="$(sed "${SLURM_ARRAY_TASK_ID}q;d" QTT/pheno_env_id_list.txt | cut -f1 -d,)"

#Rscript scripts/all_trans_weighted_plots.R
#Rscript scripts/trans_fixed.R
#Rscript scripts/trans_r_plot.R
Rscript scripts/cor_diff.R 20