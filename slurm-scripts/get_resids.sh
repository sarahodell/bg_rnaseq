#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J resids
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 10:00:00
#SBATCH --array=1-40
#SBATCH --ntasks=16
#SBATCH --mem 16G

module load R

time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f2 -d,)"


echo $time
echo $chr
#Rscript scripts/get_trans_residuals.R $time $chr 16
Rscript scripts/trans_bv2.R $time $chr 16

