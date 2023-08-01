#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J trans_resids
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 96:00:00
#SBATCH --array=91-100
#SBATCH --ntasks=4
#SBATCH --mem 4G

module load R
g
time="WD_0712" # 22
#time="WD_0718" #22
#time="WD_0720" #12
#time="WD_0727" #19 d

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_chrom_factor.txt | cut -f1 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_residuals_chrom_factor_FIXED.txt | cut -f2 -d,)"
factor="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_residuals_chrom_factor_FIXED.txt | cut -f3 -d,)"

echo $time
echo $chr
echo $factor

Rscript scripts/trans_residuals_GridLMM.R $time $factor $chr 1
