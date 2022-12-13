#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J filter_eqtl
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH --array=2-66
#SBATCH --ntasks=1
#SBATCH --mem 3G

module load R

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"


time="WD_0712"
echo $time
#f="$(sed "${SLURM_ARRAY_TASK_ID}q;d" 'eqtl/pheno_residuals_trans_eQTL_list.txt' | cut -f2 -d,)"
f="$(sed "${SLURM_ARRAY_TASK_ID}q;d" 'eqtl/pheno_trans_eQTL_list.txt' | cut -f2 -d,)"

echo $f
Rscript scripts/trans_eqtl_check.R $time $f
