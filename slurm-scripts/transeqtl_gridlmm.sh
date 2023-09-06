#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J transeqtl
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 96:00:00
#SBATCH --array=6,7,27,54
#SBATCH --ntasks=4
#SBATCH --mem 4G

module load R

#time="WD_0712" # 27
#time="WD_0718" #17
#time="WD_0720" #16
time="WD_0727" #12

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_chrom_factor.txt | cut -f1 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_chrom_factor_FIXED.txt | cut -f2 -d,)"
factor="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_chrom_factor_FIXED.txt | cut -f3 -d,)"

echo $time
echo $chr
echo $factor

Rscript scripts/transeqtl_fkeep_GridLMM.R $time $factor $chr 4
