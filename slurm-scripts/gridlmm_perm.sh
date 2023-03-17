#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression/
#SBATCH -J glm_perm
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=1-920%100
#SBATCH --ntasks=4
#SBATCH --mem=4G

module load R

time="WD_0712"
# Permuation 1000 times
factor="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_chrom_factor.txt | cut -f3 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_chrom_factor.txt | cut -f2 -d,)"

echo $time
echo $chr
echo $factor

Rscript scripts/transeqtl_permute.R $time $factor $chr 4 1000

#echo "Finished array task $SLURM_ARRAY_TASK_ID" >> report.txt
