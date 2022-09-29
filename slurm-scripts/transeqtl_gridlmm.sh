#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J transeqtl
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 96:00:00
#SBATCH --ntasks=2
#SBATCH --mem 7G

module load R/4.1.0

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"
#chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f2 -d,)"

#times=( "WD_0718" "WD_0720" "WD_0727" "WD_0712" )
#time=${times[$SLURM_ARRAY_TASK_ID]}
time="WD_0720"
echo $time
for chr in {1..10}; do
  echo $chr
  Rscript scripts/transeqtl_GridLMM.R $time $chr 2
done
