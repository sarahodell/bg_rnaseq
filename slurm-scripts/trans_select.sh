#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J tselect
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 120:00:00
#SBATCH --array=1-40%20
#SBATCH --ntasks=24
#SBATCH --mem 128G

module load R

#1,3,4,11,12,13,14,21,22,23,24,31
#Re-do 31, 32, and 38
#time="WD_0712"
#index=( 31 32 38 )
#index=( 39 40 )
#loc=${index[$SLURM_ARRAY_TASK_ID]}
time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f2 -d,)"
#factor="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_chrom_factor.txt | cut -f3 -d,)"

#echo $loc
echo $time
echo $chr

#Rscript scripts/trans_bv.R $time $chr 22
Rscript scripts/trans_select.R $time $chr 22