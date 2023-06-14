#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J eqtl
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH --array=0-3
#SBATCH --ntasks=16
#SBATCH --mem 128G

module load R

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"
#chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f2 -d,)"

times=( "WD_0718" "WD_0727" "WD_0720" "WD_0712" )
time=${times[$SLURM_ARRAY_TASK_ID]}

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"



echo $time

Rscript scripts/trans_weighted_plots.R $time
