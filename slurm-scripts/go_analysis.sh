#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J goseq
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A.txt
#SBATCH -t 96:00:00
#SBATCH --ntasks=2
#SBATCH --mem 7G

module load R

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"
#chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f2 -d,)"

#times=( "WD_0712" "WD_0718" "WD_0720" "WD_0727" )
#time=${times[$SLURM_ARRAY_TASK_ID]}
time="WD_0727"
echo $time

Rscript scripts/go_analysis.R $time
