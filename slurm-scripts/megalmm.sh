#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J megalmm
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-3
#SBATCH --ntasks=60
#SBATCH --mem=60G

module load R/4.1.0

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"
#chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f2 -d,)"

#times=( "WD_0718" "WD_0720" "WD_0727" "WD_0712" )
#time=${times[$SLURM_ARRAY_TASK_ID]}
#echo $time
#echo $chr
time="WD_0712"
run=$SLURM_ARRAY_TASK_ID
#run=4
Rscript scripts/runMegaLMM.R $time $run 59
