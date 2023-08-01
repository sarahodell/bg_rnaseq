#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J megalmm_resid
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 96:00:00
#SBATCH --array=0-11%3
#SBATCH --ntasks=60
#SBATCH --mem=60G

module load R

times=( "WD_0712" "WD_0712" "WD_0712" "WD_0718" "WD_0718" "WD_0718" "WD_0720" "WD_0720" "WD_0720" "WD_0727" "WD_0727" "WD_0727")
runs=( 1 2 3 1 2 3 1 2 3 1 2 3 )

time=${times[$SLURM_ARRAY_TASK_ID]}
run=${runs[$SLURM_ARRAY_TASK_ID]}
echo $time
echo $run

Rscript scripts/runMegaLMM_resids.R $time $run 59
