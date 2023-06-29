#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J eqtl
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-3
#SBATCH --ntasks=32
#SBATCH --mem 32G

module load R

times=( "WD_0712" "WD_0718" "WD_0720" "WD_0727" ) # "WD_0712" )
time=${times[$SLURM_ARRAY_TASK_ID]}
#time="WD_0712"
echo $time

Rscript scripts/filter_h2.R $time 32
