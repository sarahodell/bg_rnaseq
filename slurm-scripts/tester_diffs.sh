#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J genofile
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=2
#SBATCH --mem 7G

module load R
chr=$SLURM_ARRAY_TASK_ID

Rscript scripts/tester_diffs.R $chr
