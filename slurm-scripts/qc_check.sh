#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J mds
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 96:00:00
#SBATCH --array=1-4
#SBATCH --ntasks=3
#SBATCH --mem 23G

#module load java
#module load jdk
module load R
#module load samtools
#module load qualimap

#times=( 1 2 3 4 )
i=$SLURM_ARRAY_TASK_ID
echo $i
Rscript qc_check.R $i
