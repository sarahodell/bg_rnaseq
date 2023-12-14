#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J epistasis
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 12:00:00
#SBATCH --array=1-100%20
#SBATCH --ntasks=8
#SBATCH --mem 16G

module load R

rep=$SLURM_ARRAY_TASK_ID

ENV=ALL
echo $ENV
echo $rep
for chr in {1..10}; do
	echo $chr
	Rscript ../scripts/epistasis_perm.R $chr $ENV $rep 7
done