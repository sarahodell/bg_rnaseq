#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression/
#SBATCH -J pgs_sim
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=66-73
#SBATCH --ntasks=4
#SBATCH --mem=16G

module load R

Rscript scripts/sim_pgs.R $SLURM_ARRAY_TASK_ID

#chr=${bob[]}
#block=${sue[$SLURM_ARRAY_TASK_ID]}

#json="Biogemma_WGS_c${chr}_${block}.json"
#echo $json
#Rscript qtl2_array.R $SLURM_ARRAY_TASK_ID

#rep=$SLURM_ARRAY_TASK_ID
#chr=1
#chr=$SLURM_ARRAY_TASK_ID
#for chr in {1..10}; do
#  echo $chr
#  Rscript ../qtl2_array.R $chr $rep 8
#done
