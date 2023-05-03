#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J factor
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 96:00:00
#SBATCH --array=1-15
#SBATCH --ntasks=1
#SBATCH --mem 4G

module load R/4.1.0

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"
#chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f2 -d,)"

time="WD_0712"

factor="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_factors.txt | cut -f1 -d,)"

#times=( "WD_0718" "WD_0720" "WD_0727" "WD_0712" )
#time=${times[$SLURM_ARRAY_TASK_ID]}
#f=$SLURM_ARRAY_TASK_ID
echo $time
echo $factor
Rscript scripts/test_factor_eQTL.R $time $factor 1 1
Rscript scripts/test_factor_eQTL.R $time $factor 2 1
Rscript scripts/test_factor_eQTL.R $time $factor 3 1
