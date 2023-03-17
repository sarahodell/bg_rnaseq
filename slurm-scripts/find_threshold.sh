#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression/
#SBATCH -J glm_perm
#SBATCH -o /home/sodell/projects/biogemma/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=1-92
#SBATCH --ntasks=1
#SBATCH --mem=3G

module load R

time="WD_0712"
# Permuation 1000 times
factor="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/trans/${time}_factors.txt  | cut -f1 -d,)"

echo $time
echo $factor

Rscript scripts/find_threshold.R $time $factor 0.05

#echo "Finished array task $SLURM_ARRAY_TASK_ID" >> report.txt
