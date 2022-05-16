#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J deseq
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%j.txt
#SBATCH -t 5:00:00
#SBATCH --array=0-3
#SBATCH --ntasks=2
#SBATCH --mem 13G

module load java
module load jdk
module load R
module load samtools
module load qualimap

times=( "WD_0718" "WD_0720" "WD_0727" "WD_0712" )
time=${times[$SLURM_ARRAY_TASK_ID]}
echo $time
Rscript scripts/deseq.R $time
