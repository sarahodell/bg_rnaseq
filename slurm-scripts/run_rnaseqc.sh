#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J rnaseqc
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=3
#SBATCH --mem 23G

module load samtools

scripts/./align_qc.sh
