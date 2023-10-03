#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J bimbam
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-9
#SBATCH --ntasks=4
#SBATCH --mem 16G

module load R
chr=$SLURM_ARRAY_TASK_ID


Rscript scripts/founder_rarealleles.R $chr
#gunzip ../genotypes/probabilities/allele_probs/bg${chr}_wgs_alleleprobs.txt.gz | grep -f datasets/overlap/tropical_rares_markers.txt > datasets/overlap/chr${chr}_tropical_alleleprobs.txt 