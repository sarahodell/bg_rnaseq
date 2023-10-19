#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J rares
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 12:00:00
#SBATCH --array=1-40%10
#SBATCH --ntasks=8
#SBATCH --mem 16G

module load R
#chr=$SLURM_ARRAY_TASK_ID

time1="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"
chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f2 -d,)"

echo $time1
echo $chr

#Rscript scripts/rare_count.R $chr $time1 7
#Rscript scripts/rare_count2.R $chr $time1


Rscript scripts/gerp_rare_count.R $chr $time1 7

#gunzip ../genotypes/probabilities/allele_probs/bg${chr}_wgs_alleleprobs.txt.gz | grep -f datasets/overlap/tropical_rares_markers.txt > datasets/overlap/chr${chr}_tropical_alleleprobs.txt 