#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression/ASE
#SBATCH -J genofile
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=2
#SBATCH --mem 7G

module load python
module load bcftools
module load tabix
chr=$SLURM_ARRAY_TASK_ID

bcftools view -Oz -o Biogemma_WGS_Founders_transcript_variants_chr${chr}.vcf.gz Biogemma_WGS_Founders_transcript_variants.vcf.gz -r $chr

tabix -p vcf Biogemma_WGS_Founders_transcript_variants_chr${chr}.vcf.gz

python ../../scripts/old/foundergeno.py Biogemma_WGS_Founders_transcript_variants_chr${chr}.vcf.gz Biogemma_WGS_transcripts_chr${chr}.csv
