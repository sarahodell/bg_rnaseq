#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J crossmap
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH --array=1-9
#SBATCH --mem=4G
#SBATCH --ntasks=1

#module load conda
#conda activate snakemake-tutorial

chr=$SLURM_ARRAY_TASK_ID
#python /home/sodell/bin/liguowang-CrossMap-d485c9f/bin/CrossMap.py vcf datasets/AGPv3_to_B73_RefGen_v4.chain.gz datasets/ames/AmesTropical_${chr}rarealleles1.vcf /group/jrigrp/Share/assemblies/Zea_mays.B73_RefGen_v4.dna.toplevel.fa  datasets/ames/AmesTropical_${chr}rarealleles1_v4.vcf

module load bcftools
#bgzip datasets/ames/AmesTropical_${chr}rarealleles1_v4.vcf

bcftools reheader -f /group/jrigrp/Share/assemblies/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.fai datasets/ames/AmesTropical_${chr}rarealleles1_v4.vcf.gz -o datasets/ames/AmesTropical_${chr}rarealleles1_v4_reheader.vcf.gz
#tabix -p vcf datasets/ames/AmesTropical_${chr}rarealleles1_v4_reheader.vcf.gz

bcftools sort datasets/ames/AmesTropical_${chr}rarealleles1_v4_reheader.vcf.gz -Oz -o datasets/ames/AmesTropical_${chr}rarealleles1_v4_sorted.vcf.gz
tabix -p vcf datasets/ames/AmesTropical_${chr}rarealleles1_v4_sorted.vcf.gz
bcftools view datasets/ames/AmesTropical_${chr}rarealleles1_v4_sorted.vcf.gz -r $chr -Oz -o datasets/ames/AmesTropical_${chr}rarealleles_v4.vcf.gz
tabix -p vcf datasets/ames/AmesTropical_${chr}rarealleles_v4.vcf.gz

#rm datasets/ames/AmesTropical_${chr}rarealleles1_v4_reheader.vcf.gz
#rm datasets/ames/AmesTropical_${chr}rarealleles1_v4_sorted.vcf.gz
#bcftools reheader -f /group/jrigrp/Share/assemblies/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.fai datasets/ames/AmesTropical_${chr}rarealleles1_v4.vcf.gz | tabix -p vcf | bcftools sort | tabix -p vcf | bcftools view -r $chr -Oz -o datasets/ames/AmesTropical_${chr}rarealleles_v4.vcf.gz













