#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J ames
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/%A_%a.out
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/%A_%a.error
#SBATCH -t 24:00:00
#SBATCH --array=1-10
#SBATCH --mem 16G
#SBATCH --ntasks=4

module load jdk
#module load tassel
#module load vcftools

chr=$SLURM_ARRAY_TASK_ID


#/home/sodell/bin/tassel-5-standalone/run_pipeline.pl -FilterSiteBuilderPlugin -help
#-FilterTaxaPropertiesPlugin -minNotMissing 0.1 -endPlugin -export temp
#-exportType Hapmap Diploid

/home/sodell/bin/tassel-5-standalone/run_pipeline.pl -Xmx16g -importGuess datasets/ames/AmesTropical_imputed_chr${chr}.hmp.txt -FilterSiteBuilderPlugin -siteMaxAlleleFreq 0.05 -siteMinCount 10 -endPlugin -export datasets/ames/AmesTropical_${chr}rarealleles -exportType VCF

#/run_pipeline.pl -fork1 -h my_hapmap.hmp.txt -export -exportType VCF -runfork1

#Convert vcf file to plink format:

#vcftools --vcf datasets/ames/AmesTropical_${chr}rarealleles1.vcf --plink --out  datasets/ames/AmesTropical_${chr}rarealleles
