#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J common
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 12:00:00
#SBATCH --array=1-10
#SBATCH --ntasks=4
#SBATCH --mem 32G


module load bcftools

chr=$SLURM_ARRAY_TASK_ID

# B73 is major
#bcftools view /group/jrigrp/Share/genotypes/MaizeHapMapV3.2.1/uplifted_APGv4/hmp321_agpv4_chr${chr}.vcf.gz -Q 0.05[:nref] -q 0.01[:nref] -O z -o datasets/hapmap/hapmap321_${chr}_rarealleles.vcf.gz
# B73 is minor
#bcftools view /group/jrigrp/Share/genotypes/MaizeHapMapV3.2.1/uplifted_APGv4/hmp321_agpv4_chr${chr}.vcf.gz -Q 0.90[:nref] -q 0.95[:nref] -O z -o datasets/hapmap/hapmap321_${chr}_rarealleles_B73minor.vcf.gz
bcftools view /group/jrigrp/Share/genotypes/MaizeHapMapV3.2.1/uplifted_APGv4/hmp321_agpv4_chr${chr}.vcf.gz -q 0.2[:nref] -O z -o datasets/hapmap/hapmap321_${chr}_commonalleles.vcf.gz

tabix -p vcf datasets/hapmap/hapmap321_${chr}_commonalleles.vcf.gz

founderfile=../genotypes/WGS/Biogemma_WGS_Founders_filtered_snps_diploid_chr.vcf.gz

bcftools isec -n +2 datasets/hapmap/hapmap321_${chr}_commonalleles.vcf.gz $founderfile > datasets/hapmap/hapmap321_${chr}_biogemma_common_alleles.txt
bcftools query -R datasets/hapmap/hapmap321_${chr}_biogemma_common_alleles.txt -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ../genotypes/WGS/Biogemma_WGS_Founders_filtered_snps_diploid_chr.vcf.gz -o datasets/hapmap/chr${chr}_founder_common_alleles.txt

# index rare alleles vcf
#tabix -p vcf datasets/hapmap/hapmap321_${chr}_rarealleles.vcf.gz
#bcftools index -f datasets/hapmap/hapmap321_${chr}_rarealleles_B73major.vcf.gz

#bcftools merge datasets/hapmap/hapmap321_${chr}_rarealleles_B73minor.vcf.gz datasets/hapmap/hapmap321_${chr}_rarealleles_B73major.vcf.gz --use-header datasets/hapmap/hapmap321_${chr}_rarealleles_B73major.vcf.gz | bcftools index | bcftools sort -Oz -o datasets/hapmap/hapmap321_${chr}_rarealleles.vcf.gz
#bcftools index -f datasets/hapmap/hapmap321_${chr}_commonalleles.vcf.gz


#bcftools index -f datasets/hapmap/hapmap321_${chr}_rarealleles.vcf.gz

# find intersection of variants between hapmap rare alleles and Biogemma founders

#bcftools isec -n +2 datasets/hapmap/hapmap321_${chr}_rarealleles.vcf.gz $founderfile > datasets/hapmap/hapmap321_${chr}_biogemma_rare_alleles.txt

# "A632_usa","A654_inra","B73_inra","B96","C103_inra","CO255_inra","D105_inra","DK63","EP1_inra","F492", "FV252_inra", "FV2_inra","MBS847","ND245","OH43_inra","VA85","W117_inra"

# Grab founder alleles at rare variants
#bcftools query -R datasets/hapmap/hapmap321_${chr}_biogemma_rare_alleles.txt -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ../genotypes/WGS/Biogemma_WGS_Founders_filtered_snps_diploid_chr.vcf.gz -o datasets/hapmap/chr${chr}_founder_rare_alleles.txt

#bcftools query -R datasets/hapmap/hapmap321_${chr}_biogemma_rare_alleles.txt -f '%ID\t%CHROM\t%POS\n' datasets/hapmap/hapmap321_${chr}_rarealleles.vcf.gz -o datasets/hapmap/chr${chr}_hapmap_rare_alleles.txt

# Calculate allele frequencies of hapmap rare alleles
#/home/sodell/bin/./plink --threads 4 --vcf datasets/hapmap/hapmap321_${chr}_rarealleles.vcf.gz --maf 0.009 --freq counts --out datasets/hapmap/hapmap321_${chr}



# Overlap of Tropical rare and HapMap common
#tropical=datasets/ames/AmesTropical_${chr}rarealleles_v4.vcf.gz
# Hapmap common alleles that are rare in tropical
#bcftools isec -n +2 datasets/hapmap/hapmap321_${chr}_commonalleles.vcf.gz $tropical > datasets/hapmap/hapmap321_${chr}_tropical_rare_alleles.txt
# Tropical rare alleles in biogemma founders
#founderfile=../genotypes/WGS/Biogemma_WGS_Founders_filtered_snps_diploid_chr.vcf.gz
#bcftools isec -n +2 $founderfile $tropical > datasets/overlap/biogemma_${chr}_tropical_rare_alleles.txt

# Get overlap of all three
# First, use sites of overlap with biogemma to get tropical rares
#bcftools query -R datasets/overlap/biogemma_${chr}_tropical_rare_alleles.txt -f '%ID\t%CHROM\t%POS\n' $founderfile -o datasets/overlap/chr${chr}_biogemma_tropical_rare_genotypes.txt

#bcftools query -R datasets/overlap/biogemma_trop_rare_hapmap_common.txt -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ../genotypes/WGS/Biogemma_WGS_Founders_filtered_snps_diploid_chr.vcf.gz -o datasets/overlap/biogemma_founder_tropicalrare_hmpcommon_genos.txt
