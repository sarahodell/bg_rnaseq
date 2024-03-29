#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J plink
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=16
#SBATCH --mem 64G

module load plink

#/home/sodell/bin/./plink --threads 5 --bfile /home/sodell/projects/biogemma/genotypes/plink_files/600K/Biogemma_DHLines_600K_Genotypes_binary --r2 with-freqs inter-chr --ld-snp-list MegaLMM/WD_0727_Factor2_snplist.txt --ld-window-r2 0.01 --out MegaLMM/Factor2_rsquared

#/home/sodell/bin/./plink --threads 16 --bfile /home/sodell/projects/biogemma/genotypes/plink_files/600K/Biogemma_DHLines_600K_Genotypes_binary --r2 with-freqs inter-chr --ld-snp-list eqtl/data/snplist.txt --ld-window-r2 0.01 --out eqtl/data/SNP_LD

founderfile=../genotypes/WGS/Biogemma_WGS_Founders_filtered_snps_diploid_no_tester.vcf.gz

/home/sodell/bin/./plink --threads 4 --vcf $founderfile --maf 0.01 --freq counts --out datasets/Biogemma_founders


#plink --bfile Biogemma_DHLines_600K_Genotypes_binary --freqx --nonfounders --out ../ld_decay/Biogemma_DHLine_allele_freq

#plink --threads 16 --bfile genotypes/plink_files/600K/Biogemma_DHLines_600K_Genotypes_binary --r2 with-freqs inter-chr --ld-window-r2 0.9 --out stats/ld_decay/Biogemma_DHLines_rsquared_all_chroms_r2_0.9

#if [ ! -d /scratch/sodell ]; then
#  mkdir /scratch/sodell;
#fi

#plink --threads 16 --bfile /home/sodell/projects/biogemma/genotypes/plink_files/600K/Biogemma_DHLines_600K_Genotypes_binary --r2 with-freqs --ld-snp-list snplist.txt --ld-window-r2 0 --ld-window 99999 --ld-window-kb 1000000 --out Biogemma_DHLines_rsquared

#awk '{print > $1"_rsquared.ld"}' Biogemma_DHLines_rsquared.ld
#awk '$1 == 8' Biogemma_DHLines_rsquared.ld > 8_rsquared.ld
#awk '$1 == 9' test.ld > 9_test.ld

#echo "Starting c9"
#plink --threads 7 --bfile /home/sodell/projects/biogemma/genotypes/plink_files/600K/Biogemma_DHLines_600K_Genotypes_binary --r2 with-freqs --ld-snp-list c9_snplist.txt --ld-window-r2 0 --ld-window 99999 --ld-window-kb 1000000 --out 9_rsquared
#echo "Starting c10"
#plink --threads 7 --bfile /home/sodell/projects/biogemma/genotypes/plink_files/600K/Biogemma_DHLines_600K_Genotypes_binary --r2 with-freqs --ld-snp-list c10_snplist.txt --ld-window-r2 0 --ld-window 99999 --ld-window-kb 1000000 --out 10_rsquared
