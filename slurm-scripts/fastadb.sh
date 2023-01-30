#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J pseudo
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%j.txt
#SBATCH -t 48:00:00
#SBATCH --ntasks=3
#SBATCH --mem 23G

module load bcftools
module load python

python scripts/fasta_db.py MBS847 ASE/Zea_mays_B73v4_exons.fa ASE/Biogemma_WGS_Founders_transcript_variants.vcf.gz ASE/Zea_mays_MBS847_exons.fa
