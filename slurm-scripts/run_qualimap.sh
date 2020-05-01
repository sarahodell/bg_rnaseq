#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J rnaseqc
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --mem 3G

module load java
module load jdk
module load R
module load samtools
module load qualimap

sample="18048FL-06-01-01_S1_L001"
sample_path=star/results/bamfiles/${sample}.Aligned.sortedByCoord.out.bam

samtools index $sample_path

scripts/./align_qc.sh $sample
