#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J samtools
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 12:00:00
#SBATCH --array=1-198
#SBATCH --ntasks=1
#SBATCH --mem 3G


#WD_0712 78
#WD_0718 157
#WD_0720 246
#WD_0727 198

module load samtools

time="WD_0727"
path="$(sed "${SLURM_ARRAY_TASK_ID}q;d" final_bams/${time}_path_samples.txt | cut -f1 -d,)"
sample="$(sed "${SLURM_ARRAY_TASK_ID}q;d"  final_bams/${time}_path_samples.txt | cut -f2 -d,)"

echo $time
echo $path
echo $sample

outfile="final_bams/alignment_check/${time}_${sample}.bam"
sortfile="final_bams/alignment_check/${time}_${sample}.sorted.bam"
bedfile="eqtl/results/${time}_genes.bed"

samtools view $path -h -b --region-file $bedfile > $outfile
samtools sort $outfile > $sortfile
samtools index $sortfile





