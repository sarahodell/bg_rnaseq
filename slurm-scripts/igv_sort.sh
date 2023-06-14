#!/bin/bash -l
#SBATCH -D /home/sodell/bin/IGV_2.16.1
#SBATCH -J igvtools
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-78
#SBATCH --ntasks=1
#SBATCH --mem 2G


#WD_0712 78
#WD_0718 157
#WD_0720 246
#WD_0727 198
module load jdk
module load igv

# individual sample
time="WD_0712"
path="$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/sodell/projects/biogemma/expression/final_bams/${time}_path_samples.txt | cut -f1 -d,)"
sample="$(sed "${SLURM_ARRAY_TASK_ID}q;d"  /home/sodell/projects/biogemma/expression/final_bams/${time}_path_samples.txt | cut -f2 -d,)"
echo $time
echo $path
echo $sample


infile="/home/sodell/projects/biogemma/expression/final_bams/alignment_check/${time}_${sample}.bam"
outfile="/home/sodell/projects/biogemma/expression/final_bams/alignment_check/${time}_${sample}.sorted.bam"

java -Xmx1500m --module-path=lib @igv.args --module=org.igv/org.broad.igv.tools.IgvTools sort $infile $outfile
