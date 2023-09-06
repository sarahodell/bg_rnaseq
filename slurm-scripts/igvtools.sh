#!/bin/bash -l
#SBATCH -D /home/sodell/bin/IGV_2.16.1
#SBATCH -J igvtools
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH --array=1-26
#SBATCH --ntasks=1
#SBATCH --mem 2G


#WD_0712 78
#WD_0718 157
#WD_0720 246
#WD_0727 198
module load jdk
module load igv

#pwd

chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/sodell/projects/biogemma/expression/eqtl/results/all_eQTL_genes_FIXED.bed | cut -f1 -d$'\t')"
start="$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/sodell/projects/biogemma/expression/eqtl/results/all_eQTL_genes_FIXED.bed | cut -f2 -d$'\t')"
end="$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/sodell/projects/biogemma/expression/eqtl/results/all_eQTL_genes_FIXED.bed | cut -f3 -d$'\t')"
gene="$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/sodell/projects/biogemma/expression/eqtl/results/all_eQTL_genes_FIXED.bed | cut -f4 -d$'\t')"
time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/sodell/projects/biogemma/expression/eqtl/results/all_eQTL_genes_FIXED.bed | cut -f5 -d$'\t')"

genome="/group/jrigrp/Share/assemblies/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
query="${chr}:${start}-${end}"

echo $time
echo $gene
echo $query


# using a bamlist as input
founders=( "B73_inra" "A632_usa" "CO255_inra" "FV252_inra" "OH43_inra" "A654_inra" "C103_inra" "FV2_inra" "EP1_inra" "D105_inra" "W117_inra" "B96" "DK63" "F492" "ND245" "VA85" )
for founder in "${founders[@]}"; do
	echo $founder
	infile="/home/sodell/projects/biogemma/expression/final_bams/alignment_check/${time}_${gene}_${founder}.bam.list"
	outfile="/home/sodell/projects/biogemma/expression/final_bams/alignment_check/${time}_${gene}_${founder}.tdf"
	if [ -f "$infile" ]; then
		java -Xmx1500m --module-path=lib @igv.args --module=org.igv/org.broad.igv.tools.IgvTools count -w 10 --extFactor 250 --query $query $infile $outfile $genome
	fi
done

#/home/sodell/projects/biogemma/expression/final_bams/alignment_check/final_bams/batch_1/18048FL-06-02-84_S180_L002.Aligned.sortedByCoord.MKDup.Processed.out.bam

# individual sample
#path="$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/sodell/projects/biogemma/expression/final_bams/${time}_path_samples.txt | cut -f1 -d,)"
#sample="$(sed "${SLURM_ARRAY_TASK_ID}q;d"  /home/sodell/projects/biogemma/expression/final_bams/${time}_path_samples.txt | cut -f2 -d,)"
#echo $time
#echo $path
#echo $sample

#infile="/home/sodell/projects/biogemma/expression/final_bams/alignment_check/${time}_${sample}.bam"
#outfile="/home/sodell/projects/biogemma/expression/final_bams/alignment_check/${time}_${sample}.tdf"

#java -Xmx1500m --module-path=lib @igv.args --module=org.igv/org.broad.igv.tools.IgvTools count --extFactor 250 $infile $outfile $genome