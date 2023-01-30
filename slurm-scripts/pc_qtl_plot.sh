#!/bin/bash -l
#SBATCH -D /home/sodell/projects/biogemma/expression
#SBATCH -J pc_plot
#SBATCH -o /home/sodell/projects/biogemma/expression/slurm-logs/out-%A_%a.txt
#SBATCH -e /home/sodell/projects/biogemma/expression/slurm-logs/error-%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --mem 3G

module load R

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"
#chr="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f2 -d,)"

times=( "WD_0718" "WD_0720" "WD_0727" "WD_0712" )
time=${times[$SLURM_ARRAY_TASK_ID]}

#time="$(sed "${SLURM_ARRAY_TASK_ID}q;d" eqtl/eqtl_time_chrom.txt | cut -f1 -d,)"


#time="WD_0712"
echo $time
#f="$(sed "${SLURM_ARRAY_TASK_ID}q;d" 'eqtl/pheno_residuals_trans_eQTL_list.txt' | cut -f2 -d,)"
#f=$SLURM_ARRAY_TASK_ID
#echo $f
pcs=( "PC1" "PC2" "PC3")
for pc in "${pcs[@]}"; do
  Rscript scripts/pc_qtl_plots.R $time $pc
done
