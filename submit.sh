#!/bin/bash

date=$(date "+%Y_%m_%d")
# file="snakemake_logs/${date}_bodem.log"
echo $date
# echo $file
module load conda/base
conda activate /home/sodell/.conda/envs/snakemake-tutorial

#module load fastqc
#module load deprecated/multiqc/bio3 
module load star
module load samtools

#module load java
#module load jdk
#module load bcftools
#module load qualimap/2.1.1
#
#module load HTSeq/0.9.1
#module load salmon
#module load mashmap
#module load bedtools

#--unlock \
snakemake --jobs 50 --use-conda \
--rerun-incomplete \
--latency-wait 300 \
--cluster-config submit.json \
--cluster "sbatch -A jrigrp -p {cluster.p} -o {cluster.o} --mem {cluster.mem} --time {cluster.time} --job-name {cluster.name}" # -p &>> $file
