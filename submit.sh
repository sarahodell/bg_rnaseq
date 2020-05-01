#!/bin/bash

date=$(date "+%Y_%m_%d")
# file="snakemake_logs/${date}_bodem.log"
echo $date
# echo $file

module load conda3
module load fastqc
module load R
module load star
module load samtools
module load multiqc
module load java
module load jdk
module load bcftools
module load qualimap/2.1.1
module load python
module load HTSeq/0.9.1


snakemake --jobs 200 --use-conda \
--rerun-incomplete \
--latency-wait 60 \
--cluster-config submit.json \
--cluster "sbatch -A jrigrp -p {cluster.p} -o {cluster.o} --mem {cluster.mem} --time {cluster.time} --job-name {cluster.name}" # -p &>> $file
