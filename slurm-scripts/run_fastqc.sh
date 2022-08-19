#!/bin/bash

trimmomatic=/home/sodell/bin/trimmomatic
fastq_path=/home/sodell/projects/biogemma/expression/raw_reads/batch_1/
samples=("18048FL-06-01-01_S1_L001")
trim_path=/home/sodell/projects/biogemma/expression/trimmed/
qc_path=/home/sodell/projects/biogemma/expression/qc/

### Step 1: Running Trimmomatic

for s in ${samples[@]};	do
    R1=${fastq_path}${samples}_R1_001.fastq.gz
    R2=${fastq_path}${samples}_R2_001.fastq.gz
    java -jar $trimmomatic.jar PE $R1 $R2 ${trim_path}${s}_R1_001.qc.fastq.gz ${trim_path}orphans_1 ${trim_path}${s}_R2_001.qc.fastq.gz ${trim_path}orphans_2 ILLUMINACLIP:$trimmomatic/adapters/TruSeq2-PE.fa:2:40:14

done
### Step 2: Run FastQC

for s in ${samples[@]}; do
    R1=${fastq_path}${samples}_R1_001.fastq.gz
    R2=${fastq_path}${samples}_R2_001.fastq.gz
    fastqc $R1 -o $qc_path 
    fastqc $R2 -o $qc_path
done
    
