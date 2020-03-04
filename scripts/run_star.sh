#!/bin/bash

#Set Paths
threads=8
sample="18048FL-06-01-01_S1_L001"



reference=/group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa
annotation=/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.chr.gtf.gz
fastq_left=trimmed/${sample}_R1_001.qc.fastq.qz
fastq_right=trimmed/${sample}_R2_001.qc.fastq.qz

output_path=star/genome_index/pass1
output_path2=stare/genom_index/pass2
SJ_out=star/sj

### Step 1: Building the STAR index.*

STAR \
--runMode genomeGenerate \
--genomeDir $output_path \
--genomeFastaFiles $reference \
--sjdbGTFfile $annotation \
--runThreadN $threads

### Step 2: Alignment 1st Pass.

STAR \
--genomeDir $output_path \
--readFilesIn $fastq_left $fastq_right \
--runThreadN $threads \
--readFilesCommand zcat

### Step 3: Intermediate Index Generation.

STAR \
--runMode genomeGenerate \
--genomeDir $output_path2 \
--genomeFastaFiles $reference \
--runThreadN $threads \
--sjdbFileChrStartEnd $SJ_out/${sample}_sj1.tab $SJ_out/${sample}_sj2.tab

### Step 4: Alignment 2nd Pass.

STAR \
--genomeDir $output_path2 \
--readFilesIn $fastq_left $fastq_right \
--runThreadN $threads \
--readFilesCommand zcat


