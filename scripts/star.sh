#!/bin/bash

#Set Paths
threads=16
sample="18048FL-06-01-01_S1_L001"

reference=/group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa
annotation=/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.chr.gtf
fastq_left=trimmed/${sample}_R1_001.qc.fastq.gz
fastq_right=trimmed/${sample}_R2_001.qc.fastq.gz

tmp_dir=star/tmp
output_path=star/genome_index/pass1
output_path2=star/genome_index/pass2
output_path3=star/results

if [ -d $tmp_dir ]
then
    rm -r $tmp_dir
fi

echo "Step 1: Building the STAR index."

STAR \
    --runMode genomeGenerate \
    --genomeDir $output_path \
    --genomeFastaFiles $reference \
    --sjdbOverhang 100 \
    --sjdbGTFfile $annotation \
    --runThreadN $threads \
    --outTmpDir $tmp_dir


echo "Step 2: Alignment 1st Pass."
if [ -d $tmp_dir ]
then
    rm -r $tmp_dir
fi

STAR \
    --genomeDir $output_path \
    --readFilesIn $fastq_left $fastq_right \
    --outFilterMultimapScoreRange 1 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 20 \
    --alignIntronMax 500000 \
    --alignMatesGapMax 1000000 \
    --sjdbScore 2 \
    --alignSJDBoverhangMin 1 \
    --genomeLoad NoSharedMemory \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --sjdbOverhang 100 \
    --outSAMstrandField intronMotif \
    --outSAMtype None \
    --outSAMmode None \
    --readFilesCommand zcat \
    --runThreadN $threads \
    --outTmpDir $tmp_dir

echo "Step 3: Intermediate Index Generation."
if [ -d $tmp_dir ]
then
    rm -r $tmp_dir
fi

#shopt -s nullglob
#sjfiles=($output_path/*.tab)

STAR \
    --runMode genomeGenerate \
    --genomeDir $output_path2 \
    --genomeFastaFiles $reference \
    --sjdbOverhang 100 \
    --runThreadN $threads \
    --sjdbFileChrStartEnd $output_path/sjdbList.out.tab \
    --outTmpDir $tmp_dir

echo "Step 4: Alignment 2nd Pass."

if [ -d $tmp_dir ]
then
    rm -r $tmp_dir
fi

STAR \
    --genomeDir $output_path2 \
    --quantMode TranscriptomeSAM GeneCounts \
    --readFilesIn $fastq_left $fastq_right \
    --runThreadN $threads \
    --outFilterMultimapScoreRange 1 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --alignIntronMax 500000 \
    --alignMatesGapMax 1000000 \
    --sjdbScore 2 \
    --alignSJDBoverhangMin 1 \
    --genomeLoad NoSharedMemory \
    --limitBAMsortRAM 0 \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --sjdbOverhang 100 \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI NM MD AS XS \
    --outSAMunmapped Within \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMheaderHD @HD VN:1.4 \
    --outSAMattrRGline ID:$sample
    --outTmpDir $tmp_dir \
    --outFileNamePrefix $output_path3/${sample}






