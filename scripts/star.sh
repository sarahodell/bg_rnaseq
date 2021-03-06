#!/bin/bash

#Set Paths
threads=16
sample1="18048FL-06-01-01_S1_L001"
sample2="18048FL-06-01-02_S2_L001"
reference=/group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa
annotation=/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.chr.gtf
fastq_left1=/home/sodell/projects/biogemma/expression/trimmed/${sample1}_R1_001.pe.qc.fastq.gz
fastq_right1=/home/sodell/projects/biogemma/expression/trimmed/${sample1}_R2_001.pe.qc.fastq.gz
fastq_left2=/home/sodell/projects/biogemma/expression/trimmed/${sample2}_R1_001.pe.qc.fastq.gz
fastq_right2=/home/sodell/projects/biogemma/expression/trimmed/${sample2}_R2_001.pe.qc.fastq.gz
tmp_dir=/home/sodell/projects/biogemma/expression/star/tmp
output_path=/home/sodell/projects/biogemma/expression/star/results/pass1
output_path2=/home/sodell/projects/biogemma/expression/star/results/pass2
output_path3=/home/sodell/projects/biogemma/expression/star/results/bamfiles

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
    --readFilesIn $fastq_left1,$fastq_left2 $fastq_right1,$fastq_right2 \
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
    --sjdbGTFfile $annotation \
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
    --readFilesIn $fastq_left1,$fastq_left2 $fastq_right1,$fastq_right2 \
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
    --outSAMattrRGline ID:$sample \
    --outFileNamePrefix $output_path3/${sample1}. $output_path3/${sample2}.

echo "Step 5: Mark Duplicates"

STAR \
    --runMode inputAlignmentsFromBAM \
    --bamRemoveDuplicatesType UniqueIdentical \
    --readFilesCommand samtools view \
    --readFilesType SAM PE \
    --inputBAMfile $output_path3/${sample1}.Aligned.sortedByCoord.out.bam $output_path3/${sample2}.Aligned.sortedByCoord.out.bam \
    --outFileNamePrefix $output_path3/${sample1}.sortedByCoord.MkDup. $output_path3/${sample1}.sortedByCoord.MkDup.


## Clear out tmp directories
echo "Step 6: Cleaning Up..."
if [ -d $tmp_dir ]
then
    rm -r $tmp_dir
fi


rm -r $output_path3/*STAR_tmp
if [ -d STAR_tmp ]
then
    rm -r STAR_tmp

fi




