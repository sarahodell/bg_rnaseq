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

#STAR \
#    --genomeDir $output_path2 \
#    --readFilesIn $fastq_left $fastq_right \
#    --runThreadN $threads \
#    --outFilterMultimapScoreRange 1 \
#    --outFilterMultimapNmax 20 \
#    --outFilterMismatchNmax 10 \
#    --alignIntronMax 500000 \
#    --alignMatesGapMax 1000000 \
#    --sjdbScore 2 \
#    --alignSJDBoverhangMin 1 \
#    --genomeLoad NoSharedMemory \
#    --limitBAMsortRAM 0 \
#    --outFilterMatchNminOverLread 0.33 \
#    --outFilterScoreMinOverLread 0.33 \
#    --sjdbOverhang 100 \
#    --outSAMstrandField intronMotif \
#    --outSAMattributes NH HI NM MD AS XS \
#    --outSAMunmapped Within \
#    --readFilesCommand zcat \
#    --outSAMtype BAM SortedByCoordinate \
#    --outSAMheaderHD @HD VN:1.4 \
#    --outSAMattrRGline ID:$sample \
#    --outFileNamePrefix $output_path3/${sample}
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

#STAR \
#    --runMode inputAlignmentsFromBAM \
#    --bamRemoveDuplicatesType UniqueIdentical \
#    --readFilesCommand samtools view \
#    --readFilesType SAM PE \
#    --inputBAMfile $output_path3/${sample}.Aligned.sortedByCoord.out.bam \
#    --outFileNamePrefix $output_path3/${sample}.sortedByCoord.MkDup.





