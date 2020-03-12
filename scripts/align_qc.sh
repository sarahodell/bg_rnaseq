#!/bin/bash

sample=$1
#sample="18048FL-06-01-01_S1_L001"
sample_path=star/results/bamfiles/${sample}.sortedByCoord.MkDup.Processed.out.bam
qualimap=/share/apps/qualimap-2.1.1/./qualimap
reference=/group/jrigrp/Share/assembles/Zea_mays.AGPv4.dna.chr.f
gtf=/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.chr.gtf
output=/home/sodell/projects/biogemma/expression/star/qc
#threads=8

samtools index $sample_path

$qualimap rnaseq -outdir $output -pe -p strand-specific-forward -s -bam $sample_path -gtf $gtf
