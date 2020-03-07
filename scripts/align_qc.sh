#!/bin/bash

sample="18048FL-06-01-01_S1_L001"
sample_path=star/results/${sample}.bam
rnaseqc=/home/sodell/bin/rnaseqc/RNA-SeQC.jar
bamfile="${sample}|${sample_path}|testrun"
reference=/group/jrigrp/Share/assembles/Zea_mays.AGPv4.dna.chr.f
gtf=/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.chr.gtf
output=/star/qc

samtools index $sample_path

java -jar $rnaseqc -o $output -r $reference -s $bamfile -t $gtf
