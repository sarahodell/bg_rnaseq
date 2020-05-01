import csv
import os
import pandas as pd
import glob, sys

configfile: "config.yaml"

DATA_PATH="raw_reads/batch_1"

SAMPLE_TABLE=pd.read_csv('samples.tsv',sep='\t')
SAMPLES=SAMPLE_TABLE.loc[SAMPLE_TABLE['include']==True,'sample']
SAMPLES=list(SAMPLES)


#SAMPLE_FILES=glob.glob("raw_reads/batch_1/*_R1_001.fastq.gz")
#SAMPLES=[]
#for s in SAMPLE_FILES:
#    s=os.path.basename(s)
#    one=s.split('_')
#    two=('_').join(one[:3])
#    SAMPLES.append(two)

LANES = [1,2,3]

#TEST
#SAMPLES=["18048FL-06-01-01_S1_L001","18048FL-06-01-01_S1_L001","18048FL-06-01-02_S2_L001","18048FL-06-01-02_S2_L001"]
#LANES=[1]


print("Running alignment on these samples:")
print(SAMPLES)

#ALIGNMENT_TOOL="STAR"

#if ALIGNMENT_TOOL=="STAR":
#    outs=expand("final_bams/{sample}.Aligned.sortedByCoord.MKDup.bam",sample=SAMPLES)
#else if ALIGNMENT_TOOL="salmon":
#    outs=expand("{sample}_transcripts_quant/quant.sf")

rule all:
  input:
    expand("trimmed/{sample}_R1_001.pe.qc.fastq.gz",sample=SAMPLES),
    expand("trimmed/{sample}_R2_001.pe.qc.fastq.gz",sample=SAMPLES),
    expand("qc/fastqc/{sample}_R1_001.pe.qc_fastqc.zip", sample = SAMPLES),
    expand("qc/fastqc/{sample}_R2_001.pe.qc_fastqc.zip", sample = SAMPLES),
    expand("qc/bg_batch_1_L00{lane}_multiqc.html",lane=LANES),
    'pass1/sjdbList.out.tab',
    expand('{sample}_pass1/SJ.out.tab', sample=SAMPLES),
    expand('{sample}_pass2/Aligned.sortedByCoord.out.bam',sample=SAMPLES),
    expand("final_bams/{sample}.Aligned.sortedByCoord.MKDup.bam.bai",sample=SAMPLES)
#    expand("qc/rnaseqc/{sample}_stats/qualimapReport.html",sample=SAMPLES),
#    expand("qc/bamqc/{sample}_stats/qualimapReport.html",sample=SAMPLES),
#    "qc/multisampleBamQcReport.html"
#    expand('{sample}_HTSeq_union_gff3_no_gene_ID.log', sample=SAMPLES),
#    expand('{sample}_HTSeq.csv', sample=SAMPLES),


include: "rules/trimmomatic.smk"
include: "rules/fastqc.smk"
include: "rules/star.smk"
#include: "rules/salmon.smk"
#include: "rules/qualimap.smk"

#include: "rules/htseq.smk"
