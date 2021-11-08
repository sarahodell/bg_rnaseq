import csv
import os
import pandas as pd
import numpy as np
import glob, sys

configfile: "config.yaml"

DATA_PATH="raw_reads/"

SAMPLE_TABLE=pd.read_csv('metadata/BG_completed_sample_list.txt',sep='\t')
#SAMPLES=SAMPLE_TABLE["sample_name"]
#SAMPLES=list(SAMPLES)
#SAMPLES=np.unique(SAMPLES)
SAMPLE_TABLE["shortcut"] = SAMPLE_TABLE["batch"] + "/" + SAMPLE_TABLE["sample_name"]
drop_samples=["BG-P13-well-72_S170_L007","BG-P24-well-86_S182_L002"]

SAMPLE_TABLE=SAMPLE_TABLE.loc[SAMPLE_TABLE['sample_name']!=drop_samples[0]]
SAMPLE_TABLE=SAMPLE_TABLE.loc[SAMPLE_TABLE['sample_name']!=drop_samples[1]]
R1_TABLE=SAMPLE_TABLE.loc[SAMPLE_TABLE['read']==1]

SAMPLES=R1_TABLE["sample_name"]
SAMPLES=list(SAMPLES)
DIRECTORIES=list(R1_TABLE["batch"])

# 18048-85-01-03/BG-P13-well-72_S170_L007 only R2
#18048-85-lane12-16/BG-P24-well-86_S182_L002 only R1

#dir, file = glob_wildcards("raw_reads/{dir}/{file}_R*_001.fastq.qz")
#
#SAMPLE_FILES=glob.glob("raw_reads/batch_1/*_R1_001.fastq.gz")
#SAMPLES=[]
#for s in SAMPLE_FILES:
#    s=os.path.basename(s)
#    one=s.split('_')
#    two=('_').join(one[:3])
#    SAMPLES.append(two)

#LANES = [1,2,3]
#SAMPLES.remove("18048Fl-06-03-20_S212_L003")
#LANES = [7,8]
#TEST
#SAMPLES=["18048FL-06-01-01_S1_L001","18048FL-06-01-01_S1_L001","18048FL-06-01-02_S2_L001","18048FL-06-01-02_S2_L001"]
#LANES=[1]


print("Running alignment on these samples:")
print(SAMPLES)

#ALIGNMENT_TOOL="STAR"

#if ALIGNMENT_TOOL=="STAR":
#    outs=expand("finqal_bams/{sample}.Aligned.sortedByCoord.MKDup.bam",sample=SAMPLES)
#else if ALIGNMENT_TOOL="salmon":
#    outs=expand("{sample}_transcripts_quant/quant.sf")
#directories, files = glob_wildcards("raw_reads/{dir}/{file}_R*_001.fastq.gz")
#print(directories)
#files=np.unique(files)


rule all:
  input:
    #expand("raw_reads/{dir}/{sample}_R1_001.fastq.gz",zip,dir=DIRECTORIES,sample=SAMPLES),
    #expand("raw_reads/{dir}/{sample}_R2_001.fastq.gz",zip,dir=DIRECTORIES,sample=SAMPLES),
    #expand("raw_reads/{dir}/{sample}_R1_001.fastq.gz",zip, dir=DIRECTORIES, sample=SAMPLES),
    #expand("raw_reads/{dir}/{sample}_R2_001.fastq.gz",zip, dir=DIRECTORIES, sample=SAMPLES),
    expand("trimmed/{sample}_R1_001.pe.qc.fastq.gz",sample=SAMPLES),
    expand("trimmed/{sample}_R2_001.pe.qc.fastq.gz",sample=SAMPLES),
    #expand("qc/fastqc/{sample}_R1_001.pe.qc_fastqc.zip", sample = SAMPLES),
    #expand("qc/fastqc/{sample}_R2_001.pe.qc_fastqc.zip", sample = SAMPLES),
    #expand("qc/bg_{batch}_L00{lane}_multiqc.html",batch=BATCHES,lane=LANES),
    #directory("transcript_index"),
    #expand("salmon_quant/{sample}_transcripts_quant/quant.sf",sample=SAMPLES)
    #'pass1/sjdbList.out.tab',
    #expand('{sample}_pass1/SJ.out.tab', sample=SAMPLES),
    #expand('{sample}_pass2/Aligned.sortedByCoord.out.bam',sample=SAMPLES),
    #expand("final_bams/{sample}.Aligned.sortedByCoord.MKDup.bam.bai",sample=SAMPLES)
#    expand("qc/rnaseqc/{sample}_stats/qualimapReport.html",sample=SAMPLES),
#    expand("qc/bamqc/{sample}_stats/qualimapReport.html",sample=SAMPLES),
#    "qc/multisampleBamQcReport.html"
#    expand('{sample}_HTSeq_union_gff3_no_gene_ID.log', sample=SAMPLES),
#    expand('{sample}_HTSeq.csv', sample=SAMPLES),


include: "rules/trimmomatic.smk"
#include: "rules/fastqc.smk"
#include: "rules/star.smk"
#include: "rules/salmon.smk"
#include: "rules/qualimap.smk"
#include: "rules/htseq.smk"
