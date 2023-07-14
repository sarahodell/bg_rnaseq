import csv
import os
import numpy
import pandas
import glob, sys
import random

configfile: "config.yaml"

DATA_PATH="raw_reads/"

#SAMPLE_TABLE=pandas.read_csv('metadata/updated_file_genotype_info_FIXED_genotypes_fastq.csv',sep=',')

#SAMPLE_TABLE=pandas.read_csv('metadata/updated_file_genotype_info_FIXED_genotypes_star.csv')

SAMPLE_TABLE=pandas.read_csv('metadata/BG_completed_sample_list_FIXED.txt',sep='\t')
#SAMPLES=SAMPLE_TABLE["sample_name"]
#SAMPLES=list(SAMPLES)
#SAMPLES=np.unique(SAMPLES)
#SAMPLE_TABLE["shortcut"] = SAMPLE_TABLE["seq_name"]

#SAMPLE_TABLE["shortcut"] = SAMPLE_TABLE["batch"] + "/" + SAMPLE_TABLE["sample_name"]

#drop_samples=["BG-P13-well-72_S170_L007","BG-P24-well-86_S182_L002","BG-P23-well-10_S10_L001",
#"BG-P12-well-32_S33_L006","BG-P12-well-32_S33_L006","18048FL-06-01-67_S67_L001",
#"18048FL-06-01-56_S56_L001","18048Fl-06-03-18_S210_L003","18048FL-06-02-89_S185_L002",
#"18048FL-06-02-88_S184_L002","18048FL-06-02-70_S166_L002","18048Fl-06-03-67_S259_L003",
#"18048FL-06-01-31_S31_L001","18048FL-06-01-14_S14_L001","18048FL-06-01-49_S49_L001",
#"18048FL-06-02-95_S191_L002","18048FL-06-02-84_S180_L002","18048Fl-06-03-29_S221_L003"
#"18048FL-06-01-07_S7_L001","18048Fl-06-03-40_S232_L003","18048Fl-06-03-50_S242_L003",
#"18048FL-06-02-52_S148_L002","18048Fl-06-03-54_S246_L003","18048FL-06-02-85_S181_L002",
#"18048Fl-06-03-02_S194_L003","18048FL-06-02-62_S158_L002","18048FL-06-02-60_S156_L002",
#"18048FL-06-02-67_S163_L002","18048FL-06-01-26_S26_L001","18048Fl-06-03-13_S205_L003",
#"18048Fl-06-03-79_S271_L003","18048Fl-06-03-92_S284_L003","18048FL-06-01-89_S89_L001"
#"18048Fl-06-03-54_S246_L003","18048FL-06-02-60_S156_L002","18048FL-06-02-67_S163_L002",
#"18048Fl-06-03-50_S242_L003","18048Fl-06-03-18_S210_L003","18048FL-06-02-85_S181_L002",
#"18048FL-06-02-70_S166_L002","18048FL-06-02-62_S158_L002","18048Fl-06-03-13_S205_L003",
#"18048FL-06-01-14_S14_L001","18048Fl-06-03-29_S221_L003","18048FL-06-02-89_S185_L002",
#"18048Fl-06-03-67_S259_L003","18048FL-06-01-56_S56_L001","18048FL-06-01-49_S49_L001",
#"18048FL-06-02-88_S184_L002","18048FL-06-02-84_S180_L002","18048FL-06-02-95_S191_L002",
#"18048FL-06-01-67_S67_L001","18048Fl-06-03-40_S232_L003","18048Fl-06-03-02_S194_L003",
#"18048FL-06-01-07_S7_L001","18048FL-06-01-31_S31_L001","18048FL-06-02-52_S148_L002",
#"18048Fl-06-03-79_S271_L003","18048Fl-06-03-92_S284_L003","18048FL-06-01-89_S89_L001"]

#SAMPLE_TABLE=SAMPLE_TABLE.loc[SAMPLE_TABLE['sample_name']!=drop_samples[0]]
#SAMPLE_TABLE=SAMPLE_TABLE.loc[SAMPLE_TABLE['sample_name']!=drop_samples[1]]
#SAMPLE_TABLE=SAMPLE_TABLE.loc[SAMPLE_TABLE['batch']!="batch_1"]
#SAMPLE_TABLE=SAMPLE_TABLE.loc[SAMPLE_TABLE['batch']!="batch_2"]
#temp for testing (only one batch)
#SAMPLE_TABLE=SAMPLE_TABLE.loc[SAMPLE_TABLE['batch']=="batch_1"]
#SAMPLE_TABLE=SAMPLE_TABLE.loc[SAMPLE_TABLE['sample_name']!=drop_samples[5]]
#SAMPLE_TABLE = SAMPLE_TABLE[~SAMPLE_TABLE.batch.isin(drop_samples)]
#SAMPLE_TABLE=SAMPLE_TABLE.tail(n=361)

SAMPLE_TABLE=SAMPLE_TABLE.loc[SAMPLE_TABLE['read']==1]

SAMPLES=SAMPLE_TABLE["sample_name"]
SAMPLES=list(SAMPLES)
SHORT_SAMPLES=[s[:-1] for s in SAMPLES]
SAMPLE_TABLE['short_sample']=SHORT_SAMPLES
SAMPLE_TABLE['slane']=SAMPLE_TABLE.lane.apply(str)
DIRECTORIES=list(SAMPLE_TABLE["batch"])
LANES=list(SAMPLE_TABLE["lane"])
#SLANES=list(R1_TABLE["slane"])

#UNIQUE_LANES=list(np.unique(SLANES))
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
#SAMPLES=SAMPLES[:10]
#DIRECTORIES=DIRECTORIES[1:10]
print("Running alignment on this many samples:")
print(len(SAMPLES))
#SAMPLES=SAMPLE_TABLE["sample_name"]

#IND=random.choices(range(0,(len(SAMPLES)-1)),k=10)
#SAMPLES=SAMPLES[IND]
#DIRECTORIES=DIRECTORIES[IND]

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
    #expand("raw_reads/{dir}/{sample}-R1-001.fastq.gz",zip, dir=DIRECTORIES, sample=SAMPLES),
    #expand("raw_reads/{dir}/{sample}-R2-001.fastq.gz",zip, dir=DIRECTORIES, sample=SAMPLES),
    #expand("/group/runciegrp/Data/Illumina/bg/trimmed/update/{batch}/{sample}-R1-001.pe.qc.fastq.gz",zip,batch=DIRECTORIES,sample=SAMPLES),
    #expand("/group/runciegrp/Data/Illumina/bg/trimmed/update/{batch}/{sample}-R2-001.pe.qc.fastq.gz",zip,batch=DIRECTORIES,sample=SAMPLES),
    #expand("qc/fastqc/update/{batch}/{sample}-R1-001.pe.qc_fastqc.zip",zip,batch=DIRECTORIES,sample=SAMPLES),
    #expand("qc/fastqc/update/{batch}/{sample}-R2-001.pe.qc_fastqc.zip",zip,batch=DIRECTORIES,sample=SAMPLES),
    #expand("qc/fastqc/update/{batch}/{sample}-R1-001.pe.qc_fastqc.html",zip,batch=DIRECTORIES,sample=SAMPLES),
    #expand("qc/fastqc/update/{batch}/{sample}-R2-001.pe.qc_fastqc.html",zip,batch=DIRECTORIES,sample=SAMPLES),
    #expand("qc/update/bg-{batch}-extra-multiqc.html",zip,batch=DIRECTORIES)
    #"pseudos/decoys.txt",
    #"founders_transcript_index/sa.bin",
    #expand("salmon_quant/{batch}/{sample}_transcripts_quant/quant.sf",batch=DIRECTORIES,sample=SAMPLES)
    #expand("salmon_quant/{sample}_transcripts_quant/quant.sf",sample=SAMPLES)
    'star/pass1/sjdbList.out.tab',
    expand('/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass1/SJ.out.tab',zip,batch=DIRECTORIES,sample=SAMPLES),
    "star/pass2/sjdbList.out.tab",
    expand('/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass2/Aligned.sortedByCoord.out.bam',zip,batch=DIRECTORIES,sample=SAMPLES),
    expand("/group/runciegrp/Data/Illumina/bg/final_bams/{batch}/{sample}.Aligned.sortedByCoord.MKDup.Processed.out.bam.bai",zip,batch=DIRECTORIES,sample=SAMPLES),
    expand('/home/sodell/projects/biogemma/expression/htseq/{batch}/{sample}_HTSeq_union_gtf_no_gene_ID.log',zip,batch=DIRECTORIES,sample=SAMPLES)#,
    #expand('/home/sodell/projects/biogemma/expression/htseq/update/{batch}/{sample}_HTSeq.csv',zip,batch=DIRECTORIES,sample=SAMPLES)
#    expand("qc/rnaseqc/{sample}-stats/qualimapReport.html",sample=SAMPLES),
#    expand("qc/bamqc/{sample}_stats/qualimapReport.html",sample=SAMPLES),
#    "qc/multisampleBamQcReport.html"
#    expand('{sample}_HTSeq_union_gff3_no_gene_ID.log', sample=SAMPLES),
#    expand('{sample}_HTSeq.csv', sample=SAMPLES),


#include: "rules/trimmomatic.smk"
#include: "rules/fastqc.smk"
include: "rules/star.smk"
#include: "rules/salmon.smk"
#include: "rules/qualimap.smk"
include: "rules/htseq.smk"
#include: "rules/kallisto.smk"
