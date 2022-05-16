
rule trim_reads:
  input:
    lambda wc: "raw_reads/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].batch+"/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].sample_name+"_R1_001.fastq.gz",
    lambda wc: "raw_reads/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].batch+"/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].sample_name+"_R2_001.fastq.gz"

    #"raw_reads/{dir}/{sample}_R1_001.fastq.gz",
    #"raw_reads/{dir}/{sample}_R2_001.fastq.gz"
  output:
    "/group/runciegrp/Data/Illumina/bg/trimmed/{batch}/{sample}_R1_001.pe.qc.fastq.gz",
    "/group/runciegrp/Data/Illumina/bg/trimmed/{batch}/{sample}_R1_001.se.qc.fastq.gz",
    "/group/runciegrp/Data/Illumina/bg/trimmed/{batch}/{sample}_R2_001.pe.qc.fastq.gz",
    "/group/runciegrp/Data/Illumina/bg/trimmed/{batch}/{sample}_R2_001.se.qc.fastq.gz"
  params:
    trimmomatic="/home/sodell/bin/Trimmomatic-0.39/trimmomatic-0.39.jar",
    adapters="/home/sodell/bin/Trimmomatic-0.39/adapters/TruSeq2-PE.fa",
    threads=3
  run:
    shell("java -jar {params.trimmomatic} PE -threads {params.threads} {input} {output}  -phred33 LEADING:2 TRAILING:2 \
      SLIDINGWINDOW:4:15 \
      MINLEN:25 \
      ILLUMINACLIP:{params.adapters}:2:40:14")
