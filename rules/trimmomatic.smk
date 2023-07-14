
rule trim_reads:
  input:
    lambda wc: "raw_reads/"+SAMPLE_TABLE[SAMPLE_TABLE.sample_name == wc.sample].batch+"/"+SAMPLE_TABLE[SAMPLE_TABLE.sample_name == wc.sample].sample_name+"_R1_001.fastq.gz",
    lambda wc: "raw_reads/"+SAMPLE_TABLE[SAMPLE_TABLE.sample_name == wc.sample].batch+"/"+SAMPLE_TABLE[SAMPLE_TABLE.sample_name == wc.sample].sample_name+"_R2_001.fastq.gz"

    #"raw_reads/{dir}/{sample}_R1_001.fastq.gz",
    #"raw_reads/{dir}/{sample}_R2_001.fastq.gz"
  output:
    "/group/runciegrp/Data/Illumina/bg/trimmed/update/{batch}/{sample}-R1-001.pe.qc.fastq.gz",
    "/group/runciegrp/Data/Illumina/bg/trimmed/update/{batch}/{sample}-R1-001.se.qc.fastq.gz",
    "/group/runciegrp/Data/Illumina/bg/trimmed/update/{batch}/{sample}-R2-001.pe.qc.fastq.gz",
    "/group/runciegrp/Data/Illumina/bg/trimmed/update/{batch}/{sample}-R2-001.se.qc.fastq.gz"
  params:
    trimmomatic="/home/sodell/bin/Trimmomatic-0.39/trimmomatic-0.39.jar",
    adapters="/home/sodell/bin/Trimmomatic-0.39/adapters/TruSeq2-PE.fa",
    threads=3
  run:
    shell("java -jar {params.trimmomatic} PE -threads {params.threads} {input} {output}  -phred33 LEADING:2 TRAILING:2 \
      SLIDINGWINDOW:4:15 \
      MINLEN:25 \
      ILLUMINACLIP:{params.adapters}:2:40:14")
