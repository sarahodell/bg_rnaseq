rule trim_reads:
  input:
    "raw_reads/batch_1/{filename}_R1_001.fastq.gz",
    "raw_reads/batch_1/{filename}_R2_001.fastq.gz"
  output:
    "trimmed/{filename}_R1_001.pe.qc.fastq.gz",
    "trimmed/{filename}_R1_001.se.qc.fastq.gz",
    "trimmed/{filename}_R2_001.pe.qc.fastq.gz",
    "trimmed/{filename}_R2_001.se.qc.fastq.gz"
  params:
    trimmomatic="/home/sodell/bin/Trimmomatic-0.39/trimmomatic-0.39.jar",
    adapters="/home/sodell/bin/Trimmomatic-0.39/adapters/TruSeq2-PE.fa"
  run:
    shell("java -jar {params.trimmomatic} PE {input} {output}  -phred33 LEADING:2 TRAILING:2 \
      SLIDINGWINDOW:4:15 \
      MINLEN:25 \
      ILLUMINACLIP:{params.adapters}:2:40:14")
