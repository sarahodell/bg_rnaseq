rule trim_reads:
  input:
    "raw_reads/batch_1/{filename}_L001_R1_001.fastq.gz",
    "raw_reads/batch_1/{filename}_L001_R2_001.fastq.gz"
  output:
    "trimmed/{filename}_L001_R1_001.pe.qc.fastq.gz",
    "trimmed/{filename}_L001_R1_001.se.qc.fastq.gz",
    "trimmed/{filename}_L001_R2_001.pe.qc.fastq.gz",
    "trimmed/{filename}_L001_R2_001.se.qc.fastq.gz"
  params:
    trimmomatic="/home/sodell/bin/trimmomatic.jar",
    adapters="/home/sodell/bin/trimmomatic/adapters/TruSeq2-PE.fa"
  run:
    shell("java -jar {params.trimmomatic} PE {input} {output} LEADING:2 TRAILING:2 \
      SLIDINGWINDOW:4:15 \
      MINLEN:25 \
      ILLUMINACLIP:{params.adapters}:2:40:14")