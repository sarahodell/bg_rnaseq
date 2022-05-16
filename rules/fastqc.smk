
rule run_fastqc:
  input:
     lambda wc: "/group/runciegrp/Data/Illumina/bg/trimmed/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].batch+"/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].sample_name+"_R1_001.pe.qc.fastq.gz",
     lambda wc: "/group/runciegrp/Data/Illumina/bg/trimmed/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].batch+"/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].sample_name+"_R2_001.pe.qc.fastq.gz"
    #"trimmed/{batch}/{sample}.pe.qc.fastq.gz"
  output:
    "qc/fastqc/{batch}/{sample}_R1_001.pe.qc_fastqc.html",
    "qc/fastqc/{batch}/{sample}_R1_001.pe.qc_fastqc.zip",
    "qc/fastqc/{batch}/{sample}_R2_001.pe.qc_fastqc.html",
    "qc/fastqc/{batch}/{sample}_R2_001.pe.qc_fastqc.zip"
  params:
    "qc/fastqc/{batch}/"
  run:
    shell("fastqc -o {params} --noextract {input}")

rule run_multiqc:
  input:
    lambda wc: "qc/fastqc/"+R1_TABLE[R1_TABLE.batch == wc.batch].batch+"/"+R1_TABLE[R1_TABLE.batch == wc.batch].sample_name+"_R1_001.pe.qc_fastqc.zip",
    lambda wc: "qc/fastqc/"+R1_TABLE[R1_TABLE.batch == wc.batch].batch+"/"+R1_TABLE[R1_TABLE.batch == wc.batch].sample_name+"_R2_001.pe.qc_fastqc.zip"
    #lambda wc: "qc/fastqc/"+R1_TABLE[R1_TABLE.batch == wc.batch].batch+"/"+R1_TABLE[R1_TABLE.batch == wc.batch].short_sample+R1_TABLE[R1_TABLE.slane == wc.lane].slane+"_R1_001.pe.qc_fastqc.zip",
    #lambda wc: "qc/fastqc/"+R1_TABLE[R1_TABLE.batch == wc.batch].batch+"/"+R1_TABLE[R1_TABLE.batch == wc.batch].short_sample+R1_TABLE[R1_TABLE.slane == wc.lane].slane+"_R2_001.pe.qc_fastqc.zip"
    #lambda wc: "qc/fastqc/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].batch+"/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].sample_name+"_R2_001.pe.qc.fastq.gz"
    #expand("qc/fastqc/{batch}/{sample}{lane}_R1_001.pe.qc_fastqc.zip", batch=DIRECTORIES,sample=SHORT_SAMPLES,lane=LANES,allow_missing=True),
    #expand("qc/fastqc/{batch}/{sample}{lane}_R2_001.pe.qc_fastqc.zip", batch=DIRECTORIES,sample=SHORT_SAMPLES,lane=LANES,allow_missing=True)
  output:
    multiqc='qc/bg_{batch}_multiqc.html'
  params:
    #indir="qc/fastqc/{batch}",
    outdir="qc/",
    filename="bg_{batch}_multiqc.html"
  run:
    #shell("if [ -f bg_{batch}_L00{lane}_multiqc.html ]; then rm bg_{batch}_L00{lane}_multiqc.html; fi")
    #shell("if [ -d bg_{batch}_L00{lane}_multiqc ]; then rm -rf bg_{batch}_L00{lane}_multiqc; fi")
    #shell("if [ -f bg_batch_2_L008_multiqc.html ]; then rm bg_batch_2_L008_multiqc.html; fi")
    #shell("if [ -d bg_batch_2_L008_multiqc ]; then rm -rf bg_batch_2_L008_multiqc; fi")
    shell("multiqc {input} -o {params.outdir} -n {params.filename}")
    #shell("multiqc {params.indir}/*L008* -o {params.outdir} -n bg_batch_2_L008_multiqc.html")
