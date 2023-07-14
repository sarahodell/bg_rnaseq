
rule run_fastqc:
  input:
     lambda wc: "/group/runciegrp/Data/Illumina/bg/trimmed/update/"+SAMPLE_TABLE[SAMPLE_TABLE.sample_name == wc.sample].batch+"/"+SAMPLE_TABLE[SAMPLE_TABLE.sample_name == wc.sample].sample_name+"-R1-001.pe.qc.fastq.gz",
     lambda wc: "/group/runciegrp/Data/Illumina/bg/trimmed/update/"+SAMPLE_TABLE[SAMPLE_TABLE.sample_name == wc.sample].batch+"/"+SAMPLE_TABLE[SAMPLE_TABLE.sample_name == wc.sample].sample_name+"-R2-001.pe.qc.fastq.gz"
    #"trimmed/{batch}/{sample}.pe.qc.fastq.gz"
  output:
    "qc/fastqc/update/{batch}/{sample}-R1-001.pe.qc_fastqc.html",
    "qc/fastqc/update/{batch}/{sample}-R1-001.pe.qc_fastqc.zip",
    "qc/fastqc/update/{batch}/{sample}-R2-001.pe.qc_fastqc.html",
    "qc/fastqc/update/{batch}/{sample}-R2-001.pe.qc_fastqc.zip"
  params:
    "qc/fastqc/update/{batch}/"
  run:
    shell("fastqc -o {params} --noextract {input}")

rule run_multiqc:
  input:
    lambda wc: "qc/fastqc/update/"+SAMPLE_TABLE[SAMPLE_TABLE.batch == wc.batch].batch+"/"+SAMPLE_TABLE[SAMPLE_TABLE.batch == wc.batch].sample_name+"-R1-001.pe.qc_fastqc.zip",
    lambda wc: "qc/fastqc/update/"+SAMPLE_TABLE[SAMPLE_TABLE.batch == wc.batch].batch+"/"+SAMPLE_TABLE[SAMPLE_TABLE.batch == wc.batch].sample_name+"-R2-001.pe.qc_fastqc.zip"
    #lambda wc: "qc/fastqc/update/"+SAMPLE_TABLE[SAMPLE_TABLE.batch == wc.batch].batch+"/"+SAMPLE_TABLE[SAMPLE_TABLE.batch == wc.batch].short_sample+SAMPLE_TABLE[SAMPLE_TABLE.slane == wc.lane].slane+"_R1_001.pe.qc_fastqc.zip",
    #lambda wc: "qc/fastqc/update/"+SAMPLE_TABLE[SAMPLE_TABLE.batch == wc.batch].batch+"/"+SAMPLE_TABLE[SAMPLE_TABLE.batch == wc.batch].short_sample+SAMPLE_TABLE[SAMPLE_TABLE.slane == wc.lane].slane+"_R2_001.pe.qc_fastqc.zip"
    #lambda wc: "qc/fastqc/update/"+SAMPLE_TABLE[SAMPLE_TABLE.sample_name == wc.sample].batch+"/"+SAMPLE_TABLE[SAMPLE_TABLE.sample_name == wc.sample].sample_name+"_R2_001.pe.qc.fastq.gz"
    #expand("qc/fastqc/update/{batch}/{sample}{lane}_R1_001.pe.qc_fastqc.zip", batch=DIRECTORIES,sample=SHORT_SAMPLES,lane=LANES,allow_missing=True),
    #expand("qc/fastqc/update/{batch}/{sample}{lane}_R2_001.pe.qc_fastqc.zip", batch=DIRECTORIES,sample=SHORT_SAMPLES,lane=LANES,allow_missing=True)
  output:
    multiqc='qc/update/bg-{batch}-extra-multiqc.html'
  params:
    #indir="qc/fastqc/update/{batch}",
    outdir="qc/update/",
    filename="bg-{batch}-extra-multiqc.html"
  run:
    #shell("if [ -f bg_{batch}_L00{lane}_multiqc.html ]; then rm bg_{batch}_L00{lane}_multiqc.html; fi")
    #shell("if [ -d bg_{batch}_L00{lane}_multiqc ]; then rm -rf bg_{batch}_L00{lane}_multiqc; fi")
    #shell("if [ -f bg_batch_2_L008_multiqc.html ]; then rm bg_batch_2_L008_multiqc.html; fi")
    #shell("if [ -d bg_batch_2_L008_multiqc ]; then rm -rf bg_batch_2_L008_multiqc; fi")
    shell("multiqc {input} -o {params.outdir} -n {params.filename}")
    #shell("multiqc {params.indir}/*L008* -o {params.outdir} -n bg_batch_2_L008_multiqc.html")
