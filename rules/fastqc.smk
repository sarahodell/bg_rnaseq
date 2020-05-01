
rule run_fastqc:
  input:
    "trimmed/{sample}.pe.qc.fastq.gz"
  output:
    touch("qc/fastqc/{sample}.pe.qc_fastqc.html"),
    touch("qc/fastqc/{sample}.pe.qc_fastqc.zip")
  params:
    "qc/fastqc/"
  run:
    shell("fastqc -o {params} --noextract {input}")

rule run_multiqc:
  input:
    R1=expand("qc/fastqc/{sample}_R1_001.pe.qc_fastqc.zip", sample = SAMPLES),
    R2=expand("qc/fastqc/{sample}_R2_001.pe.qc_fastqc.zip", sample = SAMPLES)
  output:
    multiqc=expand('qc/bg_batch_1_L00{lane}_multiqc.html',lane=LANES)
  params:
    indir="qc/fastqc",
    outdir="qc/"
  run:
    shell("if [ -f bg_batch_1_L001_multiqc.html ]; then rm bg_batch_1_L001_multiqc.html; fi")
    shell("if [ -d bg_batch_1_L001_multiqc ]; then rm -rf bg_batch_1_L001_multiqc; fi")
    shell("if [ -f bg_batch_1_L002_multiqc.html ]; then rm bg_batch_1_L002_multiqc.html; fi")
    shell("if [ -d bg_batch_1_L002_multiqc ]; then rm -rf bg_batch_1_L002_multiqc; fi")
    shell("if [ -f bg_batch_1_L003_multiqc.html ]; then rm bg_batch_1_L003_multiqc.html; fi")
    shell("if [ -d bg_batch_1_L003_multiqc ]; then rm -rf bg_batch_1_L003_multiqc; fi")
    shell("multiqc {params.indir}/*L001* -o {params.outdir} -n bg_batch_1_L001_multiqc.html")
    shell("multiqc {params.indir}/*L002* -o {params.outdir} -n bg_batch_1_L002_multiqc.html")
    shell("multiqc {params.indir}/*L003* -o {params.outdir} -n bg_batch_1_L003_multiqc.html")
