rule run_rna_multiqc:
  input:
    lambda wc: "qc/fastqc/"+R1_TABLE[R1_TABLE.batch == wc.batch].batch+"/"+R1_TABLE[R1_TABLE.batch == wc.batch].sample_name+"_R1_001.pe.qc_fastqc.zip",
    lambda wc: "qc/fastqc/"+R1_TABLE[R1_TABLE.batch == wc.batch].batch+"/"+R1_TABLE[R1_TABLE.batch == wc.batch].sample_name+"_R2_001.pe.qc_fastqc.zip"

  output:
    multiqc='qc/star/bg_{batch}star_multiqc.html'
  params:
    #indir="qc/fastqc/{batch}",
    outdir="qc/star",
    filename="bg_{batch}_star_multiqc.html"
  run:
    #shell("if [ -f bg_{batch}_L00{lane}_multiqc.html ]; then rm bg_{batch}_L00{lane}_multiqc.html; fi")
    #shell("if [ -d bg_{batch}_L00{lane}_multiqc ]; then rm -rf bg_{batch}_L00{lane}_multiqc; fi")
    #shell("if [ -f bg_batch_2_L008_multiqc.html ]; then rm bg_batch_2_L008_multiqc.html; fi")
    #shell("if [ -d bg_batch_2_L008_multiqc ]; then rm -rf bg_batch_2_L008_multiqc; fi")
    shell("multiqc {input} -o {params.outdir} -n {params.filename}")
    #shell("multiqc {params.indir}/*L008* -o {params.outdir} -n bg_batch_2_L008_multiqc.html")


#star:
#  fn: '*Log.final.out'
#star/genecounts:
#  fn: '*ReadsPerGene.out.tab'
