import glob, sys
fastqc_input = glob.glob('test/*_L001_R?_001.fastq.gz')

fastqc_output = []
for filename in fastqc_input:
  filesplit=filename.split('/')[1]
  filesplit2=filesplit.split('.')[0]
  head=filesplit2.split('_')[0]+'_'+filesplit2.split('_')[1]
  rg=filesplit2.split('_')[3]
  new_filename = 'test/' + head + '_L001_' + rg + '_001' + '_fastqc.html'
  fastqc_output.append(new_filename)
  new_filename = 'test/' + head + '_L001_' + rg + '_001' + '.pe.qc_fastqc.html'
  fastqc_output.append(new_filename)
  
print('from these input files', fastqc_input, file=sys.stderr)
print('I constructed these output filenames', fastqc_output, file=sys.stderr)

rule all:
  input:
    "multiqc_report.html"
    
rule clean:
  shell:
    "rm -f {fastqc_output} multiqc_report.html"

rule fastqc_a_file:
  input:
    "test/{arglebarf}.fastq.gz"
  output:
    "test/{arglebarf}_fastqc.html",
    "test/{arglebarf}_fastqc.zip"
  shell:
    "fastqc {input}"

rule run_multiqc:
  input:
    fastqc_output
  output:
    "multiqc_report.html",
    directory("multiqc_data")
  shell:
    "multiqc test/"

rule trim_reads:
  input:
    "test/{filename}_L001_R1_001.fastq.gz",
    "test/{filename}_L001_R2_001.fastq.gz"
  output:
    "test/{filename}_L001_R1_001.pe.qc.fastq.gz",
    "test/{filename}_L001_R1_001.se.qc.fastq.gz",
    "test/{filename}_L001_R2_001.pe.qc.fastq.gz",
    "test/{filename}_L001_R2_001.se.qc.fastq.gz"
  shell:
    """java -jar /home/sodell/bin/trimmomatic.jar PE {input} {output} LEADING:2 TRAILING:2 \
      SLIDINGWINDOW:4:15 \
      MINLEN:25 \
      ILLUMINACLIP:/home/sodell/bin/trimmomatic/adapters/TruSeq2-PE.fa:2:40:14"""