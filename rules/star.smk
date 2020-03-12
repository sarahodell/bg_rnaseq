rule genome_index1:
    input:
        fa = "/group/jrigrp/Share/assemblies/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",
        gtf = "/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.gtf"
    output:
        "pass1/sjdbList.out.tab"
    params:
        pasdir = "pass1"
    threads: 16
    run:
        shell("if [ -d {params.pasdir} ]; then rm -rf {params.pasdir}; fi")
        shell("if [ -d _STARtmp ]; then rm -rf _STARtmp; fi")
        shell("mkdir {params.pasdir}")
        shell("STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params.pasdir} \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 100")
        shell("rm -rf _STARtmp")


rule pass1:
    input:
        R1L1 = 'trimmed/{sample}_R1_001.pe.qc.fastq.gz',
        R2L1 = 'trimmed/{sample}_R2_001.pe.qc.fastq.gz',
        sjdb = 'pass1/sjdbList.out.tab'
    output:
        '{sample}_pass1/SJ.out.tab'
    params:
        sampdir = '{sample}_pass1/',
        rmbam = '{sample}_pass1/Aligned.out.bam',
        passdir = 'pass1',
        tmpdir = '{sample}_tmp'
    threads: 16
    run:
        shell("if [ -d {params.sampdir} ]; then rm -rf {params.sampdir}; fi")
        shell("mkdir {params.sampdir}")
        shell("cd {params.sampdir}")
        shell("if [ -d {params.tmpdir} ]; then rm -rf {params.tmpdir}; fi")
        shell("STAR --runThreadN {threads} \
        --genomeDir {params.passdir} \
        --readFilesIn {input.R1L1} {input.R2L1} \
        --sjdbOverhang 100 \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --genomeLoad NoSharedMemory \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 1 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.sampdir} \
        --outTmpDir {params.tmpdir} \
        --outSAMtype BAM Unsorted")
        shell("cd ..")
        shell("if [ -f {params.rmbam} ]; then rm {params.rmbam}; fi")
        shell("rm -rf {params.tmpdir}")


rule genome_index2:
    input:
        fa = "/group/jrigrp/Share/assemblies/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",
        gtf = "/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.gtf",
        SJfiles = expand('{sample}_pass1/SJ.out.tab',sample=SAMPLES)
    output:
        "pass2/sjdbList.out.tab"
    params:
        pasdir = "pass2"
    threads: 16
    run:
        shell("if [ -d {params.pasdir} ]; then rm -rf {params.pasdir}; fi")
        shell("if [ -d _STARtmp ]; then rm -rf _STARtmp; fi")
        shell("mkdir {params.pasdir}")
        shell("STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params.pasdir} \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbFileChrStartEnd {input.SJfiles} \
        --sjdbOverhang 100")
        shell("rm -rf _STARtmp")


rule pass2:
    input:
        R1L1 = 'trimmed/{sample}_R1_001.pe.qc.fastq.gz',
        R2L1 = 'trimmed/{sample}_R2_001.pe.qc.fastq.gz',
        sjdb = 'pass2/sjdbList.out.tab',
        SJfiles = '{sample}_pass1/SJ.out.tab'
    params:
        outdir = '{sample}_pass2/',
        id = '{sample}',
        refdir = 'pass2',
        gtf = "/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.chr.gtf",
        tmpdir = '{sample}_tmp'
    output:
        '{sample}_pass2/Aligned.sortedByCoord.out.bam'
    threads: 16
    run:
        shell("if [ -d {params.outdir} ]; then rm -rf {params.outdir}; fi")
        shell("mkdir {params.outdir}")
        shell("cd {params.outdir}")
        shell("if [ -d {params.tmpdir} ]; then rm -rf {params.tmpdir}; fi")
        shell("STAR --runThreadN {threads} \
        --genomeDir {params.refdir} \
        --readFilesIn {input.R1L1} {input.R2L1} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile {params.gtf} \
        --sjdbFileChrStartEnd {input.SJfiles} \
        --outSAMattrRGline ID:{params.id}\tSM:{params.id}\tLB:None\tPL:Illumina \
        --quantMode TranscriptomeSAM GeneCounts \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 1 \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM 0 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --outSAMstrandField intronMotif \
        --outSAMattributes Standard \
        --outSAMunmapped Within \
        --outFileNamePrefix {params.outdir} \
        --outTmpDir {params.tmpdir} \
        --outSAMheaderHD @HD VN:1.4")
        shell("cd ..")
        shell("rm -rf {params.tmpdir}")
        shell("rm -rf _STARgenome")


rule mark_duplicates:
    output:
        "final_bams/{sample}.Aligned.sortedByCoord.MKDup.bam",
    input:
        "{sample}_pass2/Aligned.sortedByCoord.out.bam"
    params:
        prefix="{sample}_pass2/MkDup.",
        tmp="{sample}_pass2/MkDup.Processed.out.bam",
	    bamdir="final_bams/",
        id='{sample}'
    run:
        shell("if [ ! -d {params.bamdir}; then mkdir {params.bamdir}; fi")
        shell("if [ -d -STARtmp ]; then rm -rf _STARtmp; fi")
        shell("STAR \
	    --runMode inputAlignmentsFromBAM \
        --bamRemoveDuplicatesType UniqueIdentical \
        --readFilesCommand samtools view \
        --readFilesType SAM PE \
        --outSAMattrRGline ID:{params.id}\tSM:{params.id}\tLB:None\tPL:Illumina \
        --outSAMstrandField intronMotif \
        --outSAMattributes Standard \
        --outSAMunmapped Within \
        --inputBAMfile {input} \
        --outFileNamePrefix {params.prefix}")
        shell("mv {params.tmp} {output}")
        shell("rm -rf _STARtmp")


rule index_bam:
    output:
        "final_bams/{sample}.Aligned.sortedByCoord.MKDup.bam.bai"
    input:
        "final_bams/{sample}.Aligned.sortedByCoord.MKDup.bam"
    run:
        shell("samtools index {input}")
