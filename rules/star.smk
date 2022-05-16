rule genome_index1:
    input:
        fa = "/group/jrigrp/Share/assemblies/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",
        gtf = "/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.gtf"
    output:
        "star/pass1/sjdbList.out.tab"
    params:
        pasdir = "star/pass1"
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
        #lambda wc: "trimmed/"+R1_TABLE[R1_TABLE.batch == wc.batch].batch+"/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].sample_name+"_R1_001.pe.qc.fastq.gz",
        #lambda wc: "trimmed/"+R1_TABLE[R1_TABLE.batch == wc.batch].batch+"/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].sample_name+"_R2_001.pe.qc.fastq.gz"
        '/group/runciegrp/Data/Illumina/bg/trimmed/{batch}/{sample}_R1_001.pe.qc.fastq.gz',
        '/group/runciegrp/Data/Illumina/bg/trimmed/{batch}/{sample}_R2_001.pe.qc.fastq.gz'
        #R2L1 = 'trimmed/{batch}/{sample}_R2_001.pe.qc.fastq.gz',
    output:
        '/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass1/SJ.out.tab'
    params:
        sampdir = '/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass1/',
        rmbam = '/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass1/Aligned.out.bam',
        passdir = 'star/pass1',
        tmpdir = '/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_tmp',
        sjdb = 'star/pass1/sjdbList.out.tab',
    threads: 5
    run:
        shell("if [ -d {params.sampdir} ]; then rm -rf {params.sampdir}; fi")
        shell("mkdir {params.sampdir}")
        shell("cd {params.sampdir}")
        shell("if [ -d {params.tmpdir} ]; then rm -rf {params.tmpdir}; fi")
        shell("STAR --runThreadN {threads} \
        --genomeDir {params.passdir} \
        --readFilesIn {input} \
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
        SJfiles = expand('/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass1/SJ.out.tab',zip,batch=DIRECTORIES,sample=SAMPLES)
    output:
        "star/pass2/sjdbList.out.tab"
    params:
        pasdir = "star/pass2"
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
        #lambda wc: "trimmed/"+R1_TABLE[R1_TABLE.batch == wc.batch].batch+"/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].sample_name+"_R1_001.pe.qc.fastq.gz",
        #lambda wc: "trimmed/"+R1_TABLE[R1_TABLE.batch == wc.batch].batch+"/"+R1_TABLE[R1_TABLE.sample_name == wc.sample].sample_name+"_R2_001.pe.qc.fastq.gz"
        #'trimmed/{batch}/{sample}_R1_001.pe.qc.fastq.gz',
        #'trimmed/{batch}/{sample}_R2_001.pe.qc.fastq.gz'
        R1L1 = '/group/runciegrp/Data/Illumina/bg/trimmed/{batch}/{sample}_R1_001.pe.qc.fastq.gz',
        R2L1 = '/group/runciegrp/Data/Illumina/bg/trimmed/{batch}/{sample}_R2_001.pe.qc.fastq.gz',
        SJfiles = '/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass1/SJ.out.tab',
        pass2="star/pass2/sjdbList.out.tab"
    params:
        outdir = '/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass2/',
        id = '{sample}',
        refdir = 'star/pass2',
        gtf = "/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.chr.gtf",
        tmpdir = '/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_tmp'
    output:
        '/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass2/Aligned.sortedByCoord.out.bam'
    threads: 5
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
        "/group/runciegrp/Data/Illumina/bg/final_bams/{batch}/{sample}.Aligned.sortedByCoord.MKDup.Processed.out.bam"
    input:
        '/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass2/Aligned.sortedByCoord.out.bam'
    params:
        prefix="/group/runciegrp/Data/Illumina/bg/final_bams/{batch}/{sample}.Aligned.sortedByCoord.MKDup.",
        tmp="/group/runciegrp/Data/Illumina/bg/final_bams/{batch}/{sample}_tmp",
	    bamdir="/group/runciegrp/Data/Illumina/bg/final_bams",
        batchdir="/group/runciegrp/Data/Illumina/bg/final_bams/{batch}",
        id="{sample}"
    threads: 5
    run:
        #shell("if [ ! -d {params.bamdir}; then mkdir {params.bamdir}; fi")
        #shell("cd {params.bamdir}")
        shell("if [ ! -d {params.batchdir} ]; then mkdir {params.batchdir}; fi")
        #shell("cd {params.batchdir}")
        shell("if [ -d -STARtmp ]; then rm -rf _STARtmp; fi")
        shell("STAR --runThreadN {threads} \
	    --runMode inputAlignmentsFromBAM \
        --bamRemoveDuplicatesType UniqueIdentical \
        --limitBAMsortRAM 32000000000 \
        --readFilesCommand samtools view \
        --readFilesType SAM PE \
        --outSAMattrRGline ID:{params.id}\tSM:{params.id}\tLB:None\tPL:Illumina \
        --outSAMstrandField intronMotif \
        --outSAMattributes Standard \
        --outSAMunmapped Within \
        --outTmpDir {params.tmp} \
        --inputBAMfile {input} \
        --outFileNamePrefix {params.prefix}")
        #shell("cd ..")
        #shell("mv {params.tmp} {output}")
        #shell("rm -rf _STARtmp")


rule index_bam:
    output:
        "/group/runciegrp/Data/Illumina/bg/final_bams/{batch}/{sample}.Aligned.sortedByCoord.MKDup.Processed.out.bam.bai"
    input:
        "/group/runciegrp/Data/Illumina/bg/final_bams/{batch}/{sample}.Aligned.sortedByCoord.MKDup.Processed.out.bam"
    run:
        shell("samtools index {input}")
