##### Index Genome and Alignment Pass One #####
rule star_pass_one:
    output:
        outdir=touch(directory("star/results/pass1")),
        sjdb=expand("star/results/bamfiles/{sample}.sjdbList.out.tab",sample=SAMPLES)
    input:
        R1=expand("trimmed/{sample}_R1_001.pe.qc.fastq.gz",sample=SAMPLES),
        R2=expand("trimmed/{sample}_R2_001.pe.qc.fastq.gz",sample=SAMPLES) 
    params:
        annotation="/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.chr.gtf",
        reference="/group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa",
	tmp="star/results/bamfiles/tmp",
	prefix="star/results/bamfiles"
    threads: 16
    run:
        shell("STAR \
	--runMode genomeGenerate \
	--genomeDir {output.outdir} \
	--genomeFastaFiles {params.reference} \
	--sjdbOverhang 100 \
	--sjdbGTFfile {params.annotation} \
	--runThreadN {threads} \
	--outTmpDir {params.tmp}")
	shell("rm -r {params.tmp}")
	shell("STAR \
	--genomeDir {params.index} \
	--readFilesIn {input.R1} {input.R2} \
	--outFilterMultimapScoreRange 1 \
	--outFilterMultimapNmax 20 \
	--outFilterMismatchNmax 20 \
	--alignIntronMax 500000 \
	--alignMatesGapMax 1000000 \
	--sjdbScore 2 \
	--alignSJDBoverhangMin 1 \
	--genomeLoad NoSharedMemory \
	--outFilterMatchNminOverLread 0.33 \
	--outFilterScoreMinOverLread 0.33 \
	--sjdbOverhang 100 \
	--outSAMstrandField intronMotif \
	--outSAMtype None \
	--outFileNamePrefix {params.prefix} \
	--outSAMmode None \
	--readFilesCommand zcat \
	--runThreadN {threads} \
	--outTmpDir {params.tmp}")
	shell("rm -r {params.tmp}")

##### Second Genome Index and Alignment #####

rule star_pass_two:
    input:
        R1=expand("trimmed/{sample}_R1_001.pe.qc.fastq.gz",sample=SAMPLES),
        R2=expand("trimmed/{sample}_R2_001.pe.qc.fastq.gz",sample=SAMPLES),
        sjdb=expand("star/results/bamfiles/{sample}.sjdbList.out.tab",sample=SAMPLES)
    output:
        outdir=touch(directory("star/results/pass2")),
        bam=expand("star/results/bamfiles/{sample}.Aligned.sortedByCoord.out.bam",sample=SAMPLES)
    threads: 16
    params:
        annotation="/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.chr.gtf",
        reference="/group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa",
        index="star/results/pass2",
        tmp="star/results/bamfiles/tmp",
	prefix="star/results/bamfiles"
    run:
        shell("STAR \
	--runMode genomeGenerate \
	--genomeDir {params.index} \
	--genomeFastaFiles {params.reference} \
	--sjdbOverhang 100 \
	--runThreadN {threads} \
	--sjdbFileChrStartEnd {input.sjdb} \
	--outTmpDir {params.tmp}")
        shell("rm -r {params.tmp}")
        shell("STAR \
	--genomeDir {params.index} \
	--readFilesIn {input.R1} {input.R2} \
	--runThreadN {threads} \
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
	--sjdbOverhang 100 \
	--outSAMstrandField intronMotif \
	--outSAMattributes NH HI NM MD AS XS \
	--outSAMunmapped Within \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMheaderHD @HD VN:1.4 \
	--outFileNamePrefix {params.prefix}")
	shell("rm -r {params.tmp}")

rule mark_duplicates:
    output:
        bam="star/results/bamfiles/{sample}.sortedByCoord.MkDup.out.bam"
    input:
        "star/results/bamfiles/{sample}.Aligned.sortedByCoord.out.bam"
    params:
        tmp="star/results/bamfiles/tmp",
	prefix="star/results/bamfiles/{sample}.sortedByCoord.MkDup"
    run:
        shell("STAR \
	--runMode inputAlignmentsFromBAM \
	--bamRemoveDuplicatesType UniqueIdentical \
	--readFilesCommand samtools view \
	--readFilesType SAM PE \
	--inputBAMfile {input} \
	--outFileNamePrefix {params.prefix}")
        shell("rm -r {params.tmp}")