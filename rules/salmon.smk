rule decoy_genome:
    input:
        transcriptome="/group/jrigrp/Share/transcriptome/Zea_mays.B73_RefGen_v4.cdna.all.fa",
        genome="/group/jrigrp/Share/assemblies/Zea_mays.AGPv4.dna.chr.fa",
        annotation="eqtl/data/Zea_mays.B73_RefGen_v4.46.gtf"
    output:
        "pseudos/decoys.txt"
    params:
        outdir="pseudos"
    run:
        shell("/home/sodell/bin/SalmonTools/scripts./generateDecoyTranscriptome.sh \
        -g {input.genome} \
        -t {input.transcriptome} \
        -a {input.annotation} \
        -o {params.outdir}")

rule salmon_index:
    input:
        transcriptome="/group/jrigrp/Share/transcriptome/Zea_mays.B73_RefGen_v4.cdna.all.fa",
        decoy="pseudos/decoys.txt"
    output:
        "founders_transcript_index/sa.bin",
    params:
        outdir="founders_transcript_index"
    threads: 16
    run:
        shell("salmon index \
        -t {input.transcriptome} \
        -i {params.outdir} \
        --decoys {input.decoy} \
        --threads {threads} \
        -k 21")

rule salmon_quant:
    input:
        index = "founders_transcript_index/sa.bin",
        R1 = "trimmed/{batch}/{sample}_R1_001.pe.qc.fastq.gz",
        R2 = "trimmed/{batch}/{sample}_R2_001.pe.qc.fastq.gz"
    output:
        "salmon_quant/{batch}/{sample}_transcripts_quant/quant.sf"
    params:
        outdir = "salmon_quant/{batch}/{sample}_transcripts_quant",
        outsam = "salmon_quant/batch/{sample}_transcripts_quant/{sample}_mapping.sam",
        libtype = "IU",
        index="founders_transcript_index"
    threads: 16
    run:
        shell("salmon quant \
        -i {params.index} \
        --threads {threads} \
        --writeMappings={params.outsam} \
        -l  {params.libtype} \
        -1 {input.R1} \
        -2 {input.R2} \
        -o {params.outdir}")
