rule salmon_index:
    input:
        "/group/jrigrp/Share/transcriptome/Zea_mays.B73_RefGen_v4.cdna.all.fa"
    output:
        directory("transcript_index")
    threads: 16
    run:
        shell("salmon index \
        -t {input} \
        -i {output} \
        --threads {threads} \
        -k 21")

rule salmon_quant:
    input:
        index = "transcript_index",
        R1 = "trimmed/{sample}_R1_001.pe.qc.fastq.gz",
        R2 = "trimmed/{sample}_R2_001.pe.qc.fastq.gz"
    output:
        "salmon_quant/{sample}_transcripts_quant/quant.sf"
    params:
        outdir = "salmon_quant/{sample}_transcripts_quant",
        libtype = "A"
    threads: 16
    run:
        shell("salmon quant \
        -i {input.index} \
        --threads {threads} \
        --validateMappings \
        --minAssignedFrags 0 \
        -l  {params.libtype} \
        -1 {input.R1} \
        -2 {input.R2} \
        -o {params.outdir}")
