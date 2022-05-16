rule kallisto_index:
    input:
        "/group/jrigrp/Share/transcriptome/Zea_mays.B73_RefGen_v4.cdna.all.fa"
    output:
        "transcript_index"
    threads: 16
    run:
        shell("kallisto index \
        -t {input} \
        --threads {threads} \
        -k 31")

rule kallisto_quant:
    input:
    index = "transcript_index",
    R1 = "trimmed/{sample}_R1_001.pe.qc.fastq.gz",
    R2 = "trimmed/{sample}_R2_001.pe.qc.fastq.gz"
    output:
        "{sample}_transcripts_quant"
    params:
        libtype = "A"
    threads: 16
    run:
        shell("kallisto quant \
        -i {input.index} \
        --threads {threads} \
        -o {output}" \
        {input.R1} \
        {input.R2} \
        )
