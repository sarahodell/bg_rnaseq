rule kallisto_index:
    input:
        "/group/jrigrp/Share/transcriptome/Zea_mays.B73_RefGen_v4.cdna.all.fa"
        #"pseudos/Sample_MBS847_pseudo.fa"
    output:
        "founders_transcript_index"
    threads: 32
    run:
        shell("kallisto index \
        -t {input} \
        --threads {threads} \
        -k 31")

rule kallisto_quant:
    input:
    index = "founders_transcript_index",
    R1 = "trimmed/{batch}/{sample}_R1_001.pe.qc.fastq.gz",
    R2 = "trimmed/{batch}/{sample}_R2_001.pe.qc.fastq.gz"
    output:
        "kallisto/{sample}_transcripts_quant"
    params:
        libtype = "A"
    threads: 16
    run:
        shell("kallisto quant \
        -i {input.index} \
        --threads {threads} \
        --fr-stranded \
        --gtf eqtl/data/Zea_mays.B73_RefGen_v4.46.gtf \
        --pseudobam \
        --genomebam \
        -o {output} \
        {input.R1} \
        {input.R2}"
        )
