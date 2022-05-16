rule rnaseqc:
    input:
        bam = "final_bams/{batch}/{sample}.Aligned.sortedByCoord.MKDup.Processed.out.bam",
        gtf = "/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.gtf"
    output:
        "qc/rnaseqc/{batch}/{sample}_stats/qualimapReport.html"
    params:
        pdir = "qc/rnaseqc/{batch}",
        sampdir = "qc/rnaseqc/{batch}/{sample}_stats",
        mem = "63G",
        qualimap = "/share/apps/qualimap-2.1.1/qualimap.jar"
    run:
        shell("if [ ! -d {params.pdir} ]; then mkdir {params.pdir}; fi")
        shell("if [ -d {params.sampdir} ]; then rm -rf {params.sampdir}; fi")
        shell("qualimap rnaseq \
        -outdir {params.sampdir} \
        --java-mem-size={params.mem} \
        --paired \
        --sequencing-protocol strand-specific-forward \
        -bam {input.bam} \
        -gtf {input.gtf}")

rule bamqc:
    input:
        bam = "final_bams/{batch}/{sample}.Aligned.sortedByCoord.MKDup.Processed.out.bam",
        gff = '/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.gff3'
    output:
        "qc/bamqc/{batch}/{sample}_stats/qualimapReport.html"
    params:
        pdir = "qc/bamqc/{batch}",
        sampdir = "qc/bamqc/{batch}/{sample}_stats",
        mem = "63G",
        qualimap = "/share/apps/qualimap-2.1.1/qualimap.jar"
    threads: 8
    run:
        shell("if [ ! -d {params.pdir} ]; then mkdir {params.pdir}; fi")
        shell("if [ -d {params.sampdir} ]; then rm -rf {params.sampdir}; fi")
        shell("qualimap bamqc \
        -bam {input.bam} \
        --paint-chromosome-limits \
        --java-mem-size={params.mem} \
        --sequencing-protocol strand-specific-forward \
        -gff {input.gff} \
        -nt {threads} \
        -outdir {params.sampdir} \
        -outformat HTML \
        --skip-duplicated")

################################################################################
# Combine qualimap results
# this includes some python magic to create the params.infile txt file, which
# lists all of the individual qualimap report information
################################################################################

rule multibamqc:
    input:
        all = expand("qc/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLES)
    output:
        "qc/multisampleBamQcReport.html"
    params:
        outdir = "qc",
        infile = "qc/bamqc_list.txt",
        qualimap = "/share/apps/qualimap-2.1.1/qualimap.jar",
        mem = "20G"
    run:
        shell("find qc/bamqc -mindepth 1 -maxdepth 1 -type d | grep SamC > qc/ALL.bamqclist.txt")
        import pandas as pd
        data = pd.read_csv("qc/ALL.bamqclist.txt", sep = " ", header = None, names = ['filename'])
        data['sample'] = data['filename'].str.split('.').str[0].str.split('/').str[2].str.split('_stats').str[0]
        data = data.sort_values('sample', axis = 0, ascending = True)
        data = data[['sample','filename']]
        data.to_csv(r'qc/bamqc_list.txt', header = None, index = None, sep = ' ', mode = 'a')
        shell("qualimap multi-bamqc \
        --java-mem-size={params.mem} \
        --data {params.infile} \
        --paint-chromosome-limits \
        -outdir {params.outdir} \
        -outformat html")
