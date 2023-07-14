rule htseq:
        input:
            bam = '/group/runciegrp/Data/Illumina/bg/star/{batch}/{sample}_pass2/Aligned.sortedByCoord.out.bam',
            gtf = "/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.gtf"

        output:
            outfile='/home/sodell/projects/biogemma/expression/htseq/{batch}/{sample}_HTSeq_union_gtf_no_gene_ID.log'#,
            #csvfile='/home/sodell/projects/biogemma/expression/htseq/update/{batch}/{sample}_HTSeq.csv'
        threads: 1
        run:
            shell("htseq-count -m union -t gene -r pos {input.bam} {input.gtf} > {output.outfile}")
            #shell("grep ENS {output.outfile} | sed "s/gene://g" > {output.csvfile}")