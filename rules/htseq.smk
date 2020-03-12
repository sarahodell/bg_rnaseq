rule htseq:
        input:
            bam = '{sample}_pass2/Aligned.sortedByCoord.out.bam',
            gff = '/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.gff3'
        output:
            log = '{sample}_HTSeq_union_gff3_no_gene_ID.log',
            csv = '{sample}_HTSeq.csv'
        threads: 1
        run:
            shell("htseq-count -m union -s no -t gene -i ID -r pos -f bam {input.bam} {input.gff} &> {output.log}")
            shell("grep ENS {output.log} | sed 's/gene://g' > {output.csv}")
