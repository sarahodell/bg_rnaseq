rule bamqc:
    input:
        bam = "final_bams/update/{batch}/{sample}.Aligned.sortedByCoord.MKDup.Processed.out.bam",
        gff = '/group/jrigrp/Share/annotations/Zea_mays.B73_RefGen_v4.46.gtf'
    output:
        "qc/bamqc/update/{batch}/{sample}_stats/qualimapReport.html"
    params:
        pdir = "qc/bamqc/update/{batch}",
        sampdir = "qc/bamqc/update/{batch}/{sample}_stats",
        mem = "23G",
        qualimap = "/share/apps/qualimap-2.1.1/qualimap.jar"
    threads: 3
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


rule multibamqc:
    input:
        expand("qc/bamqc/update/{batch}/{sample}_stats/qualimapReport.html",zip,batch=DIRECTORIES,sample = SAMPLES)
    output:
        "qc/bamqc/update/multisampleBamQcReport.html"
    params:
        outdir = "qc/bamqc/update/",
        infile = "qc/bamqc/update/bamqc_list.txt",
        qualimap = "/share/apps/qualimap-2.1.1/qualimap.jar",
        mem = "23G"
    run:
        #shell("find qc/bamqc/ -mindepth 1 -maxdepth 1 -type d | grep SamC > qc/bamcqc/ALL.bamqclist.txt")
        import pandas as pd
        dirs=glob.glob("qc/bamqc/*/*_stats")
        dirs=pd.DataFrame(dirs,columns=['dirname'])
        #data = pd.read_csv("qc/bamqc/ALL.bamqclist.txt", sep = " ", header = None, names = ['filename'])
        dirs['sample'] = dirs['dirname'].str.split('/').str[3].str.split('_stats').str[0]
        dirs=dirs[['sample','dirname']]
        #data = data.sort_values('sample', axis = 0, ascending = True)
        #data = data[['sample','filename']]
        pd.DataFrame(dirs).to_csv('qc/bamqc/update/bamqc_list.txt',index=False,header=False,sep='\t') 
        #data.to_csv('qc/bamqc/bamqc_list.txt', header = None, index = None, sep = ' ', mode = 'a')
        shell("qualimap multi-bamqc \
        --java-mem-size={params.mem} \
        --data {params.infile} \
        --paint-chromosome-limits \
        -outdir {params.outdir} \
        -outformat html")
