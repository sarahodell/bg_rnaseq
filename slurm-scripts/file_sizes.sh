#!/bin/bash

cd /group/runciegrp/Data/Illumina/reads/bg/18048-69
for i in *.fastq.gz; do
    FILESIZE=$(stat -c%s "$i")
    echo "$i  $FILESIZE" >> /home/sodell/projects/biogemma/expression/raw_reads/batch_2_file_sizes.txt
done

cd /home/sodell/projects/biogemma/expression/trimmed
for i in well*pe.qc.fastq.gz; do
    FILESIZE=$(stat -c%s "$i")
    echo "$i  $FILESIZE" >> /home/sodell/projects/biogemma/expression/raw_reads/batch_2_trimmed_sizes.txt
done
