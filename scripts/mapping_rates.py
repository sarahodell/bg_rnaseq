import glob, sys
import os
import json

#SAMPLE_FILES=glob.glob("raw_reads/batch_2/*_R1_001.fastq.gz")
SAMPLE_FILES=glob.glob("raw_reads/batch_1/*_R1_001.fastq.gz")

SAMPLES=[]
for s in SAMPLE_FILES:
    s=os.path.basename(s)
    one=s.split('_')
    two=('_').join(one[:3])
    SAMPLES.append(two)


SAMPLES.remove("18048Fl-06-03-20_S212_L003")

for s in SAMPLES:
    with open('salmon_quant/{0}_transcripts_quant/aux_info/meta_info.json'.format(s)) as f:
        data = json.load(f)
    perc=data['percent_mapped']
    line="{0}\t{1}\n".format(s,perc)
    with open('salmon_quant/batch_1_mapping_rates.txt','a') as outfile:
        outfile.write(line)
