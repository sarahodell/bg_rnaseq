#!/usr/bin/env python

"""
Takes a fasta file and converts it to a database
"""
import pandas as pd
from subprocess import PIPE,Popen,STDOUT
import sys
import argparse
from Bio.Seq import Seq


#testing #####
#with open(infile) as f:
#    line = f.readline().strip('\n')


#fasta_infile="ASE/Zea_mays_B73v4_exons.fa"
#fasta_db=parse_fasta(fasta_infile)
#fasta_keys=fasta_db.keys()
#vcf_infile="ASE/Biogemma_WGS_Founders_transcript_variants.vcf.gz"
#samples,info = call_bcftools(vcf_infile)
#sample="Sample_MBS847"
#loc=samples.index(sample)
#ref_db=alt_ref_db(info,loc)
#markers=ref_db.keys()
# 0 is chrom
# 1 is start
# 2 is end
# 3 is sequence
outfile='pseudos/{0}_pseudo.fa'.format(sample)
#####
def parse_args():
    """ -h for info on arguments
    """
    parser = argparse.ArgumentParser(description="""Program description""")
    parser.add_argument("sample",type=str,help="The name of the sample")
    parser.add_argument("fasta_infile",type=str,help="""The input fasta file""")
    parser.add_argument("vcf_infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output fasta filename""")
    args=parser.parse_args()
    return args

def parse_fasta(fasta_infile):
    #input_file="pseudos/test.fa"
    f = open(fasta_infile)
    parsed_seqs = {}
    curr_seq_id = None
    curr_seq = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            linesplit=line.split(':')
            #transcript1=linesplit[0].split('.')
            transcript=linesplit[0].split('>')[1]
            #exon=transcript[1]
            #strand=linesplit[2]
            #linesplit2=linesplit[2].split(':')
            chr=linesplit[2]
            start=int(linesplit[3].split('-')[0])
            end=linesplit[3].split('-')[1]
            #end=int(end1.split('(')[0])
            #
            #gene=linesplit[3].split(':')[1]
            if curr_seq_id is not None:
                if len(curr_seq) > 4:
                    curr_seq[3]=('').join(curr_seq[3:])
                    curr_seq=curr_seq[:4]
                parsed_seqs[curr_seq_id] = curr_seq
            curr_seq_id = transcript
            curr_seq = [chr,start,end]
            continue
        curr_seq.append(line)
    if len(curr_seq) > 4:
        curr_seq[3]=('').join(curr_seq[3:])
        curr_seq=curr_seq[:4]
    parsed_seqs[curr_seq_id] = curr_seq
    return(parsed_seqs)


#Need to grab position of SNPs I want to keep for that founder
#Need to get the reference and alternate alleles (is it the same or different
#from B73 as well as from MBS847 (anything that isn't 0))
    # Make dictionary of SNPs with 0 and 1 alleles

# Need to find the correct location of the nucleotide in the fasta database and change it
def call_bcftools(vcf_infile):
    """Reads vcf file and extracts data on samples and genotypes for each marker
    Returns a list of samples and list of markers and genotypes (info)
    """
    process1 = Popen(['bcftools','query','-l',vcf_infile],stdout=PIPE,stderr=STDOUT)
    stdout,stderr=process1.communicate()
    print(stderr)
    samples = stdout.split('\n')[:-1]
    process2 = Popen(['bcftools','query','-f','%CHROM\t%POS\t%REF\t%ALT[\tGT=%GT]\n',vcf_infile],stdout=PIPE,stderr=STDOUT)
    stdout,stderr=process2.communicate()
    info = stdout.split('\n')[:-1]
    print(stderr)
    print "Read vcf file: Contains info on {0} samples and {1} markers".format(len(samples),len(info))
    return samples,info


def alt_ref_db(info,loc):
    """ Outputs a database with keys as SNP names and values as the ref and alt
    allele for the SNP. Last item in the list is the sample's genotype
    """
    ref_db = {}
    #sample_list=[]
    for i in info:
        split=i.split('\t')
        genos=split[4:]
        sample_geno=genos[loc]
        #nogood=['0','.']
        if '0' not in sample_geno and '.' not in sample_geno:
            sample_geno=int(sample_geno.split('=')[1])
            marker = "r{0}_{1}".format(split[0],split[1])
            ref=split[2]
            alt=split[3]
            if ',' in alt:
                ref_db[marker]=[ref]
                alt=alt.split(',')
                for a in alt:
                    if a=='*':
                        ref_db[marker].append('N')
                    else:
                        ref_db[marker].append(a)
            else:
                ref_db[marker]=[ref,alt]
            ref_db[marker].append(sample_geno)
    return(ref_db)


# This should be two functions
# one to make the alt,ref database for fast lookup
# one to get the sample genotype

def variant_loc(fasta_db,ref_db):
    new_fasta_db=fasta_db
    fasta_keys=new_fasta_db.keys()
    markers=ref_db.keys()
    for m in markers:
        split=m.split('_')
        chr=split[0].split('r')[1]
        pos=int(split[1])
        sample_allele=ref_db[m][-1]
        new_allele=ref_db[m][sample_allele]
        for k, v in new_fasta_db.items():
            if v[0]==chr:
                if pos < v[2] and pos > v[1]:
                    #print(k)
                    find=pos-v[1]-1
                    new_seq=v[3][:find] + new_allele + v[3][find + 1:]
                    v[3]=new_seq
                elif pos == v[1]:
                    #print(k)
                    find=0
                    new_seq=new_allele + v[3][1:]
                    v[3]=new_seq
                elif pos == v[2]:
                    #print(k)
                    find=-1
                    new_seq=v[3][:-1] + new_allele
                    v[3]=new_seq
    return(new_fasta_db)

def new_fasta(sample_fasta_db,sample):
    txt=''
    for(k,v in sample_fasta_db.items())
        header='>'+k+'::'+v[0]+':'+str(v[1])+'-'+str(v[2])+':'+sample+'\n'
        line=header+v[3]+'\n'
        txt+=line
    return(txt)

def get_pseudo_transcriptome():
    args=parse_args()
    fasta_db=parse_fasta(args.fasta_infile)
    samples,info = call_bcftools(args.vcf_infile)
    loc=samples.index(args.sample)
    ref_db=alt_ref_db(info,loc)
    sample_fasta_db=variant_loc(fasta_db,ref_db)
    txt=new_fasta(sample_fasta_db,sample)
    print "Writing out to {0}".format(args.outfile)
    with open(args.outfile,'w') as outfile:
        outfile.write(txt)


if __name__ == "__main__":
    get_pseudo_transcriptome()
