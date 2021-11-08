Scripts for running RNASeq Analysis Pipeline

FASTQ file
1) Trimmomatic
2) fastqc
3) STAR or
4) salmon

*In Progress
5) kallisto
6) htseq
7) qualimap

To make a [copy of this repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository):

```console
git clone https://github.com/sarahodell/bg_rnaseq.git
```

The full set of fastq files are located in `/group/runciegrp/Data/Illumina/reads/bg/` with each directory containing a different batch of sequencing libraries.
For ease of use, I would recommend making a symbolic link to your local directory so that we don't have to type out that path every time.

```console
cd bg_rnaseq
mkdir raw_reads
cd raw_reads
ln -s /group/runciegrp/Data/Illumina/reads/bg/ .
```

To run the snakemake pipeline, adjust the file, `Snakefile` so that it includes the sample files you'd like to run.
The path to the fastq files on farm and their associated sample information can be found in ` metadata/BG_sequencing_sample_conversion_table.txt`
Some of these files are from a different experiment, so the shortened list of files that we want to use for this project are found in `metadata/BG_completed_sample_list.txt`

To change which programs you'd like to run, edit the bottom of the file to only `include` the rules you wish to run (comment out the rest), and also in the lines after `run all:`,
comment out the output files that correspond to the programs you don't wish to run.
This involves looking through the `.smk` files in the `rules` directory and looking at the lines under `output`.

This pipeline is setup to run on the cluster, so each job is submitted to slurm separately.
We can change the partition, the amount of memory, the amount of CPU, and the amount of time requested in the file `submit.json`.
Right now, most of them are set to request from `bigmemh` or `bigmemm`.
However, if you don't have permission to access those nodes, you may have to switch it to `high` or `med`. You can also adjust the mem, cpu, and time to adjust to the requirements of those nodes.

Once this is all done, we can run the snakemake pipeline with:
```console
./submit.sh
```
