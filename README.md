# baseqDrops
A versatile pipeline for processing dataset from 10X, indrop and Drop-seq.

## Install baseqDrops
We need python3 and a package called: baseqDrops, which could be installed by:

    pip install baseqDrops

After install, you will have a runnable command `baseq-Drop`

## Config file

The pipeline need the following software or resources:

+ `star`: STAR software, for fast alignment of RNA-Seq data;
+ `samtools`: Sorting bam file;
+ `whitelistDir`: The barcode whitelist files for indrop and 10X should be placed under whitelistDir. These files could bed dowmloaded from https://github.com/basedata10/DropRNA/tree/master/whitelist;
+ `cellranger_ref_<genome>`: The key process of read alignment and tagging to genes are inspired and borrowed from the open source cellranger pipeline(https://github.com/10XGenomics/cellranger). The refernces of genome index and transcriptome can be downloaded from https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest.
In the config file, the directory of cellrange references is named as `cellranger_<genome>`.

While running command, the configures are recorded in the file called `config_drops.ini`:

    [Drops]
    samtools = /path/to/samtools
    star = /path/to/STAR
    whitelistDir = /path/to/whitelist_file_directory
    cellranger_ref_hg38 = /path/to/reference/refdata-cellranger-GRCh38-1.2.0/


## Process Steps

1. `Extract the Cell Barcode` Counting the number of each kinds of barcode; this will genrate a barcode_count.<sample>.csv;
2. `Cell Barcode correction and filtering` Correcting the cell barcode with 1bp mismatch, filtering the barcode with min number of reads;
3. `Split the reads of valid Cell Barcodes` The raw pair-end raw reads are splitted to 16 single end files for multiprocessing according to the 2bp prefix of barcode; For example, we will get: split.<sample>.<AA|AT|AC|AG...|GG>.fq
4. `Star Alignment` Fastq files runs at the same time; The bam file sorted by sequence header is generated;
5. `Reads tagging` Tagging the reads alignment position to the corresponding gene name
6. `Genrating UMI table`


## Run Command

The main config is:

+ `--config`: config file;
+ `--genome/-g`: genome version;
+ `--protocol`: [10X|indrop|dropseq]
+ `--minreads`:  Minimum reads for a barcode
+ `--name/-n` : Sample name
+ `--fq1/-1`: Read 1
+ `--fq2/-2`: Read 2
+ `--top_million_reads`: How many million reads to use, mainly for testing pipeline with fraction of reads
+ `--dir/-d`: output path

If you config the: `cellranger_ref_hg38` you can run the following:

    baseqDrops run_pipe --config ./config_drops.ini -g hg38 -p 10X --minreads 10000 -n 10X_test -1 10x_1.1.fq.gz -2 10x.2.fq.gz -d ./

