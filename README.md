# baseqDrops
A versatile pipeline for processing dataset from 10X, indrop and Drop-seq.

## Install baseqDrops
We need python3 and a package called: baseqDrops, which could be installed by:

    pip install baseqDrops

After install, you will have a runnable command `baseqDrops`

It is recommend for the computer or server to have memory >= 30Gb and CPU cores >=8 for efficient processing;

## Configuration file

The following software or resources are required:

+ `star`: STAR software, for fast alignment of RNA-Seq data to the genome;
+ `samtools`: For sorting the aligned bam file (version >=1.6);
+ `whitelistDir`: The barcode whitelist files for indrop and 10X should be placed under whitelistDir. These files could bed downloaded from https://github.com/beiseq/baseqDrops/tree/master/whitelist;
+ `cellranger_ref_<genome>`: The key process of read alignment and tagging to genes are inspired and borrowed from the open source cellranger pipeline(https://github.com/10XGenomics/cellranger). The references of genome index and transcriptome can be downloaded from https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest.
In the config file, the directory of cellranger references is named as `cellranger_<genome>`.

While running command, the configures are recorded in the file called `config_drops.ini`:

    [Drops]
    samtools = /path/to/samtools
    star = /path/to/STAR
    whitelistDir = /path/to/whitelist_file_directory
    cellranger_ref_hg38 = /path/to/reference/refdata-cellranger-GRCh38-1.2.0/

## For Help Informations
	
	baseqDrops run-pipe --help

## Process Steps

1. `Cell Barcode Counting`: Counting the existed barcodes in dataset. This will generate a file named: barcode_count_<sample>.csv;
2. `Cell Barcode Correction, Aggregating and Filtering`: Correcting the cell barcodes within 1bp mismatch and then aggregating, filtering the barcode by minimum number of reads (default 5000), this will generate a valid barcode list named: barcode_stats_<sample>.csv;
3. `Split the Reads of Valid Cell Barcodes`: The raw pair-end raw reads are splitted to 16 single-end files for multiprocessing according to the 2bp prefix of the barcode; The folder of barcode_splits contains files like: split.<sample>.<AA|AT|AC|AG...|GG>.fq;
4. `Alignment to Genome using STAR`: Several (defined by --parallel/-p) STAR programs run at the same time, the results will be at folder named as star_align; The bam files are further sorted by sequence header;
5. `Reads Tagging`: Tagging the reads alignment position to the corresponding gene name;
6. `Generating Expression Table`: Both the expression table quantified by UMI (Result.UMIs.<sample>.txt) and raw read count (Result.Reads.<sample>.txt) will be generated;

## Run Pipeline

These parameters should be provided: (or run: baseqDrops run-pipe --help for information)

+ `--outdir/-d`: Output path (default ./, the result will be stored in ./<name>);
+ `--config`: Path to the config file;
+ `--genome/-g`: Genome version [hg38/mm38/hgmm];
+ `--protocol/-p`: [10X|indrop|dropseq];
+ `--minreads`:  Minimum reads required for a barcode;
+ `--name/-n` : Name of sample, a folder of <outdir>/<name> will be created and be the main directory; 
+ `--parallel` : The number of STAR and tagging processes runs at the same time (default is 4, need more memory for larger parallel number); 
+ `--fq1/-1`: Path of Pair-end 1 sequencing file;
+ `--fq2/-2`: Path of Pair-end 2 sequencing file;
+ `--top_million_reads`: For huge dataset, you can choose to use part of the data for a quick look, the reads exceeding N million of reads will be skipped;

If your data is human origin and `cellranger_ref_hg38` has been defined in configuration file, you can run:

    baseqDrops run-pipe --config ./config_drops.ini -g hg38 -p 10X --minreads 1000 -n 10X_test -1 10x_1.1.fq.gz -2 10x.2.fq.gz -d ./

## Run by Single Steps

We also provide step-wise ways for running the pipeline, all the parameters should be provided as described above, an extra "--step" should be provided, for example:
	
	baseqDrops run-pipe --config ./config.ini -g hg38 -p dropseq --minreads 1000 -n dropseq2 --top_million_reads 20 -1 dropseq_1.1.fq.gz -2 dropseq.2.fq.gz --step count -d ./

The steps are listed:

+ `Cell Barcode Counting`:  --step count
+ `Cell Barcode Correction, Aggregating and Filtering`: --step stats
+ `Split the Reads of Valid Cell Barcodes`: --step split
+ `Alignment to Genome using STAR`: --step star
+ `Reads Tagging` : --step tagging
+ `Generating Expression Table`: --step table

## Contact

For any questions, please email to: friedpine@gmail.com
