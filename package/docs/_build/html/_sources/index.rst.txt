baseqDrops
==================
baseqDrops enables fast analysis of 10X, inDrop and DropSeq datasets. The process includes barcode detection, error correction, aggregation and spliting. The raw sequencing datas are splitted into 16 files according to their barcode prefix (first two bases). 

Install
----------
The package is host at Pypi and can be installed using pip command as:
::
   pip install baseqDrops
   
After installation, an executable command named as baseqDrops can be used for conducting analysis.

Config file
--------------
The pipeline are relied on some softwares and resources. The paths should be record in an file (named as config.drops.ini, for example). The following items shoule be provided.

* star: STAR_, for fast alignment of RNA-Seq data;
* samtools_: For manipulating bam files (recommand version 1.2);
* whitelistDir: The barcode whitelist files for indrop and 10X should be placed under whitelistDir. These files could bed downloaded from whitelist_;
* cellranger_ref: The key process of read alignment and tagging to genes are inspired and borrowed from the open source cellranger pipeline(https://github.com/10XGenomics/cellranger). The refernces of genome index and transcriptome can be downloaded from https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest. In the config file, the directory of cellrange references is named as cellranger_<genome>.

.. _STAR: https://github.com/alexdobin/STAR
.. _samtools: https://github.com/samtools/samtools/releases/
.. _whitelist: https://github.com/basedata10/DropRNA/tree/master/whitelist

Here is an example of config file (config.ini):
::
    [Drops]
    samtools = /path/to/samtools
    star = /path/to/STAR
    whitelistDir = /path/to/whitelist_file_directory
    cellranger_ref_hg38 = /path/to/reference/refdata-cellranger-GRCh38-1.2.0/

Run pipeline
-----------------
Here is a standard command:
::
    baseqDrops run_pipe --config ./config.ini -g hg38 -p 10X \
        --minreads 10000 -n 10X_sample \
        -1 10X_sample.1.fastq.gz \
        -2 10X_sample.2.fastq.gz \
        -d ./

Some key parameters:
::
    -g: genome version, could be hg38/mm38/hg_mm, the corresponding cellranger_ref should be provided in the config file;
    --config: config file;
    --genome/-g: genome version;
    --protocol: [10X|indrop|dropseq|10X_14_5]
    --minreads: Minimum reads for a barcode
    --name/-n : Sample name
    --fq1/-1: Read 1
    --fq2/-2: Read 2
    --top_million_reads: How many million reads to use, mainly for testing pipeline with fraction of reads
    --dir/-d: output path

Protocols supported
---------------------
We support the following .

.. csv-table::
    :header: "Name", "Description"
    :widths: 10, 30

    "10X", "10X"
    "10X_14_5", "10X older version, 14bp UMI and 5 bp UMI"
    "10X_14_10", "10X older version, 14bp UMI and 10 bp UMI"
    "indrop", "indrop"
    "Drop-Seq", "dropseq"

Commands
--------------
**Downsampling Reads** Get the simuated UMIs and Reads matrix on lower sequencing depth from 10%-90% of the total reads.
::
    baseqDrops run_sampling -d ./ -n indrop2 -t 8 &

APIs
------------
.. toctree::
   :maxdepth: 1

        APIs 