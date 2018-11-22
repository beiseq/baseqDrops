import os, sys
import multiprocessing as mp

def pipeline(config, genome, protocol, cells, minreads, name, fq1, fq2, dir,
             top_million_reads, step, parallel):
    """
    Run the Data Processing Pipeline...

    #. Stats and count the barcode from pair-end 1 sequences;
    #. Read the barcode counts files;
    #. Correct the barcode with 1bp mismatch;
    #. Stats the mismatch barcode reads and sequences;
    #. Determine wheather mutate on the last base (show A/T/C/G with similar ratio at the last base);
    #. Filter by whitelist;
    #. Filter by read counts (>=min_reads);
    #. Print the number of barcode and reads retained after each steps.

    Usage:
    ::
        from baseqDrops import pipeline
        pipeline("", "hg38", "10X", 1000, minreads, name, fq1, fq2, dir, top_million_reads, step, parallel)
        
        #specify the length of UMI and barcodes
        pipeline("", "hg38", "10X", 1000, minreads, name, fq1, fq2, dir, top_million_reads, step, parallel)
        
        #Run in command line
        baseqdrops  

    Protocols:

    #. For 10X: 16bp Barcode and 10 bp UMI => 10X (most commonly used)
    #. For 10X: 14bp Barcode and 5/10 bp UMIs => 10X_14_10 / 10X_14_5 (For some old version data) 
    #. For DropSeq ==> dropseq
    #. For inDrop ==> indrop

    Args:
        config: The path of configuration file;
        genome: Genome version (hg38, mm10, hg38_mm10);
        cells: Max number of cells;
        minreads: Minimum number of reads for a cell barcode (10000);
        name: Samplename;
        fq1, fq2: Path to fastq reads;
        dir: The folder for processing, a folder with samplename will be created;
        top_million_reads: Number of reads used for processing in Million;
        step: Steps to run;
        parallel: How many thread to use;

    Steps:
        count
        stats
        split
        star
        tagging
        table
    """

    from . import count_barcode
    from . import valid_barcode
    from . import split_by_barcode
    from .barcode.split_fast import split_by_barcode_faster
    from .tagging.prime3 import check_reference_files

    #Set Config Files...
    if config:
        if not os.path.exists(config):
            sys.exit("[error] The config file does not exist in: {}!".format(config))
        os.environ["BASEQCFG"] = config

    #Check Annotation Files
    print('[info] Checking Reference Files...')
    check_reference_files(genome)

    print('[info] Start Processing Your RNA-Seq Dataset ...')

    dir = os.path.abspath(os.path.join(dir, name))
    bc_counts = os.path.join(dir, "barcode_count_{}.csv".format(name))
    bc_stats = os.path.join(dir, "barcode_stats_{}.csv".format(name))
    bc_splits_dir = os.path.join(dir, "barcode_splits")
    align_dir = os.path.join(dir, "star_align")
    tagging_dir = os.path.join(dir, "read_tagging")
    tpm_table = os.path.join(dir, "Result.UMIs.{}.txt".format(name))
    reads_table = os.path.join(dir, "Result.Reads.{}.txt".format(name))

    from itertools import product
    barcode_prefix = [x[0] + x[1] for x in list(product('ATCG', repeat=2))]

    dirs = [dir, align_dir, tagging_dir, bc_splits_dir]

    for dir in dirs:
        if not os.path.exists(dir):
            os.mkdir(dir)

    #Check the existance of the files
    if not os.path.exists(fq1):
         sys.exit("[error] Fastq-1 does not exist!")

    if not os.path.exists(fq2):
         sys.exit("[error] Fastq-2 does not exist!")

    #count barcode
    if step in ["all", "count"]:
        print("[info] Counting the barcodes ...")
        count_barcode(fq1, bc_counts, protocol, min_reads=50, topreads=int(top_million_reads))

    #aggregate
    if step in ["all", "stats"]:
        print("[info] Aggregating the barcodes errors and get valid ones ...")
        valid_barcode(protocol, bc_counts, max_cell=cells, min_reads=minreads, output=bc_stats)

    #barcode split
    if step in ["all", "split"]:
        print("[info] Split the barcode for each Barcode Prefix ...")
        split_by_barcode_faster(name, protocol, bc_stats, fq1, fq2, 
            bc_splits_dir, int(top_million_reads))

    #run alignment
    if step in ["all", "star"]:
        from .star import star_align
        star_align(bc_splits_dir, align_dir, name, genome, parallel=int(parallel))

    #run reads tagging
    if step in ["all", "tagging"]:
        from .tagging.prime3 import tagging_reads
        print('[info] Tagging the reads to genes...')
        pool = mp.Pool(processes=int(parallel))
        result = []
        for bc in barcode_prefix:
            bamfile = os.path.join(align_dir, "{}/{}.sort.bam".format(bc, bc))
            outfile = os.path.join(tagging_dir, "tagging.{}.txt".format(bc))
            #tagging_reads(genome, bamfile, outfile)
            result.append(pool.apply_async(tagging_reads, (genome, bamfile, outfile,)))
        pool.close()
        pool.join()

    #run Table aggragation
    if step in ["all", "table"]:
        print('[info] Build gene expression table from the tagging files...')

        from .aggregate import read_barcode_gene_file, write_to_table
        pool = mp.Pool(processes=int(parallel))
        result = []
        for bc in barcode_prefix:
            filepath = os.path.join(tagging_dir, "tagging.{}.txt".format(bc))
            result.append(pool.apply_async(read_barcode_gene_file, (filepath, 1)))
        pool.close()
        pool.join()

        from itertools import chain
        barcodes_all = [x.get()[0] for x in result]
        barcodes_lists = list(chain(*barcodes_all))
        
        exp = {}
        UMIs_all = [x.get()[1] for x in result]
        for UMI in UMIs_all:
            for gene in UMI:
                if gene in exp:
                    exp[gene].update(UMI[gene])
                else:
                    exp[gene] = UMI[gene]
        write_to_table(barcodes_lists, exp, tpm_table, "UMI Table")  
        
        exp = {}
        Reads_all = [x.get()[2] for x in result] 
        for UMI in Reads_all:
            for gene in UMI:
                if gene in exp:
                    exp[gene].update(UMI[gene])
                else:
                    exp[gene] = UMI[gene]
        write_to_table(barcodes_lists, exp, reads_table, "Read Count Table")

if __name__ == "__main__":
    print("begin running pipeline")