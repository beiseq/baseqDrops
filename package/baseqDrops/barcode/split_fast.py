import os, sys
import pandas as pd
import numpy as np
from time import time
from ..file_reader import read_file_by_lines
from .. import extract_barcode, extract_UMI
from itertools import product

barcode_prefix = [x[0]+x[1] for x in list(product('ATCG', repeat=2))]

def open_splited_files_bacode_prefix(dir, name, barcode_prefix_used):
    files = {}
    buffers = {}
    if not os.path.exists(dir):
        try:
            os.mkdir(dir)
        except:
            sys.exit("[error] Failed to make output dir...")

    path = os.path.join(dir, "split.{}.{}.fq".format(name, barcode_prefix_used))
    files = open(path, 'w')
    buffers = []
    return (files, buffers)

def barcode_split(name, protocol, bcstats, bc_prefix_used, fq1, fq2, dir, top_million_reads):
    """ Split the barcode for Certain BC prefix.
    
    It will read the barcode stats files, generated a barcode_corrected dict which will map raw barcode to corrected ones.
    """
    
    #barcode infos
    barcode_corrected = {}
    barcode_mutate_last = []
    bc_stats = pd.read_csv(bcstats)
    bc_stats = bc_stats.replace(np.nan, '', regex=True)

    for index, row in bc_stats.iterrows():
        barcode = row['barcode']
        barcode_corrected[barcode] = barcode
        mutate_last = row['mutate_last_base']
        if mutate_last:
            barcode_mutate_last.append(barcode)
        if row['mismatch_bc']:
            barcode_mis = row['mismatch_bc'].split("_")
            for bc in barcode_mis:
                barcode_corrected[bc] = barcode

    #make buffers and file ...
    files, buffers = open_splited_files_bacode_prefix(dir, name, bc_prefix_used)
    fq1 = read_file_by_lines(fq1, top_million_reads * 1000 * 1000, 4)
    fq2 = read_file_by_lines(fq2, top_million_reads * 1000 * 1000, 4)

    counter = 0
    start = time()
    for read1 in fq1:
        read2 = next(fq2)
        seq1 = read1[1].strip()
        counter += 1
        # Each 1M reads, report the progress;
        if counter % 1000000 == 0:
            print("[info] Processed {}M lines with barcode prefix {}, {}s per Million Reads".format(counter/1000000, bc_prefix_used, round(time()-start, 2)))
            start = time()
            files.writelines("\n".join(buffers)+"\n")
            buffers = []
        
        bc = extract_barcode(protocol, seq1)

        if bc in barcode_corrected:
            bc_corrected = barcode_corrected[bc]
        else:
            continue

        bc_prefix = bc_corrected[0:2]
        if bc_prefix != bc_prefix_used:
            continue

        #get UMI
        if barcode in barcode_mutate_last:
            UMI = extract_UMI(protocol, barcode, seq1, 1)
        else:
            UMI = extract_UMI(protocol, barcode, seq1, 0)
        if UMI == "":
            continue
        
        #build new header
        header = "_".join(['@', bc_corrected, UMI, str(counter)])
        seq = read2[1].strip()
        quality = read2[3].strip()
        buffers.append("\n".join([header, seq, "+", quality]))

    files.writelines("\n".join(buffers)+"\n")
    files.close()
    return 1

def split_by_barcode_faster(name, protocol, bcstats, fq1, fq2, dir, top_million_reads):
    """
    Split the read2 files for each barcode prefix.
    ::
        from baseqdrops import split_by_barcode_fast
        split_by_barcode_faster("sample", "10X", "bc.stats.txt", "fq1", "fq2", "./")
    """

    import multiprocessing as mp
    pool = mp.Pool(processes=16)
    results = []
    for bc in barcode_prefix:
        results.append(pool.apply_async(barcode_split, [name, protocol, bcstats, bc, fq1, fq2, dir, top_million_reads]))
    pool.close()
    pool.join()
    results = [x.get() for x in results]
