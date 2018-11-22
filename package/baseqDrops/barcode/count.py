import pandas as pd
from time import time
from ..file_reader import read_file_by_lines
from .. import extract_barcode

def count_barcode(path, output, protocol, min_reads=50, topreads=1000):
    """
    Count the number of the barcodes (), write a count file.
    Usage:
    ::
        from baseqdrops import count_barcode
        count_barcodes("10X.1.fq.gz", "bc.counts.txt", "10X", min_reads=50, topreads=1000)
    
    Args:   
        min_reads (int): The barcode with reads lower than this will be discard (50).
        topreads (str): To process the top N millions of reads (1000).
    
    Returns:
        A csv file of barcode_count, by default (bc.counts.txt)
        cellbarcode/counts
    """

    bc_counts = {}
    start = time()

    print("[info] Counting the barcodes, if there contains more than {}M reads in {}, the remainings will be skipped ...".format(topreads, path))
    print("[info] Barcode with less than {} reads is discard ".format(min_reads))
    
    index = 0

    lines = read_file_by_lines(path, topreads * 1000 * 1000, 4)

    for line in lines:
        index += 1
        bc = extract_barcode(protocol, line[1])
        if index % 1000000 == 0:
            print("[info] Processed {}M lines, {}s for one million reads".format(index/1000000, round(time()-start, 2)))
            start = time()
        if bc == "":
            continue
        if bc in bc_counts:
            bc_counts[bc] += 1
        else:
            bc_counts[bc] = 1

    bc_counts_filter = []
    
    for k, v in bc_counts.items():
        if v >= min_reads:
            bc_counts_filter.append([k, v])

    print("[info] Barcode depth file: {}".format(output))
    df = pd.DataFrame(bc_counts_filter, columns=["barcode", "counts"])
    df.to_csv(output, sep=",", index=False)