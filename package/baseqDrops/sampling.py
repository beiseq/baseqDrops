import os, sys
import multiprocessing as mp
from .aggregate import read_barcode_gene_file, write_to_table

def runsampling(name, dir, ratio="10_90", parallel=8):
    """Sampling the reads to certain relative depth
    and get the UMIs and Counts table.
    
    Args:
        name: name of the sample;
        dir: path of the analysis;
        ratio: 10_100 by default (from 10 to 100% by 10%);
        
    The ratio can be:
        "10_100": 10_100 by default (from 10 to 100% by 10%);
        "integer": like 50, means sampling 50% of the raw reads;
    
    It will use the files in dir/name/read_tagging, the UMIs results at each sampling
    depth will be exported at dir/samplename/...
    """

    print('Start Aggregating Results ...')

    dir = os.path.abspath(os.path.join(dir, name))
    tagging_dir = os.path.join(dir, "read_tagging")
    
    if ratio=="10_100":
        sampling_ratios = [x*10 for x in range(1, 11)]
    else:
        try:
            sampling_ratios = [int(ratio)]
        except:
            sys.exit("Ratio should be integer rather than string...")

    for sampling_factor in sampling_ratios:

        print("Sampling Factor is {}".format(sampling_factor))

        tpm_table = os.path.join(dir, "Result.UMIs.sampling_{}.{}.txt".format(str(sampling_factor), name))
        reads_table = os.path.join(dir, "Result.Reads.sampling_{}.{}.txt".format(str(sampling_factor), name))

        from itertools import product
        barcode_prefix = [x[0] + x[1] for x in list(product('ATCG', repeat=2))]

        print('[info] Tagging the reads to genes...')

        pool = mp.Pool(processes=int(parallel))
        result = []
        for bc in barcode_prefix:
            filepath = os.path.join(tagging_dir, "tagging.{}.txt".format(bc))
            result.append(pool.apply_async(read_barcode_gene_file, (filepath, float(sampling_factor)/100)))
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
        write_to_table(barcodes_lists, exp, tpm_table)  
        
        exp = {}
        Reads_all = [x.get()[2] for x in result] 
        for UMI in Reads_all:
            for gene in UMI:
                if gene in exp:
                    exp[gene].update(UMI[gene])
                else:
                    exp[gene] = UMI[gene]
        write_to_table(barcodes_lists, exp, reads_table)
