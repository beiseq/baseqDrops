import pandas as pd
from .whitelist import read_whitelist, check_whitelist

def mutate_single_base(seq):
    mutated = []
    bases = ['A', 'T', 'C', 'G']
    for index in range(0, len(seq)):
        temp = list(seq)
        base_raw = temp[index]
        for base in bases:
            if base != base_raw:
                temp[index] = base
                mutated.append(''.join(temp))
    return mutated

def valid_barcode(protocol="", barcode_count="", max_cell=10000, min_reads=2000, output="bc.stats.txt"):
    """
    Aggregate the mismatch barcode, get the total_reads;

    Note:
        #. Read the barcode counts files;
        #. Correct the barcode with 1bp mismatch;
        #. Stats the mismatch barcode reads and sequences;
        #. Determine the last base of barcode is missed (mainly for drop-seq);
        #. Filter the barcode whitelist;
        #. Filter barcodes by read counts >= min_reads (2000);
        #. Filter barcodes by the max number of cells;
        #. Print the number of barcode and reads retained after each steps.

    Examples:
    ::
        from baseqdrops.barcode.stats import valid_barcode
        valid_barcode("10X", "bc.counts.txt",
            max_cell=10000, min_reads=2000, output="bc.stats.txt")

    This write a bc_stats.csv file (CSV) which contains:
    ::
        barcode/counts/mismatch_reads/mismatch_bc/mutate_last_base
    """

    print("[info] Correct and aggregate the barcodes in: {}".format(barcode_count))

    #Read the whitelist
    bc_white = read_whitelist(protocol) 
    
    df = pd.read_csv(barcode_count, index_col=0).sort_values("counts", ascending=False)
    
    df["mismatch_reads"] = 0
    df["mismatch_bc"] = ""
    df["mutate_last_base"] = 0
    df["total_reads"] = 0
    df["whitelist"] = 1

    #Aggregate by 1 Mismatch
    for bc in df.index.tolist():
        count = df.loc[bc, 'counts']
                
        #Count==0 means it is mutated from other barcode;
        if count == 0 or count <= 0.25 * min_reads:
            continue

        #Find 1bp mismatches of barcode;
        #The mismatch barcode should have less barcode counts;
        bc_mis = mutate_single_base(bc)

        index = df.index.isin(bc_mis)

        #The mutated barcode should contains less reads
        index_fewerreads = df['counts'] <= count
        
        df_mis = df.loc[index & index_fewerreads].sort_values("counts", ascending=False)
        print(df_mis)

        #Determine lost of last bases (Drop-Seq)
        barcode_mis = df_mis.index.tolist()
        if len(barcode_mis)>=3 and sum([1 for x in barcode_mis[0:3] if x[0:-1]==bc[0:-1]])==3:
            df.loc[bc, 'mutate_last_base'] = 1

        df.loc[bc, 'mismatch_reads'] = sum(df_mis['counts'])
        df.loc[bc, 'mismatch_bc'] = "_".join(df_mis.index.tolist())
        df.loc[bc, 'total_reads'] = sum(df_mis['counts']) + df.loc[bc, 'counts']
        
        #Set the Mutated Counts to 0 ...
        df.loc[index & index_fewerreads, 'counts'] = 0

    #Filter by whitelist
    for index, row in df.iterrows():
        valid = check_whitelist(bc_white, protocol, index)
        df.loc[index, 'whitelist'] = valid

    df_white = df[df['whitelist']==1]

    #Filter by read counts
    print("[info] Filtering the barcodes exceeds number {}".format(min_reads))
    df_white_reads = df_white.loc[df_white['total_reads'] >= min_reads].sort_values("total_reads", ascending=False)

    #Filter by max number of cells
    if df_white_reads.shape[0] >= max_cell:
        df_white_reads = df_white_reads.iloc[0:max_cell, ]

    #Print informations
    print("[info] Raw CBs count {} ==> Whitelist CBs {} ==> Abundent CBs {}".format(len(df.index), len(df_white.index), len(df_white_reads.index)))
    print("[info] Raw Reads {} ==> Whitelist Reads {} ==> Abundent Reads {}".format(sum(df['total_reads']), sum(df_white['total_reads']),
                                                                                    sum(df_white_reads['total_reads'])))
    print("[info] The valid barcode list is write to: {}".format(output))
    df_white_reads.to_csv(output)
