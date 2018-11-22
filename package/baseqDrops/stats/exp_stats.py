import pandas as pd

def stats_table(counts, umis):
    """
    Stats on the table of UMIs
    ::
        from baseq.drops import stats_table
        stats_table("counts.txt", "UMIs.txt")
    """
    df = pd.read_table(umis, index_col=0)
    #df_reads = pd.read_table(counts, index_col=0)

    #Number of cells
    cells = len(df.index)

    #Number of genes
    gene_counts = (df>0).astype(int).sum(axis=1).median()
    UMI_counts = df.sum(axis=1).median()
    #read_counts = df_reads.sum(axis=1).median()

    print("Cells Number: {}".format(cells))
    #print("Median Read Counts: {}".format(read_counts))
    print("Median UMI Counts: {}".format(UMI_counts))
    print("Median Gene Counts: {}".format(gene_counts))

    return df