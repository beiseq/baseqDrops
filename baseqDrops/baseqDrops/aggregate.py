# -*- coding: utf-8 -*-
import json, re, sys, os
import random
import heapq

def HammingDistance(seq1, seq2):
    return sum([1 for x in zip(seq1, seq2) if x[0] != x[1]])

def sampling_UMIs(UMIs, ratio):
    """ Sampling the reads of each UMI.
    """
    UMIs_sample = {}
    for UMI in UMIs:
        UMIs_sample[UMI] = sum([1 for x in range(UMIs[UMI]) if random.random()<=ratio])
    return {x:UMIs_sample[x] for x in UMIs_sample if UMIs_sample[x]>0}

def calculate_UMI_with_mismatch(UMIs):
    """ Corrected the mismatches in UMIs
    input: UMI sequences and their counts;
    return: Corrected unique UMI sequences
    """
    if len(UMIs.keys()) == 1:
        return [x for x in UMIs if UMIs[x]>0]

    UMIs = sorted(UMIs.items(), key=lambda k: k[1], reverse=True)
    UMI_info = {x[0]:x[1] for x in UMIs}
    umi_num = len(UMIs)

    if umi_num <= 10:
        for idx1 in range(0, umi_num-1):
            for idx2 in range(idx1+1, umi_num):
                umi_1 = UMIs[idx1][0]
                umi_2 = UMIs[idx2][0]
                if HammingDistance(umi_1, umi_2) <= 1:
                    UMI_info[umi_1] += UMI_info[umi_2]
                    UMI_info[umi_2] = 0

    return [x for x in UMI_info if UMI_info[x]>0]

def read_barcode_gene_file(filepath, sampling_factor=1):
    """ Return the UMIs and Counts for each tagging file.
    [barcodes list, UMIs data and Reads data]
    
    UMIs data: UMIs-->genes-->barcode ==> UMIs count
    Reads data: Reads-->genes-->barcode ==> Read count 
    """
    UMIs = {}
    Reads = {}
    barcodes = []

    with open(filepath, 'r') as file:
        for line in file:
            infos = re.split('\t', line)
            barcode = infos[0]
            barcodes.append(barcode)
            umi_counts = json.loads(infos[1])
            for gene in umi_counts:
                gene_UMIs = sampling_UMIs(umi_counts[gene], sampling_factor)
                bcs = calculate_UMI_with_mismatch(gene_UMIs)
                reads = sum(gene_UMIs.values())
                if gene in UMIs:
                    UMIs[gene][barcode] = len(bcs)
                    Reads[gene][barcode] = reads
                else:
                    UMIs[gene] = {}
                    UMIs[gene][barcode] = len(bcs)
                    Reads[gene] = {}
                    Reads[gene][barcode] = reads
    return [barcodes, UMIs, Reads]

def write_to_table(barcodes, UMIs, outfile, info=""):
    """ Count the total UMIs of all genes and Find the top 10000 genes;
    """

    header = "\t".join(["genes"]+barcodes)
    total_UMIs = {}
    genes_abundent = []

    for gene in UMIs:
        total_UMIs[gene] = sum(UMIs[gene].values())

    nlargestList = heapq.nlargest(10000, total_UMIs.values())
    for gene in UMIs:
        if total_UMIs[gene] >= nlargestList[-1]:
            genes_abundent.append(gene)

    with open(outfile, 'w') as file:
        file.write(header+"\n")
        for gene in UMIs:
            if gene not in genes_abundent:
                continue
            umis = [gene]
            for bc in barcodes:
                if bc in UMIs[gene]:
                    umis.append(UMIs[gene][bc])
                else:
                    umis.append(0)
            umis_str = "\t".join([str(x) for x in umis])
            file.write(umis_str+"\n")

    print("[Success] Write the {} to: {}".format(info, outfile))
