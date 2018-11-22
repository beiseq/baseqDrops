import os, sys
from ..config import get_config

def rev_comp(seq):
    ___tbl = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(___tbl[s] for s in seq[::-1])

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

def read_whitelist(protocol):
    """
    Read Whitelist, get whitelist from config file: Drops/whitelistDir,
    Return: A set of whitelist barcodes
    """
    whilelistdir = get_config("Drops", "whitelistDir")
    if protocol == "10X":
        bc_white = [set()]
        white_path = os.path.join(whilelistdir, "whitelist.10X.txt")
        print("[info] The white list is {}".format(white_path)) 
        with open(white_path, 'r') as infile:
            for line in infile:
               bc_white[0].add(line.strip())

    elif protocol in ["10X_14_10", "10X_14_5"]:   
        bc_white = [set()]
        white_path = os.path.join(whilelistdir, "whitelist.10X_14.txt")
        print("[info] The white list is {}".format(white_path))
        with open(white_path, 'r') as infile:
            for line in infile:
               bc_white[0].add(line.strip())

    elif protocol == "indrop":
        bc_white = [set(), set()]
        white_path = os.path.join(whilelistdir, "whitelist.indrop_1.txt")
        if not os.path.exists(white_path):
            sys.exit("WhiteList path does noe exits: {}".format(white_path))
        with open(white_path, 'r') as infile:
            for line in infile:
                bc_rev = rev_comp(line.strip())
                bc_white[0].add(bc_rev)

        white_path = os.path.join(whilelistdir, "whitelist.indrop_2.txt")
        with open(white_path, 'r') as infile:
            for line in infile:
                bc_rev = rev_comp(line.strip())
                bc_white[1].add(bc_rev)

    elif protocol == "dropseq":
        bc_white = []

    else:
        sys.exits("Not valid Protocol Type entered...")
        raise Exception("Not valid Protocol Type used.")

    return bc_white

def check_whitelist(bc_white, protocol, barcode):
    """
    Check whitelist.

    Usage:
    ::
        check_whitelist(bc_white, "10X", "ATTATATATT")
    """
    if protocol == "10X":
        if barcode in bc_white[0]:
            return 1
        else:
            return 0

    if protocol in ["10X_14_10", "10X_14_5"]:
        if barcode in bc_white[0]:
            return 1
        else:
            return 0
    
    if protocol == "dropseq":
        if barcode in ['TCAAAAGCAGTG']:
            return 0
        else:
            return 1
    
    if protocol == "indrop":
        if barcode[0:-8] in bc_white[0] and barcode[-8:] in bc_white[1]:
            return 1
        else:
            return 0