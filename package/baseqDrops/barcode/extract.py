def HammingDistance(seq1, seq2):
    return sum([1 for x in zip(seq1, seq2) if x[0] != x[1]])

def extract_barcode(protocol, seq):
    """Extract cell barcode from reads 1 sequences, if protocol not defined, return blank.
        - 10X: seq[0:16]
        - indrop: seq[0:i] + seq[i + 22 : i + 22 + 8] (i is length of barcode 1)
        - dropseq: seq[0:12]

    Usages:
    ::
        from baseqdrops import extract_barcode
        extract_barcode("10X", "ATCGATCGATCGACTAAATTTTTTT")

        #You can speficy the UMI and barcode length
        extract_barcode("10X", "ATCGATCGATCGACTAAATTTTTTT", 14, 5)

    Returns:
        barcode sequence, if no valid barcode, return ""

    """
    if protocol == "10X":
        return seq[0:16]
    elif protocol in ["10X_14_10", "10X_14_5"]: 
        return seq[0:14]
    elif protocol == "dropseq":
        return seq[0:12]
    elif protocol == "indrop":
        w1 = "GAGTGATTGCTTGTGACGCCTT"
        if w1 in seq:
            w1_pos = seq.find(w1)
            if 7 < w1_pos < 12:
                return seq[0:w1_pos] + seq[w1_pos + 22:w1_pos + 22 + 8]
        else:
            for i in range(8, 12):
                w1_mutate = seq[i:i + 22]
                if HammingDistance(w1_mutate, w1) < 2:
                    return seq[0:i] + seq[i + 22 : i + 22 + 8]
                    break
        return ""
    else:
        return ""   

def extract_UMI(protocol, barcode, seq1, missing_base=False):
    """
    Extract UMI from reads 1 sequences. If the reads length too short, return blank.
    ::
        10X: 16-26
        10X_14_10, 10X_14_5: 14:19/ 14:24
        indrop: seq1[len(barcode) + 22:len(barcode) + 22 + 6]
        dropseq: 11-19/12-20
    Args:
        barcode: Provide the barcode sequences;
        seq1: The sequence of read1;
    Usage:
    ::
        from baseqdrop import extract_UMI
        extract_UMI("", "", "AA....TTT", False)
    """
    UMI = ""
    seq1 = seq1.strip()
    if protocol == "10X" and len(seq1)>=26:
        UMI = seq1[16:26]
    if protocol == "10X_14_5" and len(seq1)>=19:
        UMI = seq1[14:19]
    if protocol == "10X_14_10" and len(seq1)>=24:
        UMI = seq1[14:24]
    if protocol == "dropseq" and len(seq1)>=20:
        if missing_base:
            UMI = seq1[11:19]
        else:
            UMI = seq1[12:20]
    if protocol == "indrop" and len(seq1)>=len(barcode) + 22 + 6:
        UMI = seq1[len(barcode) + 22:len(barcode) + 22 + 6]
    return UMI
