def get_cigar_segments(cigar, alignment_start):
    """
    Get Cigar Segments, return the regions of alignments in the read;
    """
    i = 0
    j = -1
    # get clipped blocks
    left_clip = []
    right_clip = []
    if cigar[i][0] == 5: # hardclip
        left_clip.append(cigar[i])
        i += 1
    if cigar[i][0] == 4: # softclip
        left_clip.append(cigar[i])
        i += 1
    if cigar[j][0] == 5: # hardclip
        right_clip.insert(0, cigar[j])
        j -= 1
    if cigar[j][0] == 4: # softclip
        right_clip.insert(0, cigar[j])
        j -= 1

    # get aligned blocks
    aligned_segments = []
    curr_start = alignment_start
    curr_end = curr_start
    curr_cigar = []
    for (op, length) in cigar[i:len(cigar)+j+1]:
        if op == 3: # splice junction -> start a new segment
            aligned_segments.append((curr_start, curr_end, curr_cigar))
            curr_end += length
            curr_start = curr_end
            curr_cigar = []
        else: # continue the current segment
            ref_bases = length if op != 1 else 0 # insertions don't consume ref bases
            curr_end += ref_bases
            curr_cigar.append((op, length))
    aligned_segments.append((curr_start, curr_end, curr_cigar))
    return aligned_segments