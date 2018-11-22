# # -*- coding: utf-8 -*-
# import os, bisect, sys, json
# import numpy as np
# from .bam_parser import get_cigar_segments
# from baseq.mgt.config import get_config

# def tagging_reads(genome, bam, outpath, prime5=False):
#     """
#     Tagging reads will transform the genomic position of reads to the gene name
#     Report the genes for a ...
#     ::
#         from baseq.drops import tagging_reads
#         tagging_reads("")

#     Todo:
#         Add support for 5 prime reads, in case of "dist_to_end"
#         Like: strand should be opposite...
#     """

#     #Remove the Already Exists File
#     if os.path.exists(outpath):
#         os.remove(outpath)

#     print("[info] File write to {}".format(outpath))

#     cellranger_refs = get_config("Drops", "cellranger_ref_"+genome)
#     star_index_dir = os.path.join(cellranger_refs, "star")

#     #read star index
#     chrom_name_path = os.path.join(star_index_dir, "chrName.txt")
#     chrom_start_path = os.path.join(star_index_dir, "chrStart.txt")
#     tx_path = os.path.join(star_index_dir, "transcriptInfo.tab")
#     ex_path = os.path.join(star_index_dir, "exonInfo.tab")

#     #gene transc
#     gene_transc = os.path.join(cellranger_refs, "genes_transcript.txt")

#     tx_dtype = {'names': ('tx_ids', 'tx_starts', 'tx_ends', 'tx_max_ends', 'tx_strands',
#                           'tx_num_exons', 'tx_break_idxs'),
#                 'formats': ('object', 'u8', 'u8', 'u8', 'u1', 'u2', 'u4')}
#     ex_dtype = {
#         'names': ('ex_starts', 'ex_ends', 'ex_cum_lengths'),
#         'formats': ('u8', 'u8', 'u8')
#     }
#     gene_tx_dtype = {
#         'names': ('gene', 'tx', 'genename'),
#         'formats': ('object', 'object', 'object')
#     }

#     chrom_names = np.loadtxt(chrom_name_path, dtype='object')
#     chrom_starts = np.loadtxt(chrom_start_path, dtype='u8')
#     tx_info = np.loadtxt(tx_path, tx_dtype, skiprows=1, unpack=True)
#     ex_info = np.loadtxt(ex_path, ex_dtype, skiprows=1, unpack=True)
#     ex_tx = np.loadtxt(gene_transc, gene_tx_dtype, unpack=True)

#     tx_starts = tx_info[1]
#     chrom_bins = [bisect.bisect_left(tx_starts, cs) for cs in chrom_starts]
#     chrom_info = chrom_names, chrom_starts, chrom_bins

#     ex_starts, ex_ends, ex_cum_lengths = ex_info
#     ex_breaks = np.empty((2 * ex_starts.size,), dtype=ex_starts.dtype)
#     ex_breaks[0::2] = ex_starts
#     ex_breaks[1::2] = ex_ends

#     tx2gene = {}
#     for idx in range(len(ex_tx[0])):
#         tx2gene[ex_tx[1][idx]] = ex_tx[2][idx]

#     transcript_ids, transcript_starts, transcript_ends, transcript_max_ends, transcript_strands, transcript_num_exons, transcript_exon_break_idx = tx_info

#     #This function returns the information: the genename, exon or intron, antisense?
#     #It should return a single transcript and Gene

#     def align_to_transcriptome(read_tid, read_pos, read_alen):
#         ref_offset = chrom_starts[read_tid]
#         clipped_read_start = ref_offset + read_pos
#         clipped_read_end = clipped_read_start + read_alen - 1
#         tx = bisect.bisect(transcript_starts, clipped_read_end) - 1

#         # read is at the extreme start / end of transcriptome
#         if tx == -1 or tx == transcript_starts.size:
#             return {}
#         tx_hits = {}

#         while tx >= 0 and clipped_read_start <= max(transcript_max_ends[tx], transcript_ends[tx]):
#             if clipped_read_start <= transcript_ends[tx]:
#                 tx_hits[transcript_ids[tx]] = {
#                     'start': clipped_read_start-transcript_starts[tx],
#                     'exons': transcript_num_exons[tx],
#                     'exonIdx': transcript_exon_break_idx[tx],
#                     'gene': tx2gene[transcript_ids[tx]],
#                     'strand': transcript_strands[tx]
#                 }
#             tx -= 1
#         return tx_hits


#     def align2exon(start, end, strand, ex_idx, ex_idx2):
#         ex_start = bisect.bisect(ex_breaks, start, lo=2 * ex_idx, hi=2 * (ex_idx2 + 1)) - 1
#         ex_end = bisect.bisect_left(ex_breaks, end, lo=max(ex_start, 0), hi=2 * (ex_idx2 + 1)) - 1
#         exon_len = 0
#         dist_to_end = 10000

#         if ex_start%2==0:
#             exon_len = min(ex_breaks[ex_start+1], end) - max(ex_breaks[ex_start], start)
#         elif ex_end%2==0:
#             exon_len = end - ex_breaks[ex_start+1]

#         ##Calculate The distance of ex_Start to the End of Genes...
#         ## We should Seprate by the Strand Of the Genes...

#         if strand == 1:
#             exon_to_geneend = ex_breaks[int(ex_start/2)*2 +2: 2 * (ex_idx2 + 1)]
#         else:
#             exon_to_geneend = ex_breaks[2 * ex_idx : int(ex_start / 2) * 2]

#         if len(exon_to_geneend)%2==0:
#             dist_to_end = sum(exon_to_geneend[1::2]) - sum(exon_to_geneend[::2])

#         return ex_start, ex_end, exon_len, dist_to_end

#     ##Read The Bamfile
#     ##Iterate the bam file

#     import pysam, itertools
#     from collections import Counter

#     bam_file = pysam.Samfile(bam, 'rb')
#     genome_bam_iter = itertools.groupby(bam_file, key = lambda read: read.qname)

#     res = {}
#     previous_barcode = ""

#     for qname, reads_iter in genome_bam_iter:
#         reads = list(reads_iter)
#         read_genes = {}
#         infos = qname.split('_')

#         barcode = infos[1]
#         UMI = infos[2]

#         if barcode != previous_barcode:
#             if previous_barcode == "":
#                 previous_barcode = barcode
#                 continue
#             res2 = {}
#             for gene in res:
#                 res2[gene] = Counter(res[gene])
#             res = {}
#             with open(outpath, 'a') as outfile:
#                 outfile.write(previous_barcode + '\t')
#                 json.dump(res2, outfile)
#                 outfile.write('\n')

#             previous_barcode = barcode

#         for read in reads:
#             #Get Gene For A single Alignment
#             #For each Alignment Of this Reads, we need to Get the gene Corresponed
#             #If there are multiple Gene, return None
#             gene_info = {}
#             rlen = read.rlen

#             if not read.is_unmapped:
#                 #There are multiple transcript For this Read
#                 #For each transcript, Get the Exon Coverage
#                 #If there are only One Gene, return the gene name if exon overlap ratio >50%
#                 #If there are multiple Genes

#                 tx_hits = align_to_transcriptome(read.tid, read.pos, read.alen)

#                 #Filter the tx_hits using the Strand Information
#                 #If the Read is is_reversed, the strand of transcript should be 2


#                 #For prime3
#                 required_strand = 1
#                 if read.is_reverse:
#                     required_strand = 2

#                 #For prime5
#                 required_strand = 2
#                 if read.is_reverse:
#                     required_strand = 1

#                 # tx_hits:
#                 tx_hits = {k:v for k, v in tx_hits.items() if v['strand']==required_strand}
#                 strand_info = {x['gene']:x['strand'] for x in tx_hits.values()}

#                 #gene_info: initiate the informations for each gene in the tx_hits
#                 gene_info = {x['gene']:{'overlap_len':0, 'dist_to_end':10000} for x in tx_hits.values()}

#                 for tx in tx_hits:
#                     gene = tx_hits[tx]['gene']
#                     # If the gene already have good alignment, skip
#                     if gene_info[gene]['overlap_len']>= 0.9 * rlen:
#                         continue

#                     alns = get_cigar_segments(read.cigar, tx_hits[tx]['start'])
#                     read_exon_overlap = 0

#                     for aln in alns:
#                         read_exon = align2exon(aln[0], aln[1], tx_hits[tx]['strand'],
#                                                tx_hits[tx]['exonIdx'], tx_hits[tx]['exonIdx'] + tx_hits[tx]['exons'] - 1)
#                         read_exon_overlap += read_exon[2]

#                     gene_info[gene]['overlap_len'] = max(read_exon_overlap, gene_info[gene]['overlap_len'])
#                     gene_info[gene]['dist_to_end'] = read_exon[3]

#                 gene_info = {k:v for k, v in gene_info.items() if v['overlap_len'] >= 0.5*rlen}

#             for k, v in gene_info.items():
#                 read_genes[k] = v

#         if len(list(read_genes)) == 1:
#             mapgene = list(read_genes)[0]
#         else:
#             #gene_info = {k: v for k, v in read_genes.items() if v['dist_to_end'] <= 400}

#             #For 5prime reads
#             gene_info = {k: v for k, v in read_genes.items()}
#             gene_names = [x[0] for x in sorted(list(gene_info.items()), key=lambda val: val[1]['dist_to_end'])]

#             #For 5 prime reads
#             gene_names = gene_names[::-1]
#             gene_names.append("")
#             mapgene = gene_names[0]

#         if not mapgene:
#             continue

#         if mapgene in res:
#             res[mapgene].append(UMI)
#         else:
#             res[mapgene] = [UMI]

#     bam_file.close()

#     res2 = {}
#     for gene in res:
#         res2[gene] = Counter(res[gene])

#     with open(outpath, 'a') as outfile:
#         outfile.write(barcode + '\t')
#         json.dump(res2, outfile)
#         outfile.write('\n')