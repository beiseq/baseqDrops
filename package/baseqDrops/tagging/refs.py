import os
import bisect
import numpy as np

class reference:
    """
    About the references of cellranger for counting the 10X data.
    ::
        from baseq.drops import reference
        refs = reference("XXXX")

    Todo:
        目前函数抽离的工作还没有work。
        类还无法使用。
    """
    def __init__(self, star_index_dir, cellranger_refs):
        #read star index
        self.read_index()

    def read_index(self, cellranger_refs, star_index_dir):
        """
        Read the cellranger and star index files. It has the following properities:
        ::
            chrom_names: []
            chrom_starts: [tx_ids', 'tx_starts', 'tx_ends', 'tx_max_ends', 'tx_strands', 'tx_num_exons', 'tx_break_idxs']
            tx_info: []
            ex_info: ['ex_starts', 'ex_ends', 'ex_cum_lengths']
            ex_tx: ['gene', 'tx', 'genename']
        """
        chrom_name_path = os.path.join(star_index_dir, "chrName.txt")
        chrom_start_path = os.path.join(star_index_dir, "chrStart.txt")
        tx_path = os.path.join(star_index_dir, "transcriptInfo.tab")
        ex_path = os.path.join(star_index_dir, "exonInfo.tab")

        gene_transc = os.path.join(cellranger_refs, "genes_transcript.txt")

        tx_dtype = {'names': ('tx_ids', 'tx_starts', 'tx_ends', 'tx_max_ends', 'tx_strands',
                              'tx_num_exons', 'tx_break_idxs'),
                    'formats': ('object', 'u8', 'u8', 'u8', 'u1', 'u2', 'u4')}
        ex_dtype = {
            'names': ('ex_starts', 'ex_ends', 'ex_cum_lengths'),
            'formats': ('u8', 'u8', 'u8')
        }
        gene_tx_dtype = {
            'names': ('gene', 'tx', 'genename'),
            'formats': ('object', 'object', 'object')
        }
        self.chrom_names = np.loadtxt(chrom_name_path, dtype='object')
        self.chrom_starts = np.loadtxt(chrom_start_path, dtype='u8')
        self.tx_info = np.loadtxt(tx_path, tx_dtype, skiprows=1, unpack=True)
        self.ex_info = np.loadtxt(ex_path, ex_dtype, skiprows=1, unpack=True)
        self.ex_tx = np.loadtxt(gene_transc, gene_tx_dtype, unpack=True)

        tx_starts = tx_info[1]
        chrom_bins = [bisect.bisect_left(tx_starts, cs) for cs in chrom_starts]
        chrom_info = chrom_names, chrom_starts, chrom_bins

        ex_starts, ex_ends, ex_cum_lengths = ex_info
        ex_breaks = np.empty((2 * ex_starts.size,), dtype=ex_starts.dtype)
        ex_breaks[0::2] = ex_starts
        ex_breaks[1::2] = ex_ends

    def align_to_transcriptome(self, chrID, start, alength):
        """
        Get the overlap transcription names.
        ::
            refs.get_overlap_transc("chr1", start, end)
            > ['ENST....', 'ENST...']
        """
        ref_offset = self.chrom_starts[chrID]
        clipped_read_start = ref_offset + start
        clipped_read_end = clipped_read_start + alength - 1
        tx = bisect.bisect(self.transcript_starts, clipped_read_end) - 1

        # read is at the extreme start / end of transcriptome
        if tx == -1 or tx == self.transcript_starts.size:
            return {}
        tx_hits = {}

        while tx >= 0 and clipped_read_start <= max(transcript_max_ends[tx], transcript_ends[tx]):
            if clipped_read_start <= transcript_ends[tx]:
                tx_hits[transcript_ids[tx]] = {
                    'start': clipped_read_start-transcript_starts[tx],
                    'exons': transcript_num_exons[tx],
                    'exonIdx': transcript_exon_break_idx[tx],
                    'gene': tx2gene[transcript_ids[tx]],
                    'strand': transcript_strands[tx]
                }
            tx -= 1
        return tx_hits

    def overlap_exons(self, chr, start, end):
        """
        Get the overlap exons.
        """
        pass