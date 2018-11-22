import os

star_script = """
mkdir ${workdir}
cd ${workdir}
${STAR} --genomeLoad LoadAndKeep --genomeDir ${STAR_REF} \
--readFilesIn ${file} --runThreadN 10 --outSAMunmapped Within \
--outSAMtype BAM Unsorted --outFileNamePrefix ${name}_
${samtools} sort -n -@ 10 ${name}_Aligned.out.bam -o ${name}.sort.bam
#rm ${name}_Aligned.out.bam
"""

from .command import run_cmd
from .config import get_config
from string import Template

def star_align(bc_dir, workdir, sample, genome, parallel):
    """
    Run Star Alignments in parallel.
    Uasge:
    ::
        from baseqdrops.run_star import star_align
        star_align(bc_dir, workdir, sample, genome, parallel=4)

    Results:
    ::
        Aligned.AA.bam
    """
    import multiprocessing as mp
    pool = mp.Pool(processes=int(parallel))

    STAR = get_config("Drops", "star")
    cellranger_refs = get_config("Drops", "cellranger_ref_"+genome)
    
    try:
        STAR_REF = get_config("Drops", "star_ref_" + genome)
    except:
        STAR_REF = os.path.join(cellranger_refs, "star")
    
    samtools = get_config("Drops", "samtools")

    def work_script(name, file):
        return Template(star_script).substitute(
            STAR = STAR,
            STAR_REF = STAR_REF,
            samtools = samtools,
            workdir = os.path.join(workdir, name),
            name = name,
            file = os.path.abspath(file)
        )

    from itertools import product
    barcode_prefix = [x[0] + x[1] for x in list(product('ATCG', repeat=2))]

    for bc in barcode_prefix:
        fastq = os.path.join(bc_dir, "split.{}.{}.fq".format(sample, bc))
        script = work_script(bc, fastq)
        pool.apply_async(run_cmd, ("Star Alignment", script,))

    pool.close()
    pool.join()