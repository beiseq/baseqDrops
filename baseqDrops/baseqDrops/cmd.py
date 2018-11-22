import click, sys, os

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

@cli.command(short_help="Run full pipeline for inDrop/Drop-Seq/10X")
@click.option('--config', default="", help="path of configuration file")
@click.option('--genome', '-g', help="hg38/mm10/hg38_mm10")
@click.option('--protocol', '-p', type=click.Choice(['10X', 'indrop', 'dropseq', '10X_14_5', '10X_14_10']))
@click.option('--cells', default=5000, help='Max number of cells')
@click.option('--minreads', default=5000, help='minimum reads required for a barcode')
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--fq1', '-1', default='', help='fastq1 path')
@click.option('--fq2', '-2', default='', help='fastq2 path')
@click.option('--top_million_reads', default=1000, help='Use the top X M reads (1000M for default)')
@click.option('--parallel', '-t', default='4', help='Process N splitted at the same time...')
@click.option('--outdir', '-d', default='./', help='Folder to write (./)')
@click.option('--step', default='all', help='all/count/start/split/star/tagging/table/')

def run_pipe(config, genome, protocol, cells, minreads, name, fq1, fq2, outdir, top_million_reads, step, parallel):
    from .pipeline import pipeline
    pipeline(config, genome, protocol, int(cells), int(minreads), name, fq1, fq2, outdir, top_million_reads, step, parallel)

@cli.command(short_help="Subsampling the raw reads from 0.1-0.9 and export the UMIs and Reads table")
@click.option('--dir', '-d', default='./', help='Folder (./)')
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--ratio', default='10_100', help="Ratios for sampling (default is 10_100, can also be integer like 50)")
@click.option('--parallel', '-t', default='8', help='Process N tasks at the sametime')
def run_sampling(dir, name, ratio, parallel):
    from .sampling import runsampling
    runsampling(name, dir, ratio, parallel)