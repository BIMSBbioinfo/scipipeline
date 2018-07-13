from os.path import join
import os
import pysam

from utils.assemble_pseudogenome import create_pseudo_genome
from utils.split_reads import split_reads_by_barcode
from utils.split_reads import deduplicate_reads
from utils.count_matrix import make_beds_for_intervalsize
from utils.count_matrix import sparse_count_reads_in_regions


PSGENOME_OUTDIR = join(OUT_DIR, 'barcodes')


# Trimmed reads
TRIM_PATTERN = join(OUT_DIR, '{sample}_trimmed')


def get_trim_inputs(wildcards):
    samples = config['samples'][wildcards.sample]
    return [samples[x] for x in samples]


def get_mapping_inputs(wildcards):
    samples = config['samples'][wildcards.sample]
    if not config['trim_reads']:
        # trimming is not required.
        # the raw reads will be passed to the mapper
        # directly
        return [samples[x] for x in samples]

    prefix = join(OUT_DIR, wildcards.sample + '_trimmed')

    if 'read2' in samples:
        return [prefix + '_1.fastq', prefix + '_2.fastq']
    else:
        return [prefix + '.fastq']

def is_paired(wildcards):
    if 'read2' in config['samples'][wildcards.sample]:
        return True
    else:
        return False


def bowtie_input_filetype_option(filename):
    # if the input is a fasta file, put th -f option
    if filename.endswith('.fa') or filename.endswith('.fasta') or \
      filename.endswith('.fa.gz') or filename.endswith('.fasta.gz'):
        return ' -f '
    else:
        # consider it as default
        return ''
