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


def is_paired(wildcards):
    if os.path.exists(samples[samples.Name == wildcards.sample].read2.tolist()[0]):
        return True
    else:
        return False


def get_trim_inputs(wildcards):
    samples = samples[samples.Name == wildcards.sample]
    in_ = samples[samples.Name == wildcards.sample].read1.tolist()
    if is_paired(wildcards):
        in_ += samples[samples.Name == wildcards.sample].read2.tolist()

    return in_


def get_mapping_inputs(wildcards):
    samples = samples[samples.Name == wildcards.sample]
    if not config['trim_reads']:
        # trimming is not required.
        # the raw reads will be passed to the mapper
        # directly
        return get_trim_inputs(wildcards)

    prefix = join(OUT_DIR, wildcards.sample + '_trimmed')

    if is_paired(wildcards):
        return [prefix + '_1.fastq', prefix + '_2.fastq']
    else:
        return [prefix + '.fastq']



def bowtie_input_filetype_option(filename):
    # if the input is a fasta file, put th -f option
    if filename.endswith('.fa') or filename.endswith('.fasta') or \
      filename.endswith('.fa.gz') or filename.endswith('.fasta.gz'):
        return ' -f '
    else:
        # consider it as default
        return ''
