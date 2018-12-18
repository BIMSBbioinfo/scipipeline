from os.path import join
import os
import pysam

from utils.assemble_pseudogenome import create_pseudo_genome
from utils.split_reads import split_reads_by_barcode
from utils.count_matrix import make_beds_for_intervalsize
from utils.count_matrix import sparse_count_reads_in_regions
from utils.count_matrix import make_barcode_table


PSGENOME_OUTDIR = join(OUT_DIR, 'barcodes')

def is_paired(wildcards):
    if os.path.exists(samples[samples.Name == wildcards.sample].read2.tolist()[0]):
        return True
    else:
        return False


def get_trim_inputs(wildcards):
    sam = samples[samples.Name == wildcards.sample]
    in_ = sam.read1.tolist()
    if is_paired(wildcards):
        in_ += sam.read2.tolist()
    return in_


def get_trimgalore_outputfilename(filename, extension):
    """Rename the original filename to get the
       trimgalore filename."""
    if filename.endswith('.fq') or filename.endswith('.fastq'):
        in_ = splitext(os.path.basename(filename))[0]
        in_ += extension + '.fq.gz'
    elif filename.endswith('.fq.gz') or filename.endswith('.fastq.gz'):
        filename = splitext(filename)[0] 
        in_ = splitext(os.path.basename(filename))[0]
        in_ += extension + '.fq.gz'
    return in_


def get_trimgalore_output(wildcards):
    sam = samples[samples.Name == wildcards.sample]

    if is_paired(wildcards):
        filename = sam.read1.tolist()[0]
        read1 = get_trimgalore_outputfilename(filename, '_val_1')
        filename = sam.read2.tolist()[0]
        read2 = get_trimgalore_outputfilename(filename, '_val_2')
        return [read1, read2]
    else:
        filename = sam.read1.tolist()[0]
        read1 = get_trimgalore_outputfilename(filename, '_trimmed')

    return [read1]

    
def get_mapping_inputs(wildcards):

    prefix = join(OUT_DIR, wildcards.sample, 'trimmed', 'sample' )

    if is_paired(wildcards):
        return [prefix + '_1.fastq.gz', prefix + '_2.fastq.gz']
    else:
        return [prefix + '.fastq.gz']


def bowtie_input_filetype_option(filename):

    # if the input is a fasta file, put th -f option
    if filename.endswith('.fa') or filename.endswith('.fasta') or \
      filename.endswith('.fa.gz') or filename.endswith('.fasta.gz'):
        return ' -f '
    else:
        # consider it as default
        return ''
