import os
import h5py
import numpy as np
from pysam import AlignmentFile
from pysam import view
from collections import defaultdict


def split_reads_by_barcode(barcode_bams, treatment_bam,
                           output_bam, min_mapq=None, max_mismatches=1):
    """Splits the reads by barcodes.

    This function takes a set of barcode alignments (in bam format)
    and reads (in bam format). It produces an output bam file
    containing the RG tag that hold the barcode.

    Note: The method assumes that the respective BAM files contain reads
    with corresponding read ids. Also the BAM files need to be sorted
    by read name. Furthermore, the treatment_reads must be a subset of
    the barcode reads.

    barcodes : list(str)
        List of bam-file names pointing to the barcode alignments
        towards the pseudo genomes.
    reads : str
        BAM file containing the read alignments
    output_bam : str
        Output BAM file
    max_open_files : int
        Maximum number of writers to handle in one sweep. This is limited by the
        operating system's limitation opening file handles. Default: 1000.
    min_mapq : int
        Filter for mapping quality of the barcode reads. Only reads mapq >= min_mapq
        are considered. Default: None means no filter is applied.
    max_mismatches : int
        Maximum number of mismatches. Default: 1.
    """

    aligned_cnt = 0
    unaligned_cnt = 0
    if not min_mapq:
        min_mapq = 0

    treatment_reader = AlignmentFile(treatment_bam, 'rb')

    barcode_readers = [AlignmentFile(bamfile, 'rb') for bamfile in barcode_bams]

    # start at the beginning of the bam files
    barcode_it = [reader.fetch(until_eof=True)
                  for reader in barcode_readers]
    bnames = [next(br) for br in barcode_it]
    bam_writer = AlignmentFile(output_bam + '.tmp', 'wb', template=treatment_reader)

    barcodes = set()

    for aln in treatment_reader.fetch(until_eof=True):
        # only retain the aligned reads
        if aln.is_unmapped:
            unaligned_cnt += 1
            continue
        aligned_cnt += 1

        # extract the corresponding barcodes
        for i, br_it in enumerate(barcode_it):
            while not bnames[i].query_name == aln.query_name:
                # increment barcode names until they
                # match the alignment read name.
                bnames[i] = next(br_it)

        if any([x.is_unmapped or x.is_reverse for x in bnames]):
            # one or more barcodes are abscent
            # Barcode must map to the forward strand only
            # Reverse strand matches are definitely due to sequencing errors
            continue

        if any([x.mapq < min_mapq for x in bnames]):
            # remove poor mapping quality
            continue

        if any([x.get_tag('XM') > max_mismatches for x in bnames]):
            # maximum number of mismatches exceeded
            continue

        comb_id = '_'.join([baln.reference_name for baln in bnames])
        barcodes.add(comb_id)
        
        aln.set_tag('RG', comb_id)
        bam_writer.write(aln)

    print("End batch ...")

    treatment_reader.close()
    bam_writer.close()
    for reader in barcode_readers:
        reader.close()

    f = AlignmentFile(output_bam + '.tmp', 'rb')
    header = f.header
    header['RG'] = [{'ID': combid} for combid in barcodes]
    bam_writer = AlignmentFile(output_bam, 'wb', header=header)
    for aln in f.fetch(until_eof=True):
        bam_writer.write(aln)

    f.close()
    os.remove(output_bam + '.tmp')
    bam_writer.close()


def deduplicate_reads_by_barcode(bamin, bamout):
    """This script deduplicates the original bamfile.
    
    It assumes the original bamfile to be sorted according
    to genomic positions. Furthermore, each read must be augmented
    by a corresponding RG tag indicating the barcode.

    Parameters
    ----------
    bamfile : str
        Sorted bamfile containing barcoded reads.
    output : str
        Output path to a bamfile that contains the deduplicated reads.
    """
    bamfile = AlignmentFile(bamin, 'rb')
    output = AlignmentFile(bamout, 'wb', template=bamfile)

    # grep all barcodes from the header
    barcodes = defaultdict(lambda: 0)
    
    for aln in bamfile.fetch(until_eof=True): 
        # if previous hash matches the current has
        # skip the read
        if barcodes[aln.get_tag('RG')] == hash((aln.query_name, aln.pos, aln.flag, aln.tlen)):
            continue
        else:
            barcodes[aln.get_tag('RG')] = hash((aln.query_name, aln.pos, aln.flag, aln.tlen))
            output.write(aln)
    
if __name__ == '__main__':
    barcode_bams = ['pseudo_genome_I{}_sorted.bam'.format(i) for i in [1,2,3,4]]
    treatment_bam = 'test.input.bam'
    output_dir = 'test_split'
    max_open_files=1000
    split_reads_by_barcode(barcode_bams, treatment_bam, output_dir)
