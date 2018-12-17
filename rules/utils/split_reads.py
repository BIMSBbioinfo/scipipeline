import os
import numpy as np
import pandas as pd
from pysam import AlignmentFile
from pysam import view
from collections import defaultdict


def augment_alignment_by_barcode_from_name(inbam, outbam, reftable):
    """
    This function takes a bam-file and outputs
    a bam-file with RG-tag representing the barcodes.

    The barcodes are encoded by the read name
    where the the read name is assumed to have the form:
    @<barcodename>:<number> ...

    The barcode name from the read is then used to look up
    the corresponding barcode name from a reference table
    that has the form
    
    > cat barcodesheet.tsv
    readpre1    cell1
    readpre2    cell2

    Therefore, for a the read @readpre1:12312 the corresponding
    RGID would be cell1.

    Importantly, the readprefixes that encode the barcode in the
    read names must be equally long!
    """

    # load mapping between read name prefix and barcode
    refs = pd.read_csv(reftable, sep='\t', header=[0])
    # obtain the prefix length and check if it is the same 
    # for all entries
    reflen = len(refs['Readprefix'][0])
    rmap = {}
    for row in refs.iterrows():
        rmap[row[1].Readprefix] = row[1].Name

    treatment_reader = AlignmentFile(inbam, 'rb')
    bam_writer = AlignmentFile(outbam + '.tmp', 'wb', template=treatment_reader)

    barcodes = set()
    for aln in treatment_reader.fetch(until_eof=True):
        # extract barcode between @ and the first :
        refprefix = aln.query_name[:reflen]
        barcode = rmap[refprefix]
        barcodes.add(barcode)

        aln.set_tag('RG', barcode)
        bam_writer.write(aln)

    treatment_reader.close()
    bam_writer.close()

    # update the header with the available barcodes
    f = AlignmentFile(outbam + '.tmp', 'rb')
    header = f.header.to_dict().copy()
    header['RG'] = [{'ID': bc, 'SM':bc} for bc in barcodes]
    bam_writer = AlignmentFile(outbam, 'wb', header=header)
    for aln in f.fetch(until_eof=True):
        bam_writer.write(aln)

    f.close()
    os.remove(outbam + '.tmp')
    bam_writer.close()


def split_reads_by_barcode(barcode_bams, treatment_bam,
                           output_bam, min_mapq=None, max_mismatches=1, log=None):
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
    log_content = {}
    log_content['with_barcode'] = 0
    log_content['no_barcode'] = 0
    log_content['unmapped'] = 0

    for aln in treatment_reader.fetch(until_eof=True):
        # only retain the aligned reads
        if aln.is_unmapped:
            log_content['unmapped'] += 1
            continue

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
            log_content['no_barcode'] += 1
            continue

        if any([x.mapq < min_mapq for x in bnames]):
            # remove poor mapping quality
            log_content['no_barcode'] += 1
            continue

        if any([x.get_tag('XM') > max_mismatches for x in bnames]):
            # maximum number of mismatches exceeded
            log_content['no_barcode'] += 1
            continue

        comb_id = '_'.join([baln.reference_name for baln in bnames])
        barcodes.add(comb_id)

        aln.set_tag('RG', comb_id)
        bam_writer.write(aln)
        log_content['with_barcode'] += 1

    print("End batch ...")

    treatment_reader.close()
    bam_writer.close()
    for reader in barcode_readers:
        reader.close()

    f = AlignmentFile(output_bam + '.tmp', 'rb')
    header = f.header.to_dict().copy()
    header['RG'] = [{'ID': combid, 'SM': combid} for combid in barcodes]
    bam_writer = AlignmentFile(output_bam, 'wb', header=header)
    for aln in f.fetch(until_eof=True):
        bam_writer.write(aln)

    f.close()
    os.remove(output_bam + '.tmp')
    bam_writer.close()

    #write log file
    with open(log, 'w') as f:
        f.write('Readgroup\tcounts\n')
        for icnt in log_content:
            f.write('{}\t{}\n'.format(icnt, log_content[icnt]))


def deduplicate_reads(bamin, bamout, by_rg=True):
    """This script deduplicates the original bamfile.

    Deduplication removes reads align to the same position.
    If the reads in the bamfile contain a RG tag and
    by_rg=True, deduplication is done for each group separately.

    Parameters
    ----------
    bamfile : str
        Sorted bamfile containing barcoded reads.
    output : str
        Output path to a bamfile that contains the deduplicated reads.
    by_rg : boolean
        If True, the reads will be split by group tag.
    """
    bamfile = AlignmentFile(bamin, 'rb')
    output = AlignmentFile(bamout, 'wb', template=bamfile)

    # grep all barcodes from the header
    #barcodes = set()
    last_barcode = {}

    for aln in bamfile.fetch():
        # if previous hash matches the current has
        # skip the read
        val = (aln.reference_id, aln.reference_start,
               aln.is_reverse, aln.tlen)
        if aln.has_tag('RG') and by_rg:
            rg = aln.get_tag('RG')
        else:
            rg = 'dummy'

        if rg not in last_barcode:
            output.write(aln)
            # clear dictionary
            last_barcode[rg] = val
                
        if val == last_barcode[rg]:
            continue
        else:
            output.write(aln)
            last_barcode[rg] = val

if __name__ == '__main__':
    barcode_bams = ['pseudo_genome_I{}_sorted.bam'.format(i) for i in [1,2,3,4]]
    treatment_bam = 'test.input.bam'
    output_dir = 'test_split'
    split_reads_by_barcode(barcode_bams, treatment_bam, output_dir)
