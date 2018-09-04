import os
from pysam import AlignmentFile


def make_genome_size(inbam, output):
    """
    This function derives the genome size table
    from the bam header.
    """

    treatment = AlignmentFile(inbam, 'rb')
    
    with open(output, "w") as out:
        for seq in treatment.header['SQ']:
            out.write("{}\t{}\n".format(seq['SN'], seq['LN']))


def remove_chroms(inbam, outbam, chroms):
    """
    This function takes a bam-file and outputs
    a bam-file in which the specified chromosomes
    have been removed.

    The function searches for matching chromosomes
    using regular expressions.
    For example, chroms=['chrM', '_random']
    would remove 'chrM' as well as all random chromsomes.
    E.g. chr1_KI270706v1_random.
    """

    treatment = AlignmentFile(inbam, 'rb')

    header = treatment.header
    new_chroms = []
    chrnames = []

    # tid_map is to reindex the chromosomes in the
    # new bam file.
    tid_map = [-1 for i in range(len(header['SQ']))]

    N = 0
    # make new header with valid chromsomes
    for i, seq in enumerate(header['SQ']):
        keep = True
        for chrom in chroms:
            if chrom in seq['SN']:
                keep = False
                break
        if keep:
            tid_map[i] = N
            N += 1
            new_chroms.append(seq)
            chrnames.append(seq['SN'])

    header['SQ'] = new_chroms

    bam_writer = AlignmentFile(outbam, 'wb', header=header)
#    bam_writer = AlignmentFile(outbam, 'wb', template=treatment)

    # write new bam files containing only valid chromosomes
    for aln in treatment.fetch(until_eof=True):
        if aln.is_unmapped:
            continue
        if aln.reference_name in chrnames:
            aln.tid = tid_map[aln.tid]
            bam_writer.write(aln)
        
    bam_writer.close()
    treatment.close()


def remove_low_mapq_reads(inbam, outbam, minmapq):
    """
    This function takes a bam file and it produces
    a new bam file with aligments with a minimum
    mapping quality.

    For paired-end data, only alignments where
    both mates exceed the threshold are retained.
    """

    treatment = AlignmentFile(inbam, 'rb')
    bam_writer = AlignmentFile(outbam, 'wb', template=treatment)

    waiting_for_pair = {}
    for aln in treatment.fetch(until_eof=True):
        if aln.is_paired:
            if aln.rname not in waiting_for_pair:
                # for the first mate that we encounter
                # we get here and keep the aln in the waiting list
                # if it is valid
                if aln.mapq >= minmapq and not aln.is_unmapped:
                    waiting_for_pair[aln.rname] = aln
                else:
                    waiting_for_pair[aln.rname] = None
            else:
                # for the second mate that we encounter
                # we get here
                if aln.mapq >= minmapq and not aln.is_unmapped \
                   and waiting_for_pair[aln.rname] is not None:
                    # both pairs satisfy the min mapq threshold
                    # and are mapped.
                    # write them into the output file
                    bam_writer.write(waiting_for_pair[aln.rname])
                    bam_writer.write(aln)
                
                # finally clear the waiting list to save memory
                waiting_for_pair.pop(aln.rname)
                    
        else:
            # single end
            if aln.maqp >= minmapq:
                bam_writer.write(aln)
    
    treatment.close()
    bam_writer.close()


def remove_low_cellcount_reads(inbam, outbam, mincount):
    """
    This function takes a bam file with barcodes in the 
    RG tag as input and outputs a bam file containing
    only barcodes that exceed the minimum number of aligments
    for a given barcode.

    """
    treatment = AlignmentFile(inbam, 'rb')
    header = treatment.header
    barcodecounts = {bc['ID']: 0 for bc in header['RG']}

    # first parse the file to determine the per barcode
    # alignment counts

    for aln in treatment.fetch(until_eof=True):
        if aln.is_proper_pair and aln.is_read1 or not aln.is_paired:
            rg = aln.get_tag('RG')
            barcodecounts[rg] += 1

    treatment.close()

    # make new header with the valid barcodes
    treatment = AlignmentFile(inbam, 'rb')
    header = treatment.header
    rgheader = []
    for rg in header['RG']:
       if barcodecounts[rg['ID']] >= mincount:
           rgheader.append(rg)

    header['RG'] = rgheader
    
    bam_writer = AlignmentFile(outbam, 'wb', header=header)
    for aln in treatment.fetch(until_eof=True):
        if barcodecounts[aln.get_tag('RG')] >= mincount:
            bam_writer.write(aln)

    treatment.close()
    bam_writer.close()
