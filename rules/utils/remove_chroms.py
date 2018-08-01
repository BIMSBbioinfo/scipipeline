import os
from pysam import AlignmentFile


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
    new_header = []
    chrnames = []
    # make new header with valid chromsomes
    for seq in header['SQ']:
        keep = True
        for chrom in chroms:
            if chrom in seq['SN']:
                keep = False
                break
        if keep:
            new_chroms.append(seq)
            chrnames.append(seq['SN'])

    header['SQ'] = new_chroms

    bam_writer = AlignmentFile(outbam, 'wb', header=header)

    # write new bam files containing only valid chromosomes
    for aln in treatment.fetch(until_eof=True):
        if aln.reference_name in chrnames:
            bam_writer.write(aln)
        
    treatment.close()
    bam_writer.close()

