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


def remove_chroms(inbam, outbam, rmchroms, log):
    """
    This function takes a bam-file and outputs
    a bam-file in which the specified chromosomes
    have been removed.

    The function searches for matching chromosomes
    using regular expressions.
    For example, rmchroms=['chrM', '_random']
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

    chr_to_remove_reason = {}
    # make new header with valid chromsomes
    for i, seq in enumerate(header['SQ']):
        keep = True
        for chrom in rmchroms:
            if chrom in seq['SN']:
                keep = False
                chr_to_remove_reason[seq['SN']] = chrom
                break
        if keep:
            tid_map[i] = N
            N += 1
            new_chroms.append(seq)
            chrnames.append(seq['SN'])

    new_header = {'SQ': new_chroms}

    bam_writer = AlignmentFile(outbam, 'wb', header=new_header)

    log_content = {chrom: 0 for chrom in rmchroms}
    log_content['remaining'] = 0
    log_content['unmapped'] = 0
    log_content['total'] = 0

    # write new bam files containing only valid chromosomes
    for aln in treatment.fetch(until_eof=True):
        log_content['total'] += 1
        if aln.is_unmapped:
            log_content['unmapped'] += 1
            continue
        if aln.reference_name in chrnames:
            aln.reference_id = tid_map[aln.reference_id]
            if aln.is_paired and aln.is_proper_pair:
                aln.next_reference_id = tid_map[aln.next_reference_id]
            bam_writer.write(aln)
            log_content['remaining'] += 1
        else:
            log_content[chr_to_remove_reason[aln.reference_name]] += 1

    bam_writer.close()
    treatment.close()

    #write log file
    with open(log, 'w') as f:
        f.write('Readgroup\tcounts\n')
        for icnt in log_content:
            f.write('{}\t{}\n'.format(icnt, log_content[icnt]))


def remove_low_mapq_reads(inbam, outbam, minmapq, log):
    """
    This function takes a bam file and it produces
    a new bam file with aligments with a minimum
    mapping quality.

    For paired-end data, only alignments where
    both mates exceed the threshold are retained.
    """

    treatment = AlignmentFile(inbam, 'rb')
    bam_writer = AlignmentFile(outbam, 'wb', template=treatment)

    log_content = {}
    log_content['below_mapq'] = 0
    log_content['above_mapq'] = 0
    log_content['total'] = 0

    waiting_for_pair = {}
    for aln in treatment.fetch(until_eof=True):
        log_content['total'] += 1
        if aln.is_paired:
            if aln.rname not in waiting_for_pair:
                # for the first mate that we encounter
                # we get here and keep the aln in the waiting list
                # if it is valid
                if aln.mapq >= minmapq and not aln.is_unmapped:
                    waiting_for_pair[aln.rname] = aln
                else:
                    # None marks an invalid first mate
                    waiting_for_pair[aln.rname] = None
                    log_content['below_mapq'] += 1
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
                    log_content['above_mapq'] += 2
                else:
                    # either the first mate was invalid
                    # or the second mate was below mapq
                    log_content['below_mapq'] += 1

                # finally clear the waiting list to save memory
                waiting_for_pair.pop(aln.rname)

        else:
            # single end
            if aln.maqp >= minmapq:
                bam_writer.write(aln)
                log_content['above_mapq'] += 1
            else:
                log_content['below_mapq'] += 1

    treatment.close()
    bam_writer.close()

    #write log file
    with open(log, 'w') as f:
        f.write('Readgroup\tcounts\n')
        for icnt in log_content:
            f.write('{}\t{}\n'.format(icnt, log_content[icnt]))



def remove_low_cellcount_reads(inbam, outbam, mincount, log):
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
    header = treatment.header.to_dict().copy()
    rgheader = []
    for rg in header['RG']:
       if barcodecounts[rg['ID']] >= mincount:
           rgheader.append(rg)

    header['RG'] = rgheader

    #log summary
    log_content = {}
    log_content['below_minbarcodecounts'] = 0
    log_content['above_minbarcodecounts'] = 0
    log_content['total'] = 0
    for bc in barcodecounts:
        log_content['total'] += barcodecounts[bc]
        if barcodecounts[bc] >= mincount:
            log_content['above_minbarcodecounts'] += barcodecounts[bc]
        else:
            log_content['below_minbarcodecounts'] += barcodecounts[bc]

    bam_writer = AlignmentFile(outbam, 'wb', header=header)
    for aln in treatment.fetch(until_eof=True):
        if barcodecounts[aln.get_tag('RG')] >= mincount:
            bam_writer.write(aln)

    treatment.close()
    bam_writer.close()

    #write log file
    with open(log, 'w') as f:
        f.write('Readgroup\tcounts\n')
        for icnt in log_content:
            f.write('{}\t{}\n'.format(icnt, log_content[icnt]))
