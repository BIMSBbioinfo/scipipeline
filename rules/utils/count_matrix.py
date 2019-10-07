import numpy as np
import pandas as pd
from scipy.sparse import dok_matrix, coo_matrix
from pysam import AlignmentFile
from pybedtools import BedTool
#from HTSeq import BED_Reader

def make_beds_for_intervalsize(bamfile, binsize, storage):
    """ Genome intervals for binsize.

    For a given genome and binsize,
    this function creates a bed-file containing all intervals.
    """
    # Obtain the header information
    afile = AlignmentFile(bamfile, 'rb')

    # extract genome size

    genomesize = {}
    for chrom, length in zip(afile.references, afile.lengths):
        genomesize[chrom] = length
        #offset += length//binsize + 1
    print('found {} chromosomes'.format(len(genomesize)))
    bed_content = pd.DataFrame(columns=['chr', 'start', 'end'])

    for chrom in genomesize:
        nbins = genomesize[chrom]//binsize + 1 if (genomesize[chrom] % binsize > 0) else 0
        starts = [int(i*binsize) for i in range(nbins)]
        ends = [min(int((i+1)*binsize), genomesize[chrom]) for i in range(nbins)]
        chr_ = [chrom] * nbins
        cont = {'chr': chr_, 'start': starts, 'end': ends}
        bed_entry = pd.DataFrame(cont)
        bed_content = bed_content.append(bed_entry, ignore_index=True, sort=False)
    bed_content.to_csv(storage, sep='\t', header=False, index=False,
                       columns=['chr', 'start', 'end'])


def make_barcode_table(inbam, barcodetable):
    """ This function creates a table of barcodes derived
    from the bam-file's header.

    Parameters
    ----------
    inbam : str
        Path to input bam file.
    barcodetable : str
        Output path to barcode table.
    """

    afile = AlignmentFile(bamfile, 'rb')

    if 'RG' in afile.header:
        use_group = True
    else:
        use_group = False

    if use_group:
        # extract barcodes
        barcodes = [item['ID'] for item in afile.header['RG']]
    else:
        barcodes = ['dummy']

    df = pd.DataFrame({'barcodes':barcodes})

    df.to_csv(barcodetable, sep='\t', header=True, index=True)


def sparse_count_reads_in_regions(bamfile, regions, storage, flank=0, log=None,
                                  template_length=1000,
                                  count_both_ends=False):
    """ This function obtains the counts per bins of equal size
    across the genome.

    The function automatically extracts the genome size from the
    bam file header.
    If group tags are available, they will be used to extract
    the indices from.
    Finally, the function autmatically detects whether the bam-file
    contains paired-end or single-end reads.
    Paired-end reads are counted once at the mid-point between the two
    pairs while single-end reads are counted at the 5' end.
    For paired-end reads it is optionally possible to count both read ends
    by setting count_both_ends=True.

    Parameters
    ----------
    bamfile :  str
        Path to a bamfile. The bamfile must be indexed.
    regions : str
        BED or GFF file containing the regions of interest.
    storage : str
        Path to the output hdf5 file, which contains the counts per chromsome.
    flank : int
        Extension of the regions in base pairs. Default: 0
    template_length : int
        Assumed template length. This is used when counting paired-end reads
        at the mid-point and the individual reads do not overlap with
        the given region, but the mid-point does.
    count_both_ends : bool
        Indicates whether for paired-end sequences, the ends of both mates should
        be counted separately. Default: False.
    """

    # Obtain the header information
    afile = AlignmentFile(bamfile, 'rb')

    # extract genome lengths
    if log is not None:
        f = open(log, 'w')
        fwrite = f.write
    else:
        fwrite = print
    fwrite('Make countmatrix from region\n')
    fwrite('bamfile: {}\n'.format(bamfile))
    fwrite('bedfile: {}\n\n'.format(regions))
    fwrite('get genomesize\n')
    # extract genome size
    genomesize = {}
    for chrom, length in zip(afile.references, afile.lengths):
        genomesize[chrom] = length
    fwrite('found {} chromosomes'.format(len(genomesize)))

    nreg = 0
    regfile = BedTool(regions)
    
    nreg = len(regfile)

    fwrite('number of regions to collect counts from: {}'.format(nreg))

    if 'RG' in afile.header:
        use_group = True
    else:
        use_group = False

    # get barcodes from header
    barcodes = {}
    if use_group:
        # extract barcodes
        for idx, item in enumerate(afile.header['RG']):
            barcodes[item['ID']] = idx
    else:
        barcodes['dummy'] = 0
    fwrite('found {} barcodes'.format(len(barcodes)))

    # barcode string for final table
    barcode_string = ';'.join([item['ID'] for item in afile.header['RG']])

    sdokmat = dok_matrix((nreg, len(barcodes)), dtype='int32')
    nbarcode_inregions = {key: 0 for key in barcodes}

    if count_both_ends:
        # if both ends are counted, template_length is irrelevant
        tlen = 0
    else:
        tlen = template_length

    for idx, iv in enumerate(regfile):

        iv.start -= flank
        iv.end += flank

        if iv.chrom not in genomesize:
            # skip over peaks/ regions from chromosomes
            # that are not contained in the bam file
            continue

        fetchstart = max(iv.start - tlen, 0)
        fetchend =  min(iv.end + tlen, genomesize[iv.chrom])

        for aln in afile.fetch(iv.chrom, fetchstart, fetchend):
            if aln.is_proper_pair and aln.is_read1 and not count_both_ends:

                pos = min(aln.reference_start, aln.next_reference_start)

                # count paired end reads at midpoint
                midpoint = pos + abs(aln.template_length)//2
                if midpoint >= iv.start and midpoint < iv.end:
                   sdokmat[idx, barcodes[aln.get_tag('RG') if use_group else 'dummy']] += 1
                   nbarcode_inregions[aln.get_tag('RG') if use_group else 'dummy'] += 1

            if not aln.is_paired or count_both_ends:
                # count single-end reads at 5p end
                if not aln.is_reverse:
                    if aln.reference_start >= iv.start and aln.reference_start < iv.end:
                        sdokmat[idx,
                        barcodes[aln.get_tag('RG') if use_group else 'dummy']] += 1
                else:
                    if aln.reference_start + aln.reference_length - 1 >= iv.start and \
                       aln.reference_start + aln.reference_length - 1 < iv.end:
                        sdokmat[idx,
                        barcodes[aln.get_tag('RG') if use_group else 'dummy']] += 1
                nbarcode_inregions[aln.get_tag('RG') if use_group else 'dummy'] += 1

    afile.close()

    fwrite('sparse matrix shape: {}'.format(sdokmat.shape))
    fwrite('density: {}'.format(sdokmat.nnz/np.prod(sdokmat.shape)))

    # store the results in COO sparse matrix format
    spcoo = sdokmat.tocoo()
    # sort lexicographically

    order_ = np.lexsort((spcoo.col, spcoo.row))
    indices = np.asarray([x for x in zip(spcoo.row, spcoo.col)], dtype=np.int64)[order_]
    values = spcoo.data.astype(np.float32)[order_]
    cont = {'region': indices[:,0], 'cell': indices[:, 1], 'count': values}

    df = pd.DataFrame(cont)
    with open(storage, 'w') as title:
        title.write('#  ' + barcode_string + '\n')

    df.to_csv(storage, mode = 'a', sep='\t', header=True, index=False)
        
    #main output file

    names = [key for key in barcodes]
    counts = [nbarcode_inregions[key] for key in barcodes]

    df = pd.DataFrame({'barcodes':names, 'counts':counts})

    df.to_csv(storage + '.counts', sep='\t', header=True, index=False)
    fwrite('total number of tags with barcodes: {}'.format(df.counts.sum()))
    if log is not None:
        f.close()


def get_barcode_frequency_genomewide(bamfile, storage):
    """ This function obtains the barcode frequency
    and stores it in a table.

    Parameters
    ----------
    bamfile :  str
        Path to a bamfile. The bamfile must be indexed.
    storage : str
        Path to the output hdf5 file, which contains the counts per chromsome.
    """

    # Obtain the header information
    afile = AlignmentFile(bamfile, 'rb')

    if 'RG' in afile.header:
        use_group = True
    else:
        use_group = False

    barcodes = {}
    if use_group:
        # extract barcodes
        for idx, item in enumerate(afile.header['RG']):
            barcodes[item['ID']] = 0
    else:
        barcodes['dummy'] = 0
    print('found {} barcodes'.format(len(barcodes)))

    for aln in afile.fetch(until_eof=True):
        if aln.is_proper_pair and aln.is_read1:
            barcodes[aln.get_tag('RG') if use_group else 'dummy'] += 1

        if not aln.is_paired:
            barcodes[aln.get_tag('RG') if use_group else 'dummy'] += 1

    afile.close()

    names = [key for key in barcodes]
    counts = [barcodes[key] for key in barcodes]

    df = pd.DataFrame({'barcodes':names, 'counts':counts})

    df.to_csv(storage, sep='\t', header=True, index=False)


if __name__ == '__main__':
  regions = '../../scipipe_output/hg19/countmatrix/atac_bin_1000.bed'
  in_ = '../../scipipe_output/hg19/atac.barcoded.dedup.bam'
  out = '../../scipipe_output/hg19/countmatrix/atac_bin_1000.tab'
  sparse_count_reads_in_regions(in_, regions, out, flank=0)
