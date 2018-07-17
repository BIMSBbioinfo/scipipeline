import numpy
from scipy.sparse import dok_matrix, coo_matrix
from pysam import AlignmentFile
from HTSeq import BED_Reader
import h5py

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
        nbins = genomesize[chrom]//binsize
        starts = [int(i*binsize) for i in range(nbins)]
        ends = [int((i+1)*binsize) for i in range(nbins)]
        chr_ = [chrom] * nbins
        cont = {'chr': chr_, 'start': starts, 'end': ends}
        bed_entry = pd.DataFrame(cont)
        bed_content = bed_content.append(bed_entry, ignore_index=True, sort=False)
    bed_content.to_csv(storage, sep='\t', header=False, index=False,
                       columns=['chr', 'start', 'end'])




def sparse_count_reads_in_regions(bamfile, regions, storage, flank=0,
                                  template_length=1000):
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
    """

    # Obtain the header information
    afile = AlignmentFile(bamfile, 'rb')

    # extract genome lengths
    print('get genomesize')
    # extract genome size
    genomesize = {}
    for chrom, length in zip(afile.references, afile.lengths):
        genomesize[chrom] = length
    print('found {} chromosomes'.format(len(genomesize)))

    nreg = 0
    regfile = BED_Reader(regions)
    for reg in regfile:
        nreg += 1

    if 'RG' in afile.header:
        use_group = True
    else:
        use_group = False

    barcodes = {}
    if use_group:
        # extract barcodes
        for idx, item in enumerate(afile.header['RG']):
            barcodes[item['ID']] = idx
    else:
        barcodes['dummy'] = 0
    print('found {} barcodes'.format(len(barcodes)))

    sdokmat = dok_matrix((nreg, len(barcodes)), dtype='int32')

    tlen = 200
    for idx, region in enumerate(regfile):
        iv = region.iv.copy()
        iv.start -= flank + tlen
        iv.end += flank + tlen

        if iv.start < tlen:
            iv.start = tlen
        if iv.end >= genomesize[iv.chrom] - tlen:
            iv.end = genomesize[iv.chrom] - 1 - tlen
        for aln in afile.fetch(iv.chrom, iv.start - tlen, iv.end + tlen):
            if aln.is_paired and aln.is_read1:
                # if paired end and first mate
                if aln.pos + aln.template_length//2 >= iv.start and \
                   aln.pos + aln.template_length//2 < iv.end:
                   sdokmat[idx, barcodes[aln.get_tag('RG') if use_group else 'dummy']] += 1
            if not aln.is_paired:
                if not aln.is_reverse:
                    if aln.pos >= iv.start and aln.pos < iv.end:
                        sdokmat[idx,
                        barcodes[aln.get_tag('RG') if use_group else 'dummy']] += 1
                else:
                    if aln.pos + aln.reference_length - 1 >= iv.start and \
                       aln.pos + aln.reference_length - 1 < iv.end:
                        sdokmat[idx,
                        barcodes[aln.get_tag('RG') if use_group else 'dummy']] += 1

    afile.close()

    print('sparse matrix shape: {}'.format(sdokmat.shape))
    print('density: {}'.format(sdokmat.nnz/numpy.prod(sdokmat.shape)))

    # store the results in COO sparse matrix format
    scoomat = sdokmat.tocoo()
    # sort lexicographically

    order_ = np.lexsort((spcoo.col, spcoo.row))
    indices = np.asarray([x for x in zip(spcoo.row, spcoo.col)], dtype=np.int64)[order_]
    values = spcoo.data.astype(np.float32)[order_]
#    shape = np.asarray(spcoo.shape, dtype=np.int64)
    cont = {'region': indices[:,0], 'cell': indices[:, 1], 'count': values}

    df = pd.DataFrame(cont)
    df.to_csv(storage, sep='\t', header=True, index=False)


def per_barcode_count_summary(cnt_mat, peak_counts, storage):
    df_cnt = pd.read_csv(cnt_mat, header=[0], sep='\t')
    df_peak = pd.read_csv(peak_counts, header=[0], sep='\t')

    per_bc_count = df_cnt.groupby('cell').agg('sum')['count']
    per_bc_count_in_peaks = df_peak.groupby('cell').agg('sum')['count']
    per_bc_count['count_in_peaks'] = per_bc_count_in_peaks['count']
    per_bc_count['percent_in_peaks'] = per_bc_count.apply(lambda row: row['count_in_peaks']/row['count'])

    per_bc_count.to_csv(storage, sep='\t')
