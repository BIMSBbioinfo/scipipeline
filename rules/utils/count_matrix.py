import numpy
from pysam import AlignmentFile
from HTSeq import BED_Reader
import h5py

def count_reads_in_bins(bamfile, binsize, storage, merged_storage=None):
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
    binsize : int
        Binsize used to bin the genome in consecutive bins.
    storage : str
        Path to the output hdf5 file, which contains the counts per chromsome.
    merged_storage : str
        Path to the output hdf5 file, which contains the counts 
        across the whole genome in a single dataset.
    """

    # Obtain the header information
    afile = AlignmentFile(bamfile, 'rb')

    print('get genomesize')
    # extract genome size
    genomesize = {}
    for chrom, length in zip(afile.references, afile.lengths):
        genomesize[chrom] = length
    print('found {} chromosomes'.format(len(genomesize)))

    if 'RG' in afile.header:
        use_group = True
    else:
        use_group = False

    print('get barcodes')
    barcodes = {}
    if use_group:
        # extract barcodes
        for idx, item in enumerate(afile.header['RG']):
            barcodes[item['ID']] = idx
    else:
        barcodes['dummy'] = 0
    print('found {} barcodes'.format(len(barcodes)))

    # open up a hdf5 file to dump the dataset in
    data = h5py.File(storage, 'w', driver='core')
    for chrom in genomesize:
        data.create_dataset(chrom, (genomesize[chrom]//binsize + 1, len(barcodes)),
                            dtype='int32', compression='gzip',
                            data=numpy.zeros((genomesize[chrom]//binsize + 1, len(barcodes))))

    data.create_dataset('barcode_name', data=[ numpy.string_(barcode) for barcode in barcodes])
    data.create_dataset('barcode_index', data=[ barcodes[barcode] for barcode in barcodes])

    for aln in afile:
        if aln.is_paired and aln.is_read2:
            continue
        if aln.is_paired:
            # if paired end and first mate
            idx = (aln.pos + aln.template_length//2)//binsize
        else:
            if not aln.is_reverse:
                idx = (aln.reference_start)//binsize
            else:
                idx = (aln.reference_start + aln.reference_length-1)//binsize

        data[aln.reference_name][
             idx, barcodes[aln.get_tag('RG') if use_group else 'dummy']] += 1

    afile.close()
    data.close()


def count_reads_in_regions(bamfile, regions, storage, flank=0):
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

    # open up a hdf5 file to dump the dataset in
    f = h5py.File(storage, 'w', driver='core')
    data = f.create_dataset('regions', (nreg, len(barcodes)),
                            dtype='int32', compression='gzip',
                            data=numpy.zeros((nreg, len(barcodes))))

    for barcode in barcodes:
        f.attrs[barcode] = barcodes[barcode]

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
                   data[idx, 
                   barcodes[aln.get_tag('RG') if use_group else 'dummy']] += 1
            if not aln.is_paired:
                if not aln.is_reverse:
                    if aln.pos >= iv.start and aln.pos < iv.end:
                        data[idx, 
                        barcodes[aln.get_tag('RG') if use_group else 'dummy']] += 1
                else:
                    if aln.pos + aln.reference_length - 1 >= iv.start and \
                       aln.pos + aln.reference_length - 1 < iv.end:
                        data[idx, 
                        barcodes[aln.get_tag('RG') if use_group else 'dummy']] += 1

    afile.close()
    f.close()


if __name__ == '__main__':
    #import cProfile
    from os.path import join
    bams = join('..', '..', '..', 'sciatac_data', 
                "dm6", "atac_barcode_dedup.bam")
    
    output = join('..', '..', '..', 'sciatac_data', 
                "dm6", "cnt_test.h5")
    #cProfile.run('count_reads_in_bins(bams, 10000, output)')
    count_reads_in_bins(bams, 10000, output)
