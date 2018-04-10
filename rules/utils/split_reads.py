import os
import shutil
import h5py
import numpy as np
from pysam import AlignmentFile
from pysam import view
import deeptools.countReadsPerBin as crpb


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
    bam_writer = AlignmentFile(output_bam, 'wb', template=treatment_reader)

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
    bam_writer.header['RG'] = ['ID': combid for combid in barcodes]

    treatment_reader.close()
    bam_writer.close()
    for reader in barcode_readers:
        reader.close()


def obtain_barcode_frequencies(originals, dedup, output):

    with open(output, 'w') as f:
        f.write("file\toriginal\tdeduplicated\n")
        for ofile, defile in zip(originals, dedup):
            ocnt = int(view(ofile, '-c'))
            dcnt = int(view(defile, '-c'))
            f.write("{}\t{}\t{}\n".format(os.path.basename(ofile).split('.')[0], ocnt, dcnt))


def count_reads_in_bins(bams, binsize, storage):
    """ This function obtains the counts per bin """

    # Obtain the header information
    afile = AlignmentFile(bams[0], 'rb')
    genomesize = {}
    for chrom, length in zip(afile.references, afile.lengths):
        genomesize[chrom] = length
    afile.close()

    # number of batches
    batches = len(bams)//1000 + 1 if len(bams)%1000 > 0 else 0

    # open up a hdf5 file to dump the dataset in
    data = h5py.File(storage, 'w')
    for chrom in genomesize:
        data.create_dataset(chrom, (genomesize[chrom]//binsize, len(bams)),
                            dtype='int32', compression='lzf')

    for b in range(batches):
        cr = crpb.CountReadsPerBin(bams[b*1000:(b+1)*1000], 
                                   binLength=binsize, 
                                   stepSize=binsize, 
                                   center_read=True, samFlag_exclude = 128)
        for chrom in genomesize:
            data[chrom][:,b*1000:(b+1)*1000] = cr.count_reads_in_region(chrom, 0, genomesize[chrom] -1)[0]
    data.close()


if __name__ == '__main__':
    barcode_bams = ['pseudo_genome_I{}_sorted.bam'.format(i) for i in [1,2,3,4]]
    treatment_bam = 'test.input.bam'
    output_dir = 'test_split'
    max_open_files=1000
    split_reads_by_barcode(barcode_bams, treatment_bam, output_dir, max_open_files)
