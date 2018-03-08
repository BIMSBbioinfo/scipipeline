import os
from HTSeq import BAM_Reader
from HTSeq import BAM_Writer


def split_reads_by_barcode(barcodes, reads, output_dir, max_open_files):
    """Splits the reads by barcodes.

    This function takes a set of barcode alignments (in bam format)
    and reads (in bam format). It produces an output bam file
    whose name corresponds to the barcodes (one file per barcode).

    Note: The method assumes that the respective BAM files contain reads
    with corresponding read ids. Also the BAM files need to be sorted
    by read name.

    barcodes : list(str)
        List of bam-file names pointing to the barcode alignments
        towards the pseudo genomes.
    reads : str
        BAM file containing the read alignments
    output_dir : str
        Output directory in which the split up bam files will be
        produced.
    """

    os.makedirs(output_dir, exist_ok=True)

    barcode_readers = [iter(BAM_Reader(bamfile)) for bamfile in barcodes]
    treatment_reader = BAM_Reader(reads)
    output_writers = {}
    batch_id = 0

    bnames = [next(br) for br in barcode_readers]
    aligned_cnt = 0
    unaligned_cnt = 0

    for aln in treatment_reader:
        # only retain the aligned reads
        if not aln.aligned:
            unaligned_cnt += 1
            continue
        aligned_cnt += 1
        # extract the corresponding barcodes
        for i, br_it in enumerate(barcode_readers):
            while not bnames[i].read.name == aln.read.name:
                # increment barcode names until they
                # match the alignment read name.
                bnames[i] = next(br_it)
                        
        if any([not x.aligned for x in bnames]):
           # one or more barcodes are abscent
           continue

        comb_id = '_'.join([baln.iv.chrom for baln in bnames])
        if not comb_id in output_writers:
            # instantiate a new writer for the barcode
            # if it has not already been.
            writer = BAM_Writer.from_BAM_Reader(os.path.join(output_dir, 
                                                             comb_id + '.{}.bam'.format(batch_id)), 
                                                treatment_reader)
            output_writers[comb_id] = writer

        # append the current read
        output_writers[comb_id].write(aln)

        if len(output_writers) >= max_open_files:
            # if max file size is reached, close all files
            # and start with a new batch
            for writer in output_writers:
                output_writers[writer].close()
            output_writers = dict()
            batch_id += 1

    for writer in output_writers:
        output_writers[writer].close()

    print("Split {} reads. {} unaligned reads were ignored.".format(aligned_cnt, unaligned_cnt))
