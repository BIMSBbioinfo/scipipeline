import os
import shutil
from pysam import AlignmentFile


def split_reads_by_barcode(barcode_bams, treatment_bam, 
                           output_dir, max_open_files):
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
    shutil.copy(treatment_bam, treatment_bam + '.remaining.bam')

    aligned_cnt = 0
    unaligned_cnt = 0
    current_writers = dict()
    max_reached = True

    while max_reached:

        # reset the current writers
        current_writers = dict()
        max_reached = False

        treatment_reader = AlignmentFile(treatment_bam + '.remaining.bam', 'rb')

        barcode_readers = [AlignmentFile(bamfile, 'rb') for bamfile in barcode_bams]

        # start (again) at the beginning of the bam files
        barcode_it = [reader.fetch(until_eof=True) for reader in barcode_readers]
        bnames = [next(br) for br in barcode_it]
        tmp_writer = AlignmentFile(os.path.join(output_dir, 'tmp.bam'), 'wb',
                                   template=treatment_reader)

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
                        
            if any([ x.is_unmapped for x in bnames]):
               # one or more barcodes are abscent
               continue

            comb_id = '_'.join([baln.reference_name for baln in bnames])

            if len(current_writers) >= max_open_files:
                # If max open files reached, do not open
                # any further bam_writers, but finish processing
                # of the current batch
                max_reached = True

            if not comb_id in current_writers and not max_reached:
                # instantiate a new writer for the barcode
                # if it has not already been.
                writer = AlignmentFile(os.path.join(output_dir, comb_id + '.bam'), 
                    'wb', template=treatment_reader)
                current_writers[comb_id] = writer

            if comb_id in current_writers:
                # append the current read
                current_writers[comb_id].write(aln)
            else:
                tmp_writer.write(aln)

        # close all remaining bam files
        for writer in current_writers:
            current_writers[writer].close()

        tmp_writer.close()
        treatment_reader.close()
        for reader in barcode_readers:
            reader.close()

        os.rename(os.path.join(output_dir, 'tmp.bam'),
                  treatment_bam + '.remaining.bam')

    print("Split {} reads. {} unaligned reads were ignored.".format(
        aligned_cnt, 
        unaligned_cnt))


if __name__ == '__main__':
    barcode_bams = ['pseudo_genome_I{}_sorted.bam'.format(i) for i in [1,2,3,4]]
    treatment_bam = 'test.input.bam'
    output_dir = 'test_split'
    max_open_files=1000
    split_reads_by_barcode(barcode_bams, treatment_bam, output_dir, max_open_files)