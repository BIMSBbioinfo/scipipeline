import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pysam import AlignmentFile
from pysam import view


def split_reads_by_barcode(barcode_bams, treatment_bam, 
                           output_dir, max_open_files=1000, min_mapq=None, max_mismatches=1):
    """Splits the reads by barcodes.

    This function takes a set of barcode alignments (in bam format)
    and reads (in bam format). It produces an output bam file
    whose name corresponds to the barcodes (one file per barcode).

    Note: The method assumes that the respective BAM files contain reads
    with corresponding read ids. Also the BAM files need to be sorted
    by read name. Furthermore, the treatment_reads must be a subset of
    the barcode reads.

    barcodes : list(str)
        List of bam-file names pointing to the barcode alignments
        towards the pseudo genomes.
    reads : str
        BAM file containing the read alignments
    output_dir : str
        Output directory in which the split up bam files will be
        produced.
    max_open_files : int
        Maximum number of writers to handle in one sweep. This is limited by the
        operating system's limitation opening file handles. Default: 1000.
    min_mapq : int
        Filter for mapping quality of the barcode reads. Only reads mapq >= min_mapq
        are considered. Default: None means no filter is applied.
    max_mismatches : int
        Maximum number of mismatches. Default: 1.
    """

    os.makedirs(output_dir, exist_ok=True)
    shutil.copy(treatment_bam, treatment_bam + '.remaining.bam')

    aligned_cnt = 0
    unaligned_cnt = 0
    current_writers = dict()
    max_reached = True
    if not min_mapq:
        min_mapq = 0

    while max_reached:

        # reset the current writers
        current_writers = dict()
        max_reached = False

        treatment_reader = AlignmentFile(treatment_bam + '.remaining.bam', 
                                         'rb')

        barcode_readers = [AlignmentFile(bamfile, 'rb') 
                           for bamfile in barcode_bams]

        # start (again) at the beginning of the bam files
        barcode_it = [reader.fetch(until_eof=True) 
                      for reader in barcode_readers]
        bnames = [next(br) for br in barcode_it]
        tmp_writer = AlignmentFile(os.path.join(output_dir, 'tmp.bam'), 'wb',
                                   template=treatment_reader)

        print("Start batch ... ")
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

            if len(current_writers) >= max_open_files:
                # If max open files reached, do not open
                # any further bam_writers, but finish processing
                # of the current batch
                max_reached = True

            if not comb_id in current_writers and not max_reached:
                # instantiate a new writer for the barcode
                # if it has not already been.
                writer = AlignmentFile(os.path.join(output_dir, 
                                                    comb_id + '.bam'), 
                                       'wb', template=treatment_reader)
                current_writers[comb_id] = writer

            if comb_id in current_writers:
                # append the current read
                current_writers[comb_id].write(aln)
            else:
                tmp_writer.write(aln)

        print("End batch ...")
        # close all remaining bam files
        for writer in current_writers:
            current_writers[writer].close()
        print("Closed writers")

        tmp_writer.close()
        print("Closed tmp writers")
        treatment_reader.close()
        print("Closed treatment reader")
        for reader in barcode_readers:
            reader.close()

        print("Closed barcode reader")
        os.rename(os.path.join(output_dir, 'tmp.bam'),
                  treatment_bam + '.remaining.bam')
        print("Move tmp.bam to treatment.remaining.bam")

    os.remove(treatment_bam + '.remaining.bam')
    print("Split {} reads. {} unaligned reads were ignored.".format(
        aligned_cnt, 
        unaligned_cnt))


def obtain_barcode_frequencies(originals, dedup, output):

    with open(output, 'w') as f:
        f.write("file\toriginal\tdeduplicatedi\n")
        for ofile, defile in zip(originals, dedup):
            ocnt = int(view(ofile, '-c'))
            dcnt = int(view(defile, '-c'))
            f.write("{}\t{}\t{}\n".format(os.path.basename(ofile).split('.')[0], ocnt, dcnt))


def plot_barcode_frequencies(tab_file, plotname):
    x=pd.read_csv(tab_file,sep='\t')
    f = plt.figure()
    plt.plot(np.log10(x['deduplicated'].sort_values(ascending=False).values//2))
    plt.ylabel('Log10(# pairs)')
    plt.xlabel('Barcodes')
    plt.title('Barcode frequency (deduplicated)')
    f.savefig(plotname, dpi=f.dpi)
    

def species_specificity(aln_file1, aln_file2, output_txt):
    aln1 = AlignmentFile(aln_file1, 'r').fetch(until_eof=True)
    aln2 = AlignmentFile(aln_file2, 'r').fetch(until_eof=True)

    cnt = np.zeros((2,2))
    try:
        while 1:
            m1 = not next(aln1).is_unmapped
            m2 = not next(aln2).is_unmapped
            cnt[0 if m1 else 1, 0 if m2 else 1] += 1
    finally:
        np.savetxt(output_txt, cnt, delimiter='\t')


if __name__ == '__main__':
    barcode_bams = ['pseudo_genome_I{}_sorted.bam'.format(i) for i in [1,2,3,4]]
    treatment_bam = 'test.input.bam'
    output_dir = 'test_split'
    max_open_files=1000
    split_reads_by_barcode(barcode_bams, treatment_bam, output_dir, max_open_files)
