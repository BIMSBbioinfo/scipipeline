"""This script is supposed to generate the
pseudo genome for the barcodes.
"""
import itertools
import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


def create_pseudo_genome(input_table, output_fasta):
    """ Create a fasta barcode reference from a table.

    The input table is expected to hold the reference barcodes
    in the first column of the tab-separated file. It could optionally
    hold additional columns which will be merged into the sequence id
    if present.

    For example, the input table

    > cat barcodes.tsv
    TCTCGCGCCATAGCGACCTTCGCCTCAAGATTAAG    tag1
    GAATTCGTATACCTCGACTCAGTACGATACCTAAG    tag1
    GAGATTCCTTGATTGGCGGGCTATAAAACTTATAA    tag2
    ...

    will yield the following output fasta file:

    > cat output.fasta
    >barcodes-tag1-1
    TCTCGCGCCATAGCGACCTTCGCCTCAAGATTAAG
    >barcodes-tag1-2
    GAATTCGTATACCTCGACTCAGTACGATACCTAAG
    >barcodes-tag2-3
    GAGATTCCTTGATTGGCGGGCTATAAAACTTATAA
    ...

    """

    df = pd.read_csv(input_table, header=None, sep='\t')
    record_list = []
    for i, row in df.iterrows():
        # store them in a pseudo-ref. genome in fasta format
        record_list.append(SeqRecord(Seq(row[0], IUPACAmbiguousDNA),

                                     id='{}-{}-{}'.format(
                                        os.path.basename(input_table).split('.')[0],
                                        '-'.join([str(x) for x in row[1:]]),
                                        i)))

    SeqIO.write(record_list, output_fasta, "fasta")
