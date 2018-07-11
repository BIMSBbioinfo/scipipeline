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

    df = pd.read_csv(input_table, header=None, sep='\t')
    record_list = []
    for i, item in enumerate(df[0]):

        # store them in a pseudo-ref. genome in fasta format
        record_list.append(SeqRecord(Seq(item, IUPACAmbiguousDNA),
                                     id='{}-{}'.format(os.path.basename(input_table).split('.')[0], i))) 

    SeqIO.write(record_list, output_fasta, "fasta")
