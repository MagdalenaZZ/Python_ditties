#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Counts the number of sequences in a fasta/q file',
    usage = '%(prog)s <fasta/q in>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
options = parser.parse_args()
print(fastn.count_sequences(options.infile))
