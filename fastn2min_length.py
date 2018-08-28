#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Extracts all sequences of at least a given length from a fasta/q file',
    usage = '%(prog)s <min_length> <infile> <outfile>')
parser.add_argument('min_length', type=int, help='Minimum length of sequence to take')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output fasta/q file')
options = parser.parse_args()

seq_reader = fastn.file_reader(options.infile)
f_out = utils.open_file_write(options.outfile)

for seq in seq_reader:
    if len(seq) >= options.min_length:
        print(seq, file=f_out)

utils.close(f_out)

