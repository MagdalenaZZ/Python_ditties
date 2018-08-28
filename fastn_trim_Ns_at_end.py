#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Trims any Ns off each sequence in a fasta/q file. Does nothing to gaps in the middle, just trims the ends',
    usage = '%(prog)s [options] <fasta/q in> <fasta/q out>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output fasta/q file')
options = parser.parse_args()


seq_reader = fastn.file_reader(options.infile)
fout = utils.open_file_write(options.outfile)

for seq in seq_reader:
    seq.trim_Ns()
    if len(seq):
        print(seq, file=fout)
