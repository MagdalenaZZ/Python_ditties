#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Strips /1 or /2 off the end of every read name in a fasta/q file',
    usage = '%(prog)s [options] <fasta/q in> <fasta/q out>')
parser.add_argument('infile', help='Name of input fasta/q file')
parser.add_argument('outfile', help='Name of output fasta/q file')
options = parser.parse_args()
seq_reader = fastn.file_reader(options.infile)
f = utils.open_file_write(options.outfile)

for seq in seq_reader:
    seq.strip_illumina_suffix()
    print(seq, file=f)

utils.close(f)
