#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Gets all IDs from a fasta or fastq file',
    usage = '%(prog)s <infile> <outfile>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()

seq_reader = fastn.file_reader(options.infile)
f_out = utils.open_file_write(options.outfile)

for seq in seq_reader:
    print(seq.id, file=f_out)

utils.close(f_out)
