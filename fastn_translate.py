#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Given a fasta and qual file, makes a fastq file',
    usage = '%(prog)s <in.fasta/q> <out.fasta/q>')
parser.add_argument('infile', help='Name of fasta/q file to be translated', metavar='in.fasta/q')
parser.add_argument('outfile', help='Name of output fasta/q file', metavar='out.fasta/q')
options = parser.parse_args()


seq_reader = fastn.file_reader(options.infile)
f_out = utils.open_file_write(options.outfile)

for seq in seq_reader:
    print(seq.translate(), file=f_out)

utils.close(f_out)


