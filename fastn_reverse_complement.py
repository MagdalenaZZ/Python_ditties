#!/usr/bin/env python3.3

import argparse
import fastn

parser = argparse.ArgumentParser(
    description = 'Reverse complements a fasta/q file',
    usage = '%(prog)s [options] <fasta/q in> <fasta/q out>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output fasta/q file')
options = parser.parse_args()
fastn.reverse_complement(options.infile, options.outfile)
