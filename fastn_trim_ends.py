#!/usr/bin/env python3.3

import argparse
import fastn

parser = argparse.ArgumentParser(
    description = 'Trims ends off each sequence in a fasta/q file',
    usage = '%(prog)s [options] <fasta/q in> <bases off start> <bases off end> <fasta/q out>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('start_trim', type=int, help='Number of bases to trim off start')
parser.add_argument('end_trim', type=int, help='Number of bases to trim off end')
parser.add_argument('outfile', help='Name of output fasta/q file')
options = parser.parse_args()
fastn.trim(options.infile, options.outfile, options.start_trim, options.end_trim)
