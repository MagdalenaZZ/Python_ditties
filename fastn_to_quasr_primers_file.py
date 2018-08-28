#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Converts a fasta/q file to QUASR primers format: just the sequence on each line and its reverse complement, tab separated',
    usage = '%(prog)s <fasta/q in> <outfile>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()
fastn.fastn_to_quasr_primers(options.infile, options.outfile)
