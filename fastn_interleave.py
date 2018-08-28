#!/usr/bin/env python3.3

import argparse
import fastn

parser = argparse.ArgumentParser(
    description = 'Interleaves two fasta/q files, so that reads are written alternately first/second in output file',
    usage = '%(prog)s [options] <fasta/q 1> <fasta/q 2> <outfile>')
parser.add_argument('infile_1', help='Name of first fasta/q file to be read')
parser.add_argument('infile_2', help='Name of second fasta/q file to be read')
parser.add_argument('outfile', help='Name of output fasta/q file of interleaved reads')
options = parser.parse_args()
fastn.interleave(options.infile_1, options.infile_2, options.outfile)
