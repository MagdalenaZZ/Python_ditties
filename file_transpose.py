#!/usr/bin/env python3.3

import sys
import argparse
import utils

parser = argparse.ArgumentParser(
    description = 'Transposes a file (like a matrix transpose.  Columns/rows swapped)',\
    usage = '%(prog)s [options] <infile> <outfile>')
parser.add_argument('--sep_in', help='Separator used to delimit columns of input file [all whitespace]', default=None)
parser.add_argument('--sep_out', help='Separator used to delimit columns of output file [tab]', default='\t')
parser.add_argument('infile', help='Name of input file')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()

utils.file_transpose(options.infile, options.outfile, sep_in=options.sep_in, sep_out=options.sep_out)
