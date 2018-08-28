#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Replaces all occurences of one letter with another in a fasta/q file',
    usage = '%(prog)s <fasta/q in> <outfile> <old> <new>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of output file')
parser.add_argument('old', help='Base to be replaced')
parser.add_argument('new', help='Replace with this letter')
options = parser.parse_args()
fastn.replace_bases(options.infile, options.outfile, options.old, options.new)
