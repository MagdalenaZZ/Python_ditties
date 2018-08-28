#!/usr/bin/env python3.3

import sys
import argparse
import fastn
import utils
import re

parser = argparse.ArgumentParser(
    description = 'Takes every sequence from a fasta/q file whose ID matches the given regular expression',
    usage = '%(prog)s [options] <infile> <regex> <outfile>')
parser.add_argument('-v,--invert_match', dest='invert_match', action='store_true', help='Use an ID if it doesn\'t match the regex (like grep -v)')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('regex', help='Regular expression - all IDs matching this will be kept')
parser.add_argument('outfile', help='Name of fasta/q output file')
options = parser.parse_args()

p = re.compile(options.regex)

seq_reader = fastn.file_reader(options.infile)
fout = utils.open_file_write(options.outfile)



for seq in seq_reader:
    hit = p.search(seq.id) is not None
    if hit != options.invert_match:
        print(seq, file=fout)

utils.close(fout)

