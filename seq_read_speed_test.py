#!/usr/bin/env python3.3

import sys
import argparse
import fastn

parser = argparse.ArgumentParser(
    description = 'Just reads a fasta or fastq file and print number of reads. Used for testing',\
    usage = '%(prog)s [options] <fasta/q in>')
#parser.disable_interspersed_args()
parser.add_argument('infile', help='Name of fasta/q file to be read')
options = parser.parse_args()


seq_reader = fastn.file_reader(options.infile)

counter = 0
for seq in seq_reader:
    counter += 1

print ('read', counter, 'sequences from', options.infile)

