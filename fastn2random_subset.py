#!/usr/bin/env python3.3

import sys
import argparse
import random
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Takes a random subset of reads from a fasta/q file and optionally the corresponding read ' +
                  'from a mates file.  Ouptut is interleaved if mates file given',
    usage = '%(prog)s [options] <fasta/q in> <outfile> <percent reads wanted in [0,100]>')
parser.add_argument('--mate_file', help='Name of fasta/q mates file')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('outfile', help='Name of fasta/q output file')
parser.add_argument('read_percent', type=int, help='percent of reads to take from input file')
options = parser.parse_args()

seq_reader = fastn.file_reader(options.infile)
fout = utils.open_file_write(options.outfile)
counter_in = 0
counter_out = 0

if options.mate_file:
    mate_seq_reader = fastn.file_reader(options.mate_file)

for seq in seq_reader:
    counter_in += 1
    if options.mate_file:
        try:
            mate_seq = next(mate_seq_reader)
        except StopIteration:
            print('Error! Didn\'t get mate for read', seq.id, file=sys.stderr)
            sys.exit(1)
    if random.randint(0, 100) <= options.read_percent:
        counter_out += 1
        print(seq, file=fout)
        if options.mate_file:
            print(mate_seq, file=fout)


utils.close(fout)

if options.mate_file:
    print('Used', counter_out, 'pairs from total of', counter_in)
else:
    print('Used', counter_out, 'sequences from total of', counter_in)
