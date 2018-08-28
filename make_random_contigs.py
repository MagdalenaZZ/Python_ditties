#!/usr/bin/env python3.3

import sys
import argparse
import random
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Makes a set of contigs of random bases in FASTA format',
    usage = '%(prog)s [options] <number of contigs> <contig length> <outfile>')
parser.add_argument('--first_number', type=int, help='If numbering the sequences, the first sequence gets this number [%(default)s]', default=1)
parser.add_argument('--name_by_letters', action='store_true', help='Name the contigs A,B,C,... will start at A again if you get to Z')
parser.add_argument('--prefix', help='Prefix to add to start of every sequence name', default='')
parser.add_argument('--seed', type=int, help='Seed for random number generator. Default is to use python\'s default', default=None)
parser.add_argument('contigs', type=int, help='Nunber of contigs to make')
parser.add_argument('length', type=int, help='Length of each contig')
parser.add_argument('outfile', help='Name of output file')
options = parser.parse_args()

random.seed(a=options.seed)

fout = utils.open_file_write(options.outfile)
letters = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
letters_index = 0

for i in range(options.contigs):
    if options.name_by_letters:
        name = letters[letters_index]
        letters_index += 1
        if letters_index == len(letters):
            letters_index = 0
    else:
        name = str(i + options.first_number)

    fa = fastn.Fasta(options.prefix + name, ''.join([random.choice('ACGT') for x in range(options.length)]))
    print(fa, file=fout)

utils.close(fout)

