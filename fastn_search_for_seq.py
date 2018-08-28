#!/usr/bin/env python3.3

import argparse
import fastn
import utils

parser = argparse.ArgumentParser(
    description = 'Searches for an exact match on a given string, and its reverese complement, in a fasta/q file. Case insensitive.',
    usage = '%(prog)s [options] <fasta/q in> <search_string>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('search_string', help='String to search for in the file')
options = parser.parse_args()


seq_reader = fastn.file_reader(options.infile)

for seq in seq_reader:
    hits = seq.search(options.search_string)
    for hit in hits:
        print(seq.id, hit[0]+1, hit[1], sep='\t')
